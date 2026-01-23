import sqlite3
import argparse
import sys
import os
from collections import defaultdict
import time

K = 10

def parse_fasta_stream(file_path):
    """
    Yields chunks of sequence from a FASTA file, handling line breaks.
    """
    with open(file_path, 'r') as f:
        buffer = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # New record, reset buffer
                buffer = ""
                continue
            
            # It's sequence data
            # Check for Ns? We can filter later.
            seq_part = line.upper()
            
            # If we have a buffer from previous line, combine
            if buffer:
                combined = buffer + seq_part
                yield combined
                # Keep the last K-1 bases for the next overlap
                buffer = combined[-(K-1):] if len(combined) >= K-1 else combined
            else:
                if len(seq_part) < K:
                    buffer = seq_part
                else:
                    yield seq_part
                    buffer = seq_part[-(K-1):]

def build_database(fasta_path, db_path):
    print(f"Building database from {fasta_path}...")
    start_time = time.time()
    
    # In-memory counter
    # 4^10 is ~1M keys. This easily fits in RAM.
    counts = defaultdict(int)
    
    total_kmers_processed = 0
    
    try:
        # We need a robust reader that handles cross-line k-mers correctly.
        # The stream generator approach above yields chunks, but we need to be careful not to double count.
        # Modified approach: 
        # Iterate file. Keep a running string buffer.
        # Process k-mers in buffer. 
        # Keep overlapping tail.
        
        with open(fasta_path, 'r') as f:
            overlap = ""
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    overlap = "" # Reset on new sequence
                    continue
                
                seq = line.upper()
                
                # Filter invalid chars if necessary, but assuming standard FASTA
                # If N is present, we might want to split by N to avoid N-containing kmers
                
                full_seq = overlap + seq
                
                # Split by 'N' so we don't count invalid k-mers
                subsequences = full_seq.replace('N', ' ').split()
                
                for sub in subsequences:
                    if len(sub) < K:
                        continue
                    
                    # Sliding window
                    for i in range(len(sub) - K + 1):
                        kmer = sub[i : i+K]
                        counts[kmer] += 1
                        total_kmers_processed += 1
                
                # Prepare overlap for next line
                # We only need the last K-1 bases of the RAW line (before N splitting logic?)
                # Actually, the biological sequence continues across lines.
                # If the line ends with N, the overlap is effectively useless for the next line's start 
                # (unless we skip Ns but treat them as gaps). 
                # Standard approach: N breaks the chain.
                # So we just take the last K-1 chars of the full_seq to be the next overlap.
                if len(full_seq) >= K-1:
                    overlap = full_seq[-(K-1):]
                else:
                    overlap = full_seq
        
        print(f"Counting complete. Processed {total_kmers_processed} k-mers.")
        print(f"Unique k-mers found: {len(counts)}")
        
        # Save to SQLite
        print(f"Saving to SQLite database: {db_path}...")
        conn = sqlite3.connect(db_path)
        c = conn.cursor()
        c.execute('DROP TABLE IF EXISTS kmer_counts')
        c.execute('CREATE TABLE kmer_counts (kmer TEXT PRIMARY KEY, count INTEGER)')
        
        # Batch insert
        c.executemany('INSERT INTO kmer_counts VALUES (?, ?)', counts.items())
        
        conn.commit()
        
        # Create an index on count for filtering?
        # Not strictly necessary for PK lookups, but good if we want "find all kmers with frequency > N"
        c.execute('CREATE INDEX IF NOT EXISTS idx_count ON kmer_counts (count)')
        
        conn.close()
        print(f"Database built successfully in {time.time() - start_time:.2f} seconds.")

    except FileNotFoundError:
        print(f"Error: File {fasta_path} not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def generate_neighbors(kmer):
    """Generate all k-mers with hamming distance 1."""
    bases = ['A', 'C', 'G', 'T']
    neighbors = []
    for i in range(len(kmer)):
        original_base = kmer[i]
        for base in bases:
            if base != original_base:
                neighbor = kmer[:i] + base + kmer[i+1:]
                neighbors.append(neighbor)
    return neighbors

def map_locations(fasta_path, db_path, min_repeats=2):
    print(f"Mapping locations from {fasta_path} to {db_path}...")
    print(f"Only indexing k-mers appearing >= {min_repeats} times.")
    
    if not os.path.exists(db_path):
        print("Error: Database does not exist. Run 'build' first.")
        return

    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    # Create locations table
    c.execute('''CREATE TABLE IF NOT EXISTS kmer_locations 
                 (kmer TEXT, chrom TEXT, pos INTEGER)''')
    c.execute('CREATE INDEX IF NOT EXISTS idx_kmer_loc ON kmer_locations (kmer)')
    
    # Pre-load counts to determine which k-mers to index?
    # Loading 1M+ keys into a set is fine for RAM.
    print("Loading target k-mer set (this might take a moment)...")
    c.execute("SELECT kmer FROM kmer_counts WHERE count >= ?", (min_repeats,))
    # We use a set for O(1) lookups
    target_kmers = set(row[0] for row in c.fetchall())
    print(f"Indexing locations for {len(target_kmers)} unique k-mers.")
    
    batch_data = []
    batch_size = 100000
    total_indexed = 0
    
    with open(fasta_path, 'r') as f:
        current_chrom = "unknown"
        overlap = ""
        # We need to track absolute position within the chromosome.
        # Line-based reading makes this slightly tricky because of newlines.
        # We'll maintain a 1-based index `chrom_pos_start` for the start of the current line.
        
        chunk_idx = 1 # 1-based genomic coordinate
        
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                current_chrom = line[1:].split()[0] # Take first word as ID
                overlap = "" 
                chunk_idx = 1 # Reset position for new chromosome
                continue
            
            seq = line.upper()
            full_seq = overlap + seq
            
            # Position calculation:
            # The 'full_seq' starts at (chunk_idx - len(overlap)) effectively?
            # No, 'overlap' is the tail of the PREVIOUS line.
            # So the new 'seq' starts at chunk_idx.
            # The k-mers start from index 0 of full_seq.
            # If `i` is index in full_seq:
            #   If i < len(overlap): It started in previous line (but we only process it now if we missed it? 
            #   actually we process *sliding window*.
            #   Correct logic:
            #   The first base of `full_seq` corresponds to coordinate: `chunk_idx - len(overlap)`
            
            start_coord = chunk_idx - len(overlap)
            
            # We want to iterate k-mers in full_seq
            # But we must avoid re-counting k-mers that were fully in 'overlap' (already processed).
            # Actually, in the 'build' phase we just processed everything.
            # Here we must be precise.
            # A k-mer starts at index `i` in `full_seq`.
            # Its genomic position is `start_coord + i`.
            
            for i in range(len(full_seq) - K + 1):
                kmer = full_seq[i : i+K]
                
                # Filter N
                if 'N' in kmer:
                    continue
                
                if kmer in target_kmers:
                    # Store (kmer, chrom, pos)
                    pos = start_coord + i
                    batch_data.append((kmer, current_chrom, pos))
            
            if len(batch_data) >= batch_size:
                c.executemany('INSERT INTO kmer_locations VALUES (?, ?, ?)', batch_data)
                conn.commit()
                total_indexed += len(batch_data)
                batch_data = []
                print(f"Indexed {total_indexed} locations...", end='\r')

            # Prepare next overlap
            if len(full_seq) >= K-1:
                overlap = full_seq[-(K-1):]
            else:
                overlap = full_seq
                
            chunk_idx += len(seq)
            
    if batch_data:
        c.executemany('INSERT INTO kmer_locations VALUES (?, ?, ?)', batch_data)
        conn.commit()
        total_indexed += len(batch_data)
        
    print(f"\nMapping complete. Total locations stored: {total_indexed}")
    conn.close()

def query_probe(db_path, probe_seq, check_mismatches=False):
    if not os.path.exists(db_path):
        print(f"Error: Database {db_path} does not exist.")
        return

    probe_seq = probe_seq.upper()
    print(f"\nAnalyzing Probe: {probe_seq}")
    
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    kmers_in_probe = []
    if len(probe_seq) < K:
        print(f"Probe is shorter than K={K}. Cannot extract k-mers.")
        return

    print(f"\n{'Position':<10} {'K-mer':<15} {'Exact Count':<15} {'Mismatch Info / Locations'}")
    print("-" * 100)
    
    # Check if locations table exists
    has_locations = False
    try:
        c.execute('SELECT 1 FROM kmer_locations LIMIT 1')
        has_locations = True
    except sqlite3.OperationalError:
        pass
    
    max_exact_freq = 0
    
    for i in range(len(probe_seq) - K + 1):
        kmer = probe_seq[i : i+K]
        c.execute('SELECT count FROM kmer_counts WHERE kmer = ?', (kmer,))
        row = c.fetchone()
        count = row[0] if row else 0
        
        info_parts = []
        
        # 1. Location info (if available)
        if has_locations and count > 0:
            c.execute('SELECT chrom, pos FROM kmer_locations WHERE kmer = ? LIMIT 5', (kmer,))
            locs = c.fetchall()
            if locs:
                loc_str = ", ".join([f"{l[0]}:{l[1]}" for l in locs])
                if count > 5:
                    loc_str += f" ..."
                info_parts.append(f"Locs: {loc_str}")
        
        # 2. Mismatch info
        if check_mismatches:
            neighbors = generate_neighbors(kmer)
            neighbor_hits = []
            for n in neighbors:
                c.execute('SELECT count FROM kmer_counts WHERE kmer = ?', (n,))
                nrow = c.fetchone()
                if nrow:
                    # neighbor_hits.append(f"{n}({nrow[0]})")
                     neighbor_hits.append(f"{n}({nrow[0]})")
            if neighbor_hits:
                info_parts.append(f"Neighbors: {', '.join(neighbor_hits[:3])}")

        info_str = " | ".join(info_parts)
        print(f"{i+1:<10} {kmer:<15} {count:<15} {info_str}")
        
        if count > max_exact_freq:
            max_exact_freq = count
        kmers_in_probe.append((kmer, count))

    conn.close()
    print("-" * 100)
    if max_exact_freq > 1000:
        print("\n[WARNING] Highly repetitive sequences detected.")

def get_probe_metrics(db_path, probe_seq):
    """
    Analyzes a probe sequence against the k-mer database and returns metrics.
    Returns:
        dict: {
            'max_kmer_count': int,
            'mean_kmer_count': float,
            'high_freq_kmers': list of (kmer, count)
        }
    """
    if not os.path.exists(db_path):
        return None

    probe_seq = probe_seq.upper()
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    
    max_count = 0
    total_count = 0
    kmer_counts = []
    
    num_kmers = len(probe_seq) - K + 1
    if num_kmers < 1:
        conn.close()
        return {'max_kmer_count': 0, 'mean_kmer_count': 0, 'high_freq_kmers': []}

    for i in range(num_kmers):
        kmer = probe_seq[i : i+K]
        c.execute('SELECT count FROM kmer_counts WHERE kmer = ?', (kmer,))
        row = c.fetchone()
        count = row[0] if row else 0
        
        if count > max_count:
            max_count = count
        total_count += count
        
        if count > 100: # Threshold for "high frequency" tracking
            kmer_counts.append((kmer, count))
            
    conn.close()
    
    return {
        'max_kmer_count': max_count,
        'mean_kmer_count': total_count / num_kmers if num_kmers > 0 else 0,
        'high_freq_kmers': kmer_counts
    }


def main():
    parser = argparse.ArgumentParser(description="K-mer Frequency Database Builder & Query Tool")
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # Build
    parser_build = subparsers.add_parser("build", help="Count k-mer frequencies")
    parser_build.add_argument("--fasta", required=True, help="Input FASTA")
    parser_build.add_argument("--db", required=True, help="Output DB")
    
    # Map (New)
    parser_map = subparsers.add_parser("map", help="Index k-mer locations (Recommended for transcriptomes/small genomes)")
    parser_map.add_argument("--fasta", required=True, help="Input FASTA")
    parser_map.add_argument("--db", required=True, help="Existing DB from 'build' step")
    parser_map.add_argument("--min-repeats", type=int, default=2, help="Only index k-mers appearing >= N times (default: 2)")

    # Query
    parser_query = subparsers.add_parser("query", help="Query a probe")
    parser_query.add_argument("--db", required=True, help="SQLite DB")
    parser_query.add_argument("--probe", required=True, help="Probe sequence")
    parser_query.add_argument("--check-mismatches", action="store_true", help="Check for 1-mismatch off-targets")
    
    args = parser.parse_args()
    
    if args.command == "build":
        build_database(args.fasta, args.db)
    elif args.command == "map":
        map_locations(args.fasta, args.db, args.min_repeats)
    elif args.command == "query":
        query_probe(args.db, args.probe, args.check_mismatches)
if __name__ == "__main__":
    main()
