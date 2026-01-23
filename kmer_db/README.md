# K-mer Frequency Database Builder

This tool allows you to count 10-mer frequencies in a genome (e.g., Human Genome) and store them in a local SQLite database. You can then use this database to analyze the off-target potential of designed probes by checking the frequency of their constituent k-mers.

## Features

- **Efficient Counting**: Uses a sliding window approach to count 10-mers.
- **SQLite Storage**: Stores frequencies in a portable `.db` file.
- **Probe Analysis**: Analyzing a probe sequence to see if it contains highly repetitive subsequences.

## Usage

### 1. Build the Database

```bash
python kmer_tool.py build --fasta input_genome.fa --db human_kmer_counts.db
```

### 2. Map Locations (Optional)
To query *where* repetitive sequences are located, you must map the genome locations. 
*Warning: This can result in a very large database for large genomes.*

```bash
# Only map locations for k-mers that appear at least twice
python kmer_tool.py map --fasta input_genome.fa --db output.db --min-repeats 2
```

### 3. Analyze a Probe
Run a basic query. If you have run `map`, it will also show you the genomic coordinates (Chrom:Position) of the hits.
```bash
python kmer_tool.py query --db output.db --probe ATGCATGCATGC
```

Run an off-target analysis that checks for 1-mismatch neighbors:
```bash
python kmer_tool.py query --db output.db --probe ATGCATGCATGC --check-mismatches
```

## Requirements

- Python 3.x
- No external dependencies required (uses standard `sqlite3` and `argparse`).
