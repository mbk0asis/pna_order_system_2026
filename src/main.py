import argparse
import sys
import os
import pandas as pd
from Bio import SeqIO
from . import designer
from . import config

# Mock function for email notification
def notify_production(best_probes, gene):
    print("\n--- [AUTOMATION] NOTIFYING PRODUCTION TEAM ---")
    print(f"Subject: New PNA Order for {gene}")
    print("Body: Please synthesize the following sequences:")
    for idx, row in best_probes.iterrows():
        print(f" - Probe_{idx}: {row['probe_sequence']} (Tm: {row['tm_pna_giesen']}C)")
    print("----------------------------------------------")

def notify_finance(best_probes, gene):
    print("\n--- [AUTOMATION] NOTIFYING FINANCE TEAM ---")
    print(f"Subject: Invoice Generation for {gene}")
    print(f"Body: {len(best_probes)} oligos ordered. Please prepare invoice.")
    print("----------------------------------------------")

def main():
    parser = argparse.ArgumentParser(description="PNA Oligo Design Automation System")
    parser.add_argument("--sequence", help="Raw target sequence string (mRNA)")
    parser.add_argument("--file", help="Path to FASTA file containing target sequence")
    parser.add_argument("--gene", help="Gene Symbol (Mock lookup if sequence not provided)", default="TargetGene")
    parser.add_argument("--output", help="Output CSV file for results", default="pna_results.csv")
    
    args = parser.parse_args()
    
    targets = [] # List of tuples: (gene_name, sequence)
    
    # 1. Input Source Processing
    if args.file:
        if os.path.exists(args.file):
            print(f"Reading batch from {args.file}...")
            for record in SeqIO.parse(args.file, "fasta"):
                targets.append((record.id, str(record.seq)))
        else:
            print(f"Error: File {args.file} not found.")
            sys.exit(1)
            
    elif args.sequence:
        targets.append((args.gene, args.sequence))
        
    elif args.gene:
        # Handle comma-separated list
        gene_list = [g.strip() for g in args.gene.split(',')]
        from . import ncbi
        
        for g_sym in gene_list:
            try:
                print(f"Fetching {g_sym}...")
                seq, acc_id, desc = ncbi.fetch_mrna_by_symbol(g_sym)
                # Use Symbol as ID, but maybe store accession in description? 
                # For simplicity, we use Symbol as gene_name
                targets.append((g_sym, seq))
                print(f" - Found {acc_id}")
            except Exception as e:
                print(f" - Error fetching {g_sym}: {e}")
                
    else:
        print("Please provide input via --file, --sequence, or --gene")
        sys.exit(1)
        
    if not targets:
        print("No valid targets found.")
        sys.exit(1)
        
    print(f"\nProcessing {len(targets)} target(s)...")
    
    all_results = []
    
    # 2. Batch Processing
    for gene_name, seq in targets:
        print(f"Designing for {gene_name} ({len(seq)} bp)...")
        # Initialize Designer
        pna_ds = designer.PNADesigner(seq, gene_name)
        # Run (using defaults in CLI)
        candidates = pna_ds.design_probes()
        
        # Add Gene Name column to distinguish in batch file
        candidates['Gene'] = gene_name
        
        all_results.append(candidates)
        
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        
        # Save
        final_df.to_csv(args.output, index=False)
        print(f"\nBatch analysis saved to {args.output}")
        
        # Filter Best for Notifications
        # We need to filter 'valid' ones from the big dataframe
        valid_only = final_df[final_df['valid'] == True]
        
        if not valid_only.empty:
            # We can pick top 3 per gene
            best_list = []
            for g, group in valid_only.groupby('Gene'):
                 # Sort by simple Tm distance from center (approx)
                 target_median = (config.TARGET_TM_MIN + config.TARGET_TM_MAX) / 2
                 group['tm_diff'] = abs(group['tm_pna_giesen'] - target_median)
                 best_list.append(group.sort_values('tm_diff').head(3))
            
            best_df = pd.concat(best_list)
            
            print("\nTop Candidates Per Gene:")
            print(best_df[['Gene', 'probe_sequence', 'tm_pna_giesen', 'gc_content']].to_string(index=False))
            
            notify_production(best_df, "Batch_Order")
            notify_finance(best_df, "Batch_Order")
        else:
             print("\nNo candidates passed strict filtering for any gene.")
    else:
        print("No results generated.")

if __name__ == "__main__":
    main()
