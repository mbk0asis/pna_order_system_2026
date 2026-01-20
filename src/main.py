import argparse
import sys
import os
from Bio import SeqIO
from . import designer

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
    
    target_seq = ""
    gene_name = args.gene
    
    if args.file:
        if os.path.exists(args.file):
            record = next(SeqIO.parse(args.file, "fasta"))
            target_seq = str(record.seq)
            gene_name = record.id
        else:
            print(f"Error: File {args.file} not found.")
            sys.exit(1)
    elif args.sequence:
        target_seq = args.sequence
    elif args.gene:
        # Fetch from NCBI
        try:
            from . import ncbi
            target_seq, acc_id, desc = ncbi.fetch_mrna_by_symbol(args.gene)
            print(f"Successfully downloaded {acc_id}: {desc}")
            gene_name = args.gene # Keep the symbol as the main identifier
        except Exception as e:
            print(f"Error fetching gene '{args.gene}': {e}")
            sys.exit(1)
    else:
        print("Please provide a sequence using --sequence, --file, or a gene symbol using --gene")
        sys.exit(1)
        
    print(f"Designing PNA probes for: {gene_name}")
    print(f"Sequence Length: {len(target_seq)} bp")
    
    # Run Designer
    pna_designer = designer.PNADesigner(target_seq, gene_name)
    all_candidates = pna_designer.design_probes()
    
    # Save all results (including failed ones for debugging)
    all_candidates.to_csv(args.output, index=False)
    print(f"\nFull analysis saved to {args.output}")
    
    # Get Best Candidates
    best = pna_designer.get_best_candidates(top_n=3)
    
    if not best.empty:
        print("\nTop Recommended Probes:")
        print(best[['probe_sequence', 'length', 'tm_pna_giesen', 'gc_content', 'purine_content']].to_string())
        
        # Trigger Automation
        notify_production(best, gene_name)
        notify_finance(best, gene_name)
    else:
        print("\nNo probes met the strict design criteria. Check the output CSV for failure reasons.")

if __name__ == "__main__":
    main()
