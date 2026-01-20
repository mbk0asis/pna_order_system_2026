from Bio import Entrez, SeqIO
from . import config
import sys

# Set email for NCBI Entrez
Entrez.email = config.ENTREZ_EMAIL

def fetch_mrna_by_symbol(gene_symbol: str, organism: str = "Homo sapiens") -> tuple[str, str, str]:
    """
    Searches NCBI for the canonical mRNA (RefSeq) of a given gene symbol.
    
    Args:
        gene_symbol (str): The gene symbol (e.g., "GAPDH").
        organism (str): The organism to filter by (default "Homo sapiens").
        
    Returns:
        tuple: (sequence_string, accession_id, description)
    """
    # 1. Search for the ID
    # Use strict filtering for RefSeq mRNA
    term = f"{gene_symbol}[Gene Name] AND {organism}[Organism] AND biomol_mrna[PROP] AND srcdb_refseq[PROP]"
    
    print(f"Searching NCBI for: {term}...")
    
    try:
        handle = Entrez.esearch(db=config.ENTREZ_DB, term=term, retmax=5, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        raise Exception(f"Failed to connect to NCBI Entrez: {e}")
    
    id_list = record["IdList"]
    
    if not id_list:
        raise ValueError(f"No mRNA Reference Sequences found for gene '{gene_symbol}' in {organism}.")
        
    # We take the first hit, assuming relevance sort puts the best NM_ match first.
    # Often there are multiple variants (NM_001, NM_002). We just pick the first one which is usually isoform 1 or most common.
    # A more advanced version would let the user pick.
    best_id = id_list[0]
    
    # 2. Fetch the Sequence
    print(f"Fetching sequence for Accession ID: {best_id}...")
    try:
        handle = Entrez.efetch(db=config.ENTREZ_DB, id=best_id, rettype="fasta", retmode="text")
        seq_record = next(SeqIO.parse(handle, "fasta"))
        handle.close()
    except Exception as e:
        raise Exception(f"Failed to fetch sequence for ID {best_id}: {e}")
        
    return str(seq_record.seq), seq_record.id, seq_record.description
