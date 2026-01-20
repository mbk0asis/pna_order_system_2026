from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import math
from . import config

def calculate_dna_tm(sequence: str) -> float:
    """
    Calculates the Nearest-Neighbor Tm for the DNA/DNA duplex version of the sequence.
    Uses SantaLucia parameters as per Giesen et al. requirements.
    
    Args:
        sequence (str): The nucleotide sequence (5'->3').
        
    Returns:
        float: The calculated DNA Tm.
    """
    # Create Seq object
    seq_obj = Seq(sequence)
    
    # Calculate DNA Tm using Biopython
    # Giesen standard conditions: 100nM DNA (using standard), 100mM NaCl
    # We use standard DNA values appropriate for the Giesen formula context.
    tm_dna = mt.Tm_NN(
        seq_obj, 
        dnac1=50, # 50nM primer
        dnac2=50, # 50nM target
        saltcorr=7, # Sodium correction
        Na=config.DNA_TM_SALT_CONC
    )
    return tm_dna

def calculate_pna_tm_giesen(sequence: str, dna_tm: float = None) -> float:
    """
    Calculates PNA/DNA Tm using the Giesen et al. (1998) empirical formula.
    
    Formula:
    Tm(PNA/DNA) = 20.79 + 0.83*Tm(DNA) - 26.13*F_pyr + 0.44*Length
    
    Args:
        sequence (str): The PNA sequence.
        dna_tm (float, optional): Pre-calculated DNA Tm. If None, it will be calculated.
        
    Returns:
        float: The predicted PNA/DNA Tm.
    """
    seq_upper = sequence.upper()
    length = len(seq_upper)
    
    if length == 0:
        return 0.0
        
    if dna_tm is None:
        dna_tm = calculate_dna_tm(sequence)
        
    # Calculate Fraction of Pyrimidines (C + T)
    # Note: PNA typically uses T, but sometimes U. We assume standard notation.
    pyr_count = seq_upper.count('C') + seq_upper.count('T') + seq_upper.count('U')
    f_pyr = pyr_count / length
    
    # Giesen Formula
    tm_pna = 20.79 + (0.83 * dna_tm) - (26.13 * f_pyr) + (0.44 * length)
    
    return round(tm_pna, 2)
