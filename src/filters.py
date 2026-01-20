from Bio.Seq import Seq
from . import config

def calculate_gc_content(sequence: str) -> float:
    seq_upper = sequence.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    if len(sequence) == 0:
        return 0.0
    return (gc_count / len(sequence)) * 100.0

def calculate_purine_content(sequence: str) -> float:
    seq_upper = sequence.upper()
    purine_count = seq_upper.count('A') + seq_upper.count('G')
    if len(sequence) == 0:
        return 0.0
    return (purine_count / len(sequence)) * 100.0

def check_g_runs(sequence: str, max_consecutive: int = config.MAX_CONSECUTIVE_G) -> bool:
    """
    Returns False if G-run exceeds limit (fail), True if OK (pass).
    """
    g_run = 'G' * (max_consecutive + 1)
    return g_run not in sequence.upper()

def check_purine_window(sequence: str, window_size: int = config.PURINE_WINDOW_SIZE, max_purines: int = config.PURINE_WINDOW_MAX) -> bool:
    """
    Returns False if any window has too many purines.
    """
    seq_upper = sequence.upper()
    for i in range(len(seq_upper) - window_size + 1):
        window = seq_upper[i : i + window_size]
        purines = window.count('A') + window.count('G')
        if purines > max_purines:
            return False
    return True

def check_self_complementarity(sequence: str, max_bp: int = config.MAX_SELF_COMPLEMENTARITY_BP) -> bool:
    """
    Checks for self-complementarity stretches greater than max_bp.
    Simple heuristic: check if any substring of length (max_bp + 1) matches the reverse complement of any other part.
    
    Returns True if safe (no significant self-comp), False if fails.
    """
    seq_len = len(sequence)
    threshold = max_bp + 1
    seq_obj = Seq(sequence)
    rc_seq = str(seq_obj.reverse_complement())
    
    # We look for common substrings between sequence and its reverse complement
    # Usually we care about hairpin stems.
    # A simple approach: Slide a window of 'threshold' length across seq
    # and see if it exists in rc_seq.
    
    # Note: This is an approximation. A rigorous check involves folding energy.
    # But strictly speaking, "self-complementary stretches > 4 bp" implies simple base pairing.
    
    for i in range(seq_len - threshold + 1):
        kmer = sequence[i : i + threshold]
        # For a hairpin, the complementary part must be elsewhere. 
        # But simply checking existence in RC checks for potential dimerization or hairpinning.
        if kmer in rc_seq:
            # We must ensure it's not just matching itself in the reverse direction if palindromic
            # but for probe design, any capable structure is usually bad.
            return False
            
    return True

def validate_candidate(sequence: str) -> dict:
    """
    Run all filters on a candidate sequence.
    Returns a dict with 'valid' (bool) and 'errors' (list of strings).
    """
    errors = []
    
    gc = calculate_gc_content(sequence)
    if not (config.GC_CONTENT_MIN_STRICT <= gc <= config.GC_CONTENT_MAX_STRICT):
        errors.append(f"GC Content {gc:.1f}% out of range ({config.GC_CONTENT_MIN_STRICT}-{config.GC_CONTENT_MAX_STRICT}%)")
        
    purine = calculate_purine_content(sequence)
    if purine > config.PURINE_CONTENT_MAX:
        errors.append(f"Purine Content {purine:.1f}% > {config.PURINE_CONTENT_MAX}%")
        
    g_perc = (sequence.upper().count('G') / len(sequence)) * 100
    if g_perc > config.G_CONTENT_MAX:
        errors.append(f"G Content {g_perc:.1f}% > {config.G_CONTENT_MAX}%")
        
    if not check_g_runs(sequence):
        errors.append(f"Contains >{config.MAX_CONSECUTIVE_G} consecutive Gs")
        
    if not check_purine_window(sequence):
        errors.append(f"Fails Purine Window Check (>{config.PURINE_WINDOW_MAX} purines in {config.PURINE_WINDOW_SIZE} bp)")
        
    if not check_self_complementarity(sequence):
        errors.append(f"Self-complementarity detected > {config.MAX_SELF_COMPLEMENTARITY_BP} bp")
        
    return {
        "valid": len(errors) == 0,
        "errors": errors,
        "metrics": {
            "gc": gc,
            "purine": purine,
            "g_percent": g_perc
        }
    }
