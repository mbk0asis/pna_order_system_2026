from . import config
from . import filters
from . import thermo
import pandas as pd

class PNADesigner:
    def __init__(self, target_sequence: str, gene_name: str = "Unknown"):
        self.target_sequence = target_sequence
        self.gene_name = gene_name
        self.candidates = []
        
    def design_probes(self, lengths=None, tm_min=config.TARGET_TM_MIN, tm_max=config.TARGET_TM_MAX):
        if lengths is None:
            lengths = range(config.MIN_LENGTH, config.MAX_LENGTH + 1)
            
        seq = self.target_sequence.upper()
        # Ensure we are designing against the target.
        # usually ISH probes are COMPLEMENTARY to the mRNA.
        # So we take tiles from the mRNA (Target) and then find the COMPLEMENT.
        
        # NOTE: The "sequence" passed to analysis tools usually refers to the PROBE sequence.
        # 1. Tile the mRNA.
        # 2. Complement the tile to get the PROBE.
        # 3. Analyze the PROBE.
        
        from Bio.Seq import Seq
        seq_obj = Seq(seq)
        
        results = []
        
        for L in lengths:
            for i in range(len(seq) - L + 1):
                target_site = seq[i : i + L] # Sequence on mRNA
                
                # Probe is Reverse Complement of the target site
                probe_seq_obj = Seq(target_site).reverse_complement()
                probe_seq = str(probe_seq_obj)
                
                # 1. Run Filters
                validation = filters.validate_candidate(probe_seq)
                
                # 2. Calculate Tm
                dna_tm = thermo.calculate_dna_tm(probe_seq)
                pna_tm = thermo.calculate_pna_tm_giesen(probe_seq, dna_tm)
                
                # Check Tm Range
                tm_valid = (tm_min <= pna_tm <= tm_max)
                if not tm_valid:
                    validation['errors'].append(f"Tm {pna_tm}C out of range ({tm_min}-{tm_max})")
                    validation['valid'] = False
                
                results.append({
                    "start": i + 1,
                    "end": i + L,
                    "length": L,
                    "target_site": target_site,
                    "probe_sequence": probe_seq,
                    "gc_content": validation['metrics']['gc'],
                    "purine_content": validation['metrics']['purine'],
                    "tm_dna_nn": round(dna_tm, 2),
                    "tm_pna_giesen": pna_tm,
                    "valid": validation['valid'],
                    "errors": "; ".join(validation['errors'])
                })
                
        self.candidates = pd.DataFrame(results)
        return self.candidates

    def get_best_candidates(self, top_n=5):
        if self.candidates.empty:
            return pd.DataFrame()
        
        # Filter for valid only
        valid_df = self.candidates[self.candidates['valid'] == True].copy()
        
        if valid_df.empty:
            return pd.DataFrame() # No valid candidates
            
        # Sorting Criteria: 
        # 1. GC closest to 50%?
        # 2. Higher Tm (specificity)?
        # 3. Lower Purine?
        
        # Let's sort by how close Tm is to center of range (75C)
        target_tm = (config.TARGET_TM_MIN + config.TARGET_TM_MAX) / 2
        valid_df['tm_diff'] = abs(valid_df['tm_pna_giesen'] - target_tm)
        
        return valid_df.sort_values('tm_diff').head(top_n)
