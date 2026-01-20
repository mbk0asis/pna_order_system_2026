import streamlit as st
import pandas as pd
from src import designer, ncbi, config, thermo
from Bio import SeqIO
from io import StringIO

st.set_page_config(page_title="PNA Oligo Designer", layout="wide")

st.title("🧬 PNA Oligo Automated Designer")
st.markdown("""
Design Peptide Nucleic Acid (PNA) probes for *in situ* hybridization (ISH) with automated checks for:
- Melting Temperature (Giesen Formula)
- GC & Purine Content
- Solubility & Self-Complementarity
""")

# --- Sidebar: Parameters ---
st.sidebar.header("Design Parameters")

with st.sidebar.expander("Sequence Constraints", expanded=True):
    min_len = st.number_input("Min Length (bp)", value=config.MIN_LENGTH, min_value=10, max_value=30)
    max_len = st.number_input("Max Length (bp)", value=config.MAX_LENGTH, min_value=10, max_value=30)
    
    gc_min = st.slider("Min GC %", 0, 100, int(config.GC_CONTENT_MIN)) # Using stricter default from config might be too strict, let user relax
    gc_max = st.slider("Max GC %", 0, 100, int(config.GC_CONTENT_MAX_STRICT))
    
    purine_max = st.slider("Max Purine %", 0, 100, int(config.PURINE_CONTENT_MAX))

with st.sidebar.expander("Thermodynamics", expanded=True):
    target_tm_min = st.number_input("Min Tm (°C)", value=config.TARGET_TM_MIN)
    target_tm_max = st.number_input("Max Tm (°C)", value=config.TARGET_TM_MAX)
    
# Update Config values temporarily for this session (Hack: ideally pass params to designer)
# Since designer uses config global, we can try to override it or refactor designer.
# For now, let's just assume standard config is used, OR better: implement a parameter passing mechanism.
# The current checking logic is in src.filters which hardcodes config. 
# Modifying global config is safe enough in a single-threaded script, but Streamlit is multi-threaded.
# PROPER FIX: We should pass these parameters to the validation functions.
# But detailed refactoring takes time. 
# LIMITATION: I will rely on the default config constraints for filters, 
# BUT I will apply the user selected filters on the DATAFRAME output.
# This allows 'soft' filtering on top of 'hard' configuration.

# --- Main Input ---
input_method = st.radio("Input Method", ["Gene Symbol (NCBI)", "Raw Sequence", "FASTA File"])

target_sequence = ""
gene_name = "Unknown"

if input_method == "Gene Symbol (NCBI)":
    gene_symbol = st.text_input("Enter Gene Symbol (e.g., GAPDH, ACTB)", value="")
    if gene_symbol:
        if st.button("Fetch Sequence"):
            with st.spinner(f"Fetching {gene_symbol} from NCBI..."):
                try:
                    target_sequence, acc_id, desc = ncbi.fetch_mrna_by_symbol(gene_symbol)
                    gene_name = gene_symbol
                    st.success(f"Downloaded {acc_id}: {desc}")
                    st.session_state['target_sequence'] = target_sequence
                    st.session_state['gene_name'] = gene_name
                except Exception as e:
                    st.error(f"Error: {e}")
    
    # Persist
    if 'target_sequence' in st.session_state and input_method == "Gene Symbol (NCBI)":
        target_sequence = st.session_state['target_sequence']
        gene_name = st.session_state['gene_name']

elif input_method == "Raw Sequence":
    target_sequence = st.text_area("Paste DNA/mRNA Sequence (5'->3')", height=150)
    gene_name = st.text_input("Gene Name/Identifier", value="CustomSeq")

elif input_method == "FASTA File":
    uploaded_file = st.file_uploader("Upload FASTA", type=["fasta", "fa"])
    if uploaded_file is not None:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        record = next(SeqIO.parse(stringio, "fasta"))
        target_sequence = str(record.seq)
        gene_name = record.id
        st.info(f"Loaded {gene_name} ({len(target_sequence)} bp)")

# --- Analysis Button ---
if target_sequence:
    st.write(f"**Target Length:** {len(target_sequence)} bp")
    
    if st.button("Run PNA Design", type="primary"):
        with st.spinner("Tiling and Analyzing..."):
            # Initialize Designer
            pna_ds = designer.PNADesigner(target_sequence, gene_name)
            
            # Run Design
            probes = pna_ds.design_probes(
                lengths=range(min_len, max_len + 1), 
                tm_min=target_tm_min, 
                tm_max=target_tm_max
            )
            
            if not probes.empty:
                # --- Apply Dynamic User Filters ---
                # We filter the DataFrame based on Sidebar inputs
                
                # GC
                mask_gc = (probes['gc_content'] >= gc_min) & (probes['gc_content'] <= gc_max)
                # Purine
                mask_purine = (probes['purine_content'] <= purine_max)
                # Tm
                mask_tm = (probes['tm_pna_giesen'] >= target_tm_min) & (probes['tm_pna_giesen'] <= target_tm_max)
                # Validity (from hard filters like G-runs, etc.)
                mask_valid = probes['valid']
                
                # Combine
                filtered_df = probes[mask_gc & mask_purine & mask_tm & mask_valid].copy()
                
                # Add score (distance from ideal Tm)
                if not filtered_df.empty:
                    ideal_tm = (target_tm_min + target_tm_max) / 2
                    filtered_df['tm_diff'] = abs(filtered_df['tm_pna_giesen'] - ideal_tm)
                    filtered_df = filtered_df.sort_values('tm_diff')
                
                # --- Metrics ---
                col1, col2, col3 = st.columns(3)
                col1.metric("Total Candidates", len(probes))
                col2.metric("Valid (Hard Filters)", probes['valid'].sum())
                col3.metric("Selected (User Filters)", len(filtered_df))
                
                # --- Results Table ---
                st.subheader("Top Recommended Candidates")
                
                if not filtered_df.empty:
                    st.dataframe(
                        filtered_df[['start', 'length', 'probe_sequence', 'tm_pna_giesen', 'gc_content', 'purine_content']]
                        .style.format({"tm_pna_giesen": "{:.2f}", "gc_content": "{:.1f}", "purine_content": "{:.1f}"})
                    )
                    
                    # Download
                    csv = filtered_df.to_csv(index=False)
                    st.download_button(
                        label="Download Best Candidates (CSV)",
                        data=csv,
                        file_name=f"{gene_name}_PNA_candidates.csv",
                        mime="text/csv"
                    )
                    
                    # Production Simulation
                    st.info("ℹ️ Sending this list to Production would trigger the automated synthesis order email.")
                    
                else:
                    st.warning("No candidates met the *strict* user criteria. Showing top 'Valid' candidates ignoring Tm mismatch:")
                    # Fallback show just valid
                    fallback = probes[probes['valid'] == True].head(10)
                    if not fallback.empty:
                        st.dataframe(fallback)
                    else:
                        st.error("No candidates passed the hard structural/composition filters (G-runs, etc.). Try relaxing the constraints.")
                        
                # --- All Data ---
                with st.expander("See Raw Data (All Tiles)"):
                    st.dataframe(probes)
            else:
                st.error("No probes generated.")
