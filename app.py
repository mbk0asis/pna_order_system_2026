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

targets = [] # List of tuples (gene_name, sequence)

if input_method == "Gene Symbol (NCBI)":
    gene_input = st.text_input("Enter Gene Symbol(s) (comma separated)", value="")
    if gene_input:
        if st.button("Fetch Sequences"):
            gene_list = [g.strip() for g in gene_input.split(',') if g.strip()]
            
            with st.spinner(f"Fetching {len(gene_list)} genes from NCBI..."):
                fetched_targets = []
                for sym in gene_list:
                    try:
                        seq, acc_id, desc = ncbi.fetch_mrna_by_symbol(sym)
                        fetched_targets.append((sym, seq))
                        st.success(f"Downloaded {sym} ({acc_id})")
                    except Exception as e:
                        st.error(f"Error fetching {sym}: {e}")
                
                if fetched_targets:
                    st.session_state['targets'] = fetched_targets

    # Persist
    if 'targets' in st.session_state and input_method == "Gene Symbol (NCBI)":
        targets = st.session_state['targets']
        st.info(f"Loaded {len(targets)} targets: {', '.join([t[0] for t in targets])}")

elif input_method == "Raw Sequence":
    target_sequence = st.text_area("Paste DNA/mRNA Sequence (5'->3')", height=150)
    gene_name = st.text_input("Gene Name/Identifier", value="CustomSeq")
    if target_sequence:
        targets = [(gene_name, target_sequence)]

elif input_method == "FASTA File":
    uploaded_file = st.file_uploader("Upload FASTA (Single or Multi-record)", type=["fasta", "fa"])
    if uploaded_file is not None:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        # Parse all records
        for record in SeqIO.parse(stringio, "fasta"):
            targets.append((record.id, str(record.seq)))
        st.info(f"Loaded {len(targets)} sequences from file.")

# --- Analysis ---
if targets:
    if st.button("Run PNA Design", type="primary"):
        all_probes_list = []
        
        progress_bar = st.progress(0)
        
        for idx, (g_name, t_seq) in enumerate(targets):
            # Update Progress
            progress_bar.progress((idx + 1) / len(targets))
            
            # Run Design
            pna_ds = designer.PNADesigner(t_seq, g_name)
            probes = pna_ds.design_probes(
                lengths=range(min_len, max_len + 1), 
                tm_min=target_tm_min, 
                tm_max=target_tm_max
            )
            
            if not probes.empty:
                probes['Gene'] = g_name
                all_probes_list.append(probes)
        
        progress_bar.empty()
        
        if all_probes_list:
            full_df = pd.concat(all_probes_list, ignore_index=True)
            
            # --- Apply Dynamic User Filters ---
            mask_gc = (full_df['gc_content'] >= gc_min) & (full_df['gc_content'] <= gc_max)
            mask_purine = (full_df['purine_content'] <= purine_max)
            mask_tm = (full_df['tm_pna_giesen'] >= target_tm_min) & (full_df['tm_pna_giesen'] <= target_tm_max)
            mask_valid = full_df['valid'] # Only basic structure validity
            
            filtered_df = full_df[mask_gc & mask_purine & mask_tm & mask_valid].copy()
            
            # Metrics
            st.divider()
            col1, col2 = st.columns(2)
            col1.metric("Total Candidates Generated", len(full_df))
            col2.metric("Selected Candidates (Filtered)", len(filtered_df))
            
            if not filtered_df.empty:
                st.subheader("Top Candidates (Grouped by Gene)")
                
                # Show top 3 per gene
                display_df = filtered_df.groupby('Gene').apply(lambda x: x.head(3)).reset_index(drop=True)
                
                st.dataframe(
                    display_df[['Gene', 'start', 'probe_sequence', 'tm_pna_giesen', 'gc_content']]
                    .style.format({"tm_pna_giesen": "{:.2f}"})
                )
                
                # Download Full Batch
                csv = filtered_df.to_csv(index=False)
                st.download_button(
                    label="Download All Selected Candidates (CSV)",
                    data=csv,
                    file_name="batch_pna_results.csv",
                    mime="text/csv"
                )
            else:
                st.warning("No candidates met the filtered criteria.")
                
        else:
            st.error("No valid probes generated for any target.")
