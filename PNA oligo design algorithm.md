# PNA Oligo Design Algorithm for In Situ Hybridization (ISH)

Peptide Nucleic Acid (PNA) oligos are synthetic DNA mimics with a neutral **polyamide backbone**. The lack of electrostatic repulsion leads to exceptionally high binding affinity and specificity, requiring a distinct algorithmic approach compared to standard DNA/RNA probes.

---

## 1. Core Design Parameters

Unlike DNA-FISH, PNA-ISH favors shorter, high-affinity probes to ensure optimal cell penetration and specificity.

| Parameter | Recommended Range | Rationale |
| :--- | :--- | :--- |
| **Length** | 12 – 18 bases | Optimal balance of specificity, kinetics, and cell penetration. |
| **GC Content** | 35% – 60% | Prevents aggregation while maintaining stable melting temperatures. |
| **Purine Content** | < 60% (Ideal < 50%) | High purine content (A, G) causes PNA to become insoluble and aggregate. |
| **G-Clusters** | Max 3 consecutive Gs | Prevents G-quadruplex-like aggregation and "sticky" non-specific binding. |
| **Mismatches** | 0 (Perfect Match) | PNA is highly sensitive; a single mismatch can drop $T_m$ by 10–15°C. |
| **Self-Complement**| < 5 bp | PNA-PNA interactions are very strong; avoid stable hairpins or dimers. |

---

## 2. Thermodynamic Prediction

PNA/DNA duplexes are significantly more stable than DNA/DNA duplexes. A 15-mer PNA often has a $T_m$ equivalent to a 25-mer DNA.

### $T_m$ Empirical Formula
A commonly used model for predicting $T_m$ for PNA/DNA duplexes (at 5 μM PNA) is:

$$T_{m,PNA} = c_0 + c_1 \cdot T_{m,DNA} + c_2 \cdot f_{GC} + c_3 \cdot L$$

*Where:*
*   $f_{GC}$: Fraction of GC content.
*   $L$: Length of the PNA.
*   $T_{m,DNA}$: Melting temperature of the equivalent DNA/DNA duplex.

> [!TIP]
> **Low-Salt Buffers:** Because PNA backbones are neutral, $T_m$ is largely independent of salt concentration. This allows ISH to be performed in low-salt buffers, destabilizing target RNA/DNA secondary structures and making them more accessible.

---

## 3. Automated Design Workflow

An optimized PNA design algorithm typically follows these six logical steps:

1.  **Sequence Tiling:** The target region (mRNA or genomic DNA) is "tiled" into candidate windows (12–18 nt).
2.  **Thermodynamic Filtering:** Candidates are filtered based on predicted $T_m$ (Target: 65°C – 75°C for ISH).
3.  **Composition Filtering:**
    *   Flag any sequences with $>3$ consecutive Gs.
    *   Discard sequences with $>60\%$ purine content.
    *   Maintain GC content between 35% and 60%.
4.  **Secondary Structure Analysis:** Identify and remove probes forming stable hairpins or self-dimers ($>4$ bp).
5.  **Specificity (Homology) Search:**
    *   Align candidates against the entire genome/transcriptome (BLAST/Bowtie2).
    *   **Specificity Target:** $>99\%$. Discard candidates with $>70\%$ homology to off-target sites.
6.  **Solubility Assessment:** Suggest modifications (e.g., O-linkers or Lysine tags) for hydrophobic sequences.

---

## 4. Synthesis & Optimization for ISH

When the algorithm identifies a high-quality sequence, it is prepared for synthesis with the following considerations:

*   **Solubility Tags:** If the sequence is hydrophobic, add **O-linkers** or **Lysine residues** to the N-terminus.
*   **Labeling:** Fluorophores (e.g., Cy3, Alexa Fluor) are added to the N-terminus via a linker to prevent steric hindrance during hybridization.
*   **Formamide Usage:** In experimental setups, each 1% formamide typically reduces $T_m$ by ~1°C, allowing for tunable stringency.

---

## 5. Recommended Software & Databases

### Databases
*   **GenBank:** For general mRNA and genomic sequences.
*   **SILVA / RDP-II:** For ribosomal RNA (16S/23S/28S) target identification.

### Design Tools
*   **PNA Tool (PNA Bio):** Specialist calculator for PNA $T_m$ and solubility.
*   **OligoMiner / PaintSHOP:** Genome-scale design tools (requires adjusting length to 12–18 nt).
*   **Probe Match (RDP-II):** For calculating specificity and sensitivity against rRNA databases.

---

### Comparison Example
| Metric | Target Values |
| :--- | :--- |
| **Specificity** | $(nL_s / T_{nL}) \times 100 > 99\%$ |
| **Sensitivity** | $(L_s / T_{Ls}) \times 100 > 95\%$ |