# PNA Oligo Production System

This is an automated system for designing Peptide Nucleic Acid (PNA) probes for *In Situ* Hybridization (ISH), based on the algorithmic principles defined in this workspace.

## Features
- **Giesen Tm Calculation**: Implements the standard empirical formula for PNA/DNA duplexes.
- **Biophysical Filtering**: Checks for GC content, Purine content (overall + windowed), G-runs, and Self-complementarity.
- **Automation Triggers**: Simulates notification of Production and Finance teams.

## Installation

1. Install Python 3.9+.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

Run the system from the root directory:

```bash
# Analyze a raw sequence string
python -m src.main --sequence "AGCTAGCTAGCTGATCGTAGCTGATCGTAGCTAG" --gene "TestGene"

# Analyze a FASTA file
python -m src.main --file path/to/sequence.fasta
```

## System Structure

- `src/config.py`: Defines all design thresholds (Tm: 70-80C, GC: 40-60%, etc.).
- `src/filters.py`: Logic for rejecting invalid sequences.
- `src/thermo.py`: Giesen formula implementation.
- `src/designer.py`: Core tiling and probe generation engine.
