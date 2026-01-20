# Biopython is required for this config to be fully utilized in the main app, 
# but this file defines the constants.

# --- NCBI Configuration ---
ENTREZ_EMAIL = "pna_designer_automation@example.com" # PLEASE UPDATE THIS
ENTREZ_DB = "nucleotide"


# --- Global Design Constraints ---

# Length
MIN_LENGTH = 12
MAX_LENGTH = 15
# For aqueous applications, max is ~23, but sticking to optimal 15-18 per docs.

# Melting Temperature (Tm)
TARGET_TM_MIN = 30.0
TARGET_TM_MAX = 60.0
# Note: Experimental Tm is often 5-10C lower than calculated.

# Composition
GC_CONTENT_MIN = 30.0 # From order system.md (range 30-70), Algorithmic principles says 40-60.
# We will use the stricter "Algorithmic principles" range as default but allow relaxation.
GC_CONTENT_MIN_STRICT = 40.0
GC_CONTENT_MAX_STRICT = 60.0

G_CONTENT_MAX = 35.0
MAX_CONSECUTIVE_G = 3

# Purine Content (A + G)
PURINE_CONTENT_MAX = 50.0 # Aqueous limit
PURINE_WINDOW_SIZE = 10
PURINE_WINDOW_MAX = 6 # Max purines in any 10-base window

# Structural
MAX_SELF_COMPLEMENTARITY_BP = 4 # >4bp is rejected

# Salt Correction for DNA Tm Calculation (Giesen formula standard conditions)
DNA_TM_SALT_CONC = 100.0 # mM NaCl
DNA_TM_PH = 7.0
