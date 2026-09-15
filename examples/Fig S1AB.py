"""
Fig S1AB: CSV-only decile analysis (unified example), ATAC log2FC.

This example demonstrates the DecileScatterPipeline which works entirely
from CSV data without extracting signal arrays from BigWig files. It shows
two variants, selected by which `CHIP_COL`/`CHIP_YLABEL` pair is active below:
- Sorting by ChIP log2 fold-change
- Sorting by ChIP baseMean (mean expression)

Optionally filters data by category (e.g., FACT category).

Workflow:
1. Load CSV with ChIP and ATAC columns
2. Optionally filter by category
3. Sort by selected column (log2FC or baseMean) descending
4. Create decile groups
5. Generate dual-axis boxplot figures

"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import DecileScatterPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")

CHIP_YLABEL = "ChIP (log2 fold change)"
N_DECILES = 10


###############################################
################## FIGURE S1A ##################
###############################################

ATAC_COL_A = "ssrp1 ATAC log2 Fold Change"
CHIP_COL_A = "ssrp1 ChIP log2 fold change"
ATAC_YLABEL_A = "ATAC (log2 fold change)"
OUTPUT_PATH_A = "examples/outputs/figs/figS1A.svg" # Output file (set to None to skip saving)
TITLE_A = "SSRP1"


###############################################
################## FIGURE S1B ##################
###############################################

ATAC_COL_B = "supt16h ATAC log2 Fold Change"
CHIP_COL_B = "suptt16h ChIP log2 fold change"
ATAC_YLABEL_B = "ATAC (log2 fold change)"
OUTPUT_PATH_B = "examples/outputs/figs/figS1B.svg" # Output file (set to None to skip saving)
TITLE_B = "SUPT16H"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print(f"Starting CSV-based decile analysis...")

    pipeA = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_A,
        rna_col=ATAC_COL_A,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=ATAC_YLABEL_A,
        label_fontsize=25,
        title= TITLE_A,
    )

    fig = pipeA.run()

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        fig.savefig(OUTPUT_PATH_A, dpi=300)
        print(f"  - Figure A saved to: {OUTPUT_PATH_A}")

    pipeB = DecileScatterPipeline(
        csv_path=CSV_PATH,
        chip_col=CHIP_COL_B,
        rna_col=ATAC_COL_B,
        n_deciles=N_DECILES,
        chip_ylabel=CHIP_YLABEL,
        rna_ylabel=ATAC_YLABEL_B,
        label_fontsize=25,
        title= TITLE_B,
    )

    fig = pipeB.run()

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        fig.savefig(OUTPUT_PATH_B, dpi=300)
        print(f"  - Figure B saved to: {OUTPUT_PATH_B}")

    plt.show()

