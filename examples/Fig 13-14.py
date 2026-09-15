"""
Fig 13-14: Circos plot comparing multi-signal log2 fold-change across conditions.

This example demonstrates the CircPlotPipeline which generates circos plots
comparing ChIP-seq, ATAC-seq, and RNA-seq log2 fold-change values for specified
genes across two conditions (e.g., SSRP1 and SUPT16H vs control).

The circular visualization shows three concentric rings for each signal type,
with connections highlighting genes that match user-specified combinations
of high/low signal changes.

Workflow:
1. Specify a list of genes of interest
2. Load expression data from CSV
3. Extract log2 fold-change values for ChIP, ATAC, and RNA for each condition
4. Generate two circos plots (one per condition)
5. Color genes by magnitude of change and show connections between selected combinations
6. (Optional) Export gene data to TSV files for further analysis
"""

import os

import matplotlib.pyplot as plt
import numpy as np
from Pyntegrate.pipelines import CircPlotPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")

# List of genes to include (Ensembl IDs or Gene Symbols)
GENES = [
    "SERPINE1",
    "F3",
    "CYR61",
    "TERT",
    "RPSAP52",
    "NEAT1"
]

# Column names for SSRP1 condition
SSRP1_CHIP_COL = "ssrp1 ChIP log2 fold change"
SSRP1_RNA_COL = "ssrp1 RNA log2 fold change"
SSRP1_ATAC_COL = "ssrp1 ATAC log2 Fold Change"
SSRP1_CHIP_BASEMEAN_COL = "ssrp1 ChIP baseMean"
SSRP1_RNA_BASEMEAN_COL = "ssrp1 RNA baseMean"
SSRP1_ATAC_BASEMEAN_COL = "ssrp1 ATAC baseMean"


# Column names for SUPT16H condition
SUPT16H_CHIP_COL = "suptt16h ChIP log2 fold change"
SUPT16H_RNA_COL = "supt16h RNA log2 fold change"
SUPT16H_ATAC_COL = "supt16h ATAC log2 Fold Change"
SUPT16H_CHIP_BASEMEAN_COL = "supt16h ChIP baseMean"
SUPT16H_RNA_BASEMEAN_COL = "supt16h RNA baseMean"
SUPT16H_ATAC_BASEMEAN_COL = "supt16h ATAC baseMean"


# ============================================================================
# Configuration: Threshold
# ============================================================================
# Threshold: values >= threshold are "high", values <= -threshold are "low"
THRESHOLD = 1.0 # Log2 fold-change threshold for categorising signal as high/low


# Allowed combinations: list of (chip_state, rna_state, atac_state) tuples
# Each state can be "high" or "low"
# Examples:
# - [("high", "low", "low"), ("low", "high", "high")]:  Current default behavior
# - [("high", "high", "high"), ("low", "low", "low")]:   Only all-high or all-low
# - [("high", "high", "high"), ("high", "high", "low"), ("high", "low", "high"), ("low", "high", "high")]:  All high-count combinations

ALLOWED_COMBINATIONS = [
    ("low","low","low"),
    ("low","low","high"),
    ("low","high","low"),
    ("low","high","high"),
    ("high","low","low"),
    ("high","low","high"),
    ("high","high","low"),
    ("high","high","high")

]

# Output files (set to None to skip saving)
OUTPUT_PATH_SSRP1 = "examples/outputs/figs/fig13.svg"
OUTPUT_PATH_SUPT16H = "examples/outputs/figs/fig14.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting circos plot generation...")
    print(f"  - Genes to include: {len(GENES)}")
    print(f"  - Threshold: {THRESHOLD}")
    print(f"  - Allowed combinations: {ALLOWED_COMBINATIONS}")

    pipe_with_tsv = CircPlotPipeline(
        csv_path=CSV_PATH,
        delimiter=",",
        gene_list=GENES,
        ssrp1_chip_col=SSRP1_CHIP_COL,
        ssrp1_rna_col=SSRP1_RNA_COL,
        ssrp1_atac_col=SSRP1_ATAC_COL,
        ssrp1_chip_basemean_col=SSRP1_CHIP_BASEMEAN_COL,
        ssrp1_rna_basemean_col=SSRP1_RNA_BASEMEAN_COL,
        ssrp1_atac_basemean_col=SSRP1_ATAC_BASEMEAN_COL,
        supt16h_chip_col=SUPT16H_CHIP_COL,
        supt16h_rna_col=SUPT16H_RNA_COL,
        supt16h_atac_col=SUPT16H_ATAC_COL,
        supt16h_chip_basemean_col=SUPT16H_CHIP_BASEMEAN_COL,
        supt16h_rna_basemean_col=SUPT16H_RNA_BASEMEAN_COL,
        supt16h_atac_basemean_col=SUPT16H_ATAC_BASEMEAN_COL,
        save_tsv=True,  # Enable TSV export
        tsv_output_dir="./examples/outputs/circplot_data/fig_13_14",  # Optional: specify output directory
        threshold=THRESHOLD,
        allowed_combinations=ALLOWED_COMBINATIONS,
        identifier_col="gene name",
        value_label="log2 Fold Change"  # Optional: specify label for transformed values
    )

    print("  - Loading CSV data...")
    print("  - Extracting log2 fold-change values...")
    print("  - Generating circos plots and exporting data...")
    results = pipe_with_tsv.run()

    # Unpack results (now includes TSV paths)
    (fig_ssrp1_data, tsv_ssrp1_path) = results[0]
    (fig_supt16h_data, tsv_supt16h_path) = results[1]

    print("  - Figures and data exported:")
    print(f"    • SSRP1 vs Control circos plot + TSV: {tsv_ssrp1_path}")
    print(f"    • SUPT16H vs Control circos plot + TSV: {tsv_supt16h_path}")

    if OUTPUT_PATH_SSRP1:
        os.makedirs(os.path.dirname(OUTPUT_PATH_SSRP1), exist_ok=True)
        fig_ssrp1_data.savefig(OUTPUT_PATH_SSRP1, dpi=300)
        print(f"    • SSRP1 figure saved to: {OUTPUT_PATH_SSRP1}")
    if OUTPUT_PATH_SUPT16H:
        os.makedirs(os.path.dirname(OUTPUT_PATH_SUPT16H), exist_ok=True)
        fig_supt16h_data.savefig(OUTPUT_PATH_SUPT16H, dpi=300)
        print(f"    • SUPT16H figure saved to: {OUTPUT_PATH_SUPT16H}")

    print("\nDisplaying plots...")
    plt.show()
