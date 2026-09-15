"""
Extra: Multi-signal comparison (3-way scatter plot).

This example demonstrates the MultiSignalComparisonPipeline which compares
log2 fold-change values across multiple modalities (e.g., ChIP-seq, ATAC-seq,
RNA-seq) and produces 2D/3D scatter plots showing the relationships between
signal changes.

Workflow:
1. Load expression results CSV with multiple log2 fold-change columns
2. Extract ChIP, ATAC, and RNA log2FC values for each gene
3. Generate scatter plots comparing all three modalities
4. Visualise correlations between signal and expression changes
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import MultiSignalComparisonPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data")

# ============================================================================
# Configuration (adjust these paths and column names to your data)
# ============================================================================
EXPRESSION_TABLE_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")
CHIP_LFC_COL = "ssrp1 ChIP log2 fold change"
ATAC_LFC_COL = "ssrp1 ATAC log2 Fold Change"
RNA_LFC_COL = "ssrp1 RNA log2 fold change"
TITLE = "Multi-Signal Comparison: ChIP vs ATAC vs RNA"

# Output file (set to None to skip saving)
OUTPUT_PATH = "examples/outputs/figs/fig_extra.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting multi-signal comparison analysis...")

    pipe = MultiSignalComparisonPipeline(
        expression_table_path=EXPRESSION_TABLE_PATH,
        lfc_col_1=CHIP_LFC_COL,
        lfc_col_2=ATAC_LFC_COL,
        lfc_col_rna=RNA_LFC_COL,
        title=TITLE,
    )

    print("  - Loading expression results...")
    print("  - Extracting log2 fold-change columns...")
    print("  - Generating 3-way comparison scatter plots...")
    fig = pipe.run()

    print("  - Figure generated: Multi-signal scatter plot")

    if OUTPUT_PATH:
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH}")

    print("\nDisplaying plot...")
    plt.show()
