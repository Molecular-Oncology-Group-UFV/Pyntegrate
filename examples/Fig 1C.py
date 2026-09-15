"""
Fig 1C: RNA-seq knockdown log2FC correlation for significant genes.

This example demonstrates RNAKnockdownCorrelationPipeline, which produces a hexbin
density scatter plot correlating the RNA-seq log2 fold-changes of two knockdown
conditions (SSRP1 KD and SUPT16H KD) for genes that are significantly changed in
at least one condition.

Workflow:
1. Load a merged table (TSV) containing RNA log2FC and padj columns for both conditions
2. Retain genes where both padj values are non-missing AND at least one padj < threshold
3. Compute the Spearman correlation between the two RNA log2FC columns
4. Generate a hexbin density scatter plot annotated with R and sample size

"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import RNAKnockdownCorrelationPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths and column names to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table.tsv")

# X axis: SSRP1 knockdown RNA-seq log2 fold change
XCOL = "ssrp1 RNA log2 fold change"

# Y axis: SUPT16H knockdown RNA-seq log2 fold change
YCOL = "supt16h RNA log2 fold change"

# Adjusted p-value columns for each knockdown condition
PADJ_COL_1 = "ssrp1 RNA padj"
PADJ_COL_2 = "supt16h RNA padj"

# Significance threshold: keep genes where at least one padj < threshold
PADJ_THRESHOLD = 0.1

# Hexbin gridsize
BINS = 50

# Axis labels and title
XLABEL = r"$\log_2$ [SSRP1 RNA]"
YLABEL = r"$\log_2$ [SUPT16H RNA]"
TITLE = "RNA-seq upon FACT depletion"

# Output file (set to None to skip saving)
OUTPUT_PATH = "examples/outputs/figs/fig1C.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting RNA-seq knockdown correlation analysis (Figure 1C)...")

    pipe = RNAKnockdownCorrelationPipeline(
        csv_path=CSV_PATH,
        xcol=XCOL,
        ycol=YCOL,
        padj_col_1=PADJ_COL_1,
        padj_col_2=PADJ_COL_2,
        padj_threshold=PADJ_THRESHOLD,
        bins=BINS,
        xlabel=XLABEL,
        ylabel=YLABEL,
        title=TITLE,
        hex_low_color="#B3B3B3",   # grey70 — matches R scale_fill_gradient low
        hex_high_color="black",
        delimiter=None,
        figsize=(7.08, 5.91),      # 12 cm × 10 cm at 300 dpi
        fontsize=22,              # base font size for labels and annotations
    )

    print("  - Loading table and filtering to significant genes...")
    print(f"    (padj < {PADJ_THRESHOLD} in at least one condition)")
    print("  - Computing Spearman correlation...")
    print("  - Generating hexbin scatter plot...")
    fig = pipe.run()

    if OUTPUT_PATH:
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH}")

    print("\nDisplaying plot...")
    plt.show()
