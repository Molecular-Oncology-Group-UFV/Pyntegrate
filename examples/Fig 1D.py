"""
Fig 1D: ATAC-seq knockdown log2FC correlation for significant peaks.

This example demonstrates ATACKnockdownCorrelationPipeline, which produces a hexbin
density scatter plot correlating the ATAC-seq log2 fold-changes of two knockdown
conditions (SUPT16H KD vs SSRP1 KD) for ATAC-seq peaks significantly differential
in at least one condition.

Workflow:
1. Load a merged table (TSV) containing ATAC log2FC and padj columns for both conditions
2. Remove rows with any missing value in the LFC or padj columns
3. Retain peaks where at least one padj < threshold
4. Compute the Spearman correlation between the two ATAC log2FC columns
5. Generate a hexbin density scatter plot annotated with R and sample size

"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import ATACKnockdownCorrelationPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths and column names to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "merged_deseq2_results.tsv")

# X axis: SUPT16H knockdown ATAC-seq log2 fold change  (note: SUPT16H on X in Fig 6D)
XCOL = "supt16h ATAC log2 Fold Change"

# Y axis: SSRP1 knockdown ATAC-seq log2 fold change
YCOL = "ssrp1 ATAC log2 Fold Change"

# Adjusted p-value columns for each knockdown condition
PADJ_COL_1 = "ssrp1 ATAC padj"
PADJ_COL_2 = "supt16h ATAC padj"

# Significance threshold: keep peaks where at least one padj < threshold
PADJ_THRESHOLD = 0.05

# Hexbin gridsize
BINS = 35

# Axis labels and title
XLABEL = r"$\log_2$ [SUPT16H ATAC]"
YLABEL = r"$\log_2$ [SSRP1 ATAC]"
TITLE = "ATAC-seq upon FACT depletion"

# Output file (set to None to skip saving)
OUTPUT_PATH = "examples/outputs/figs/fig1D.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting ATAC-seq knockdown correlation analysis (Figure 1D)...")

    pipe = ATACKnockdownCorrelationPipeline(
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
        hex_low_color="#B3B3B3",  
        hex_high_color="black",
        delimiter=None,
        figsize=(7.08, 5.91),      
        fontsize=22,              # base font size for labels and annotations
    )

    print("  - Loading table and filtering to significant peaks...")
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
