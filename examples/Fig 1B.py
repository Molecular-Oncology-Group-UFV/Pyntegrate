"""
Fig 1B: ChIP/input log2 ratio correlation filtered by FACT category.

This example demonstrates ChIPCategoryCorrelationPipeline, which produces a hexbin
density scatter plot correlating the ChIP vs input log2 ratios of two FACT components
(SSRP1 and SUPT16H) at TSS windows classified as high for at least one component.

Workflow:
1. Load a merged table (TSV) containing ChIP log2FC and FACT category columns
2. Retain rows where the SSRP1 OR SUPT16H FACT category equals "high"
3. Compute the Spearman correlation between the two ChIP log2FC columns
4. Generate a hexbin density scatter plot annotated with R and sample size

"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import ChIPCategoryCorrelationPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths and column names to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table.tsv")

# X axis: SSRP1 ChIP vs input log2 fold change
XCOL = "ssrp1 ChIP log2 fold change"

# Y axis: SUPT16H ChIP vs input log2 fold change
YCOL = "suptt16h ChIP log2 fold change"

# Columns containing FACT category labels (e.g. "high" / "low")
CATEGORY_COL_1 = "SSRP1 FACT category"
CATEGORY_COL_2 = "SUPT16H FACT category"

# Keep rows where at least one FACT category equals this value
CATEGORY_VALUE = "high"

# Hexbin gridsize (number of hexagons along each axis)
BINS = 60

# Axis labels and title
XLABEL = r"$\log_2$ [SSRP1 ChIP/Input]"
YLABEL = r"$\log_2$ [SUPT16H ChIP/Input]"
TITLE = "SUPT16H / SSRP1\nChIP correlation"

# Output file (set to None to skip saving)
OUTPUT_PATH = "examples/outputs/figs/fig1B.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting ChIP category correlation analysis (Figure 1B)...")

    pipe = ChIPCategoryCorrelationPipeline(
        csv_path=CSV_PATH,
        xcol=XCOL,
        ycol=YCOL,
        category_col_1=CATEGORY_COL_1,
        category_col_2=CATEGORY_COL_2,
        category_value=CATEGORY_VALUE,
        bins=BINS,
        xlabel=XLABEL,
        ylabel=YLABEL,
        title=TITLE,
        hex_low_color="#D9D9D9",   # grey85 — matches R scale_fill_gradient low
        hex_high_color="black",
        delimiter=None,
        figsize=(7.08, 5.91),      # 12 cm × 10 cm at 300 dpi
        fontsize=22,              # base font size for labels and annotations
    )

    print("  - Loading table and filtering to FACT-high TSS windows...")
    print("  - Computing Spearman correlation...")
    print("  - Generating hexbin scatter plot...")
    fig = pipe.run()

    if OUTPUT_PATH:
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH}")

    print("\nDisplaying plot...")
    plt.show()
