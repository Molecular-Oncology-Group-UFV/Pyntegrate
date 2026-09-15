"""
Fig 6: RNA violin panels in a single horizontal figure.

This example only configures the panel list and calls the unified
`RNAViolinPipeline`. The number of violins is controlled by the length
of `VIOLIN_SPECS`:

* one entry -> one violin
* four entries -> four violins
* any other number -> that many violins on a single axes

All violins are drawn horizontally in one figure, with no subplot
borders between them.
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import RNAViolinPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (edit the list to control how many violins appear)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")

VIOLIN_SPECS = [
    {
        "title": "SSRP1 High",
        "rna_lfc_col": "ssrp1 RNA log2 fold change",
        "rna_padj_col": "ssrp1 RNA padj",
        "category_col": "SSRP1 FACT category",
        "category_filter": "high",
        "padj_threshold": 0.1,
    },
    {
        "title": "SUPT16H High",
        "rna_lfc_col": "supt16h RNA log2 fold change",
        "rna_padj_col": "supt16h RNA padj",
        "category_col": "SUPT16H FACT category",
        "category_filter": "high",
        "padj_threshold": 0.1,
    },
    {
        "title": "SSRP1 Low",
        "rna_lfc_col": "ssrp1 RNA log2 fold change",
        "rna_padj_col": "ssrp1 RNA padj",
        "category_col": "SSRP1 FACT category",
        "category_filter": "low",
        "padj_threshold": 0.1,
    },
    {
        "title": "SUPT16H Low",
        "rna_lfc_col": "supt16h RNA log2 fold change",
        "rna_padj_col": "supt16h RNA padj",
        "category_col": "SUPT16H FACT category",
        "category_filter": "low",
        "padj_threshold": 0.1,
    },
]

# Output file (set to None to skip saving)
OUTPUT_PATH = "examples/outputs/figs/fig6.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting RNA violin analysis...")
    print(f"  - Violin count: {len(VIOLIN_SPECS)}")

    pipe = RNAViolinPipeline(
        csv_path=CSV_PATH,
        panel_specs=VIOLIN_SPECS,
        figure_title="RNA expression distribution panels",
        ylabel="RNA log2 fold change",
        title="RNA Expression Distribution",
        label_fontsize=20,
        figsize_per_panel = (4, 7)
    )

    print("  - Loading data and building a single horizontal figure...")
    fig = pipe.run()

    print("  - Figure generated")

    if OUTPUT_PATH:
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        fig.savefig(OUTPUT_PATH, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH}")

    print("\nDisplaying plot...")
    plt.show()
