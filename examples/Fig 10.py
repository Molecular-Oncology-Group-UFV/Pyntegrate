"""
Fig 10: Volcano plot from expression results.

This example demonstrates the VolcanoPipeline which loads a CSV file
with expression results and produces volcano plots with optional intensity-based
alpha transparency and gene labeling.

The pipeline can optionally split data by a categorical column (e.g., FACT category)
to generate separate plots for each category plus a combined plot.
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import VolcanoPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")
GENE_ID_COL = "EnsemblGeneID"
GENE_NAME_COL = "gene name"

# Thresholds for significance
LOG2FC_UPPER = 1.5
LOG2FC_LOWER = -1.5
NLOG10_UPPER = 2.0

# Thresholds for gene labelling
TEXT_NLOG10 = None  # Set to None to skip labelling by threshold
TEXT_LOG2FC = None  # Set to None to skip labelling by threshold

HIGHLIGHT_GENES = ["SERPINE1",
    "F3",
    "CYR61",
    "TERT",
    "RPSAP52",
    "NEAT1"]


PADJ_COL_A = "ssrp1 ChIP padj"
LFC_COL_A = "ssrp1 ChIP log2 fold change"
INTENSITY_COL_A = "ssrp1 ChIP baseMean"
TITLE_A = "SSRP1 ChIP"
OUTPUT_PATH_A = "examples/outputs/figs/fig10A.svg"


PADJ_COL_B = "supt16h ChIP padj"
LFC_COL_B = "suptt16h ChIP log2 fold change"
INTENSITY_COL_B = "supt16h ChIP baseMean"
TITLE_B = "SUPT16H ChIP"
OUTPUT_PATH_B = "examples/outputs/figs/fig10B.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting volcano plot generation...")

    pipeA = VolcanoPipeline(
        csv_path=CSV_PATH,
        gene_id_col=GENE_ID_COL,
        gene_name_col=GENE_NAME_COL,
        padj_col=PADJ_COL_A,
        lfc_col=LFC_COL_A,
        intensity_col=INTENSITY_COL_A,
        log2fc_upper=LOG2FC_UPPER,
        log2fc_lower=LOG2FC_LOWER,
        nlog10_upper=NLOG10_UPPER,
        text_nlog10=TEXT_NLOG10,
        text_log2fc=TEXT_LOG2FC,
        highlight_genes=HIGHLIGHT_GENES,
        fontsize=20,
        title=TITLE_A,
    )

    print("  - Loading expression results...")
    print("  - Generating volcano plots...")
    figs = pipeA.run()

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        if isinstance(figs, dict):
            for name, fig in figs.items():
                path = OUTPUT_PATH_A.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")
        else:
            path = OUTPUT_PATH_A.format(name="all")
            figs.savefig(path, dpi=300)
            print(f"  - Figure saved to: {path}")



    pipeB = VolcanoPipeline(
        csv_path=CSV_PATH,
        gene_id_col=GENE_ID_COL,
        gene_name_col=GENE_NAME_COL,
        padj_col=PADJ_COL_B,
        lfc_col=LFC_COL_B,
        intensity_col=INTENSITY_COL_B,
        log2fc_upper=LOG2FC_UPPER,
        log2fc_lower=LOG2FC_LOWER,
        nlog10_upper=NLOG10_UPPER,
        text_nlog10=TEXT_NLOG10,
        text_log2fc=TEXT_LOG2FC,
        highlight_genes=HIGHLIGHT_GENES,
        fontsize=20,
        title=TITLE_B,
    )

    print("  - Loading expression results...")
    print("  - Generating volcano plots...")
    figs = pipeB.run()

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        if isinstance(figs, dict):
            for name, fig in figs.items():
                path = OUTPUT_PATH_B.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")
        else:
            path = OUTPUT_PATH_B.format(name="all")
            figs.savefig(path, dpi=300)
            print(f"  - Figure saved to: {path}")

    print("\nDisplaying plots...")
    plt.show()
