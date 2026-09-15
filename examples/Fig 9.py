"""
Fig 9: ATAC-seq analysis with category heatmaps.

This example demonstrates the CategoryHeatmapPipeline which loads a single signal file
(ATAC-seq, ChIP-seq, BAM, BED, etc.), extracts signal around TSS, maps genes to category
labels (high/low), and produces heatmaps split by category.

Workflow:
1. Load signal file (specify file type in signal_extension)
2. Extract signal around TSS
3. Build gene → category mapping from CSV
4. Split data into category groups
5. Generate three heatmap views (unsorted, sorted by mean, sorted by peak position)
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import CategoryHeatmapPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
DB_PATH = os.path.join(DATA_DIR, "gencodehg19.gtf.db")
CATEGORY_TABLE_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")
GENOME = "hg19"
BINS = 100
WINDOW = (-1000, 2000)

CATEGORIES = ["high", "low"]  # Categories to display

## FIG 9A ##
SIGNAL_PATH_A = os.path.join(DATA_DIR, "ATAC-seq", "CONTROL.mRp.clN.bigWig")
SIGNAL_EXTENSION_A = "bigwig"     # File type: "bigwig", "bam", "bed", etc.
CATEGORY_COL_A = "CONTROL FACT category"
# Output files (set to None to skip saving). {name} is filled in with each
# figure's key ("heatmap_no_sorted", "heatmap_sorted_mean", "heatmap_sorted_max").
OUTPUT_PATH_A = "examples/outputs/figs/fig9A_{name}.svg"

## FIG 9B ##
SIGNAL_PATH_B = os.path.join(DATA_DIR, "ATAC-seq", "SSRP1.mRp.clN.bigWig")
SIGNAL_EXTENSION_B = "bigwig"     # File type: "bigwig", "bam", "bed", etc.
CATEGORY_COL_B = "SSRP1 FACT category"
# Output files (set to None to skip saving). {name} is filled in with each
# figure's key ("heatmap_no_sorted", "heatmap_sorted_mean", "heatmap_sorted_max").
OUTPUT_PATH_B = "examples/outputs/figs/fig9B_{name}.svg"

## FIG 9C ##
SIGNAL_PATH_C = os.path.join(DATA_DIR, "ATAC-seq", "SUPT16H.mRp.clN.bigWig")
SIGNAL_EXTENSION_C = "bigwig"     # File type: "bigwig", "bam", "bed", etc.
CATEGORY_COL_C = "SUPT16H FACT category"
# Output files (set to None to skip saving). {name} is filled in with each
# figure's key ("heatmap_no_sorted", "heatmap_sorted_mean", "heatmap_sorted_max").
OUTPUT_PATH_C = "examples/outputs/figs/fig9C_{name}.svg"
# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":

    pipeA = CategoryHeatmapPipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_A,
        category_table_path=CATEGORY_TABLE_PATH,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_A,
        category_col=CATEGORY_COL_A,
        bins=BINS,
        window=WINDOW,
        categories=CATEGORIES,
        figsize=(4.75, 8),
    )

    figs = pipeA.run()

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        for name, fig in figs.items():
            if name == "heatmap_sorted_mean":
                path = OUTPUT_PATH_A.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")


    pipeB = CategoryHeatmapPipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_B,
        category_table_path=CATEGORY_TABLE_PATH,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_B,
        category_col=CATEGORY_COL_B,
        bins=BINS,
        window=WINDOW,
        categories=CATEGORIES,
        figsize=(4.75, 8),
    )

    figs = pipeB.run()

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        for name, fig in figs.items():
            if name == "heatmap_sorted_mean":
                path = OUTPUT_PATH_B.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")


    pipeC = CategoryHeatmapPipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_C,
        category_table_path=CATEGORY_TABLE_PATH,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_C,
        category_col=CATEGORY_COL_C,
        bins=BINS,
        window=WINDOW,
        categories=CATEGORIES,
        figsize=(4.75, 8),
    )

    figs = pipeC.run()

    if OUTPUT_PATH_C:
        os.makedirs(os.path.dirname(OUTPUT_PATH_C), exist_ok=True)
        for name, fig in figs.items():
            if name == "heatmap_sorted_mean":
                path = OUTPUT_PATH_C.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")
    plt.show()
