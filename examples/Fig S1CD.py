"""
Fig S1CD: ChIP-seq profile split by ATAC-seq/RNA-seq log2FC quadrants.

This example demonstrates the GeneGroupProfilePipeline, which classifies genes
into four quadrants based on ATAC-seq and RNA-seq log2 fold-change thresholds
read from final_table.csv, then plots the average ChIP-seq signal around the
TSS for each quadrant as a separate line on the same axis.

The condition being analyzed (e.g. SUPT16H, SSRP1, CONTROL knockdown, ...) is
not hardcoded: it's selected by pointing SIGNAL_PATH, ATAC_COL, and RNA_COL at
that condition's ChIP-seq file and log2FC columns in the Configuration section
below.

Quadrants (condition's ATAC-seq/RNA-seq log2FC, threshold = +/-THRESHOLD):
All genes need to have padj < PADJ_THRESHOLD in both ATAC and RNA to be included in a quadrant.:
1. ATAC high (lfc > THRESHOLD)  + RNA high (lfc > THRESHOLD)
2. ATAC high (lfc > THRESHOLD)  + RNA low  (lfc < -THRESHOLD)
3. ATAC low  (lfc < -THRESHOLD) + RNA low  (lfc < -THRESHOLD)
4. ATAC low  (lfc < -THRESHOLD) + RNA high (lfc > THRESHOLD)

Workflow:
1. Load the configured condition's ChIP-seq signal (BigWig) around TSS windows
2. Load final_table.csv and assign each gene to a quadrant based on its
   configured condition's ATAC-seq and RNA-seq log2 fold change
3. Average the ChIP-seq signal per quadrant
4. Plot one line per quadrant (quadrants with no matching genes are skipped)
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import GeneGroupProfilePipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
DB_PATH = os.path.join(DATA_DIR, "gencodehg19.gtf.db")
GENOME = "hg19"
BINS = 100
WINDOW = (-1000, 1000)
CATEGORY_CSV_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")
ID_COL = "EnsemblGeneID"
THRESHOLD = 1.0  # high = lfc > THRESHOLD, low = lfc < -THRESHOLD
PADJ_THRESHOLD = 0.05  # a gene joins a quadrant only if both ATAC and RNA padj < this


### FIG S1C ###

SIGNAL_PATH_C = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_ssrp1ChIP.bw")
SIGNAL_EXTENSION_C = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
# CSV with per-gene ATAC-seq/RNA-seq log2 fold-change columns
ATAC_COL_C = "ssrp1 ATAC log2 Fold Change"
RNA_COL_C = "ssrp1 RNA log2 fold change"
# Columns with per-gene ATAC-seq/RNA-seq adjusted p-values
ATAC_PADJ_COL_C = "ssrp1 ATAC padj"
RNA_PADJ_COL_C = "ssrp1 RNA padj"
# Output file (set to None to skip saving)
OUTPUT_PATH_C = "examples/outputs/figs/figS1C.svg"


### FIG S1D ###

SIGNAL_PATH_D = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_supt16hChIP.bw")
SIGNAL_EXTENSION_D = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
# CSV with per-gene ATAC-seq/RNA-seq log2 fold-change columns
ATAC_COL_D = "supt16h ATAC log2 Fold Change"
RNA_COL_D = "supt16h RNA log2 fold change"
# Columns with per-gene ATAC-seq/RNA-seq adjusted p-values
ATAC_PADJ_COL_D = "supt16h ATAC padj"
RNA_PADJ_COL_D = "supt16h RNA padj"
# Output file (set to None to skip saving)
OUTPUT_PATH_D = "examples/outputs/figs/figS1D.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":

    pipeC = GeneGroupProfilePipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_C,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_C,
        bins=BINS,
        window=WINDOW,
        category_table_path=CATEGORY_CSV_PATH,
        atac_col=ATAC_COL_C,
        rna_col=RNA_COL_C,
        atac_padj_col=ATAC_PADJ_COL_C,
        rna_padj_col=RNA_PADJ_COL_C,
        id_col=ID_COL,
        threshold=THRESHOLD,
        padj_threshold=PADJ_THRESHOLD,
        name="SSRP1",
        fontsize=17,
    )

    fig = pipeC.run()

    if OUTPUT_PATH_C:
        os.makedirs(os.path.dirname(OUTPUT_PATH_C), exist_ok=True)
        fig.savefig(OUTPUT_PATH_C, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_C}")


    pipeD = GeneGroupProfilePipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_D,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_D,
        bins=BINS,
        window=WINDOW,
        category_table_path=CATEGORY_CSV_PATH,
        atac_col=ATAC_COL_D,
        rna_col=RNA_COL_D,
        atac_padj_col=ATAC_PADJ_COL_D,
        rna_padj_col=RNA_PADJ_COL_D,
        id_col=ID_COL,
        threshold=THRESHOLD,
        padj_threshold=PADJ_THRESHOLD,
        name="SUPT16H",
        fontsize=17,
    )

    fig = pipeD.run()

    if OUTPUT_PATH_D:
        os.makedirs(os.path.dirname(OUTPUT_PATH_D), exist_ok=True)
        fig.savefig(OUTPUT_PATH_D, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_D}")


    plt.show()
