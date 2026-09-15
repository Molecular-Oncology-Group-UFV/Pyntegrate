"""
Fig 5: Signal heatmap integrated with RNA-seq expression status, FACT-filtered.

This example demonstrates the IntegratedHeatmapPipeline which performs
a complete multi-modal integration: loads signal and RNA-seq data,
synchronises gene identifiers, optionally filters by category, and
produces an integrated heatmap with RNA up/down/unchanged overlays.

Workflow:
1. Load signal file (specify file type in signal_extension)
2. Load RNA-seq results (DEseq2-style table)
3. Normalise and synchronise gene identifiers between datasets
4. Optionally build category mapping and filter to selected groups
5. Generate integrated heatmap with RNA expression categories overlaid
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import IntegratedHeatmapPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
DB_PATH = os.path.join(DATA_DIR, "gencodehg19.gtf.db")
EXPRESSION_TABLE_PATH = os.path.join(DATA_DIR, "final_table_filtered.tsv")
GENOME = "hg19"
GENE_ID_COL = "EnsemblGeneID"
BINS = 300
WINDOW = (-1000, 2000)
PADJ_THRESHOLD = 0.05  # padjusted threshold for classifying genes as up/down-regulated


###############################################
################## FIGURE 5A ##################
###############################################

SIGNAL_PATH_A = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_ssrp1ChIP.bw")
SIGNAL_EXTENSION_A = "bigwig"     # File type: "bigwig", "bam", "bed", etc.

LFC_COL_A = "ssrp1 RNA log2 fold change"
PADJ_COL_A = "ssrp1 RNA padj"
DIVIDED_BY_A = "RNA"
SIGNAL_NAME_A = "SSRP1 - ChIP-seq"
SHOW_HEATMAPS_A = None  # True = Only shows heatmap and signal overlay. False = Only shows signal overlay divided by lfc. None = shows everything.
CATEGORY_COL_A = "SSRP1 FACT category"
CATEGORIES_A = ["high", "low"]  # Categories to display

# Output file (set to None to skip saving)
OUTPUT_PATH_A = "examples/outputs/figs/fig5A.svg"
# Up/down/unchanged gene counts as .txt (set to None to skip saving)
COUNTS_OUTPUT_PATH_A = (
    os.path.splitext(OUTPUT_PATH_A)[0] + "_counts.txt" if OUTPUT_PATH_A else None
)



###############################################
################## FIGURE 5B ##################
###############################################

SIGNAL_PATH_B = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_supt16hChIP.bw")
SIGNAL_EXTENSION_B = "bigwig"     # File type: "bigwig", "bam", "bed", etc.

LFC_COL_B = "supt16h RNA log2 fold change"
PADJ_COL_B = "supt16h RNA padj"
DIVIDED_BY_B = "RNA"
SIGNAL_NAME_B = "SUPT16H - ChIP-seq"
SHOW_HEATMAPS_B = None  # True = Only shows heatmap and signal overlay. False = Only shows signal overlay divided by lfc. None = shows everything.
CATEGORY_COL_B = "SUPT16H FACT category"
CATEGORIES_B = ["high", "low"]  # Categories to display

# Output file (set to None to skip saving)
OUTPUT_PATH_B = "examples/outputs/figs/fig5B.svg"
# Up/down/unchanged gene counts as .txt (set to None to skip saving)
COUNTS_OUTPUT_PATH_B = (
    os.path.splitext(OUTPUT_PATH_B)[0] + "_counts.txt" if OUTPUT_PATH_B else None
)

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":
    print("Starting signal + RNA integration analysis...")
    
    pipeA = IntegratedHeatmapPipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_A,
        expression_table_path=EXPRESSION_TABLE_PATH,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_A,
        gene_id_col=GENE_ID_COL,
        lfc_col=LFC_COL_A,
        divided_by = DIVIDED_BY_A,
        category_table_path=EXPRESSION_TABLE_PATH,
        category_col = CATEGORY_COL_A,
        categories=CATEGORIES_A,  
        bins=BINS,
        window=WINDOW,
        signal_name=SIGNAL_NAME_A,
        threshold=1,  # Log2 fold-change threshold for categorising RNA as up/down/unchanged
        padj_col=PADJ_COL_A,  # Column with RNA adjusted p-values
        padj_threshold=PADJ_THRESHOLD,  # Genes are up/down only if padj < this threshold
        show_heatmap_panels= SHOW_HEATMAPS_A,  # Show heatmap panels for each category
        counts_output_path=COUNTS_OUTPUT_PATH_A,
    )

    print("  - Loading ChIP signal...")
    print("  - Loading RNA-seq results...")
    print("  - Synchronising gene identifiers...")
    print("  - Building FACT category map...")
    print("  - Generating integrated heatmap...")
    fig = pipeA.run()

    print("  - Figure A generated: ChIP heatmap + RNA overlay")

    if COUNTS_OUTPUT_PATH_A:
        print(f"  - Gene counts saved to: {COUNTS_OUTPUT_PATH_A}")

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        fig.savefig(OUTPUT_PATH_A, dpi=300)
        print(f"  - Figure A saved to: {OUTPUT_PATH_A}")



    print("Starting signal + RNA integration analysis...")
    
    pipeB = IntegratedHeatmapPipeline(
        db_path=DB_PATH,
        signal_path=SIGNAL_PATH_B,
        expression_table_path=EXPRESSION_TABLE_PATH,
        genome=GENOME,
        signal_extension=SIGNAL_EXTENSION_B,
        gene_id_col=GENE_ID_COL,
        lfc_col=LFC_COL_B,
        divided_by = DIVIDED_BY_B,
        category_table_path=EXPRESSION_TABLE_PATH,
        category_col = CATEGORY_COL_B,
        categories=CATEGORIES_B,  
        bins=BINS,
        window=WINDOW,
        signal_name=SIGNAL_NAME_B,
        threshold=1,  # Log2 fold-change threshold for categorising RNA as up/down/unchanged
        padj_col=PADJ_COL_B,  # Column with RNA adjusted p-values
        padj_threshold=PADJ_THRESHOLD,  # Genes are up/down only if padj < this threshold
        show_heatmap_panels= SHOW_HEATMAPS_B,  # Show heatmap panels for each category
        counts_output_path=COUNTS_OUTPUT_PATH_B,
    )

    print("  - Loading ChIP signal...")
    print("  - Loading RNA-seq results...")
    print("  - Synchronising gene identifiers...")
    print("  - Building FACT category map...")
    print("  - Generating integrated heatmap...")
    fig = pipeB.run()

    print("  - Figure B generated: ChIP heatmap + RNA overlay")

    if COUNTS_OUTPUT_PATH_B:
        print(f"  - Gene counts saved to: {COUNTS_OUTPUT_PATH_B}")

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        fig.savefig(OUTPUT_PATH_B, dpi=300)
        print(f"  - Figure B saved to: {OUTPUT_PATH_B}")

    print("\nDisplaying plot...")
    plt.show()
