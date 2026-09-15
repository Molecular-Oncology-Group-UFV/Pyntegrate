"""
Fig S2: Signal overlay profiles, FACT high/low split.

This example demonstrates the OverlayProfilePipeline which loads two signal files
(BigWig, BAM, BED, etc.), retains only genes present in both, and produces dual-y-axis
profile figures showing both signals.

Optionally filters by genomic categories (e.g., FACT high/low) to visualize category-specific
profiles side by side.

Workflow:
1. Load two signals independently (specify file types in signal_1_extension and signal_2_extension)
2. Find common genes between the two datasets
3. Generate signal arrays for common genes
4. Optionally filter by category (FACT high/low)
5. Plot average profiles on twin y-axes for easy comparison
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import OverlayProfilePipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
DB_PATH = os.path.join(DATA_DIR, "gencodehg19.gtf.db")
SIGNAL_EXTRA_1_PATH = os.path.join(DATA_DIR, "ATAC-seq", "CONTROL.mRp.clN.bigWig")
SIGNAL_EXTRA_1_EXTENSION = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
GENOME = "hg19"
BINS = 100
WINDOW = (-1000, 1000)

SIGNAL_1_PATH_A = os.path.join(DATA_DIR, "ATAC-seq", "SSRP1.mRp.clN.bigWig")
SIGNAL_2_PATH_A = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_ssrp1ChIP.bw")
SIGNAL_1_EXTENSION_A = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
SIGNAL_2_EXTENSION_A = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
OUTPUT_PATH_A = "examples/outputs/figs/figS2A.svg" # Output file (set to None to skip saving)
TITLE_A = "SSRP1 ChIP-seq-ATAC-seq HIGH FACT"
CATEGORIES_A = ["high"]  # Categories to display


SIGNAL_1_PATH_B = os.path.join(DATA_DIR, "ATAC-seq", "SSRP1.mRp.clN.bigWig")
SIGNAL_2_PATH_B = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_ssrp1ChIP.bw")
SIGNAL_1_EXTENSION_B = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
SIGNAL_2_EXTENSION_B = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
OUTPUT_PATH_B = "examples/outputs/figs/figS2B.svg" # Output file (set to None to skip saving)
TITLE_B = "SSRP1 ChIP-seq-ATAC-seq LOW FACT"
CATEGORIES_B = ["low"]  # Categories to display


SIGNAL_1_PATH_C = os.path.join(DATA_DIR, "ATAC-seq", "SUPT16H.mRp.clN.bigWig")
SIGNAL_2_PATH_C = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_supt16hChIP.bw")
SIGNAL_1_EXTENSION_C = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
SIGNAL_2_EXTENSION_C = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
OUTPUT_PATH_C = "examples/outputs/figs/figS2C.svg" # Output file (set to None to skip saving)
TITLE_C = "SUPT16H ChIP-seq-ATAC-seq HIGH FACT"
CATEGORIES_C = ["high"]  # Categories to display


SIGNAL_1_PATH_D = os.path.join(DATA_DIR, "ATAC-seq", "SUPT16H.mRp.clN.bigWig")
SIGNAL_2_PATH_D = os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_supt16hChIP.bw")
SIGNAL_1_EXTENSION_D = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
SIGNAL_2_EXTENSION_D = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
OUTPUT_PATH_D = "examples/outputs/figs/figS2D.svg" # Output file (set to None to skip saving)
TITLE_D = "SUPT16H ChIP-seq-ATAC-seq LOW FACT"
CATEGORIES_D = ["low"]  # Categories to display


# ============================================================================
# Run Pipeline - Scenario 1: All genes (no category filtering)
# ============================================================================
if __name__ == "__main__":

    pipeA = OverlayProfilePipeline(
        db_path=DB_PATH,
        signal_1_path=SIGNAL_1_PATH_A,
        signal_2_path=SIGNAL_2_PATH_A,
        genome=GENOME,
        signal_1_extension=SIGNAL_1_EXTENSION_A,
        signal_2_extension=SIGNAL_2_EXTENSION_A,
        bins=BINS,
        window=WINDOW,
        name=TITLE_A,
        fontsize=15,
        extra_signal_1_path=SIGNAL_EXTRA_1_PATH,
        extra_signal_1_extension=SIGNAL_EXTRA_1_EXTENSION,
        label_extra_sig1="CONTROL ATAC",
        category_table_path=os.path.join(DATA_DIR, "final_table_filtered.tsv"),
        category_col="SSRP1 FACT category",
        categories=CATEGORIES_A,
    )

    fig = pipeA.run()

    if OUTPUT_PATH_A:
        os.makedirs(os.path.dirname(OUTPUT_PATH_A), exist_ok=True)
        fig.savefig(OUTPUT_PATH_A, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_A}")



    pipeB = OverlayProfilePipeline(
        db_path=DB_PATH,
        signal_1_path=SIGNAL_1_PATH_B,
        signal_2_path=SIGNAL_2_PATH_B,
        genome=GENOME,
        signal_1_extension=SIGNAL_1_EXTENSION_B,
        signal_2_extension=SIGNAL_2_EXTENSION_B,
        bins=BINS,
        window=WINDOW,
        name=TITLE_B,
        fontsize=15,
        extra_signal_1_path=SIGNAL_EXTRA_1_PATH,
        extra_signal_1_extension=SIGNAL_EXTRA_1_EXTENSION,
        label_extra_sig1="CONTROL ATAC",
        category_table_path=os.path.join(DATA_DIR, "final_table_filtered.tsv"),
        category_col="SSRP1 FACT category",
        categories=CATEGORIES_B,

    )

    fig = pipeB.run()

    if OUTPUT_PATH_B:
        os.makedirs(os.path.dirname(OUTPUT_PATH_B), exist_ok=True)
        fig.savefig(OUTPUT_PATH_B, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_B}")


    pipeC = OverlayProfilePipeline(
        db_path=DB_PATH,
        signal_1_path=SIGNAL_1_PATH_C,
        signal_2_path=SIGNAL_2_PATH_C,
        genome=GENOME,
        signal_1_extension=SIGNAL_1_EXTENSION_C,
        signal_2_extension=SIGNAL_2_EXTENSION_C,
        bins=BINS,
        window=WINDOW,
        name=TITLE_C,
        fontsize=15,
        extra_signal_1_path=SIGNAL_EXTRA_1_PATH,
        extra_signal_1_extension=SIGNAL_EXTRA_1_EXTENSION,
        label_extra_sig1="CONTROL ATAC",
        category_table_path=os.path.join(DATA_DIR, "final_table_filtered.tsv"),
        category_col="SUPT16H FACT category",
        categories=CATEGORIES_C,
    )

    fig = pipeC.run()

    if OUTPUT_PATH_C:
        os.makedirs(os.path.dirname(OUTPUT_PATH_C), exist_ok=True)
        fig.savefig(OUTPUT_PATH_C, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_C}")



    pipeD = OverlayProfilePipeline(
        db_path=DB_PATH,
        signal_1_path=SIGNAL_1_PATH_D,
        signal_2_path=SIGNAL_2_PATH_D,
        genome=GENOME,
        signal_1_extension=SIGNAL_1_EXTENSION_D,
        signal_2_extension=SIGNAL_2_EXTENSION_D,
        bins=BINS,
        window=WINDOW,
        name=TITLE_D,
        fontsize=15,
        extra_signal_1_path=SIGNAL_EXTRA_1_PATH,
        extra_signal_1_extension=SIGNAL_EXTRA_1_EXTENSION,
        label_extra_sig1="CONTROL ATAC",
        category_table_path=os.path.join(DATA_DIR, "final_table_filtered.tsv"),
        category_col="SUPT16H FACT category",
        categories=CATEGORIES_D,
    )

    fig = pipeD.run()

    if OUTPUT_PATH_D:
        os.makedirs(os.path.dirname(OUTPUT_PATH_D), exist_ok=True)
        fig.savefig(OUTPUT_PATH_D, dpi=300)
        print(f"  - Figure saved to: {OUTPUT_PATH_D}")

    plt.show()
