"""
Fig 1E: ChIP-seq IP vs INPUT analysis.

This example demonstrates the PeakHeatmapPipeline which loads IP and INPUT signal files
(BigWig, BAM, BED, etc.), generates signal arrays around TSS, computes IP - INPUT subtraction,
and produces profile plots and sorted heatmaps.

Workflow:
1. Load ChIP IP and INPUT signal files (specify file type in ip_extension and input_extension)
2. Extract signal around transcription start sites (TSS)
3. Compute normalized profiles
4. Calculate peaks (IP - INPUT)
5. Generate three heatmap views: unsorted, sorted by mean, sorted by peak position
"""

import os

import matplotlib.pyplot as plt
from Pyntegrate.pipelines import PeakHeatmapPipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

# ============================================================================
# Configuration (adjust these paths to your data)
# ============================================================================
DB_PATH = os.path.join(DATA_DIR, "gencode.mm10.annotation.gtf.db")
IP_BW = os.path.join(DATA_DIR, "bam", "SRR1204544_sort_nondup.bw")
INPUT_BW = os.path.join(DATA_DIR, "bam", "SRR1204546_sort_nondup.bw")
GENOME = "mm10"
IP_EXTENSION = "bigwig"      # File type: "bigwig", "bam", "bed", etc.
INPUT_EXTENSION = "bigwig"   # File type: "bigwig", "bam", "bed", etc.
BINS = 100
WINDOW = (-1000, 1000)

# Output files (set to None to skip saving). {name} is filled in with each
# figure's key ("profile", "heatmap_no_sorted", "heatmap_sorted_mean", "heatmap_sorted_max").
OUTPUT_PATH = "examples/outputs/figs/fig1E_{name}.svg"

# ============================================================================
# Run Pipeline
# ============================================================================
if __name__ == "__main__":

    pipe = PeakHeatmapPipeline(
        db_path=DB_PATH,
        ip_path=IP_BW,
        input_path=INPUT_BW,
        genome=GENOME,
        ip_extension=IP_EXTENSION,
        input_extension=INPUT_EXTENSION,
        bins=BINS,
        window=WINDOW,
    )

    figs = pipe.run()

    if OUTPUT_PATH:
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        for name, fig in figs.items():
            if name == "profile":
                path = OUTPUT_PATH.format(name=name)
                fig.savefig(path, dpi=300)
                print(f"  - Figure saved to: {path}")

    plt.show()
