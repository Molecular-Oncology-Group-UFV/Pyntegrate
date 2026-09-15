"""
Fig 1A: Gene meta-profile with scaled gene bodies.

This example demonstrates the GeneScaleRegionsProfilePipeline, which produces
a deepTools `plotProfile --scaleRegions`-style meta-gene figure: each gene is
scaled to a common body length, flanked by fixed-length upstream/downstream
regions, and the average signal across all annotated genes is plotted for
each input signal file (e.g. SSRP1 and SUPT16H ChIP-seq).

Workflow:
1. Load an annotated gene BED file and one or more BigWig/BAM signal files
2. Scale each gene body to `region_body_length` bins, with fixed-length
   flanking regions before/after
3. Average the per-gene signal across all genes for each signal file
4. Plot one line per signal file, optionally saving the figure to disk
"""
import os

from Pyntegrate.pipelines import GeneScaleRegionsProfilePipeline

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")

GENES_BED_PATH = os.path.join(DATA_DIR, "genes.bed.gz")
SIGNAL_PATHS = [
    os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_ssrp1ChIP.bw"),
    os.path.join(DATA_DIR, "ChIP-seq", "GSE98758_supt16hChIP.bw"),
]
SAMPLE_LABELS = ["SSRP1", "SUPT16H"]
COLORS = ["#0b1b60", "#7cc576"]
OUTPUT_PATH = "examples/outputs/figs/fig1A.svg"
REGIONS_LABEL = "Meta analysis of annotated genes. n={n}"

GeneScaleRegionsProfilePipeline(
    genes_bed_path=GENES_BED_PATH,
    signal_paths=SIGNAL_PATHS,
    sample_labels=SAMPLE_LABELS,
    colors=COLORS,
    output_path=OUTPUT_PATH,
    regions_label=REGIONS_LABEL,
    before_region_start_length=500,
    region_body_length=1000,
    after_region_start_length=500,
    bin_size=10,
    chunk_size=5000,
    skip_zeros=True,
    processes=4,
    figsize=(7.7244, 3.1496),
    ylabel="ChIP-seq Signal",
    title="",
    legend_location="upper right",
    show=True,
).run()
