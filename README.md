# Pyntegrate

Pyntegrate is a library for genomic data analysis, built on the foundation of another library called [metaseq](https://github.com/daler/metaseq)  created by Ryan Dale. The update of this library is driven by the need to analyze a wider variety of data types. The types of data supported by Pyntegrate are:

1.  ChIP-seq
2.  RNA-seq
3.  DNase-seq/ATAC-seq

With these types of data accepted by the library, various functions have been created to enable the generation of graphs with them and to allow for cross-analysis.

This library also utilizes software called [HOMER](http://homer.ucsd.edu/homer/).It is especially used for identifying ChIP-seq peak enrichments and distinguishing between different types of DNA regions based on intervals.


## Installation

Requires Python 3.9+ and Git.

```console
sudo apt update
sudo apt install -y git build-essential
```

Install directly from GitHub:

```console
pip install git+https://github.com/SalvaaGomez/Pyntegrate.git
```

Or from a local clone, in editable mode (recommended for development or for running the `examples/` scripts):

```console
git clone https://github.com/SalvaaGomez/Pyntegrate.git
cd Pyntegrate
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Or with Conda, installing the runtime dependencies explicitly:

```console
conda create --name pyntegrateTest python=3.10 -y
conda activate pyntegrate_final

conda install -c conda-forge numpy=1.26.3 -y
conda install -c bioconda bedtools -y
conda install pandas -y
conda install openpyxl -y
conda install bioconda::samtools -y
conda install conda-forge::pycirclize -y

cd ~/pyntegrate
git clone git@github.com:SalvaaGomez/Pyntegrate.git
cd Pyntegrate/
pip install -e .
```

For the full functionality, `bedtools` is also required:
- Using conda: `conda install -c bioconda bedtools`
- Without conda: `sudo apt-get install bedtools`

See [docs/installation.md](./docs/installation.md) for a Conda-based walkthrough and troubleshooting.

## HOMER Installation

To use the HOMER-dependent functions (optional), follow the installation instructions at [this link](http://homer.ucsd.edu/homer/introduction/install.html).

## Getting Started

Pyntegrate accepts genomic data with the following extensions:

- `bigwig`, `bam`: full functionality (ChIP-seq, ATAC-seq)
- `bed`, `bigbed`, `gff`, `gtf`, `vcf`: most functionality, but no read counting
- `csv`: full functionality (RNA-seq)

The recommended entry point is the `Pyntegrate.pipelines` module: each pipeline class wraps loading, processing, and plotting for one analysis type into a single call. For example, comparing ChIP-seq IP vs INPUT signal around TSS:

```python
from Pyntegrate.pipelines import PeakHeatmapPipeline

pipe = PeakHeatmapPipeline(
    db_path="/path/to/annotations.gtf.db",   # gffutils FeatureDB
    ip_path="/path/to/ip_signal.bw",
    input_path="/path/to/input_signal.bw",
    genome="mm10",
    ip_extension="bigwig",
    input_extension="bigwig",
    bins=100,
    window=(-1000, 1000),
)
figs = pipe.run()          # dict of matplotlib Figures
figs["profile"].savefig("profile.png", dpi=300)
```

For RNA-seq differential-expression tables:

```python
from Pyntegrate.results_table import DEseq2Results

data = DEseq2Results(
    "/path/to/rna.csv",
    id_column="GeneID", pval_column="padj",
    lfc_column="log2FoldChange", mean_column="baseMean",
)
```

The full set of `SeqSignalAnalysis` functions used internally by the pipelines is documented in [docs/api/seqsignalanalysis.md](./docs/api/seqsignalanalysis.md).

## Usage Examples

The [`examples/`](./examples) directory contains runnable, self-contained scripts reproducing a FACT (SSRP1/SUPT16H) knockdown ChIP-seq/ATAC-seq/RNA-seq analysis, one script per figure (or figure panel), each with a configuration section at the top and a docstring explaining the workflow. See [`docs/examples.md`](./docs/examples.md) for the full gallery and a quick-reference table.

| Figure | Script | Pipeline | Output |
|--------|--------|----------|--------|
| 1A | [Fig 1A.py](<./examples/Fig 1A.py>) | `GeneScaleRegionsProfilePipeline` | Scaled-gene-body meta-profile |
| 2A/2B | [Fig 2.py](<./examples/Fig 2.py>) | `IntegratedHeatmapPipeline` | Signal heatmap with RNA up/down overlay |
| 10A/10B | [Fig 10.py](<./examples/Fig 10.py>) | `VolcanoPipeline` | ChIP-seq volcano plots |
| 13/14 | [Fig 13-14.py](<./examples/Fig 13-14.py>) | `CircPlotPipeline` | Multi-signal circos comparison |

Sample output from these pipelines:

<p>
  <img src="./docs/images_examples/fig1a-genescale-profile.png" alt="Meta-gene profile" width="490">
  <img src="./docs/images_examples/fig2a-integrated-heatmap-rna-ssrp1.png" alt="Integrated ChIP + RNA heatmap" width="330">
  <img src="./docs/images_examples/fig10a-volcano-ssrp1-chip.png" alt="Volcano plot" height="320">
  <img src="./docs/images_examples/fig13-circplot-ssrp1.png" alt="Circos plot" height="320">
</p>

More output samples for every example are in [`docs/images_examples/`](./docs/images_examples), and the full walkthrough for each figure is in [`docs/examples.md`](./docs/examples.md).
