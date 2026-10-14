# Pyntegrate 🧬

[![Python 3.9-3.12](https://img.shields.io/badge/python-3.9--3.12-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python library for the integrative analysis and visualization of high-throughput sequencing data: ChIP-seq, RNA-seq, and DNase-seq/ATAC-seq.

Pyntegrate is built on the foundation of [metaseq](https://github.com/daler/metaseq), created by Ryan Dale, and extends it to analyze and cross-compare a wider variety of data types. Each common analysis is wrapped in a one-call pipeline that returns publication-ready `matplotlib` figures. Peak calling and peak annotation optionally rely on [HOMER](http://homer.ucsd.edu/homer/).

## 🚀 Quick Start

### Installation

Requires Python 3.9-3.12 (3.10 recommended) and Git:

```bash
sudo apt update
sudo apt install -y git build-essential
```

#### Option 1: Main environment (Conda)

Conda installs only what pip cannot (`bedtools`, `samtools`, and `pybedtools`, which often fails to compile under pip); `pip install -e .` then installs the remaining Python dependencies at the versions pinned in `requirements.txt`:

```bash
# Create and activate the main environment
conda create --name pyntegrate -c conda-forge -c bioconda \
    python=3.10 bedtools samtools pybedtools=0.9.1 -y
conda activate pyntegrate

# Clone and install Pyntegrate in this environment
git clone https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate.git
cd Pyntegrate
pip install -e .
```

`conda-forge` must take priority over `bioconda`. If you see errors like `excluded by strict repo priority`, fix the channel order with `conda config --add channels bioconda && conda config --add channels conda-forge`.

#### Option 2: One-command script (venv)

From a local clone, `setup_env.sh` creates a `venv`, installs Pyntegrate in editable mode, and installs `bedtools` automatically (via `conda`, `apt-get`, or `brew`, whichever is available). This is handy for development or for running the `examples/` scripts:

```bash
git clone https://github.com/Molecular-Oncology-Group-UFV/Pyntegrate.git
cd Pyntegrate
./setup_env.sh
source .venv/bin/activate
```

#### Optional: HOMER

HOMER is only needed for the peak-calling and peak-annotation functions. Follow the [HOMER installation instructions](http://homer.ucsd.edu/homer/introduction/install.html) and make sure `makeTagDirectory`, `findPeaks`, and `annotatePeaks.pl` are on your `PATH`.

#### Verify the installation

```bash
python -c "import Pyntegrate; print(Pyntegrate.__version__)"
bedtools --version
```

### Basic usage

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

## 📚 Documentation

- **[Complete Documentation](https://molecular-oncology-group-ufv.github.io/Pyntegrate/)**: overview, supported data types, and the list of pipelines
- **[Installation Guide](https://molecular-oncology-group-ufv.github.io/Pyntegrate/installation.html)**: detailed installation instructions and troubleshooting
- **[Getting Started](https://molecular-oncology-group-ufv.github.io/Pyntegrate/getting_started.html)**: step-by-step walkthrough of the main pipelines
- **[API Reference](https://molecular-oncology-group-ufv.github.io/Pyntegrate/api/pipelines.html)**: complete API documentation
- **[Examples](https://molecular-oncology-group-ufv.github.io/Pyntegrate/examples.html)**: real-world usage examples, one script per figure

## 🖼️ Usage Examples

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

More output samples for every example are in [`docs/images_examples/`](./docs/images_examples).

## 📋 Requirements

Command-line tools (not installable through pip):

| Tool | Required for |
| ---- | ------------ |
| `bedtools` | Interval operations (required for full functionality) |
| `samtools` | BAM-based workflows |
| HOMER *(optional)* | Peak calling and peak annotation |

Python dependencies, installed automatically by pip:

| Package | Version pinned in `requirements.txt` |
| ------- | ------------------------------------ |
| adjustText | `==1.1.1` |
| biopython | `==1.83` |
| bx-python | `==0.10.0` |
| Cython | `==3.0.8` |
| deeptools | any |
| fisher | `==0.1.14` |
| gffutils | `==0.12` |
| matplotlib | `==3.8.2` |
| numpy | `==1.26.3` |
| openpyxl | any |
| pandas | `==2.1.4` |
| pybedtools | `==0.9.1` |
| pyBigWig | `==0.3.22` |
| pycirclize | `==1.10.1` |
| PyYAML | `==6.0.1` |
| pysam | `==0.22.0` |
| scikit-learn | `==1.4.1.post1` |
| scipy | `==1.11.4` |
| seaborn | `==0.13.2` |

## 📄 License

This project is licensed under the **MIT License**. See the [LICENSE.txt](LICENSE.txt) file for details.
