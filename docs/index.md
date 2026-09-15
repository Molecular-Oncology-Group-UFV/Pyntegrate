# Pyntegrate

**Integrative analysis of high-throughput sequencing data.**

Pyntegrate is a Python library for genomic data analysis, built on the foundation of [metaseq](https://github.com/daler/metaseq), created by Ryan Dale. It extends that foundation to analyze and cross-compare a wider variety of data types — ChIP-seq, RNA-seq, and ATAC-seq — in a single, consistent API.

[:material-download: Install Pyntegrate](installation.md){ .md-button .md-button--primary }
[:material-rocket-launch: Getting Started](getting_started.md){ .md-button }

---

## Supported Data Types

| Data type | Typical input | What Pyntegrate does with it |
|---|---|---|
| **ChIP-seq** | BigWig / BAM | Signal extraction around TSS, IP − INPUT peak calling, heatmaps |
| **RNA-seq** | CSV (DESeq2-style) | Differential expression parsing, volcano plots, decile/violin analysis |
| **DNase-seq / ATAC-seq** | BigWig / BAM | Accessibility signal extraction, category heatmaps |

## Key Features

- :material-set-merge: **Multi-omics integration** — analyze ChIP-seq, RNA-seq, and ATAC-seq data together in one pipeline
- :material-chart-bell-curve: **Peak analysis** — integrated with [HOMER](http://homer.ucsd.edu/homer/) for ChIP-seq peak enrichment
- :material-chart-areaspline: **Visualization** — heatmaps, volcano plots, decile plots, circos plots, and more, out of the box
- :material-dna: **Gene annotation** — classify genomic regions (promoter, intergenic, etc.) from interval data

## Quick Start

To get started with Pyntegrate, you'll need:

1. Genomic signal data (BigWig or BAM files)
2. Genomic feature annotations (BED, GTF, or GFF files)
3. Gene expression data (CSV files for RNA-seq)

With that in place, Pyntegrate lets you:

- Load and analyze ChIP-seq and ATAC-seq signal
- Integrate RNA-seq expression data
- Generate comparative, publication-ready visualizations
- Perform cross-omics correlation analysis

!!! tip "Next Steps"
    - [:material-download: Install Pyntegrate](installation.md)
    - [:material-book-open-page-variant: Follow the Getting Started guide](getting_started.md)
    - [:material-image-multiple: Explore the Examples gallery](examples.md)

## License

Pyntegrate is distributed under the MIT License. See [LICENSE.txt](https://github.com/SalvaaGomez/Pyntegrate/blob/main/LICENSE.txt) for details.

## Support

For issues, feature requests, and contributions, visit the [GitHub repository](https://github.com/SalvaaGomez/Pyntegrate).
