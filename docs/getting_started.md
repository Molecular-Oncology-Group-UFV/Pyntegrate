# Getting Started with Pyntegrate

This guide walks through the basic Pyntegrate workflow: pointing a pipeline at your files and getting back ready-to-save `matplotlib` figures.

!!! tip "Prefer runnable code?"
    Every pattern shown here also exists as a complete, runnable script in [`examples/`](https://github.com/SalvaaGomez/Pyntegrate/tree/main/examples) — see the [Examples gallery](examples.md).

## Basic Workflow

The recommended entry point is the `Pyntegrate.pipelines` module: each pipeline class wraps loading, processing, and plotting for one analysis type into a single call. You instantiate a pipeline with file paths and parameters, then call `.run()` to get back one or more `matplotlib.figure.Figure` objects — no separate "loading" step is required, the pipeline reads the files itself.

### 1. ChIP-seq IP vs. INPUT (`PeakHeatmapPipeline`)

```python
from Pyntegrate.pipelines import PeakHeatmapPipeline

pipeline = PeakHeatmapPipeline(
    db_path="/path/to/annotations.gtf.db",  # gffutils FeatureDB
    ip_path="/path/to/ip_signal.bw",
    input_path="/path/to/input_signal.bw",
    genome="mm10",
    ip_extension="bigwig",
    input_extension="bigwig",
    bins=100,
    window=(-1000, 1000),
)
figs = pipeline.run()  # dict of Figures
figs["profile"].savefig("profile.png", dpi=300)
figs["heatmap_sorted_mean"].savefig("heatmap.png", dpi=300)
```

| Parameter | Description |
|---|---|
| `db_path` | Path to a gffutils `FeatureDB` (`.db` file) built from a GTF/GFF |
| `ip_path` / `input_path` | Path to the IP / INPUT signal file (BigWig, BAM, BED, ...) |
| `genome` | Reference genome (e.g. `"hg19"`, `"mm10"`) |
| `bins` | Number of signal bins across the TSS window |
| `window` | `(upstream, downstream)` bp around each TSS |

**Returns:** a `dict` with keys `"profile"`, `"heatmap_no_sorted"`, `"heatmap_sorted_mean"`, `"heatmap_sorted_max"`.

### 2. Loading RNA-seq Data

The signal-based pipelines above take a signal file directly. Pipelines that integrate RNA-seq (or ATAC-seq) results take a **CSV/TSV path** directly too — no separate loading step needed, as shown in the next section. `Pyntegrate.results_table` also provides a lower-level `DEseq2Results` object model, useful if you want to filter/inspect a differential-expression table yourself before or instead of using a pipeline:

```python
from Pyntegrate.results_table import DEseq2Results

rna = DEseq2Results(
    "/path/to/rna-seq.csv",
    id_column="GeneID",
    pval_column="padj",
    lfc_column="log2FoldChange",
    mean_column="baseMean",
)
up = rna.upregulated(thresh=0.05, idx=False)   # subsetted DEseq2Results
```

The CSV should contain at least a gene-ID column plus `log2FoldChange` and `padj` (or equivalently named columns you pass to `lfc_column`/`pval_column`).

### 3. Creating Visualizations

=== "Signal integrated with RNA-seq status"

    ```python
    from Pyntegrate.pipelines import IntegratedHeatmapPipeline

    pipeline = IntegratedHeatmapPipeline(
        db_path="/path/to/annotations.gtf.db",
        signal_path="/path/to/chip_signal.bw",
        expression_table_path="/path/to/rna-seq.csv",
        genome="hg19",
        gene_id_col="GeneID",
        lfc_col="log2FoldChange",
        divided_by="RNA",
        threshold=1.0,
        signal_name="ChIP-seq",
    )
    fig = pipeline.run()
    fig.savefig("integrated_output.png", dpi=300)
    ```

=== "Volcano Plot"

    ```python
    from Pyntegrate.pipelines import VolcanoPipeline

    pipeline = VolcanoPipeline(
        csv_path="/path/to/rna-seq.csv",
        gene_id_col="GeneID",
        padj_col="padj",
        lfc_col="log2FoldChange",
    )
    fig = pipeline.run()
    fig.savefig("volcano.png", dpi=300)
    ```

Every pipeline follows this same shape: construct it with paths/parameters, call `.run()`, then `.savefig(...)` on whatever it returns (a single `Figure`, or a `dict`/`tuple` of them — see each pipeline's docstring in the [API reference](api/pipelines.md) for its exact return type).

## Data File Formats

### Supported Signal Files

| Format | Notes |
|---|---|
| **BigWig** (`.bw`) | Recommended for large-scale, pre-computed signal data |
| **BAM** (`.bam`) | Aligned reads; signal is computed on the fly |
| **BED** (`.bed`) | Genomic intervals |

Pass the matching type via `*_extension="bigwig"` / `"bam"` / `"bed"` on the pipeline.

### Supported Feature Files

| Format | Notes |
|---|---|
| **GTF / GFF** | Gene annotations — pre-built into a gffutils `.db` file, passed as `db_path` |
| **BED** | Genomic intervals, e.g. for `GeneScaleRegionsProfilePipeline` |

### RNA-seq / ATAC-seq Table Format

A CSV or TSV with at least:

| Column | Description |
|---|---|
| A gene/peak ID column | Passed as `gene_id_col` (or similar) — must match IDs in `db_path` |
| A log2 fold-change column | Passed as `lfc_col` |
| An adjusted p-value column | Passed as `padj_col` (where applicable) |

## Example: Complete Analysis Workflow

```python
from Pyntegrate.pipelines import IntegratedHeatmapPipeline

pipeline = IntegratedHeatmapPipeline(
    db_path="genes.gtf.db",
    signal_path="chipseq.bw",
    expression_table_path="rna_data.csv",
    genome="hg38",
    gene_id_col="GeneID",
    lfc_col="log2FoldChange",
    divided_by="RNA",
    bins=100,
    threshold=1.0,
    signal_name="ChIP-seq",
)

fig = pipeline.run()
fig.savefig("analysis_output.png", dpi=300)
```

## Common Tasks

### Filter by Category (e.g. FACT high/low)

Most pipelines accept `category_table_path` + `category_col` + `categories` to subset genes by a label column in a CSV, rather than requiring you to pre-filter the signal yourself:

```python
pipeline = IntegratedHeatmapPipeline(
    db_path="genes.gtf.db",
    signal_path="chipseq.bw",
    expression_table_path="rna_data.csv",
    genome="hg38",
    gene_id_col="GeneID",
    lfc_col="log2FoldChange",
    divided_by="RNA",
    category_table_path="fact_categories.csv",
    category_col="fact_category",
    categories=["high", "low"],
)
```

### Create Higher-Resolution Heatmaps

Increase the `bins` parameter for finer resolution:

```python
pipeline = PeakHeatmapPipeline(
    db_path="genes.gtf.db",
    ip_path="chipseq.bw",
    input_path="input.bw",
    genome="hg38",
    bins=200,  # higher resolution than the default 100
)
```

## Next Steps

- [:material-api: Explore the API Reference](api/pipelines.md)
- [:material-image-multiple: View Complete Examples](examples.md)
- [:material-github: Check out the GitHub repository](https://github.com/SalvaaGomez/Pyntegrate)
