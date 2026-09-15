"""
High-level pipeline classes for common Pyntegrate analysis workflows.

Each pipeline encapsulates a complete, reproducible analysis. Instantiate
with file paths and parameters, then call `.run()` to obtain
ready-to-display `matplotlib` figures.

## Available pipelines

| Class | Description |
|-------|-------------|
| `PeakHeatmapPipeline` | ChIP-seq IP − INPUT heatmaps |
| `CategoryHeatmapPipeline` | Signal heatmaps split by genomic category |
| `OverlayProfilePipeline` | Dual-axis average-profile overlay |
| `IntegratedHeatmapPipeline` | Signal heatmap with RNA up/down overlay |
| `VolcanoPipeline` | Volcano plot from differential expression results |
| `DecilePipeline` | Signal decile boxplot with expression scatter overlay (BigWig/BAM signal + gene DB) |
| `DecileScatterPipeline` | CSV-only decile analysis with flexible column selection and optional category filter |
| `MultiSignalComparisonPipeline` | Multi-modal log2FC scatter plots |
| `RNAViolinPipeline` | RNA expression violin with significance overlay; supports single or multi-panel horizontal layout |
| `CircPlotPipeline` | Circos multi-signal comparison plot |
| `ChIPCategoryCorrelationPipeline` | Hexbin ChIP correlation filtered by FACT category |
| `RNAKnockdownCorrelationPipeline` | Hexbin RNA-seq correlation for significant KD genes |
| `ATACKnockdownCorrelationPipeline` | Hexbin ATAC-seq correlation for significant KD peaks |
| `GeneScaleRegionsProfilePipeline` | Meta-profile over annotated genes using scaled regions |
| `GeneGroupProfilePipeline` | Average signal profile split by ATAC/RNA log2FC quadrants |

## Quick start

```python
from Pyntegrate import PeakHeatmapPipeline

pipe = PeakHeatmapPipeline(
    db_path="gencode.mm10.db",
    ip_path="IP.bw",
    input_path="INPUT.bw",
    genome="mm10",
)
figures = pipe.run()
figures["heatmap_sorted_mean"].savefig("heatmap.png", dpi=150)
```
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .SeqSignalAnalysis import (
    SignalArrayGenerator,
    SignalArrays,
    PeakAnalysis,
    Plotting,
    FactGeneUtils,
    Filtering,
    ChipRnaUtils,
)


class PeakHeatmapPipeline:
    """
    Generate peak heatmaps from ChIP-seq IP vs INPUT subtraction.

    Loads IP and INPUT signal files, generates signal arrays around
    the TSS of each annotated gene, computes the IP − INPUT subtracted
    signal, and produces four ready-to-display figures.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `ip_path` | `str` | — | Path to the IP signal file (BigWig, BAM, BED, …). |
    | `input_path` | `str` | — | Path to the INPUT (control) signal file. |
    | `genome` | `str` | — | Genome identifier (e.g. `"mm10"`, `"hg19"`). |
    | `ip_extension` | `str` | `"bigwig"` | File type of the IP signal. One of `"bigwig"`, `"bam"`, `"bed"`. |
    | `input_extension` | `str` | `"bigwig"` | File type of the INPUT signal. |
    | `bins` | `int` | `100` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 1000)` | `(upstream, downstream)` bp around each TSS. |

    ## Examples

    ```python
    from Pyntegrate import PeakHeatmapPipeline
    import matplotlib.pyplot as plt

    pipe = PeakHeatmapPipeline(
        db_path="gencode.mm10.db",
        ip_path="IP.bw",
        input_path="INPUT.bw",
        genome="mm10",
    )
    figs = pipe.run()

    # Save individual figures
    figs["heatmap_sorted_mean"].savefig("heatmap_mean.png", dpi=150)
    plt.show()
    ```
    """

    def __init__(
        self,
        db_path: str,
        ip_path: str,
        input_path: str,
        genome: str,
        ip_extension: str = "bigwig",
        input_extension: str = "bigwig",
        bins: int = 100,
        window: Tuple[int, int] = (-1000, 1000),
    ) -> None:
        self.db_path = db_path
        self.ip_path = ip_path
        self.input_path = input_path
        self.genome = genome
        self.ip_extension = ip_extension
        self.input_extension = input_extension
        self.bins = bins
        self.window = window

    def run(self) -> Dict[str, plt.Figure]:
        """
        Execute the pipeline and return all figures.

        ## Returns

        `dict[str, matplotlib.figure.Figure]` — Dictionary with keys:

        | Key | Description |
        |-----|-------------|
        | `"profile"` | Average IP vs INPUT signal profile around TSS. |
        | `"heatmap_no_sorted"` | Unsorted per-gene heatmap of IP − INPUT signal. |
        | `"heatmap_sorted_mean"` | Heatmap sorted by per-gene mean signal (descending). |
        | `"heatmap_sorted_max"` | Heatmap sorted by the position of maximum signal. |
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)

        _, arrays, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path,
            filePath=self.ip_path,
            extensionFile=self.ip_extension,
            inputSignalPath=self.input_path,
            extensionInput=self.input_extension,
            genome=self.genome,
            bins=self.bins,
            x=x,
            normalize_to_cpm=True,
        )

        arrays_ip, arrays_input = SignalArrays.values_array(arrays["ip"], arrays["input"])
        fig_profile = Plotting.distance_from_tss_chipSeq(arrays_ip, arrays_input, xAxes=x)

        peaks = PeakAnalysis.calculate_peaks(arrays_ip, arrays_input)
        fig_hm = Plotting.heatmap_no_sorted(peaks, x)
        fig_hm_mean = Plotting.heatmap_sorted_by_meanValues(peaks, x)
        fig_hm_max = Plotting.heatmap_sorted_by_maxValueIndex(peaks, x)

        return {
            "profile": fig_profile,
            "heatmap_no_sorted": fig_hm,
            "heatmap_sorted_mean": fig_hm_mean,
            "heatmap_sorted_max": fig_hm_max,
        }

class CategoryHeatmapPipeline:
    """
    Generate signal heatmaps split by genomic category (e.g. FACT high/low).

    Loads a signal file, maps annotated genes to category labels from a
    CSV table, and produces three heatmap figures — unsorted,
    sorted by mean, and sorted by max-value index — with each category
    shown as a distinct coloured track.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `signal_path` | `str` | — | Path to the signal file (BigWig, BAM, BED, …). |
    | `category_table_path` | `str` | — | Path to CSV containing gene-to-category mappings. |
    | `genome` | `str` | — | Genome identifier (e.g. `"hg19"`, `"mm10"`). |
    | `category_col` | `str` | — | Column name in the CSV that holds the category labels. |
    | `signal_extension` | `str` | `"bigwig"` | File type of the signal. One of `"bigwig"`, `"bam"`, `"bed"`. |
    | `bins` | `int` | `100` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 2000)` | `(upstream, downstream)` bp around each TSS. |
    | `categories` | `list[str]` | `["high", "low"]` | Category values to include. Genes not in this list are excluded. |
    | `figsize` | `tuple[float, float] \| None` | `None` | Override the auto-calculated figure size `(width, height)` in inches. |

    ## Examples

    ```python
    from Pyntegrate import CategoryHeatmapPipeline

    pipe = CategoryHeatmapPipeline(
        db_path="gencodehg19.db",
        signal_path="ATAC.bw",
        category_table_path="fact_categories.csv",
        genome="hg19",
        category_col="fact_category",
        categories=["high", "low"],
    )
    figs = pipe.run()
    figs["heatmap_sorted_mean"].savefig("category_heatmap.png", dpi=150)
    ```
    """

    def __init__(
        self,
        db_path: str,
        signal_path: str,
        category_table_path: str,
        genome: str,
        category_col: str,
        signal_extension: str = "bigwig",
        bins: int = 100,
        window: Tuple[int, int] = (-1000, 2000),
        categories: Optional[List[str]] = None,
        figsize: Optional[Tuple[float, float]] = None,
    ) -> None:
        self.db_path = db_path
        self.signal_path = signal_path
        self.category_table_path = category_table_path
        self.genome = genome
        self.category_col = category_col
        self.signal_extension = signal_extension
        self.bins = bins
        self.window = window
        self.categories = categories if categories is not None else ["high", "low"]
        self.figsize = figsize

    def run(self) -> Dict[str, plt.Figure]:
        """
        Execute the pipeline and return all figures.

        ## Returns

        `dict[str, matplotlib.figure.Figure]` — Dictionary with keys:

        | Key | Description |
        |-----|-------------|
        | `"heatmap_no_sorted"` | Unsorted per-gene heatmap, colour-coded by category. |
        | `"heatmap_sorted_mean"` | Heatmap sorted by per-gene mean signal (descending). |
        | `"heatmap_sorted_max"` | Heatmap sorted by the position of maximum signal. |
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)

        features, signal_array, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path,
            filePath=self.signal_path,
            extensionFile=self.signal_extension,
            genome=self.genome,
            bins=self.bins,
            x=x,
            normalize_to_cpm=True,
        )

        signal_values = SignalArrays.values_array(signal_array)

        gene_ids = [FactGeneUtils.get_gene_id_from_feature(f) for f in features]
        fact_map = FactGeneUtils.build_fact_map(self.category_table_path, fact_col=self.category_col)
        subset_by = np.array(
            [fact_map.get(gid, "other") if gid is not None else "other" for gid in gene_ids],
            dtype=object,
        )

        if self.categories:
            mask = np.isin(subset_by, self.categories)
            arr = signal_values[mask, :]
            sub = subset_by[mask]
        else:
            arr = signal_values
            sub = subset_by

        n_cat = len(self.categories)
        subset_arg = None if n_cat < 2 else sub
        order_arg = self.categories if n_cat >= 2 else None
        fill = [dict(alpha=0.25) for _ in self.categories] if n_cat >= 2 else dict(alpha=0.25)

        fig_hm = Plotting.heatmap_no_sorted(arr, x, subset_by=subset_arg, subset_order=order_arg, fill_kwargs=fill, figsize=self.figsize)
        fig_mean = Plotting.heatmap_sorted_by_meanValues(arr, x, subset_by=subset_arg, subset_order=order_arg, fill_kwargs=fill, figsize=self.figsize)
        fig_max = Plotting.heatmap_sorted_by_maxValueIndex(arr, x, subset_by=subset_arg, subset_order=order_arg, fill_kwargs=fill, figsize=self.figsize)

        return {
            "heatmap_no_sorted": fig_hm,
            "heatmap_sorted_mean": fig_mean,
            "heatmap_sorted_max": fig_max,
        }

class OverlayProfilePipeline:
    """
    Compare two signals via a dual-axis average-profile overlay.

    Loads two signal files, retains genes present in both, and
    produces a dual-y-axis figure where each y-axis corresponds to
    one signal. Optionally splits profiles by genomic category
    (e.g. FACT high/low), drawing one overlaid profile per category.

    Optionally overlays a single extra DNase/ATAC signal and/or a single
    extra ChIP signal (e.g. an input/control track) on top of the main
    profiles — each plotted once, on its own axis and with its own color
    and legend entry, independent of category filtering.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `signal_1_path` | `str` | — | Path to the first signal file (BigWig, BAM, BED, …). Plotted on the left y-axis. |
    | `signal_2_path` | `str` | — | Path to the second signal file. Plotted on the right y-axis. |
    | `genome` | `str` | — | Genome identifier (e.g. `"hg19"`, `"mm10"`). |
    | `signal_1_extension` | `str` | `"bigwig"` | File type of signal 1. |
    | `signal_2_extension` | `str` | `"bigwig"` | File type of signal 2. |
    | `bins` | `int` | `100` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 1000)` | `(upstream, downstream)` bp around each TSS. |
    | `name` | `str` | `""` | Figure title. |
    | `category_table_path` | `str \| None` | `None` | Path to CSV with gene category labels. Required for category filtering. |
    | `category_col` | `str \| None` | `None` | Column in the category CSV holding the category labels for `signal_1`/`signal_2`. |
    | `category_col_extra` | `str \| None` | `None` | Column in the category CSV holding the category labels for the extra signals. Falls back to `category_col` when `None`. |
    | `categories` | `list[str] \| None` | `None` | Categories to display. `[]` or `None` shows all genes combined. If one category is `"*"`, it means all genes (e.g. `["*", "high"]` compares all genes vs high genes). |
    | `fontsize` | `int` | `12` | Font size for axis labels and title. |
    | `extra_signal_1_path` | `str \| None` | `None` | Optional path to a single extra DNase/ATAC signal file (e.g. input/control), plotted on the left y-axis. |
    | `extra_signal_2_path` | `str \| None` | `None` | Optional path to a single extra ChIP signal file, plotted on the right y-axis. |
    | `extra_signal_1_extension` | `str` | `"bigwig"` | File type of the extra DNase/ATAC signal. |
    | `extra_signal_2_extension` | `str` | `"bigwig"` | File type of the extra ChIP signal. |
    | `label_extra_sig1` | `str` | `"Extra Signal 1"` | Legend label for the extra DNase/ATAC signal. |
    | `label_extra_sig2` | `str` | `"Extra Signal 2"` | Legend label for the extra ChIP signal. |

    ## Examples

    ```python
    from Pyntegrate import OverlayProfilePipeline

    # All genes (no category filtering)
    pipe = OverlayProfilePipeline(
        db_path="gencodehg19.db",
        signal_1_path="chip.bw",
        signal_2_path="atac.bw",
        genome="hg19",
        name="ChIP vs ATAC",
    )
    fig = pipe.run()
    fig.savefig("overlay_profile.png", dpi=150)
    ```

    ```python
    # Split by FACT category
    pipe = OverlayProfilePipeline(
        db_path="gencodehg19.db",
        signal_1_path="chip.bw",
        signal_2_path="atac.bw",
        genome="hg19",
        category_table_path="fact_categories.csv",
        category_col="supt16h FACT category",
        categories=["high", "low"],
        name="ChIP vs ATAC (by FACT)",
    )
    fig = pipe.run()
    ```

    ```python
    # With an extra DNase/ATAC control signal overlaid
    pipe = OverlayProfilePipeline(
        db_path="gencodehg19.db",
        signal_1_path="chip.bw",
        signal_2_path="atac.bw",
        genome="hg19",
        name="ChIP vs ATAC",
        extra_signal_1_path="control_atac.bw",
        label_extra_sig1="Control ATAC",
    )
    fig = pipe.run()
    ```
    """

    def __init__(
        self,
        db_path: str,
        signal_1_path: str,
        signal_2_path: str,
        genome: str,
        signal_1_extension: str = "bigwig",
        signal_2_extension: str = "bigwig",
        bins: int = 100,
        window: Tuple[int, int] = (-1000, 1000),
        name: str = "",
        category_table_path: Optional[str] = None,
        category_col: Optional[str] = None,
        category_col_extra: Optional[str] = None,
        categories: Optional[List[str]] = None,
        fontsize: int = 12,
        extra_signal_1_path: Optional[str] = None,
        extra_signal_2_path: Optional[str] = None,
        extra_signal_1_extension: str = "bigwig",
        extra_signal_2_extension: str = "bigwig",
        label_extra_sig1: str = "Extra Signal 1",
        label_extra_sig2: str = "Extra Signal 2",
    ) -> None:
        self.db_path = db_path
        self.signal_1_path = signal_1_path
        self.signal_2_path = signal_2_path
        self.genome = genome
        self.signal_1_extension = signal_1_extension
        self.signal_2_extension = signal_2_extension
        self.bins = bins
        self.window = window
        self.name = name
        self.category_table_path = category_table_path
        self.category_col = category_col
        self.categories = categories if categories is not None else []
        self.fontsize = fontsize
        self.extra_signal_1_path = extra_signal_1_path
        self.extra_signal_2_path = extra_signal_2_path
        self.extra_signal_1_extension = extra_signal_1_extension
        self.extra_signal_2_extension = extra_signal_2_extension
        self.label_extra_sig1 = label_extra_sig1
        self.label_extra_sig2 = label_extra_sig2
        self.category_col_extra = category_col_extra


    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        When category filtering is enabled, one profile per category is
        drawn on the same axes. When disabled, all genes are combined
        into a single profile.

        ## Returns

        `matplotlib.figure.Figure` — Dual-axis overlay profile figure.
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)
        _, signal1, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path, filePath=self.signal_1_path, extensionFile=self.signal_1_extension,
            genome=self.genome, bins=self.bins, x=x,
        )
        _, signal2, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path, filePath=self.signal_2_path, extensionFile=self.signal_2_extension,
            genome=self.genome, bins=self.bins, x=x,
        )
        extra_sig_1 = None
        if self.extra_signal_1_path:
            _, extra_sig_1, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
                dbPath=self.db_path, filePath=self.extra_signal_1_path, extensionFile=self.extra_signal_1_extension,
                genome=self.genome, bins=self.bins, x=x,
            )
        extra_sig_2 = None
        if self.extra_signal_2_path:
            _, extra_sig_2, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
                dbPath=self.db_path, filePath=self.extra_signal_2_path, extensionFile=self.extra_signal_2_extension,
                genome=self.genome, bins=self.bins, x=x,
            )
        # filter_shared can't accept None, so only pass signals that were actually loaded
        signals_to_filter = [signal1, signal2]
        if extra_sig_1 is not None:
            signals_to_filter.append(extra_sig_1)
        if extra_sig_2 is not None:
            signals_to_filter.append(extra_sig_2)
        filtered = Filtering.filter_shared("id", *signals_to_filter)
        sig1_good, sig2_good = filtered[0], filtered[1]

        next_idx = 2
        extra_sig_1_good = None
        if extra_sig_1 is not None:
            extra_sig_1_good = filtered[next_idx]
            next_idx += 1
        extra_sig_2_good = None
        if extra_sig_2 is not None:
            extra_sig_2_good = filtered[next_idx]

        # Determine which categories to use
        use_categories = (
            self.category_table_path is not None
            and self.category_col is not None
            and (self.categories is not None and len(self.categories) > 0)
        )

        if use_categories:
            # Build FACT map and extract gene IDs
            fact_map = FactGeneUtils.build_fact_map(self.category_table_path, fact_col=self.category_col)
            if self.category_col_extra is not None:
                fact_map_extra = FactGeneUtils.build_fact_map(self.category_table_path, fact_col=self.category_col_extra)
            else:
                fact_map_extra = fact_map
            # Extract gene IDs from both signals
            gene_ids_1 = [d.get("id_norm") or str(d.get("id", "")).split(".")[0] for d in sig1_good]
            gene_ids_2 = [d.get("id_norm") or str(d.get("id", "")).split(".")[0] for d in sig2_good]
            if extra_sig_1_good is not None:
                gene_ids_extra_1 = [d.get("id_norm") or str(d.get("id", "")).split(".")[0] for d in extra_sig_1_good]
            if extra_sig_2_good is not None:
                gene_ids_extra_2 = [d.get("id_norm") or str(d.get("id", "")).split(".")[0] for d in extra_sig_2_good]

            # Normalize IDs if not done yet
            for i, d in enumerate(sig1_good):
                if "id_norm" not in d:
                    d["id_norm"] = str(d.get("id", "")).split(".")[0]
                    gene_ids_1[i] = d["id_norm"]

            for i, d in enumerate(sig2_good):
                if "id_norm" not in d:
                    d["id_norm"] = str(d.get("id", "")).split(".")[0]
                    gene_ids_2[i] = d["id_norm"]
            
            if extra_sig_1_good is not None:
                for i, d in enumerate(extra_sig_1_good):
                    if "id_norm" not in d:
                        d["id_norm"] = str(d.get("id", "")).split(".")[0]
                        gene_ids_extra_1[i] = d["id_norm"]
            
            if extra_sig_2_good is not None:
                for i, d in enumerate(extra_sig_2_good):
                    if "id_norm" not in d:
                        d["id_norm"] = str(d.get("id", "")).split(".")[0]
                        gene_ids_extra_2[i] = d["id_norm"]

            # Map categories
            subset_by_1 = np.array(
                [fact_map.get(gid, "other") for gid in gene_ids_1],
                dtype=object,
            )
            subset_by_2 = np.array(
                [fact_map.get(gid, "other") for gid in gene_ids_2],
                dtype=object,
            )

            if extra_sig_1_good is not None:
                subset_by_extra_1 = np.array(
                    [fact_map_extra.get(gid, "other") for gid in gene_ids_extra_1],
                    dtype=object,
                )
            if extra_sig_2_good is not None:
                subset_by_extra_2 = np.array(
                    [fact_map_extra.get(gid, "other") for gid in gene_ids_extra_2],
                    dtype=object,
                )

            # Creates lists of arrays, one per category
            chip_arrays = []
            dnase_arrays = []
            sig1_extra_arrays = []
            sig2_extra_arrays = []
            category_names = []

            for category in self.categories:
                if str(category).strip() == "*":
                    dnase_arrays.append(SignalArrays.values_array(sig1_good))
                    chip_arrays.append(SignalArrays.values_array(sig2_good))
                    if extra_sig_1_good is not None:
                        sig1_extra_arrays.append(SignalArrays.values_array(extra_sig_1_good))
                    if extra_sig_2_good is not None:
                        sig2_extra_arrays.append(SignalArrays.values_array(extra_sig_2_good))
                    category_names.append("All")
                    continue

                cat_lower = str(category).lower()

                # Filter signal 1
                mask_1 = np.array([s.lower() == cat_lower for s in subset_by_1])
                if mask_1.any():
                    sig1_cat = [d for d, m in zip(sig1_good, mask_1) if m]
                    dnase_arrays.append(SignalArrays.values_array(sig1_cat))
                    # Filter signal 2 - aligned indices
                    mask_2 = np.array([s.lower() == cat_lower for s in subset_by_2])
                    sig2_cat = [d for d, m in zip(sig2_good, mask_2) if m]
                    chip_arrays.append(SignalArrays.values_array(sig2_cat))
                    if extra_sig_1_good is not None:
                        mask_extra_1 = np.array([s.lower() == cat_lower for s in subset_by_extra_1])
                        sig1_extra_cat = [d for d, m in zip(extra_sig_1_good, mask_extra_1) if m]
                        sig1_extra_arrays.append(SignalArrays.values_array(sig1_extra_cat))
                    if extra_sig_2_good is not None:
                        mask_extra_2 = np.array([s.lower() == cat_lower for s in subset_by_extra_2])
                        sig2_extra_cat = [d for d, m in zip(extra_sig_2_good, mask_extra_2) if m]
                        sig2_extra_arrays.append(SignalArrays.values_array(sig2_extra_cat))
                    category_names.append(category)

            # If no arrays were created, fall back to all genes
            if not chip_arrays:
                sig1_vals = [SignalArrays.values_array(sig1_good)]
                sig2_vals = [SignalArrays.values_array(sig2_good)]
                sig1_extra_vals = [SignalArrays.values_array(extra_sig_1_good)] if extra_sig_1_good is not None else None
                sig2_extra_vals = [SignalArrays.values_array(extra_sig_2_good)] if extra_sig_2_good is not None else None
                category_names = ["All"]
            else:
                sig1_vals = dnase_arrays
                sig2_vals = chip_arrays
                sig1_extra_vals = sig1_extra_arrays if extra_sig_1_good is not None else None
                sig2_extra_vals = sig2_extra_arrays if extra_sig_2_good is not None else None
        else:
            # No category filtering: wrap single arrays in lists
            sig1_vals = [SignalArrays.values_array(sig1_good)]
            sig2_vals = [SignalArrays.values_array(sig2_good)]
            sig1_extra_vals = [SignalArrays.values_array(extra_sig_1_good)] if extra_sig_1_good is not None else None
            sig2_extra_vals = [SignalArrays.values_array(extra_sig_2_good)] if extra_sig_2_good is not None else None
            category_names = ["All"]

        return Plotting.twoSignalsOverlay(sig1_vals, 
                                          sig2_vals,
                                          xAxes=x, 
                                          name=self.name, 
                                          category_names=category_names, 
                                          fontsize=self.fontsize, 
                                          extra_dnase_array= sig1_extra_vals,
                                          extra_dnase_label = self.label_extra_sig1,
                                          extra_dnase_color = "#1f77b4",
                                          extra_chip_array = sig2_extra_vals,
                                          extra_chip_label= self.label_extra_sig2,
                                          extra_chip_color = "#ff7f0e")


class GeneGroupProfilePipeline:
    """
    Average signal profile split by ATAC-seq/RNA-seq log2FC quadrants.

    Loads a single signal (typically ChIP-seq) around TSS windows and a CSV
    table (e.g. ``final_table.csv``) with per-gene ATAC-seq and RNA-seq log2
    fold-change columns. Genes are assigned to one of four quadrants based on
    high/low thresholds applied to both columns, and the average signal is
    plotted as one line per quadrant on a single axis.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `signal_path` | `str` | — | Path to the signal file (BigWig, BAM, BED, …). |
    | `genome` | `str` | — | Genome identifier (e.g. `"hg19"`, `"mm10"`). |
    | `category_table_path` | `str` | — | Path to CSV with a gene-ID column plus the ATAC/RNA log2FC columns. |
    | `atac_col` | `str` | — | Column holding the ATAC-seq log2 fold change. |
    | `rna_col` | `str` | — | Column holding the RNA-seq log2 fold change. |
    | `atac_padj_col` | `str` | `None` | Column holding the ATAC-seq adjusted p-values. |
    | `rna_padj_col` | `str` | `None` | Column holding the RNA-seq adjusted p-values. |
    | `padj_threshold` | `float` | `0` | padjusted threshold for classifying genes based on ATAC and RNA. Applied only when both `atac_padj_col` and `rna_padj_col` are set **and** this is non-zero; then a gene joins a quadrant only if both its ATAC and RNA padj are `< padj_threshold`. Must be a float > 0. |
    | `id_col` | `str` | — | Gene-ID column name in `category_table_path` (e.g. `"EnsemblGeneID"`). |
    | `signal_extension` | `str` | `"bigwig"` | File type of the signal. |
    | `bins` | `int` | `100` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 1000)` | `(upstream, downstream)` bp around each TSS. |
    | `threshold` | `float` | `1.0` | log2FC strictly above `threshold` counts as "high"; strictly below `-threshold` counts as "low". |
    | `name` | `str` | `""` | Figure title. |
    | `colors` | `list[str] \| None` | `None` | Line color per quadrant, in `"ATAC high/RNA high", "ATAC high/RNA low", "ATAC low/RNA low", "ATAC low/RNA high"` order. |
    | `fontsize` | `int` | `12` | Base font size for labels/legend. |

    ## Examples

    ```python
    from Pyntegrate.pipelines import GeneGroupProfilePipeline

    pipe = GeneGroupProfilePipeline(
        db_path="gencodehg19.db",
        signal_path="chip.bw",
        genome="hg19",
        category_table_path="final_table.csv",
        atac_col="supt16h ATAC log2 Fold Change",
        rna_col="supt16h RNA log2 fold change",
        id_col="EnsemblGeneID",
        name="ChIP-seq by ATAC/RNA quadrant",
    )
    fig = pipe.run()
    ```
    """

    def __init__(
        self,
        db_path: str,
        signal_path: str,
        genome: str,
        category_table_path: str,
        atac_col: str,
        rna_col: str,
        id_col: str,
        signal_extension: str = "bigwig",
        bins: int = 100,
        window: Tuple[int, int] = (-1000, 1000),
        threshold: float = 1.0,
        name: str = "",
        colors: Optional[List[str]] = None,
        fontsize: int = 12,
        atac_padj_col: Optional[str] = None,
        rna_padj_col: Optional[str] = None,
        padj_threshold: float = 0,
    ) -> None:
        self.db_path = db_path
        self.signal_path = signal_path
        self.genome = genome
        self.category_table_path = category_table_path
        self.atac_col = atac_col
        self.rna_col = rna_col
        self.id_col = id_col
        self.signal_extension = signal_extension
        self.bins = bins
        self.window = window
        self.threshold = threshold
        self.name = name
        self.colors = colors
        self.fontsize = fontsize
        self.atac_padj_col = atac_padj_col
        self.rna_padj_col = rna_padj_col
        self.padj_threshold = padj_threshold

    def _build_quadrant_map(self) -> Dict[str, str]:
        df = pd.read_csv(self.category_table_path, sep=None, engine="python")
        df["ENSG_stripped"] = df[self.id_col].astype(str).apply(FactGeneUtils.strip_version)
        atac = pd.to_numeric(df[self.atac_col], errors="coerce")
        rna = pd.to_numeric(df[self.rna_col], errors="coerce")
        high, low = self.threshold, -self.threshold

        # Optional adjusted p-value filtering. Applied only when both padj columns
        # are provided and `padj_threshold` is non-zero; then a gene joins a
        # quadrant only if both its ATAC and RNA padj are below the threshold.
        # Empty/non-numeric padj values are coerced to NaN, and `NaN < threshold`
        # is False, so such genes never qualify. When off, `significant` is all-True
        # and the quadrant assignment is identical to the lfc-only logic.
        apply_padj = (
            bool(self.atac_padj_col)
            and bool(self.rna_padj_col)
            and self.padj_threshold != 0
        )
        if apply_padj:
            if isinstance(self.padj_threshold, bool) or not isinstance(self.padj_threshold, (int, float)) or self.padj_threshold <= 0:
                raise ValueError(
                    f"padj_threshold must be a float greater than 0, got {self.padj_threshold!r}."
                )
            atac_padj = pd.to_numeric(df[self.atac_padj_col], errors="coerce")
            rna_padj = pd.to_numeric(df[self.rna_padj_col], errors="coerce")
            significant = (atac_padj < self.padj_threshold) & (rna_padj < self.padj_threshold)
        else:
            significant = pd.Series(True, index=df.index)

        labels = pd.Series(np.nan, index=df.index, dtype=object)
        labels[(atac > high) & (rna > high) & significant] = "ATAC high / RNA high"
        labels[(atac > high) & (rna < low) & significant] = "ATAC high / RNA low"
        labels[(atac < low) & (rna < low) & significant] = "ATAC low / RNA low"
        labels[(atac < low) & (rna > high) & significant] = "ATAC low / RNA high"

        return dict(zip(df["ENSG_stripped"], labels))

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Single-axis profile with one line per
        ATAC/RNA quadrant (quadrants with no matching genes are skipped).
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)

        _, signal_array, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path,
            filePath=self.signal_path,
            extensionFile=self.signal_extension,
            genome=self.genome,
            bins=self.bins,
            x=x,
        )

        quadrants = (
            "ATAC high / RNA high",
            "ATAC high / RNA low",
            "ATAC low / RNA low",
            "ATAC low / RNA high",
        )

        quadrant_map = self._build_quadrant_map()
        gene_ids = [
            d.get("id_norm") or FactGeneUtils.strip_version(str(d.get("id", "")))
            for d in signal_array
        ]

        arrays = []
        for quadrant in quadrants:
            members = [d for d, gid in zip(signal_array, gene_ids) if quadrant_map.get(gid) == quadrant]
            arrays.append(SignalArrays.values_array(members) if members else np.empty((0, self.bins)))

        return Plotting.quadrant_signal_profile(
            arrays,
            xAxes=x,
            group_labels=list(quadrants),
            colors=self.colors,
            name=self.name,
            fontsize=self.fontsize,
        )


class IntegratedHeatmapPipeline:
    """
    Signal heatmap integrated with RNA-seq expression status.

    Loads a signal file and a differential expression CSV, synchronises
    gene IDs, and produces an integrated heatmap where genes are
    colour-coded by RNA expression status and padj values (up, down, unchanged). An
    optional category filter (e.g. FACT high/low) can further subset
    the gene list.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `signal_path` | `str` | — | Path to the signal file (BigWig, BAM, BED, …). |
    | `expression_table_path` | `str` | — | Path to CSV containing differential expression results. |
    | `genome` | `str` | — | Genome identifier (e.g. `"hg19"`, `"mm10"`). |
    | `gene_id_col` | `str` | — | Column in the CSV containing gene IDs. |
    | `lfc_col` | `str` | — | Column in the CSV containing log2 fold-change values. |
    | `divided_by` | `str` | — | Label shown on the filtered axis. |
    | `signal_extension` | `str` | `"bigwig"` | File type of the signal. |
    | `category_table_path` | `str \| None` | `None` | Path to CSV with gene category labels. Required for category filtering. |
    | `category_col` | `str \| None` | `None` | Column in the category CSV holding the labels. |
    | `bins` | `int` | `300` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 2000)` | `(upstream, downstream)` bp around each TSS. |
    | `categories` | `list[str]` | `[]` | Category values to display. Empty list includes all genes. |
    | `threshold` | `float` | `1.0` | Absolute log2FC threshold for classifying genes as up/down-regulated. |
    | `padj_col` | `str` | `None` | Column in the CSV containing padjusted values. Padj filtering is applied only when this is set **and** `padj_threshold` is non-zero. |
    | `padj_threshold` | `float` | `0` | padjusted threshold for classifying genes as up/down-regulated. When non-zero (with `padj_col` set), a gene is up/down only if its lfc passes `threshold` **and** its padj `< padj_threshold`; everything else is unchanged. Must be a float > 0. |
    | `signal_name` | `str` | `"Signal"` | Label shown on the signal (heatmap) axis. |
    | `show_heatmap_panels` | `bool \| None` | `None` | Whether to render the per-status heatmap panels. |
    | `counts_output_path` | `str \| None` | `None` | If given, write the up/down/unchanged gene counts to this `.txt` file. |

    ## Examples

    ```python
    from Pyntegrate import IntegratedHeatmapPipeline

    pipe = IntegratedHeatmapPipeline(
        db_path="gencodehg19.db",
        signal_path="chip.bw",
        expression_table_path="deseq2_results.csv",
        genome="hg19",
        gene_id_col="GeneID",
        lfc_col="log2FoldChange",
        divided_by="RNA",
        threshold=1.5,
        signal_name="ChIP-seq",
    )
    fig = pipe.run()
    fig.savefig("integrated_heatmap.png", dpi=150)
    ```
    """

    def __init__(
        self,
        db_path: str,
        signal_path: str,
        expression_table_path: str,
        genome: str,
        gene_id_col: str,
        lfc_col: str,
        divided_by: str,
        signal_extension: str = "bigwig",
        category_table_path: Optional[str] = None,
        category_col: Optional[str] = None,
        bins: int = 300,
        window: Tuple[int, int] = (-1000, 2000),
        categories: Optional[List[str]] = [],
        threshold: float = 1.0,
        signal_name: str = "Signal",
        show_heatmap_panels: bool = None,
        counts_output_path: Optional[str] = None,
        padj_col: Optional[str] = None,
        padj_threshold: float = 0,
    ) -> None:
        self.db_path = db_path
        self.signal_path = signal_path
        self.expression_table_path = expression_table_path
        self.genome = genome
        self.gene_id_col = gene_id_col
        self.lfc_col = lfc_col
        self.signal_extension = signal_extension
        self.category_table_path = category_table_path
        self.category_col = category_col
        self.bins = bins
        self.window = window
        self.categories = categories
        self.threshold = threshold
        self.signal_name = signal_name
        self.divided_by = divided_by
        self.show_heatmap_panels = show_heatmap_panels
        self.counts_output_path = counts_output_path
        self.padj_col = padj_col
        self.padj_threshold = padj_threshold

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Integrated heatmap with RNA
        expression status overlaid on the signal track.
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)

        features, signal_arr, tsses, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path, filePath=self.signal_path, extensionFile=self.signal_extension,
            genome=self.genome, bins=self.bins, x=x,
        )

        # Load expression data as DataFrame
        df = pd.read_csv(self.expression_table_path, sep=None, engine="python")

        # Extract TSS gene IDs
        tsses_gene_ids = []
        for e in tsses:
            gid = None
            if hasattr(e, "attrs"):
                raw = e.attrs.get("gene_id")
                if raw:
                    gid = str(raw).split(".")[0]
            if gid:
                tsses_gene_ids.append(gid)

        # Normalise signal IDs
        for row in signal_arr:
            if row.get("id"):
                row["id"] = str(row["id"]).split(".")[0]
        # Normalise expression IDs
        df["gene_id"] = df[self.gene_id_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
        df = df.drop_duplicates(subset="gene_id", keep="first").copy()

        df[self.gene_id_col] = df["gene_id"]
        df = df.set_index("gene_id", drop=False)
        filtered_df = df.reindex([g for g in tsses_gene_ids if g])
        # Wrap in a simple object to pass to Filtering
        class SimpleExpr:
            def __init__(self, data, lfc_col):
                self.data = data
                self.lfc_column = lfc_col

        expression = SimpleExpr(filtered_df, self.lfc_col)

        signal_good, expression = Filtering.chip_or_atac_genes_not_used_with_rna(
            expression, signal_arr, id="id", rna_column_name=self.gene_id_col, delete_no_peaks=False,
        )
        expression.data = expression.data.dropna(subset=[self.lfc_col])

        # Reindex to signal order, dropping genes without an lfc value from both sides
        valid_ids = set(expression.data.index)
        signal_good = [row for row in signal_good if row["id"] in valid_ids]
        final_ids = [row["id"] for row in signal_good]
        expression.data = expression.data.loc[expression.data.index.intersection(final_ids)].copy()
        expression.data = expression.data.reindex(final_ids)
        expression.data = expression.data.loc[~expression.data.index.isna()].copy()
        final_ids = expression.data.index.tolist()

        signal_final = SignalArrays.values_array(signal_good)

        # Optional category filtering
        subset_plot = None
        order_arg = None
        fill = dict(alpha=0.25)

        if self.category_table_path and self.category_col:
            category_map = FactGeneUtils.build_fact_map(self.category_table_path, fact_col=self.category_col)
            subset_by = np.array([category_map.get(gid, "other") for gid in final_ids], dtype=object)

            if self.categories:
                mask = np.isin(subset_by, self.categories)
                signal_final = signal_final[mask, :]
                subset_plot = subset_by[mask]
                expression.data = expression.data.iloc[mask].copy()
                order_arg = self.categories if len(self.categories) >= 2 else None
                fill = [dict(alpha=0.25) for _ in self.categories] if len(self.categories) >= 2 else dict(alpha=0.25)

        return Plotting.atac_or_chip_with_rna(
            rna=expression, signal_values=signal_final, xAxis=x,
            signalName=self.signal_name,
            divided_by=self.divided_by,
            threshold=self.threshold,
            padj_col=self.padj_col,
            padj_threshold=self.padj_threshold,
            subset_by=subset_plot, subset_order=order_arg,
            line_kwargs=dict(), fill_kwargs=fill,
            show_heatmap_panels=self.show_heatmap_panels,
            counts_output_path=self.counts_output_path,
        )

class VolcanoPipeline:
    """
    Volcano plot from differential expression results with optional category split.

    Loads a CSV with expression data and generates volcano plots using
    log2 fold-change and adjusted p-value thresholds. Supports splitting
    by a categorical column (e.g. FACT high/low) and optional
    intensity-based alpha transparency.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to CSV with differential expression results. |
    | `gene_id_col` | `str` | `"GeneID"` | Column with gene identifiers (used for deduplication). |
    | `gene_name_col` | `str \| None` | `None` | Column with human-readable gene names for labels. |
    | `padj_col` | `str` | `"padj"` | Column with adjusted p-values. |
    | `lfc_col` | `str` | `"log2FoldChange"` | Column with log2 fold-change values. |
    | `fact_col` | `str \| None` | `None` | Column to split data by category. If `None`, one combined plot is produced. |
    | `intensity_col` | `str \| None` | `None` | Column for intensity values driving point alpha transparency. |
    | `log2fc_upper` | `float` | `1.5` | Positive LFC threshold for significance colouring. |
    | `log2fc_lower` | `float` | `-1.5` | Negative LFC threshold for significance colouring. |
    | `nlog10_upper` | `float` | `2.0` | `-log10(padj)` threshold for significance colouring. |
    | `text_nlog10` | `float \| None` | `None` | `-log10(padj)` threshold for gene-label annotation; a gene is labeled only when this **and** `text_log2fc` are both exceeded. Both default to `None`, so auto-labeling is off by default (only `highlight_genes` is labeled) unless you set both. |
    | `text_log2fc` | `float \| None` | `None` | `\|log2FC\|` threshold for gene-label annotation; a gene is labeled only when this **and** `text_nlog10` are both exceeded. Both default to `None`, so auto-labeling is off by default (only `highlight_genes` is labeled) unless you set both. |
    | `highlight_genes` | `list[str] \| None` | `None` | Gene names to highlight regardless of thresholds. |
    | `categories` | `list[str] \| None` | `None` | Categories from `fact_col` to plot. `None` or `[]` → single combined plot. |
    | `fontsize` | `int` | `12` | Font size for axis labels and gene annotations. |
    | `title` | `str` | `""` | Title for the plot(s). |

    ## Examples

    ```python
    from Pyntegrate import VolcanoPipeline

    # Single plot (no category split)
    pipe = VolcanoPipeline(
        csv_path="deseq2_results.csv",
        padj_col="padj",
        lfc_col="log2FoldChange",
    )
    fig = pipe.run()
    fig.savefig("volcano.png", dpi=150)
    ```

    ```python
    # Per-category plots
    pipe = VolcanoPipeline(
        csv_path="deseq2_results.csv",
        padj_col="padj",
        lfc_col="log2FoldChange",
        fact_col="fact_category",
        categories=["high", "low"],
    )
    figs = pipe.run()   # dict: {"high": Figure, "low": Figure}
    figs["high"].savefig("volcano_high.png", dpi=150)
    ```
    """

    def __init__(
        self,
        csv_path: str,
        gene_id_col: str = "GeneID",
        gene_name_col: str | None = None,
        padj_col: str = "padj",
        lfc_col: str = "log2FoldChange",
        fact_col: str | None = None,
        intensity_col: str | None = None,
        log2fc_upper: float = 1.5,
        log2fc_lower: float = -1.5,
        nlog10_upper: float = 2.0,
        text_nlog10: float | None = None,
        text_log2fc: float | None = None,
        highlight_genes: list[str] | None = None,
        categories: list[str] | None = None,
        fontsize: int = 12,
        title: str = "",
    ) -> None:
        self.csv_path = csv_path
        self.gene_id_col = gene_id_col
        self.gene_name_col = gene_name_col
        self.padj_col = padj_col
        self.lfc_col = lfc_col
        self.fact_col = fact_col
        self.intensity_col = intensity_col
        self.log2fc_upper = log2fc_upper
        self.log2fc_lower = log2fc_lower
        self.nlog10_upper = nlog10_upper
        self.text_nlog10 = text_nlog10
        self.text_log2fc = text_log2fc
        self.highlight_genes = highlight_genes
        self.categories = categories if categories is not None else []
        self.fontsize = fontsize
        self.title = title

    def run(self) -> plt.Figure | tuple[plt.Figure, ...] | dict:
        """
        Execute the pipeline and return the figure(s).

        ## Returns

        `matplotlib.figure.Figure | dict[str, matplotlib.figure.Figure]` —
        A single figure when `fact_col` is `None` or `categories` is
        empty; a dictionary keyed by lower-cased category name when
        `fact_col` and `categories` are both specified.
        """
        # Load and prepare data
        df = pd.read_csv(self.csv_path, sep=None, engine="python")
        work = df.copy()

        # Prepare fact column if specified
        if self.fact_col is not None:
            work[self.fact_col] = work[self.fact_col].fillna("other")

        # Prepare p-value column
        work[self.padj_col] = pd.to_numeric(work[self.padj_col], errors="coerce")
        work = work.sort_values(by=self.padj_col, ascending=True, na_position="last")
        work = work.drop_duplicates(subset=[self.gene_id_col], keep="first").reset_index(drop=True)

        return Plotting.volcano_with_intensity_alpha(
            df=work,
            padj_col=self.padj_col,
            lfc_col=self.lfc_col,
            gene_id_col=self.gene_id_col,
            gene_name_col=self.gene_name_col,
            intensity_col=self.intensity_col,
            log2FoldChange_upper=self.log2fc_upper,
            log2FoldChange_lower=self.log2fc_lower,
            nlog10_upper=self.nlog10_upper,
            text_nlog10=self.text_nlog10,
            text_log2FoldChange=self.text_log2fc,
            highlight_genes=self.highlight_genes,
            fontsize=self.fontsize,
            title=self.title,
            fact_col=self.fact_col,
            categories=self.categories,
        )

class DecilePipeline:
    """
    Signal decile boxplot with expression scatter overlay (signal from BigWig/BAM).

    Sorts annotated genes by their mean signal value, divides them into
    `n_groups` quantile bins (deciles), and produces a dual-axis
    figure: a signal boxplot on the left y-axis and an expression
    log2FC strip-plot on the right y-axis.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `db_path` | `str` | — | Path to the gffutils `FeatureDB` (`.db` file). |
    | `signal_path` | `str` | — | Path to the signal file (BigWig, BAM, BED, …). |
    | `expression_table_path` | `str` | — | Path to CSV with differential expression results. |
    | `genome` | `str` | — | Genome identifier (e.g. `"hg19"`, `"mm10"`). |
    | `gene_id_col` | `str` | — | Column in the expression CSV with gene IDs. |
    | `padj_col` | `str` | — | Column with adjusted p-values. |
    | `lfc_col` | `str` | — | Column with log2 fold-change values. |
    | `signal_extension` | `str` | `"bigwig"` | File type of the signal. |
    | `category_table_path` | `str \| None` | `None` | Path to CSV with gene category labels for optional filtering. |
    | `category_col` | `str \| None` | `None` | Column in the category CSV holding labels. |
    | `category_filter` | `str \| None` | `None` | If set, only genes whose category matches this string are included. |
    | `bins` | `int` | `100` | Number of signal bins across the TSS window. |
    | `window` | `tuple[int, int]` | `(-1000, 1000)` | `(upstream, downstream)` bp around each TSS. |
    | `n_groups` | `int` | `10` | Number of quantile groups (deciles). |
    | `signal_ylabel` | `str` | `"Signal (mean)"` | Left y-axis label. |
    | `expression_ylabel` | `str` | `"Expression (log2FC)"` | Right y-axis label. |

    ## Examples

    ```python
    from Pyntegrate import DecilePipeline

    pipe = DecilePipeline(
        db_path="gencodehg19.db",
        signal_path="chip.bw",
        expression_table_path="deseq2.csv",
        genome="hg19",
        gene_id_col="GeneID",
        padj_col="padj",
        lfc_col="log2FoldChange",
        n_groups=10,
        signal_ylabel="ChIP-seq (mean CPM)",
    )
    fig = pipe.run()
    fig.savefig("decile_scatter.png", dpi=150)
    ```
    """

    def __init__(
        self,
        db_path: str,
        signal_path: str,
        expression_table_path: str,
        genome: str,
        gene_id_col: str,
        padj_col: str,
        lfc_col: str,
        signal_extension: str = "bigwig",
        category_table_path: Optional[str] = None,
        category_col: Optional[str] = None,
        category_filter: Optional[str] = None,
        bins: int = 100,
        window: Tuple[int, int] = (-1000, 1000),
        n_groups: int = 10,
        signal_ylabel: str = "Signal (mean)",
        expression_ylabel: str = "Expression (log2FC)",
    ) -> None:
        self.db_path = db_path
        self.signal_path = signal_path
        self.expression_table_path = expression_table_path
        self.genome = genome
        self.gene_id_col = gene_id_col
        self.padj_col = padj_col
        self.lfc_col = lfc_col
        self.signal_extension = signal_extension
        self.category_table_path = category_table_path
        self.category_col = category_col
        self.category_filter = category_filter
        self.bins = bins
        self.window = window
        self.n_groups = n_groups
        self.signal_ylabel = signal_ylabel
        self.expression_ylabel = expression_ylabel

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Dual-axis figure with signal
        decile boxplot (left y) and expression scatter (right y).
        """
        x = np.linspace(self.window[0], self.window[1], self.bins)

        _, signal1, _, _ = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
            dbPath=self.db_path, filePath=self.signal_path, extensionFile=self.signal_extension,
            genome=self.genome, bins=self.bins, x=x,
        )

        df_expr = pd.read_csv(self.expression_table_path, low_memory=False)

        # Normalise expression IDs
        df_expr["gene_id"] = df_expr[self.gene_id_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
        rna_ids = set(df_expr["gene_id"].dropna().astype(str).values)

        # Normalise signal IDs and filter
        for d in signal1:
            d["id_norm"] = str(d["id"]).split(".")[0]

        signal1_good = []
        seen = set()
        for d in signal1:
            gid = d["id_norm"]
            if gid in rna_ids and gid not in seen:
                signal1_good.append(d)
                seen.add(gid)

        # Optional category filter
        if self.category_table_path and self.category_col and self.category_filter:
            cat_map = FactGeneUtils.build_fact_map(self.category_table_path, fact_col=self.category_col)
            signal1_good = [
                d for d in signal1_good
                if cat_map.get(d["id_norm"], "other") == self.category_filter
            ]

        # Sort by mean
        signal1_good = ChipRnaUtils.sort_chip_by_mean_desc(signal1_good)
        signal1_ids = [str(d["id"]).split(".")[0] for d in signal1_good]

        # Deduplicate and reindex RNA
        df_expr_unique = ChipRnaUtils.deduplicate_rna_by_id(
            df_expr, id_col="gene_id", strategy="best_padj", padj_col=self.padj_col,
        )
        df_expr_ordered = (
            df_expr_unique.set_index("gene_id").reindex(signal1_ids).reset_index(drop=False)
        )

        if len(df_expr_ordered) == 0:
            raise ValueError("Empty intersection after ID normalisation.")

        # Decile binning
        signal_vals = SignalArrays.values_array(signal1_good)
        sig_means = np.nanmean(signal_vals, axis=1)
        expr_vals = pd.to_numeric(df_expr_ordered[self.lfc_col], errors="coerce").to_numpy(dtype=float)
        expr_vals = ChipRnaUtils.replace_inf_with_clip(pd.Series(expr_vals)).to_numpy(dtype=float)

        sig_data = pd.DataFrame({"Signal": sig_means})
        try:
            sig_data["Group"] = pd.qcut(sig_data["Signal"], q=self.n_groups,
                                        labels=np.arange(1, self.n_groups + 1)).astype(int)
        except ValueError:
            ranks = sig_data["Signal"].rank(method="average")
            sig_data["Group"] = pd.qcut(ranks, q=self.n_groups,
                                        labels=np.arange(1, self.n_groups + 1)).astype(int)

        scatter_data = pd.DataFrame({
            "Group": sig_data["Group"].values,
            "Expression": expr_vals,
        }).dropna(subset=["Expression"])

        return Plotting.decile_boxplot_rna_scatter(
            sig_data, scatter_data,
            chip_col="Signal", rna_col="Expression", group_col="Group",
            chip_ylabel=self.signal_ylabel, rna_ylabel=self.expression_ylabel,
        )

class DecileScatterPipeline:
    """
    CSV-only decile analysis with flexible column selection and optional category filter.

    Reads a pre-merged CSV, sorts genes by any numeric column, divides
    them into `n_deciles` quantile bins, and produces a dual-axis
    figure: a signal boxplot (left y) and an RNA log2FC scatter (right y).
    No genome database or signal file is required — all data come from the CSV.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to CSV with at least the chip and RNA columns. |
    | `chip_col` | `str` | — | Column used to sort genes into deciles (any numeric column). |
    | `rna_col` | `str` | — | Column with RNA log2 fold-change for the scatter overlay. |
    | `n_deciles` | `int` | `10` | Number of quantile groups. |
    | `chip_range` | `tuple[float, float] \| None` | `None` | Optional fixed `(min, max)` range for the ChIP boxplot y-axis. |
    | `chip_ylabel` | `str` | `"ChIP"` | Left y-axis label. |
    | `rna_ylabel` | `str` | `"RNA (log2 fold change)"` | Right y-axis label. |
    | `title` | `str \| None` | `"Box Plot (ChIP) + Scatter (RNA)"` | Figure title. |
    | `category_col` | `str \| None` | `None` | Column to filter rows by category. If `None`, no filtering is applied. |
    | `categories` | `str` | `""` | Single category string to keep (case-insensitive). Empty string `""` disables filtering. |
    | `label_fontsize` | `int` | `25` | Font size for axis labels. |

    ## Examples

    ```python
    from Pyntegrate import DecileScatterPipeline

    # Sort by LFC, all genes
    pipe = DecileScatterPipeline(
        csv_path="merged.csv",
        chip_col="chip_log2FC",
        rna_col="rna_log2FC",
        chip_ylabel="ChIP log2FC",
        title="ChIP LFC deciles",
    )
    fig = pipe.run()
    fig.savefig("decile_lfc.png", dpi=150)
    ```

    ```python
    # Filter to FACT-high genes only, sort by baseMean
    pipe = DecileScatterPipeline(
        csv_path="merged.csv",
        chip_col="chip_basemean",
        rna_col="rna_log2FC",
        chip_ylabel="ChIP baseMean",
        category_col="fact_category",
        categories="high",
        title="FACT-high deciles",
    )
    fig = pipe.run()
    ```
    """

    def __init__(
        self,
        csv_path: str,
        chip_col: str,
        rna_col: str,
        n_deciles: int = 10,
        chip_range: tuple[float, float] | None = None,
        chip_ylabel: str = "ChIP",
        rna_ylabel: str = "RNA (log2 fold change)",
        title: str| None = "Box Plot (ChIP) + Scatter (RNA)",
        category_col: str | None = None,
        categories: str = "",
        label_fontsize: int = 25,
    ) -> None:
        self.csv_path = csv_path
        self.chip_col = chip_col
        self.rna_col = rna_col
        self.n_deciles = n_deciles
        self.chip_ylabel = chip_ylabel
        self.chip_range = chip_range
        self.rna_ylabel = rna_ylabel
        self.title = title
        self.category_col = category_col
        self.categories = categories
        self.label_fontsize = label_fontsize

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        When `category_col` and `categories` are both set, only the
        matching rows are included before decile binning.

        ## Returns

        `matplotlib.figure.Figure` — Dual-axis decile figure.
        """
        df = pd.read_csv(self.csv_path, sep=None, engine="python")

        # Convert to numeric
        df[self.chip_col] = pd.to_numeric(df[self.chip_col], errors="coerce")
        df[self.rna_col] = pd.to_numeric(df[self.rna_col], errors="coerce")

        # Filter by category if specified and not empty string
        if self.category_col is not None and self.categories:
            df[self.category_col] = df[self.category_col].fillna("other")
            mask = df[self.category_col].str.lower() == self.categories.lower()
            df = df[mask].copy()

        # Remove rows with NaN ChIP values
        df = df.dropna(subset=[self.chip_col]).copy()

        if len(df) == 0:
            raise ValueError("No valid ChIP data after filtering.")

        # Sort by ChIP column descending
        df = df.sort_values(by=self.chip_col, ascending=False).reset_index(drop=True)

        # Create decile groups
        try:
            df["Decile"] = pd.qcut(
                df[self.chip_col],
                q=self.n_deciles,
                labels=np.arange(1, self.n_deciles + 1),
            ).astype(int)
        except ValueError:
            # Fallback if there are many duplicate values
            ranks = df[self.chip_col].rank(method="average")
            df["Decile"] = pd.qcut(
                ranks,
                q=self.n_deciles,
                labels=np.arange(1, self.n_deciles + 1),
            ).astype(int)

        # Prepare data for plotting
        chip_data = df[["Decile", self.chip_col]].rename(
            columns={self.chip_col: "Chip"}
        )
        scatter_data = df[["Decile", self.rna_col]].rename(
            columns={self.rna_col: "RNA"}
        ).dropna(subset=["RNA"])

        return Plotting.decile_boxplot_rna_scatter(
            chip_data,
            scatter_data,
            chip_col="Chip",
            rna_col="RNA",
            group_col="Decile",
            chip_ylabel=self.chip_ylabel,
            chip_range=self.chip_range,
            rna_ylabel=self.rna_ylabel,
            title=self.title,
            label_fontsize=self.label_fontsize
        )

class MultiSignalComparisonPipeline:
    """
    Multi-modal scatter plots comparing signal and RNA log2 fold-changes.

    Loads a pre-merged CSV that contains LFC columns from multiple
    experimental modalities (e.g. ChIP-seq, ATAC-seq, RNA-seq) and
    produces a set of 2-D scatter plots comparing pairs of LFC values,
    coloured by RNA expression status.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `expression_table_path` | `str` | — | Path to CSV containing the merged results with LFC columns from all modalities. |
    | `lfc_col_1` | `str` | — | Column name for the first signal log2 fold-change (e.g. ChIP LFC). |
    | `lfc_col_2` | `str` | — | Column name for the second signal log2 fold-change (e.g. ATAC LFC). |
    | `lfc_col_rna` | `str` | — | Column name for the RNA log2 fold-change (used for colouring). |
    | `title` | `str` | `""` | Figure suptitle. |

    ## Examples

    ```python
    from Pyntegrate import MultiSignalComparisonPipeline

    pipe = MultiSignalComparisonPipeline(
        expression_table_path="merged_results.csv",
        lfc_col_1="chip_log2FC",
        lfc_col_2="atac_log2FC",
        lfc_col_rna="rna_log2FC",
        title="ChIP vs ATAC vs RNA",
    )
    fig = pipe.run()
    fig.savefig("multi_signal.png", dpi=150)
    ```
    """

    def __init__(
        self,
        expression_table_path: str,
        lfc_col_1: str,
        lfc_col_2: str,
        lfc_col_rna: str,
        title: str = "",
    ) -> None:
        self.expression_table_path = expression_table_path
        self.lfc_col_1 = lfc_col_1
        self.lfc_col_2 = lfc_col_2
        self.lfc_col_rna = lfc_col_rna
        self.title = title

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Multi-panel scatter comparison figure.
        """
        from .results_table import ResultsTable

        df = pd.read_csv(self.expression_table_path, sep=None, engine="python")
        data = ResultsTable(df)
        data.lfc_column = self.lfc_col_rna

        return Plotting.all_signal_together_new(
            rna=data,
            chip_lfc_column=self.lfc_col_1,
            atac_lfc_column=self.lfc_col_2,
            title=self.title,
        )

class RNAViolinPipeline:
    """
    RNA expression distribution violin with significance gene overlay.

    This pipeline supports two modes:

    - **Single-panel mode**: provide `rna_lfc_col` and `rna_padj_col`
        as before, and the pipeline produces one violin plot.
    - **Multi-panel mode**: provide `panel_specs`, a list of panel
        dictionaries. The pipeline then renders all violins side by side
        on a single horizontal axes, without subplot borders or panel
        divisions.

    Each panel dictionary must include:

    - ``title``
    - ``rna_lfc_col``
    - ``rna_padj_col``

    Optional keys:

    - ``category_col``
    - ``category_filter``
    - ``padj_threshold``

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to CSV with differential expression results. |
    | `rna_lfc_col` | `str \| None` | `None` | Column with RNA log2 fold-change values (single-panel mode). |
    | `rna_padj_col` | `str \| None` | `None` | Column with adjusted p-values (single-panel mode). |
    | `category_col` | `str \| None` | `None` | Column to filter rows by category (single-panel mode). If `None`, all rows are used. |
    | `category_filter` | `str \| None` | `None` | Category value to keep (exact match, single-panel mode). Only applied when `category_col` is also set. |
    | `padj_threshold` | `float` | `0.1` | Significance threshold. Genes with `padj < threshold` are highlighted. |
    | `ylabel` | `str` | `"RNA log2 fold change"` | Y-axis label. |
    | `title` | `str` | `"RNA Expression Distribution"` | Figure title (single-panel mode). |
    | `panel_specs` | `list[dict] \| None` | `None` | Optional list of panel definitions for multi-violin mode. |
    | `figure_title` | `str \| None` | `None` | Title for the multi-panel figure. Defaults to `title` when omitted. |
    | `figsize_per_panel` | `tuple[float, float]` | `(4, 9.0)` | Approximate size per panel in multi-panel mode. |
    | `zero_pad_fraction` | `float` | `0.10` | Fractional padding applied around zero when computing multi-panel y-limits. |
    | `zero_minimum_pad` | `float` | `1.0` | Minimum absolute padding applied around zero when computing multi-panel y-limits. |
    | `label_fontsize` | `int` | `15` | Font size for axis/panel labels. |

    ## Examples

    ```python
    from Pyntegrate import RNAViolinPipeline

    pipe = RNAViolinPipeline(
        csv_path="deseq2_results.csv",
        rna_lfc_col="log2FoldChange",
        rna_padj_col="padj",
        category_col="fact_category",
        category_filter="high",
        padj_threshold=0.05,
        title="FACT-high genes: RNA expression",
    )
    fig = pipe.run()
    fig.savefig("violin_high.png", dpi=150)
    ```

    ```python
    # Multi-panel mode: one figure, many violins, no subplot divisions
    pipe = RNAViolinPipeline(
        csv_path="final_table.csv",
        panel_specs=[
            {"title": "SSRP1 high", "rna_lfc_col": "ssrp1 RNA log2 fold change", "rna_padj_col": "ssrp1 RNA padj", "category_col": "SSRP1 FACT category", "category_filter": "high"},
            {"title": "SSRP1 low", "rna_lfc_col": "ssrp1 RNA log2 fold change", "rna_padj_col": "ssrp1 RNA padj", "category_col": "SSRP1 FACT category", "category_filter": "low"},
        ],
        figure_title="SSRP1 RNA distribution",
    )
    fig = pipe.run()
    ```
    """

    def __init__(
        self,
        csv_path: str,
        rna_lfc_col: str | None = None,
        rna_padj_col: str | None = None,
        category_col: Optional[str] = None,
        category_filter: Optional[str] = None,
        padj_threshold: float = 0.1,
        ylabel: str = "RNA log2 fold change",
        title: str = "RNA Expression Distribution",
        panel_specs: Optional[List[Dict]] = None,
        figure_title: Optional[str] = None,
        figsize_per_panel: Tuple[float, float] = (4, 9.0),
        zero_pad_fraction: float = 0.10,
        zero_minimum_pad: float = 1.0,
        label_fontsize: int = 15,
    ) -> None:
        self.csv_path = csv_path
        self.rna_lfc_col = rna_lfc_col
        self.rna_padj_col = rna_padj_col
        self.category_col = category_col
        self.category_filter = category_filter
        self.padj_threshold = padj_threshold
        self.ylabel = ylabel
        self.title = title
        self.panel_specs = panel_specs
        self.figure_title = figure_title or title
        self.figsize_per_panel = figsize_per_panel
        self.zero_pad_fraction = zero_pad_fraction
        self.zero_minimum_pad = zero_minimum_pad
        self.label_fontsize = label_fontsize

    @staticmethod
    def _centered_ylim(values: np.ndarray, pad_fraction: float, minimum_pad: float) -> Tuple[float, float]:
        """
        Compute a symmetric y-axis range centred on zero.

        The limit is set to the largest absolute value in `values` plus a
        padding term, so the returned range is always `(-limit, limit)`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `values` | `np.ndarray` | — | Values to derive the range from. Non-finite entries are ignored. |
        | `pad_fraction` | `float` | — | Fractional padding applied to the largest absolute value. |
        | `minimum_pad` | `float` | — | Minimum absolute padding, also used as a fallback scale when all values are zero. |

        ## Returns

        `tuple[float, float]` — `(-limit, limit)`. Returns `(-1.0, 1.0)` when
        `values` has no finite entries.
        """
        values = np.asarray(values, dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            return -1.0, 1.0

        max_abs = float(np.max(np.abs(values)))
        if max_abs == 0:
            max_abs = minimum_pad

        pad = max(max_abs * pad_fraction, minimum_pad)
        limit = max_abs + pad
        return -limit, limit

    @staticmethod
    def _prepare_panel_data(df: pd.DataFrame, spec: Dict) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Filter and clean the source DataFrame for a single multi-panel violin.

        Applies the panel's optional category filter, coerces the LFC and
        padj columns to numeric, drops rows with missing LFC values, and
        derives the significance mask from the padj threshold.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `df` | `pd.DataFrame` | — | Full source DataFrame (shared across all panels). |
        | `spec` | `dict` | — | Panel definition. Required keys: `rna_lfc_col`, `rna_padj_col`. Optional keys: `category_col`, `category_filter`, `padj_threshold` (default `0.1`), `title`. |

        ## Returns

        `tuple[pd.DataFrame, pd.Series]` — The filtered DataFrame and a
        boolean mask (aligned to it) marking rows with `padj < threshold`.

        Raises `ValueError` if no rows remain after filtering.
        """
        work = df.copy()

        category_col = spec.get("category_col")
        category_filter = spec.get("category_filter")
        if category_col and category_filter:
            work = work[work[category_col] == category_filter].copy()

        lfc_col = spec["rna_lfc_col"]
        padj_col = spec["rna_padj_col"]
        threshold = float(spec.get("padj_threshold", 0.1))

        work[lfc_col] = pd.to_numeric(work[lfc_col], errors="coerce")
        work[padj_col] = pd.to_numeric(work[padj_col], errors="coerce")
        work = work.dropna(subset=[lfc_col]).copy()

        if work.empty:
            title = spec.get("title", lfc_col)
            raise ValueError(f"No valid RNA data for panel '{title}'.")

        sig_mask = work[padj_col] < threshold
        return work, sig_mask

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Single violin figure or a single-
        axes multi-violin figure depending on the chosen mode.
        """
        df = pd.read_csv(self.csv_path, sep=None, engine="python")

        if self.panel_specs:
            panel_frames = []
            all_values: list[np.ndarray] = []
            panel_titles: list[str] = []

            for spec in self.panel_specs:
                panel_df, sig_mask = self._prepare_panel_data(df, spec)
                title = spec.get("title", spec["rna_lfc_col"])
                panel_titles.append(title)

                panel_frames.append(
                    pd.DataFrame(
                        {
                            "Panel": title,
                            "RNA": panel_df[spec["rna_lfc_col"]].to_numpy(dtype=float),
                            "Significant": sig_mask.to_numpy(dtype=bool),
                        }
                    )
                )
                all_values.append(panel_df[spec["rna_lfc_col"]].to_numpy(dtype=float))

            combined = pd.concat(panel_frames, ignore_index=True)
            y_limits = self._centered_ylim(
                np.concatenate(all_values),
                pad_fraction=self.zero_pad_fraction,
                minimum_pad=self.zero_minimum_pad,
            )

            return Plotting.rna_violin_multi_panel(
                combined,
                panel_titles,
                y_limits,
                ylabel=self.ylabel,
                figure_title=self.figure_title,
                figsize_per_panel=self.figsize_per_panel,
                label_fontsize=self.label_fontsize,
            )

        if self.rna_lfc_col is None or self.rna_padj_col is None:
            raise ValueError("rna_lfc_col and rna_padj_col are required when panel_specs is not provided.")

        # Single-panel mode
        if self.category_col and self.category_filter:
            df = df[df[self.category_col] == self.category_filter].copy()

        df[self.rna_lfc_col] = pd.to_numeric(df[self.rna_lfc_col], errors="coerce")
        df[self.rna_padj_col] = pd.to_numeric(df[self.rna_padj_col], errors="coerce")
        df = df.dropna(subset=[self.rna_lfc_col]).copy()

        if len(df) == 0:
            raise ValueError("No valid RNA data after filtering.")

        sig_mask = df[self.rna_padj_col] < self.padj_threshold

        return Plotting.rna_violin_with_significance(
            rna_values=df[self.rna_lfc_col].values,
            significant_mask=sig_mask.values,
            title=self.title,
            ylabel=self.ylabel,
            fontsize=self.label_fontsize,
        )

class CircPlotPipeline:
    """
    Circos plot comparing multi-signal changes across two FACT conditions.

    Produces two circos plots — one per condition (e.g. SSRP1 KD and
    SUPT16H KD) — where each arc segment represents one gene and the
    arc colour encodes the combined log2FC from three modalities:
    ChIP-seq, ATAC-seq, and RNA-seq. Optionally exports gene data
    to TSV files for downstream analysis.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to CSV with per-gene log2FC values for both conditions and all three modalities. |
    | `gene_list` | `list[str]` | — | Gene IDs or symbols to include in the circos plot. |
    | `ssrp1_chip_col` | `str` | — | Column with SSRP1 ChIP log2 fold-change. |
    | `ssrp1_rna_col` | `str` | — | Column with SSRP1 RNA log2 fold-change. |
    | `ssrp1_atac_col` | `str` | — | Column with SSRP1 ATAC log2 fold-change. |
    | `supt16h_chip_col` | `str` | — | Column with SUPT16H ChIP log2 fold-change. |
    | `supt16h_rna_col` | `str` | — | Column with SUPT16H RNA log2 fold-change. |
    | `supt16h_atac_col` | `str` | — | Column with SUPT16H ATAC log2 fold-change. |
    | `ssrp1_chip_basemean_col` | `str \| None` | `None` | Optional SSRP1 ChIP baseMean column exported to TSV. |
    | `ssrp1_rna_basemean_col` | `str \| None` | `None` | Optional SSRP1 RNA baseMean column exported to TSV. |
    | `ssrp1_atac_basemean_col` | `str \| None` | `None` | Optional SSRP1 ATAC baseMean column exported to TSV. |
    | `supt16h_chip_basemean_col` | `str \| None` | `None` | Optional SUPT16H ChIP baseMean column exported to TSV. |
    | `supt16h_rna_basemean_col` | `str \| None` | `None` | Optional SUPT16H RNA baseMean column exported to TSV. |
    | `supt16h_atac_basemean_col` | `str \| None` | `None` | Optional SUPT16H ATAC baseMean column exported to TSV. |
    | `identifier_col` | `str` | — | Required. Column used to identify genes: both to match `gene_list` entries against the input table, and to extract gene IDs during threshold-based auto-selection when `gene_list` is empty. |
    | `save_tsv` | `bool` | `False` | If `True`, export gene data to TSV files. |
    | `tsv_output_dir` | `str` | `None` | Directory for TSV output. Falls back to `"."` when `tsv_output_dir` is `None`. |
    | `threshold` | `float` | `1.0` | Threshold for classifying signals as high/low. Compared against the *transformed* value when `value_transform` is given (in `Plotting.circplot`; not in the `gene_list`-empty auto-selection branch below, which always compares raw values). |
    | `allowed_combinations` | `list[tuple[str, str, str]] \| None` | `None` | Allowed (chip, rna, atac) combinations. Each should be `"high"` or `"low"`. If `None`, defaults to `[("high", "low", "low"), ("low", "high", "high")]`. |
    | `value_transform` | `Callable[[float], float] \| None` | `None` | Optional function applied to every ChIP/RNA/ATAC value (both conditions) right after reading it from the CSV, before it is used anywhere in `Plotting.circplot` — see that method for details (e.g. `lambda v: -np.log10(v) if v > 0 else float("nan")` to plot `-log10(padj)`). If `None` (default), values are used as read — identical to the previous behavior. **Not** applied in the `gene_list`-empty auto-selection branch (see `run()`), which always filters on raw column values. |
    | `value_label` | `str` | `"log2 Fold Change"` | Text describing what's plotted, shown in each figure's title as `f"{condition} vs Control - {value_label}"` (same label used for both the SSRP1 and SUPT16H plots). Defaults to `"log2 Fold Change"` — identical to the previous hardcoded title. Set this to match `value_transform` (e.g. `"-log10(padj)"`). |
    | `delimiter` | `str` | `","` | Column delimiter used when reading `csv_path`. |

    ## Examples

    ```python
    from Pyntegrate import CircPlotPipeline

    pipe = CircPlotPipeline(
        csv_path="merged_results.csv",
        gene_list=["ENSG00000139618", "ENSG00000157764"],
        ssrp1_chip_col="ssrp1_chip_lfc",
        ssrp1_rna_col="ssrp1_rna_lfc",
        ssrp1_atac_col="ssrp1_atac_lfc",
        supt16h_chip_col="supt16h_chip_lfc",
        supt16h_rna_col="supt16h_rna_lfc",
        supt16h_atac_col="supt16h_atac_lfc",
        identifier_col="GeneID",
        save_tsv=True,
        threshold=1.5,
        allowed_combinations=[("high", "high", "high"), ("low", "low", "low")],
    )
    results = pipe.run()
    # With save_tsv=False:
    #   results = (fig_ssrp1, fig_supt16h)
    # With save_tsv=True:
    #   results = ((fig_ssrp1, tsv_ssrp1), (fig_supt16h, tsv_supt16h))
    
    fig_ssrp1, tsv_ssrp1 = results[0] if isinstance(results[0], tuple) else (results[0], None)
    fig_supt16h, tsv_supt16h = results[1] if isinstance(results[1], tuple) else (results[1], None)
    
    fig_ssrp1.savefig("circos_ssrp1.png", dpi=150)
    fig_supt16h.savefig("circos_supt16h.png", dpi=150)
    
    if tsv_ssrp1:
        print(f"Data exported to: {tsv_ssrp1}")
    ```
    """

    def __init__(
        self,
        csv_path: str,
        gene_list: List[str],
        ssrp1_chip_col: str,
        ssrp1_rna_col: str,
        ssrp1_atac_col: str,
        supt16h_chip_col: str,
        supt16h_rna_col: str,
        supt16h_atac_col: str,
        identifier_col: str,
        *,
        ssrp1_chip_basemean_col: Optional[str] = None,
        ssrp1_rna_basemean_col: Optional[str] = None,
        ssrp1_atac_basemean_col: Optional[str] = None,
        supt16h_chip_basemean_col: Optional[str] = None,
        supt16h_rna_basemean_col: Optional[str] = None,
        supt16h_atac_basemean_col: Optional[str] = None,
        save_tsv: bool = False,
        tsv_output_dir: Optional[str] = None,
        threshold: float = 1.0,
        allowed_combinations: Optional[List[Tuple[str, str, str]]] = None,
        value_transform: Optional[Callable[[float], float]] = None,
        value_label: str = "log2 Fold Change",
        delimiter: str = ",",
    ) -> None:
        self.csv_path = csv_path
        self.delimiter = delimiter
        self.gene_list = gene_list
        self.ssrp1_chip_col = ssrp1_chip_col
        self.ssrp1_rna_col = ssrp1_rna_col
        self.ssrp1_atac_col = ssrp1_atac_col
        self.ssrp1_chip_basemean_col = ssrp1_chip_basemean_col
        self.ssrp1_rna_basemean_col = ssrp1_rna_basemean_col
        self.ssrp1_atac_basemean_col = ssrp1_atac_basemean_col
        self.supt16h_chip_col = supt16h_chip_col
        self.supt16h_rna_col = supt16h_rna_col
        self.supt16h_atac_col = supt16h_atac_col
        self.supt16h_chip_basemean_col = supt16h_chip_basemean_col
        self.supt16h_rna_basemean_col = supt16h_rna_basemean_col
        self.supt16h_atac_basemean_col = supt16h_atac_basemean_col
        self.identifier_col = identifier_col
        self.save_tsv = save_tsv
        self.tsv_output_dir = tsv_output_dir or "."
        self.threshold = threshold
        self.allowed_combinations = allowed_combinations
        self.value_transform = value_transform
        self.value_label = value_label

    def run(self):
        """
        Execute the pipeline and return both circos figures and optionally TSV paths.

        ## Returns

        If `save_tsv=False`:
        `tuple[matplotlib.figure.Figure, matplotlib.figure.Figure]` —
        `(fig_ssrp1, fig_supt16h)` — one circos plot per condition.

        If `save_tsv=True`:
        `tuple[tuple[matplotlib.figure.Figure, str], tuple[matplotlib.figure.Figure, str]]` —
        `((fig_ssrp1, tsv_path_ssrp1), (fig_supt16h, tsv_path_supt16h))` — figures and their data paths.
        """
        if self.gene_list is None or len(self.gene_list) == 0:
            df = pd.read_csv(self.csv_path, delimiter=self.delimiter)

            genes_ssrp1 = set(
                df.loc[
                    (df[self.ssrp1_chip_col].abs() > self.threshold) &
                    (df[self.ssrp1_atac_col].abs() > self.threshold) &
                    (df[self.ssrp1_rna_col].abs() > self.threshold),
                    self.identifier_col   
                ]
            )
            genes_supt16h = set(
                df.loc[
                    (df[self.supt16h_chip_col].abs() > self.threshold) &
                    (df[self.supt16h_atac_col].abs() > self.threshold) &
                    (df[self.supt16h_rna_col].abs() > self.threshold),
                    self.identifier_col   
                ]
            )
            gene_list_input = list(genes_ssrp1.union(genes_supt16h))
        elif len(self.gene_list) == 1:
            raise ValueError("Circos plot requires at least 2 genes. Please provide a gene list with 2 or more entries.")
        else:
            gene_list_input = self.gene_list

        # Create TSV filenames if saving
        if self.save_tsv:
            tsv_ssrp1_path = os.path.join(self.tsv_output_dir, "circplot_ssrp1_data.tsv")
            tsv_supt16h_path = os.path.join(self.tsv_output_dir, "circplot_supt16h_data.tsv")
        else:
            tsv_ssrp1_path = None
            tsv_supt16h_path = None
        

        result_ssrp1 = Plotting.circplot(
            self.csv_path,
            gene_list_input,
            self.ssrp1_chip_col,
            self.ssrp1_rna_col,
            self.ssrp1_atac_col,
            self.ssrp1_chip_basemean_col,
            self.ssrp1_rna_basemean_col,
            self.ssrp1_atac_basemean_col,
            identifier_col=self.identifier_col,
            save_tsv=self.save_tsv,
            tsv_output_path=tsv_ssrp1_path,
            threshold=self.threshold,
            allowed_combinations=self.allowed_combinations,
            value_transform=self.value_transform,
            value_label=self.value_label,
        )
        result_supt16h = Plotting.circplot(
            self.csv_path,
            gene_list_input,
            self.supt16h_chip_col,
            self.supt16h_rna_col,
            self.supt16h_atac_col,
            self.supt16h_chip_basemean_col,
            self.supt16h_rna_basemean_col,
            self.supt16h_atac_basemean_col,
            identifier_col=self.identifier_col,
            save_tsv=self.save_tsv,
            tsv_output_path=tsv_supt16h_path,
            threshold=self.threshold,
            allowed_combinations=self.allowed_combinations,
            value_transform=self.value_transform,
            value_label=self.value_label,
        )

        return result_ssrp1, result_supt16h

class ChIPCategoryCorrelationPipeline:
    """
    Hexbin density scatter of ChIP/input log2 ratios filtered by FACT category.

    Reads a table with ChIP log2 fold-change columns for two FACT components
    (e.g. SSRP1 and SUPT16H), keeps rows where at least one component is
    classified as a given category (default ``"high"``), and returns a hexbin
    density scatter plot annotated with Spearman R and sample size.

    Corresponds to Figure 1B: correlation of SSRP1 and SUPT16H ChIP vs input
    log2 ratios at TSS windows classified as high for at least one FACT component.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to the input table (TSV or CSV). |
    | `xcol` | `str` | — | Column with x-axis log2FC values (e.g. SSRP1 ChIP log2FC). |
    | `ycol` | `str` | — | Column with y-axis log2FC values (e.g. SUPT16H ChIP log2FC). |
    | `category_col_1` | `str` | — | First category column (e.g. SSRP1 FACT category). |
    | `category_col_2` | `str` | — | Second category column (e.g. SUPT16H FACT category). |
    | `category_value` | `str` | `"high"` | Value to filter by. Rows where either column equals this value are retained. |
    | `bins` | `int` | `60` | Hexbin gridsize. |
    | `xlabel` | `str \| None` | `None` | X-axis label. Defaults to `xcol`. |
    | `ylabel` | `str \| None` | `None` | Y-axis label. Defaults to `ycol`. |
    | `title` | `str` | `"ChIP correlation"` | Plot title. |
    | `hex_low_color` | `str` | `"#D9D9D9"` | Hexbin colour for low counts (≈ grey85). |
    | `hex_high_color` | `str` | `"black"` | Hexbin colour for high counts. |
    | `delimiter` | `str \| None` | `None` | Column delimiter of the input file. `None` lets pandas infer the separator. |
    | `figsize` | `tuple[float, float]` | `(5, 4)` | Figure size in inches. |
    | `fontsize` | `int` | `12` | Base font size for labels and annotations. |

    ## Examples

    ```python
    from Pyntegrate import ChIPCategoryCorrelationPipeline

    pipe = ChIPCategoryCorrelationPipeline(
        csv_path="final_table.csv",
        xcol="ssrp1_chip_log2_fold_change",
        ycol="supt16h_chip_log2_fold_change",
        category_col_1="ssrp1 fact category",
        category_col_2="supt16h fact category",
        category_value="high",
        bins=60,
        xlabel=r"$\\log_2$ [SSRP1 ChIP/Input]",
        ylabel=r"$\\log_2$ [SUPT16H ChIP/Input]",
        title="SUPT16H / SSRP1\\nChIP correlation",
        delimiter=",",
    )
    fig = pipe.run()
    fig.savefig("fig1B.svg")
    ```
    """

    def __init__(
        self,
        csv_path: str,
        xcol: str,
        ycol: str,
        category_col_1: str,
        category_col_2: str,
        category_value: str = "high",
        bins: int = 60,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: str = "ChIP correlation",
        hex_low_color: str = "#D9D9D9",
        hex_high_color: str = "black",
        delimiter: Optional[str] = None,
        figsize: Tuple[float, float] = (5, 4),
        fontsize: int = 12,
    ) -> None:
        self.csv_path = csv_path
        self.xcol = xcol
        self.ycol = ycol
        self.category_col_1 = category_col_1
        self.category_col_2 = category_col_2
        self.category_value = category_value
        self.bins = bins
        self.xlabel = xlabel if xlabel is not None else xcol
        self.ylabel = ylabel if ylabel is not None else ycol
        self.title = title
        self.hex_low_color = hex_low_color
        self.hex_high_color = hex_high_color
        self.delimiter = delimiter
        self.figsize = figsize
        self.fontsize = fontsize

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Hexbin density scatter with Spearman R annotation.
        """
        from .plotutils import hexbin_correlation_scatter

        if self.delimiter is None:
            # Let pandas infer the separator (supports CSV or TSV)
            df = pd.read_csv(self.csv_path, sep=None, engine="python")
        else:
            df = pd.read_csv(self.csv_path, sep=self.delimiter, low_memory=False)

        cat1, cat2 = self.category_col_1, self.category_col_2
        val = self.category_value

        mask = (
            (df[cat1].notna() & (df[cat1].astype(str).str.lower() == val.lower())) |
            (df[cat2].notna() & (df[cat2].astype(str).str.lower() == val.lower()))
        )
        df = df[mask].copy()

        df[self.xcol] = pd.to_numeric(df[self.xcol], errors="coerce")
        df[self.ycol] = pd.to_numeric(df[self.ycol], errors="coerce")
        df = df.dropna(subset=[self.xcol, self.ycol])

        return hexbin_correlation_scatter(
            df[self.xcol].to_numpy(dtype=float),
            df[self.ycol].to_numpy(dtype=float),
            bins=self.bins,
            xlabel=self.xlabel,
            ylabel=self.ylabel,
            title=self.title,
            hex_low_color=self.hex_low_color,
            hex_high_color=self.hex_high_color,
            figsize=self.figsize,
            fontsize=self.fontsize,
            annotation_fontsize_factor=0.8,
        )

class RNAKnockdownCorrelationPipeline:
    """
    Hexbin density scatter of RNA-seq knockdown log2 fold-changes for significant genes.

    Reads a table with RNA-seq log2FC and adjusted p-value columns for two knockdown
    conditions (e.g. SSRP1 KD and SUPT16H KD), retains genes where both padj values
    are non-missing and at least one satisfies the significance threshold, and returns
    a hexbin density scatter plot annotated with Spearman R and sample size.

    Corresponds to Figure 1C: correlation of SSRP1 and SUPT16H knockdown RNA-seq
    log2 fold-changes for genes significantly changed in at least one condition.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to the input table (TSV or CSV). |
    | `xcol` | `str` | — | Column with x-axis RNA log2FC (e.g. SSRP1 RNA log2FC). |
    | `ycol` | `str` | — | Column with y-axis RNA log2FC (e.g. SUPT16H RNA log2FC). |
    | `padj_col_1` | `str` | — | Adjusted p-value column for condition 1 (e.g. SSRP1 RNA padj). |
    | `padj_col_2` | `str` | — | Adjusted p-value column for condition 2 (e.g. SUPT16H RNA padj). |
    | `padj_threshold` | `float` | `0.1` | Significance threshold; rows where either padj < threshold are retained. |
    | `bins` | `int` | `50` | Hexbin gridsize. |
    | `xlabel` | `str \| None` | `None` | X-axis label. Defaults to `xcol`. |
    | `ylabel` | `str \| None` | `None` | Y-axis label. Defaults to `ycol`. |
    | `title` | `str` | `"RNA-seq upon FACT depletion"` | Plot title. |
    | `hex_low_color` | `str` | `"#B3B3B3"` | Hexbin colour for low counts (≈ grey70). |
    | `hex_high_color` | `str` | `"black"` | Hexbin colour for high counts. |
    | `delimiter` | `str \| None` | `None` | Column delimiter of the input file. `None` lets pandas infer the separator. |
    | `figsize` | `tuple[float, float]` | `(5, 4)` | Figure size in inches. |
    | `fontsize` | `int` | `12` | Base font size for labels and annotations. |

    ## Examples

    ```python
    from Pyntegrate import RNAKnockdownCorrelationPipeline

    pipe = RNAKnockdownCorrelationPipeline(
        csv_path="final_table.csv",
        xcol="ssrp1_rna_log2_fold_change",
        ycol="supt16h_rna_log2_fold_change",
        padj_col_1="ssrp1 rna padj",
        padj_col_2="supt16h rna padj",
        padj_threshold=0.1,
        bins=50,
        xlabel=r"$\\log_2$ [SSRP1 RNA]",
        ylabel=r"$\\log_2$ [SUPT16H RNA]",
        title="RNA-seq upon FACT depletion",
        delimiter=",",
    )
    fig = pipe.run()
    fig.savefig("fig1C.svg")
    ```
    """

    def __init__(
        self,
        csv_path: str,
        xcol: str,
        ycol: str,
        padj_col_1: str,
        padj_col_2: str,
        padj_threshold: float = 0.1,
        bins: int = 50,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: str = "RNA-seq upon FACT depletion",
        hex_low_color: str = "#B3B3B3",
        hex_high_color: str = "black",
        delimiter: Optional[str] = None,
        figsize: Tuple[float, float] = (5, 4),
        fontsize: int = 12,
    ) -> None:
        self.csv_path = csv_path
        self.xcol = xcol
        self.ycol = ycol
        self.padj_col_1 = padj_col_1
        self.padj_col_2 = padj_col_2
        self.padj_threshold = padj_threshold
        self.bins = bins
        self.xlabel = xlabel if xlabel is not None else xcol
        self.ylabel = ylabel if ylabel is not None else ycol
        self.title = title
        self.hex_low_color = hex_low_color
        self.hex_high_color = hex_high_color
        self.delimiter = delimiter
        self.figsize = figsize
        self.fontsize = fontsize

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Hexbin density scatter with Spearman R annotation.
        """
        from .plotutils import hexbin_correlation_scatter

        if self.delimiter is None:
            df = pd.read_csv(self.csv_path, sep=None, engine="python")
        else:
            df = pd.read_csv(self.csv_path, sep=self.delimiter, low_memory=False)

        df[self.padj_col_1] = pd.to_numeric(df[self.padj_col_1], errors="coerce")
        df[self.padj_col_2] = pd.to_numeric(df[self.padj_col_2], errors="coerce")

        mask = (
            df[self.padj_col_1].notna() &
            df[self.padj_col_2].notna() &
            (
                (df[self.padj_col_1] < self.padj_threshold) |
                (df[self.padj_col_2] < self.padj_threshold)
            )
        )
        df = df[mask].copy()

        df[self.xcol] = pd.to_numeric(df[self.xcol], errors="coerce")
        df[self.ycol] = pd.to_numeric(df[self.ycol], errors="coerce")
        df = df.dropna(subset=[self.xcol, self.ycol])

        return hexbin_correlation_scatter(
            df[self.xcol].to_numpy(dtype=float),
            df[self.ycol].to_numpy(dtype=float),
            bins=self.bins,
            xlabel=self.xlabel,
            ylabel=self.ylabel,
            title=self.title,
            hex_low_color=self.hex_low_color,
            hex_high_color=self.hex_high_color,
            figsize=self.figsize,
            fontsize=self.fontsize,
            annotation_fontsize_factor=0.9,
        )

class ATACKnockdownCorrelationPipeline:
    """
    Hexbin density scatter of ATAC-seq knockdown log2 fold-changes for significant peaks.

    Reads a table with ATAC-seq log2FC and adjusted p-value columns for two knockdown
    conditions (e.g. SUPT16H KD and SSRP1 KD), removes rows with missing LFC or padj
    values in any of the four columns, keeps peaks where at least one padj satisfies the
    significance threshold, and returns a hexbin density scatter plot annotated with
    Spearman R and sample size.

    Corresponds to Figure 1D: hexbin density scatter of SUPT16H vs SSRP1 knockdown
    ATAC-seq log2 fold-changes for significantly differential peaks.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `csv_path` | `str` | — | Path to the input table (TSV or CSV). |
    | `xcol` | `str` | — | Column with x-axis ATAC log2FC (e.g. SUPT16H ATAC log2FC). |
    | `ycol` | `str` | — | Column with y-axis ATAC log2FC (e.g. SSRP1 ATAC log2FC). |
    | `padj_col_1` | `str` | — | Adjusted p-value column for condition 1 (e.g. SSRP1 ATAC padj). |
    | `padj_col_2` | `str` | — | Adjusted p-value column for condition 2 (e.g. SUPT16H ATAC padj). |
    | `padj_threshold` | `float` | `0.05` | Significance threshold; rows where either padj < threshold are retained. |
    | `bins` | `int` | `35` | Hexbin gridsize. |
    | `xlabel` | `str \| None` | `None` | X-axis label. Defaults to `xcol`. |
    | `ylabel` | `str \| None` | `None` | Y-axis label. Defaults to `ycol`. |
    | `title` | `str` | `"ATAC-seq upon FACT depletion"` | Plot title. |
    | `hex_low_color` | `str` | `"#B3B3B3"` | Hexbin colour for low counts (≈ grey70). |
    | `hex_high_color` | `str` | `"black"` | Hexbin colour for high counts. |
    | `delimiter` | `str \| None` | `None` | Column delimiter of the input file. `None` lets pandas infer the separator. |
    | `figsize` | `tuple[float, float]` | `(5, 4)` | Figure size in inches. |
    | `fontsize` | `int` | `12` | Base font size for labels and annotations. |

    ## Examples

    ```python
    from Pyntegrate import ATACKnockdownCorrelationPipeline

    pipe = ATACKnockdownCorrelationPipeline(
        csv_path="merged_deseq2_results.tsv",
        xcol="supt16h_atac_log2_fold_change",
        ycol="ssrp1_atac_log2_fold_change",
        padj_col_1="ssrp1 atac padj",
        padj_col_2="supt16h atac padj",
        padj_threshold=0.05,
        bins=35,
        xlabel=r"$\\log_2$ [SUPT16H ATAC]",
        ylabel=r"$\\log_2$ [SSRP1 ATAC]",
        title="ATAC-seq upon FACT depletion",
    )
    fig = pipe.run()
    fig.savefig("fig1D.svg")
    ```
    """

    def __init__(
        self,
        csv_path: str,
        xcol: str,
        ycol: str,
        padj_col_1: str,
        padj_col_2: str,
        padj_threshold: float = 0.05,
        bins: int = 35,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: str = "ATAC-seq upon FACT depletion",
        hex_low_color: str = "#B3B3B3",
        hex_high_color: str = "black",
        delimiter: Optional[str] = None,
        figsize: Tuple[float, float] = (5, 4),
        fontsize: int = 12,
    ) -> None:
        self.csv_path = csv_path
        self.xcol = xcol
        self.ycol = ycol
        self.padj_col_1 = padj_col_1
        self.padj_col_2 = padj_col_2
        self.padj_threshold = padj_threshold
        self.bins = bins
        self.xlabel = xlabel if xlabel is not None else xcol
        self.ylabel = ylabel if ylabel is not None else ycol
        self.title = title
        self.hex_low_color = hex_low_color
        self.hex_high_color = hex_high_color
        self.delimiter = delimiter
        self.figsize = figsize
        self.fontsize = fontsize

    def run(self) -> plt.Figure:
        """
        Execute the pipeline and return the figure.

        ## Returns

        `matplotlib.figure.Figure` — Hexbin density scatter with Spearman R annotation.
        """
        from .plotutils import hexbin_correlation_scatter

        if self.delimiter is None:
            df = pd.read_csv(self.csv_path, sep=None, engine="python")
        else:
            df = pd.read_csv(self.csv_path, sep=self.delimiter, low_memory=False)

        df[self.padj_col_1] = pd.to_numeric(df[self.padj_col_1], errors="coerce")
        df[self.padj_col_2] = pd.to_numeric(df[self.padj_col_2], errors="coerce")
        df[self.xcol] = pd.to_numeric(df[self.xcol], errors="coerce")
        df[self.ycol] = pd.to_numeric(df[self.ycol], errors="coerce")

        df = df.dropna(subset=[self.xcol, self.ycol, self.padj_col_1, self.padj_col_2])

        mask = (
            (df[self.padj_col_1] < self.padj_threshold) |
            (df[self.padj_col_2] < self.padj_threshold)
        )
        df = df[mask].copy()

        return hexbin_correlation_scatter(
            df[self.xcol].to_numpy(dtype=float),
            df[self.ycol].to_numpy(dtype=float),
            bins=self.bins,
            xlabel=self.xlabel,
            ylabel=self.ylabel,
            title=self.title,
            hex_low_color=self.hex_low_color,
            hex_high_color=self.hex_high_color,
            figsize=self.figsize,
            fontsize=self.fontsize,
            annotation_fontsize_factor=0.8,
        )

class GeneScaleRegionsProfilePipeline:
    """
    Meta-profile over annotated genes using a scale-regions layout.

    Delegates to deepTools ``computeMatrix scale-regions`` +
    ``plotProfile``.  The BED or BED.GZ file is decompressed to a
    temporary directory if needed (``computeMatrix`` does not accept
    ``.bed.gz``), and the matrix is written to a temporary file that is
    cleaned up after the figure is saved.

    The expected input BED has one gene per row and includes strand in
    the sixth column, which matches the repository's ``data/genes.bed.gz``.

    ``deepTools`` must be installed and available on ``$PATH`` (e.g. via
    ``conda activate pyntegrate``).

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `genes_bed_path` | `str` | — | Path to a BED or BED.GZ file with annotated genes. |
    | `signal_paths` | `list[str]` | — | BigWig files to summarise, one line per sample. |
    | `sample_labels` | `list[str]` | — | Legend labels for each BigWig track. |
    | `colors` | `list[str]` | — | Line colours for each sample. |
    | `output_path` | `str \| None` | `None` | SVG/PNG file to save. If `None`, nothing is saved. |
    | `regions_label` | `str \| None` | `None` | Passed as ``--regionsLabel`` to ``plotProfile``. |
    | `before_region_start_length` | `int` | `500` | Upstream flank length in bp (``--beforeRegionStartLength``). |
    | `region_body_length` | `int` | `1000` | Scaled gene-body length in bp (``--regionBodyLength``). |
    | `after_region_start_length` | `int` | `500` | Downstream flank length in bp (``--afterRegionStartLength``). |
    | `bin_size` | `int` | `10` | Bin size in bp (``--binSize``). |
    | `chunk_size` | `int` | `500` | Accepted for compatibility; not used by deepTools. Clamped to a minimum of `1`. |
    | `skip_zeros` | `bool` | `True` | Pass ``--skipZeros`` to ``computeMatrix``. |
    | `processes` | `int \| None` | `None` | ``--numberOfProcessors`` for ``computeMatrix``. |
    | `figsize` | `tuple[float, float]` | `(4.7244, 3.1496)` | Figure size in inches; converted to cm for ``plotProfile``. |
    | `dpi` | `int` | `300` | ``--dpi`` for ``plotProfile``. |
    | `ylabel` | `str` | `"ChIP-seq Signal"` | ``--yAxisLabel`` for ``plotProfile``. |
    | `title` | `str` | `""` | ``--plotTitle`` for ``plotProfile``. |
    | `legend_location` | `str` | `"upper right"` | ``--legendLocation`` for ``plotProfile``. Hyphens from the input are normalised to spaces. |
    | `show` | `bool` | `True` | Open the saved figure with ``xdg-open`` after saving. |
    """

    def __init__(
        self,
        genes_bed_path: str,
        signal_paths: List[str],
        sample_labels: List[str],
        colors: List[str],
        output_path: Optional[str] = None,
        regions_label: Optional[str] = None,
        before_region_start_length: int = 500,
        region_body_length: int = 1000,
        after_region_start_length: int = 500,
        bin_size: int = 10,
        chunk_size: int = 500,
        skip_zeros: bool = True,
        processes: Optional[int] = None,
        figsize: Tuple[float, float] = (4.7244, 3.1496),
        dpi: int = 300,
        ylabel: str = "ChIP-seq Signal",
        title: str = "",
        legend_location: str = "upper right",
        show: bool = True,
    ) -> None:
        if len(signal_paths) != len(sample_labels):
            raise ValueError("signal_paths and sample_labels must have the same length")
        if len(signal_paths) != len(colors):
            raise ValueError("signal_paths and colors must have the same length")
        if before_region_start_length % bin_size != 0:
            raise ValueError("before_region_start_length must be divisible by bin_size")
        if after_region_start_length % bin_size != 0:
            raise ValueError("after_region_start_length must be divisible by bin_size")
        if region_body_length % bin_size != 0:
            raise ValueError("region_body_length must be divisible by bin_size")

        self.genes_bed_path = genes_bed_path
        self.signal_paths = signal_paths
        self.sample_labels = sample_labels
        self.colors = colors
        self.output_path = output_path
        self.regions_label = regions_label
        self.before_region_start_length = before_region_start_length
        self.region_body_length = region_body_length
        self.after_region_start_length = after_region_start_length
        self.bin_size = bin_size
        self.chunk_size = max(1, chunk_size)
        self.skip_zeros = skip_zeros
        self.processes = processes
        self.figsize = figsize
        self.dpi = dpi
        self.ylabel = ylabel
        self.title = title
        self.legend_location = legend_location.replace("-", " ")
        self.show = show


    def run(self) -> Optional[str]:
        """
        Execute the pipeline via deepTools ``computeMatrix scale-regions`` +
        ``plotProfile`` and save the result to ``output_path``.

        ## Returns

        `str | None` — Absolute path of the saved figure, or ``None`` when
        ``output_path`` was not provided.
        """
        return Plotting.generate_gene_scale_regions_profile(
            genes_bed_path=self.genes_bed_path,
            signal_paths=self.signal_paths,
            sample_labels=self.sample_labels,
            colors=self.colors,
            output_path=self.output_path,
            regions_label=self.regions_label,
            before_region_start_length=self.before_region_start_length,
            region_body_length=self.region_body_length,
            after_region_start_length=self.after_region_start_length,
            bin_size=self.bin_size,
            skip_zeros=self.skip_zeros,
            processes=self.processes,
            figsize=self.figsize,
            dpi=self.dpi,
            ylabel=self.ylabel,
            title=self.title,
            legend_location=self.legend_location,
            show=self.show,
        )
