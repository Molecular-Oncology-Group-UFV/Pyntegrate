"""
High-level analysis classes for genomic signal extraction, visualization, and multi-omics integration.

Provides utilities for TSS interval generation, windowed signal array extraction, peak enrichment
calculation, and a comprehensive set of plotting methods combining ChIP-seq, ATAC-seq, and RNA-seq data.

## Classes

| Class | Description |
|-------|-------------|
| `TSSUtils` | Generate and transform TSS intervals from gffutils databases |
| `SignalArrayGenerator` | Extract windowed signal arrays around TSS from BigWig/BAM files |
| `SignalArrays` | Convert dict-type signal arrays to NumPy matrices |
| `PeakAnalysis` | Compute IP vs INPUT peak enrichment scores |
| `Plotting` | High-level heatmap, profile, and multi-omics visualization |
| `Homer` | HOMER CLI wrappers for peak calling and annotation |
| `Filtering` | Intersect and align gene sets across modalities |
| `BigWigUtils` | Convert BigWig files to BED-like text format |
| `AnnotationAnalysis` | Integrate peak data with HOMER and FACT annotations |
| `FactGeneUtils` | Gene ID normalization and FACT category mapping |
| `ChipRnaUtils` | Match, deduplicate, and reshape chip signal / RNA-seq tables for decile pipelines |
"""
from __future__ import annotations


from ._genomic_signal import *
from .array_helpers import * 
from .plotutils import * 
from .results_table import * 


import os
import re
import sys
import copy
import subprocess
import multiprocessing
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union


import numpy as np
import pandas as pd
import gffutils
import pybedtools
from pybedtools.featurefuncs import TSS
from gffutils.helpers import asinterval
from matplotlib import pyplot as plt

try:
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
except (ImportError, AttributeError):
    # Matplotlib version compatibility issue - 3D plotting not critical for core functionality
    pass
import inspect

__all__ = [
    name for name, obj in inspect.getmembers(sys.modules[__name__])
    if inspect.isclass(obj) and obj.__module__ == __name__
]
class TSSUtils:
   """
   Utilities for generating and transforming TSS (Transcription Start Site) intervals.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `tss_generator()` | Yields one TSS interval for each annotated transcript |
   | `_fix_chr_to_match_bigwig()` | Normalizes chromosome naming to match BigWig style |
   """


   @staticmethod
   def tss_generator(db: gffutils.FeatureDB) -> Iterable[pybedtools.Interval]:
       """
       Generator that yields one TSS interval for each annotated transcript.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `db` | `gffutils.FeatureDB` | — | A gffutils database containing transcript features. |

       ## Yields

       `pybedtools.Interval` — TSS interval (1 bp upstream, 0 bp downstream) for each transcript.

       ## Examples

       ```python
       db = gffutils.FeatureDB("annotations.db")
       tsses = list(TSSUtils.tss_generator(db))
       assert len(tsses) > 0
       ```
       """
       for transcript in db.features_of_type("transcript"):
           yield TSS(asinterval(transcript), upstream=1, downstream=0)


   @staticmethod
   def _fix_chr_to_match_bigwig(
       tsses_bedtool: pybedtools.BedTool,
       bigwig_path: str,
   ) -> pybedtools.BedTool:
       """
       Adjusts the chromosome naming of each interval to match a BigWig naming style.

       Inspects chromosome names in the BigWig file: adds `'chr'` prefix when missing if BigWig
       uses `'chr1'`-style, or strips it when BigWig uses bare `'1'`-style.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `tsses_bedtool` | `pybedtools.BedTool` | — | BedTool whose chromosome field may need normalization. |
       | `bigwig_path` | `str` | — | Path to the BigWig file used as naming reference. |

       ## Returns

       `pybedtools.BedTool` — New BedTool with adjusted chromosome naming.

       ## Raises

       - `ImportError` — If `pyBigWig` is not installed.
       - `RuntimeError` — If the BigWig file cannot be opened.

       ## Examples

       ```python
       fixed = TSSUtils._fix_chr_to_match_bigwig(tsses_window, "signal.bigWig")
       assert isinstance(fixed, pybedtools.BedTool)
       ```
       """
       import pyBigWig  

       with pyBigWig.open(bigwig_path) as bw:
           bw_chroms = list(bw.chroms().keys())


       bigwig_uses_chr = any(c.startswith("chr") for c in bw_chroms)


       fixed_intervals: List[pybedtools.Interval] = []
       for g in tsses_bedtool:
           fields = list(g.fields)
           chrom = fields[0]


           if bigwig_uses_chr and not chrom.startswith("chr"):
               chrom = "chr" + chrom
           elif (not bigwig_uses_chr) and chrom.startswith("chr"):
               chrom = chrom[3:]


           fields[0] = chrom
           fixed_intervals.append(pybedtools.create_interval_from_list(fields))


       return pybedtools.BedTool(fixed_intervals)

class SignalArrayGenerator:
   """
   Utilities to generate signal arrays around TSS windows from a gffutils database.

   Extracts windowed signal data (ChIP-seq, ATAC-seq) around transcription start sites with
   optional INPUT normalization and CPM (counts-per-million) scaling.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `generate_arrays_from_tsses_from_db()` | Generates TSS features and windowed signal arrays with optional INPUT subtraction. |
   """


   @staticmethod
   def generate_arrays_from_tsses_from_db(
       dbPath: str,
       filePath: str,
       extensionFile: str,
       genome: str,
       bins: int,
       x: np.ndarray | List[float] | List[int],
       inputSignalPath: Optional[str] = None,
       extensionInput: Optional[str] = None,
       normalize_to_cpm: bool = True,
   ):
       """
       Generate TSS features and windowed signal arrays (IP, optional INPUT).

       IP signal is always required via `filePath`/`extensionFile`. INPUT is optional via
       `inputSignalPath`/`extensionInput`. Returns a 4-tuple whose second element depends on
       whether INPUT is provided:

       - **Without INPUT:** `(features, ip_array_good, tsses, tsses_window)`
       - **With INPUT:** `(features, {"ip": ip_array_good, "input": input_array_good}, tsses, tsses_window)`

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `dbPath` | `str` | — | Path to gffutils FeatureDB database. |
       | `filePath` | `str` | — | Path to IP signal file. |
       | `extensionFile` | `str` | — | File type label for IP (e.g., `"bigwig"`, `"bam"`). |
       | `genome` | `str` | — | Genome ID for `pybedtools.slop` (e.g., `"mm10"`, `"hg19"`). |
       | `bins` | `int` | — | Number of bins for signal aggregation across the window. |
       | `x` | `array-like` | — | X-axis coordinates describing the window around each TSS. |
       | `inputSignalPath` | `str` | `None` | Path to INPUT signal file (optional). |
       | `extensionInput` | `str` | `None` | File type for INPUT; required when `inputSignalPath` is set. |
       | `normalize_to_cpm` | `bool` | `True` | Divide values by `mapped_read_count / 1e6` for CPM normalization. |

       ## Raises

       - `ValueError` — If only one of `inputSignalPath`/`extensionInput` is provided.
       - `ValueError` — If IP and INPUT BigWig files require different chromosome naming styles.

       ## Examples

       ```python
       # IP only
       features, ip_arr, tsses, tsses_win = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
           dbPath="annotations.db",
           filePath="chip.bigWig",
           extensionFile="bigwig",
           genome="mm10",
           bins=200,
           x=np.linspace(-1000, 1000, 200),
       )

       # With INPUT
       features, arrays, tsses, tsses_win = SignalArrayGenerator.generate_arrays_from_tsses_from_db(
           dbPath="annotations.db",
           filePath="chip.bigWig",
           extensionFile="bigwig",
           genome="mm10",
           bins=200,
           x=np.linspace(-1000, 1000, 200),
           inputSignalPath="input.bigWig",
           extensionInput="bigwig",
       )
       assert list(arrays.keys()) == ["ip", "input"]
       ```
       """

       # --- Validate optional INPUT pair ---
       if (inputSignalPath is None) ^ (extensionInput is None):
           raise ValueError(
               "If you pass inputSignalPath you must also pass extensionInput (and vice versa)."
           )


       # --- Base TSS ---
       db = gffutils.FeatureDB(dbPath)
       tsses = pybedtools.BedTool(TSSUtils.tss_generator(db)).saveas()


       # remove_duplicates comes from your local imports (.persistence or .array_helpers)
       remove_duplicates(tsses.fn)
       tsses = pybedtools.BedTool(tsses.fn)


       # --- Window derived from x ---
       x = np.asarray(x)
       upstream = int(abs(x.min()))
       downstream = int(abs(x.max()))


       tsses_window_raw = tsses.slop(
           l=upstream,
           r=downstream,
           genome=genome,
           s=True,
       )


       # --- Match chr naming when using BigWig(s) ---
       def _maybe_fix_chr(
           bedtool: pybedtools.BedTool,
           ext: str,
           path: str,
       ) -> pybedtools.BedTool:
           if str(ext).lower() == "bigwig":
               return TSSUtils._fix_chr_to_match_bigwig(bedtool, path)
           return bedtool


       tsses_for_ip = _maybe_fix_chr(tsses_window_raw, extensionFile, filePath)


       if inputSignalPath is not None:
           tsses_for_input = _maybe_fix_chr(
               tsses_window_raw, extensionInput, inputSignalPath
           )


           # If different, IP and INPUT do not share naming (chr vs no-chr). Fail early.
           if str(tsses_for_ip) != str(tsses_for_input):
               raise ValueError(
                   "IP and INPUT appear to use different chromosome naming (chr vs no-chr). "
                   "A single interval set cannot be used for both. Normalize both files to the "
                   "same naming style (or convert one)."
               )
           tsses_window = tsses_for_ip
       else:
           tsses_window = tsses_for_ip


       # --- Signals and arrays ---
       processes = multiprocessing.cpu_count()


       ip_signal = genomic_signal(filePath, extensionFile)
       ip_arr = ip_signal.array(tsses_window, bins=bins, processes=processes)


       ip_array_good: List[Dict[str, Any]] = []
       for item in ip_arr:
           if not item:
               continue
           item0 = item[0]
           if normalize_to_cpm:
               item0["values"] /= ip_signal.mapped_read_count() / 1e6
           ip_array_good.append(item0)


       if inputSignalPath is None:
           features = tsses
           return features, ip_array_good, tsses, tsses_window


       input_signal = genomic_signal(inputSignalPath, extensionInput)  # type: ignore[arg-type]
       in_arr = input_signal.array(tsses_window, bins=bins, processes=processes)


       input_array_good: List[Dict[str, Any]] = []
       for item in in_arr:
           if not item:
               continue
           item0 = item[0]
           if normalize_to_cpm:
               item0["values"] /= input_signal.mapped_read_count() / 1e6
           input_array_good.append(item0)


       features = tsses
       arrays = {"ip": ip_array_good, "input": input_array_good}
       return features, arrays, tsses, tsses_window

class SignalArrays:
   """
   Helper methods for transforming arrays returned by `genomic_signal().array(...)`.

   Converts nested dictionary arrays (output from genomic signal extraction) into flat NumPy
   matrices suitable for visualization and statistical analysis.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `values_array()` | Converts list-of-dicts signal arrays into NumPy matrices. |
   """


   @staticmethod
   def values_array(
       array_ip: List[Dict[str, Any]],
       array_input: Optional[List[Dict[str, Any]]] = None,
   ):
       """
       Convert dict-type signal arrays into pure NumPy arrays of values.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `array_ip` | `list[dict]` | — | IP signal array; each entry must contain a `'values'` key with a 1-D numeric array. |
       | `array_input` | `list[dict]` | `None` | Optional INPUT array with the same structure as `array_ip`. |

       ## Returns

       `np.ndarray` — Shape `(N, B)` of IP values when `array_input` is `None`.
       `tuple[np.ndarray, np.ndarray]` — `(IP_values, INPUT_values)` both with shape `(N, B)` when INPUT is provided.

       ## Raises

       - `KeyError` — If any dict entry lacks the `'values'` key.

       ## Examples

       ```python
       ip_vals = SignalArrays.values_array(ip_array_good)
       assert ip_vals.shape[0] == len(ip_array_good)

       ip_vals, in_vals = SignalArrays.values_array(ip_array_good, input_array_good)
       assert ip_vals.shape == in_vals.shape
       ```
       """
       arrays_ip = [x["values"] for x in array_ip]
       arrays_ip = np.array(arrays_ip)


       if array_input is None:
           return arrays_ip


       arrays_input = [y["values"] for y in array_input]
       arrays_input = np.array(arrays_input)
       return arrays_ip, arrays_input

class PeakAnalysis:
   """
   Functions for computing peak-like enrichment between IP and INPUT signal arrays.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `calculate_peaks_with_gene_name()` | Computes `IP − INPUT` per gene, preserving gene metadata. |
   | `calculate_peaks()` | Computes `IP − INPUT` for two numeric matrices. |
   """


   @staticmethod
   def calculate_peaks_with_gene_name(
       arrays_ip: List[Dict[str, Any]],
       arrays_input: List[Dict[str, Any]],
   ) -> np.ndarray:
       """
       Compute per-gene `IP − INPUT` difference while preserving gene metadata.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `arrays_ip` | `list[dict]` | — | IP signal as list of dicts. Must include keys: `'values'`, `'gene_name'`, `'chr'`, `'start'`, `'end'`, `'strand'`, `'id'`. |
       | `arrays_input` | `list[dict]` | — | INPUT signal as list of dicts. Must include `'values'` and align 1-to-1 with `arrays_ip`. |

       ## Returns

       `np.ndarray` — 1-D array of dicts; each entry has `'values'` set to `IP − INPUT` plus metadata copied from `arrays_ip`.

       ## Examples

       ```python
       out = PeakAnalysis.calculate_peaks_with_gene_name(ip_array_good, input_array_good)
       assert len(out) == len(ip_array_good)
       ```
       """
       normalized_subtracted = np.array([])
       for index, _ in enumerate(arrays_ip):
           normalized_subtracted = np.append(
               normalized_subtracted,
               {
                   "values": arrays_ip[index]["values"] - arrays_input[index]["values"],
                   "gene_name": arrays_ip[index]["gene_name"],
                   "chr": arrays_ip[index]["chr"],
                   "start": arrays_ip[index]["start"],
                   "end": arrays_ip[index]["end"],
                   "strand": arrays_ip[index]["strand"],
                   "id": arrays_ip[index]["id"],
               },
           )
       return normalized_subtracted
   @staticmethod
   def calculate_peaks(arrays_ip: np.ndarray, arrays_input: np.ndarray) -> np.ndarray:
       """
       Compute per-row difference `IP − INPUT` for two numeric matrices.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `arrays_ip` | `np.ndarray` | — | IP signal values with shape `(N, B)`. |
       | `arrays_input` | `np.ndarray` | — | INPUT signal values with shape `(N, B)`. |

       ## Returns

       `np.ndarray` — Difference matrix with shape `(N, B)`.

       ## Examples

       ```python
       peaks = PeakAnalysis.calculate_peaks(ip_vals, in_vals)
       assert peaks.shape == ip_vals.shape
       ```
       """
       if len(arrays_ip) != len(arrays_input):
           raise ValueError(
               "arrays_ip and arrays_input have different lengths (%d vs %d)"
               % (len(arrays_ip), len(arrays_input))
           )
       peaks = arrays_ip - arrays_input
       return peaks

class Plotting:
    """
    High-level plotting methods for genomic signal visualization across multi-omics data.

    All methods are static and return a `matplotlib.figure.Figure` (or the figure-like object
    produced by `imshow_panels`). They can be called without instantiating the class.

    ## Methods

    | Method | Description |
    |--------|-------------|
    | `heatmap_no_sorted()` | Unsorted heatmap of signal values. |
    | `heatmap_sorted_by_meanValues()` | Heatmap sorted by row mean signal intensity. |
    | `heatmap_sorted_by_maxValueIndex()` | Heatmap sorted by position of maximum value per row. |
    | `distance_from_tss_chipSeq()` | Line profiles showing average IP/INPUT signal vs TSS distance. |
    | `twoSignalsOverlay()` | Dual-axis overlay of multiple ChIP-seq and ATAC/DNase average profiles. |
    | `atac_or_chip_with_rna()` | Heatmap with RNA expression context (up/down/unchanged genes). |
    | `all_signal_together_new()` | Three 2-D scatter plots comparing RNA, ChIP, and ATAC log2FCs. |
    | `chip_atac_rna_diverse()` | Calls `atac_or_chip_with_rna` for both ChIP and ATAC signals. |
    | `decile_boxplot_rna_scatter()` | Dual-axis boxplot (ChIP by deciles) + RNA scatter points. |
    | `rna_violin_with_significance()` | Violin plot of RNA log2FC with significant genes overlaid. |
    | `rna_violin_multi_panel()` | Multiple RNA log2FC violins side by side sharing one axes. |
    | `generate_gene_scale_regions_profile()` | Wrapper for the deepTools scale-regions meta-profile figure. |
    | `volcano_with_intensity_alpha()` | Volcano plot(s) with optional intensity-based alpha transparency. |
    | `circplot()` | Circos plot for multi-omics comparative visualization. |
    """

    @staticmethod
    def _ensure_list_kwargs(
        kw: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]],
        n: int,
    ) -> Optional[List[Dict[str, Any]]]:
        """
        For categorical mode where imshow_panels will call _pick_kwargs(..., col):
        - None -> None (let imshow_panels default it)
        - dict -> broadcast to list length n
        - list/tuple -> validate length n
        """
        if kw is None:
            return None
        if isinstance(kw, dict):
            return [dict(kw) for _ in range(n)]
        kw_list = list(kw)
        if len(kw_list) != n:
            raise ValueError(f"Expected {n} kwargs entries, got {len(kw_list)}.")
        return [dict(x) for x in kw_list]

    @staticmethod
    def heatmap_no_sorted(
        signal_values: np.ndarray,
        xAxis: np.ndarray,
        subset_by: Optional[np.ndarray] = None,
        subset_order: Optional[Sequence[str]] = None,
        title: str = "",
        line_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        fill_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        **kwargs,
    ):
        """
        Display signal values as an unsorted heatmap with average profile above.

        When `subset_by` is provided, data are split into one panel per category
        (e.g. FACT high / FACT low) without additional sorting.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `signal_values` | `np.ndarray` | — | 2-D array of shape `(N_genes, N_bins)` with signal. |
        | `xAxis` | `np.ndarray` | — | 1-D array of x-axis coordinates (e.g. distance from TSS). |
        | `subset_by` | `np.ndarray` | `None` | 1-D array of category labels per gene. |
        | `subset_order` | `sequence[str]` | `None` | Ordered list controlling panel order and visible categories. |
        | `title` | `str` | `""` | Figure title (used in single-panel mode). |
        | `line_kwargs` | `dict` or `list[dict]` | `None` | Kwargs forwarded to the average-profile line. |
        | `fill_kwargs` | `dict` or `list[dict]` | `None` | Kwargs forwarded to the confidence-interval fill. |
        | `**kwargs` | — | — | Additional keyword arguments passed to `imshow_panels`. |

        ## Returns

        `matplotlib.figure.Figure` — Figure produced by `imshow_panels`.

        ## Examples

        ```python
        # Single panel (no categories)
        fig = Plotting.heatmap_no_sorted(signal_matrix, x_axis, title="ChIP signal")

        # Two panels (split by FACT category)
        fig = Plotting.heatmap_no_sorted(
            signal_matrix, x_axis,
            subset_by=fact_labels, subset_order=["high", "low"],
            fill_kwargs=[dict(alpha=0.25), dict(alpha=0.25)],
        )
        ```
        """
        # ---- 1 panel / no categories: pass through (None is OK) ----
        if subset_by is None:
            return imshow_panels(
                signal_values,
                x=xAxis,
                vmin=5,
                vmax=99,
                percentile=True,
                subset_by=None,
                subset_order=None,
                single_category_title=title,
                line_kwargs=line_kwargs,   # can be None or dict
                fill_kwargs=fill_kwargs,   # can be None or dict
                **kwargs,
            )

        # ---- categorical mode: need per-column kwargs list or None ----
        # number of groups that will actually be plotted
        subset_by_arr = np.asarray(subset_by, dtype=object)
        if subset_order is None:
            seen = []
            for c in subset_by_arr.tolist():
                if c not in seen:
                    seen.append(c)
            groups = seen
        else:
            groups = [g for g in subset_order if g in subset_by_arr]

        n_groups = len(groups)
        # imshow_panels uses _pick_kwargs(..., col) where col = i*n_groups + j,
        # so total columns per FACT = n_groups; for 1 FACT it's n_groups.
        # If user passes dict, broadcast to n_groups is enough for single FACT.
        # If arr is dict of FACTs, imshow_panels expects list length n_groups*n_facts.
        # Easiest: broadcast to n_groups * n_facts (imshow_panels itself knows n_facts after it normalizes).
        # We'll approximate n_facts from signal_values: if passed as dict later, user should pass list already.
        # Here signal_values is np.ndarray (your Plotting methods pass arrays), so n_facts=1.
        n_cols = n_groups

        line_list = Plotting._ensure_list_kwargs(line_kwargs, n_cols)
        fill_list = Plotting._ensure_list_kwargs(fill_kwargs, n_cols)

        return imshow_panels(
            signal_values,
            x=xAxis,
            vmin=5,
            vmax=99,
            percentile=True,
            subset_by=subset_by_arr,
            subset_order=list(subset_order) if subset_order is not None else None,
            line_kwargs=line_list,
            fill_kwargs=fill_list,
            **kwargs,
        )

    @staticmethod
    def heatmap_sorted_by_meanValues(
        signal_values: np.ndarray,
        xAxis: np.ndarray,
        subset_by: Optional[np.ndarray] = None,
        subset_order: Optional[Sequence[str]] = None,
        title: str = "",
        line_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        fill_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        **kwargs,
    ):
        """
        Display signal values as a heatmap sorted by each row's mean signal.

        Rows with the highest mean signal appear at the top. When `subset_by` is provided,
        each category panel is sorted independently.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `signal_values` | `np.ndarray` | — | 2-D array of shape `(N_genes, N_bins)` with signal data. |
        | `xAxis` | `np.ndarray` | — | 1-D array of x-axis coordinates, shape `(N_bins,)`. |
        | `subset_by` | `np.ndarray` | `None` | Category label per gene; splits into one panel per category. |
        | `subset_order` | `sequence[str]` | `None` | Ordered list of category labels for panel ordering. |
        | `title` | `str` | `""` | Figure title for single-panel mode. |
        | `line_kwargs` | `dict` or `list[dict]` | `None` | Forwarded to average-profile line plot. |
        | `fill_kwargs` | `dict` or `list[dict]` | `None` | Forwarded to confidence-interval fill. |
        | `**kwargs` | — | — | Additional keyword arguments passed to `imshow_panels`. |

        ## Returns

        `matplotlib.figure.Figure` — Figure produced by `imshow_panels`.

        ## Examples

        ```python
        # Sort all genes by mean signal
        fig = Plotting.heatmap_sorted_by_meanValues(signal_matrix, x_axis)

        # Sort within categories
        fig = Plotting.heatmap_sorted_by_meanValues(
            signal_matrix, x_axis,
            subset_by=categories, subset_order=["high", "low"],
        )
        ```
        """
        sort_by = signal_values.mean(axis=1)

        if subset_by is None:
            return imshow_panels(
                signal_values,
                x=xAxis,
                vmin=5,
                vmax=99,
                percentile=True,
                subset_by=None,
                subset_order=None,
                sort_by=sort_by,
                single_category_title=title,
                line_kwargs=line_kwargs,
                fill_kwargs=fill_kwargs,
                **kwargs,
            )

        subset_by_arr = np.asarray(subset_by, dtype=object)
        if subset_order is None:
            seen = []
            for c in subset_by_arr.tolist():
                if c not in seen:
                    seen.append(c)
            groups = seen
        else:
            groups = [g for g in subset_order if g in subset_by_arr]

        n_groups = len(groups)
        n_cols = n_groups

        line_list = Plotting._ensure_list_kwargs(line_kwargs, n_cols)
        fill_list = Plotting._ensure_list_kwargs(fill_kwargs, n_cols)

        return imshow_panels(
            signal_values,
            x=xAxis,
            vmin=5,
            vmax=99,
            percentile=True,
            subset_by=subset_by_arr,
            subset_order=list(subset_order) if subset_order is not None else None,
            sort_by=sort_by,
            line_kwargs=line_list,
            fill_kwargs=fill_list,
            **kwargs,
        )

    @staticmethod
    def heatmap_sorted_by_maxValueIndex(
        signal_values: np.ndarray,
        xAxis: np.ndarray,
        subset_by: Optional[np.ndarray] = None,
        subset_order: Optional[Sequence[str]] = None,
        title: str = "",
        line_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        fill_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        **kwargs,
    ):
        """
        Display signal values as a heatmap sorted by the column index of each row's maximum value.

        Genes whose peak signal occurs earliest (leftmost bin) appear at the top, creating a
        characteristic diagonal pattern useful for visualising sharp binding profiles.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `signal_values` | `np.ndarray` | — | 2-D array of shape `(N_genes, N_bins)` with signal data. |
        | `xAxis` | `np.ndarray` | — | 1-D array of x-axis coordinates, shape `(N_bins,)`. |
        | `subset_by` | `np.ndarray` | `None` | Category label per gene; splits into one panel per category. |
        | `subset_order` | `sequence[str]` | `None` | Ordered list of category labels for panel ordering. |
        | `title` | `str` | `""` | Figure title for single-panel mode. |
        | `line_kwargs` | `dict` or `list[dict]` | `None` | Forwarded to average-profile line plot. |
        | `fill_kwargs` | `dict` or `list[dict]` | `None` | Forwarded to confidence-interval fill. |
        | `**kwargs` | — | — | Additional keyword arguments passed to `imshow_panels`. |

        ## Returns

        `matplotlib.figure.Figure` — Figure produced by `imshow_panels`.

        ## Examples

        ```python
        # Show diagonal pattern of peak positions
        fig = Plotting.heatmap_sorted_by_maxValueIndex(signal_matrix, x_axis)

        # Within categories with different peak positions
        fig = Plotting.heatmap_sorted_by_maxValueIndex(
            signal_matrix, x_axis,
            subset_by=categories, subset_order=["early", "late"],
        )
        ```
        """
        sort_by = np.argmax(signal_values, axis=1)

        if subset_by is None:
            return imshow_panels(
                signal_values,
                x=xAxis,
                vmin=5,
                vmax=99,
                percentile=True,
                subset_by=None,
                subset_order=None,
                sort_by=sort_by,
                single_category_title=title,
                line_kwargs=line_kwargs,
                fill_kwargs=fill_kwargs,
                **kwargs,
            )

        subset_by_arr = np.asarray(subset_by, dtype=object)
        if subset_order is None:
            seen = []
            for c in subset_by_arr.tolist():
                if c not in seen:
                    seen.append(c)
            groups = seen
        else:
            groups = [g for g in subset_order if g in subset_by_arr]

        n_groups = len(groups)
        n_cols = n_groups

        line_list = Plotting._ensure_list_kwargs(line_kwargs, n_cols)
        fill_list = Plotting._ensure_list_kwargs(fill_kwargs, n_cols)

        return imshow_panels(
            signal_values,
            x=xAxis,
            vmin=5,
            vmax=99,
            percentile=True,
            subset_by=subset_by_arr,
            subset_order=list(subset_order) if subset_order is not None else None,
            sort_by=sort_by,
            line_kwargs=line_list,
            fill_kwargs=fill_list,
            **kwargs,
        )


    @staticmethod
    def distance_from_tss_chipSeq(
        arrays_ip: np.ndarray,
        arrays_input: np.ndarray,
        xAxes: np.ndarray,
        name: str = "",
    ):
        """
        Plot the average IP and INPUT profile vs distance from TSS.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `arrays_ip` | `np.ndarray` | — | IP values with shape `(N, B)`. |
        | `arrays_input` | `np.ndarray` | — | INPUT values with shape `(N, B)`. |
        | `xAxes` | `np.ndarray` | — | X-axis positions (typically distance from TSS in bp). |
        | `name` | `str` | `""` | Figure title. |

        ## Returns

        `matplotlib.figure.Figure` — Created matplotlib figure.

        ## Examples

        ```python
        fig = Plotting.distance_from_tss_chipSeq(ip_vals, in_vals, x_axis, name="My sample")
        assert fig is not None
        ```
        """
        x = xAxes
        fig = plt.figure(figsize=(7.68, 5.91))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.18, right=0.98, top=0.9, bottom=0.15)

        ax.plot(x, arrays_ip.mean(axis=0), color="r", label="IP")
        ax.plot(x, arrays_input.mean(axis=0), color="k", label="input")


        ax.axvline(0, linestyle=":", color="k")
        ax.set_xlabel("Distance from TSS (bp)", fontsize=22)
        ax.set_ylabel("Average read coverage\n(per million mapped reads)", fontsize=22)
        ax.legend(loc="best", fontsize=20)
        ax.tick_params(labelsize=17.6)
        ax.set_title(name, fontsize=22)
        return fig


    @staticmethod
    def twoSignalsOverlay(
        dnase_arrays: List[np.ndarray],
        chip_arrays: List[np.ndarray],
        xAxes: np.ndarray,
        name: str = "",
        category_names: Optional[List[str]] = None,
        fontsize: int = 12,
        extra_dnase_array: Optional[List[np.ndarray]] = None,
        extra_dnase_label: str = "Extra DNase",
        extra_dnase_color: str = "#1f77b4",
        extra_chip_array: Optional[List[np.ndarray]] = None,
        extra_chip_label: str = "Extra ChIP",
        extra_chip_color: str = "#ff7f0e",
    ):
        """
        Plot multiple DNase/ATAC-seq and ChIP-seq profiles on two Y-axes sharing the same X-axis.

        ATAC/DNase profiles use the left Y-axis; ChIP-seq uses the right Y-axis. Categories are
        differentiated by line styles (`-`, `--`, `:`, `-.`). A legend appears only when multiple
        categories are present or an extra signal is supplied.

        Optionally, an extra DNase/ATAC signal and/or an extra ChIP signal can be overlaid (e.g.
        an input/control track), each split into one array per category just like `dnase_arrays`/
        `chip_arrays`. Each extra array is plotted on its own axis (`extra_dnase_array` on the
        left axis, `extra_chip_array` on the right axis) so the curves keep the correct scale,
        using its own color, sharing the same per-category line styles, with a single legend
        entry for the whole extra signal.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `chip_arrays` | `list[np.ndarray]` | — | List of ChIP-seq arrays, each with shape `(N, B)`. |
        | `dnase_arrays` | `list[np.ndarray]` | — | List of DNase/ATAC-seq arrays, same length as `chip_arrays`. |
        | `xAxes` | `np.ndarray` | — | X-axis positions. |
        | `name` | `str` | `""` | Figure title. |
        | `category_names` | `list[str]` | `None` | Names for each category (e.g., `["high", "low"]`). |
        | `fontsize` | `int` | `12` | Base font size for axis labels, title, and legend. |
        | `extra_dnase_array` | `list[np.ndarray]` | `None` | Optional extra DNase/ATAC arrays, one per category (same length/order as `category_names`), plotted on the left axis. |
        | `extra_dnase_label` | `str` | `"Extra DNase"` | Legend label for the extra DNase signal. |
        | `extra_dnase_color` | `str` | `"#1f77b4"` | Color for the extra DNase signal. |
        | `extra_chip_array` | `list[np.ndarray]` | `None` | Optional extra ChIP arrays, one per category (same length/order as `category_names`), plotted on the right axis. |
        | `extra_chip_label` | `str` | `"Extra ChIP"` | Legend label for the extra ChIP signal. |
        | `extra_chip_color` | `str` | `"#ff7f0e"` | Color for the extra ChIP signal. |

        ## Returns

        `matplotlib.figure.Figure` — Created matplotlib figure.

        ## Examples

        ```python
        fig = Plotting.twoSignalsOverlay(
            [dnase_high, dnase_low],
            [chip_high, chip_low],
            x_axis,
            name="ChIP vs ATAC",
            category_names=["high", "low"],
            extra_dnase_array=[dnase_input_high, dnase_input_low],
            extra_dnase_label="Input DNase",
            extra_chip_array=[chip_input_high, chip_input_low],
            extra_chip_label="Input ChIP",
        )
        ```
        """
        import matplotlib.lines as mlines
        from matplotlib.ticker import ScalarFormatter

        x = xAxes
        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.182, right=0.888, top=0.9, bottom=0.119)
        atac_color = "#000000"  
        chip_color = "#d62728" 

        # Line styles for categories
        line_styles = ["-", "--", ":", "-."]

        # Use category names for legend, or generate default names
        if category_names is None:
            category_names = [f"Category {i}" for i in range(len(dnase_arrays))]


        # Plot ATAC/DNase profiles (left y-axis)
        for i, (dnase_array, name_cat) in enumerate(zip(dnase_arrays, category_names)):
            if dnase_array.shape[0] > 0:  # Only plot if data exists
                dnase_mean = dnase_array.mean(axis=0)
                linestyle = line_styles[i % len(line_styles)]
                ax.plot(x, dnase_mean, color=atac_color, linestyle=linestyle)

        # Optional single extra DNase/ATAC signal (left y-axis, own color)
        if extra_dnase_array is not None:
            for i, (extra_dnase_array, name_cat) in enumerate(zip(extra_dnase_array, category_names)):
                if extra_dnase_array.shape[0] > 0:  # Only plot if data exists
                    extra_dnase_mean = extra_dnase_array.mean(axis=0)
                    linestyle = line_styles[i % len(line_styles)]
                    ax.plot(x, extra_dnase_mean, color=extra_dnase_color, linestyle=linestyle, linewidth=2)

        ax.set_xlabel("Distance from the TSS", fontsize=fontsize*1.1)
        ax.set_ylabel("ATAC-seq signal intensity", color=atac_color, fontsize=fontsize*1.2)
        ax.axvline(0, linestyle=":", color="k", alpha=0.3)
        ax.set_title(name, fontsize=fontsize*1.3)

        # Secondary axis for ChIP-seq (right y-axis)
        ax2 = ax.twinx()
        for i, (chip_array, name_cat) in enumerate(zip(chip_arrays, category_names)):
            if chip_array.shape[0] > 0:  # Only plot if data exists
                chip_mean = chip_array.mean(axis=0)
                linestyle = line_styles[i % len(line_styles)]
                ax2.plot(x, chip_mean, color=chip_color, linestyle=linestyle)

        # Optional single extra ChIP signal (right y-axis, own color)
        if extra_chip_array is not None:
            for i, (extra_chip_array, name_cat) in enumerate(zip(extra_chip_array, category_names)):
                if extra_chip_array.shape[0] > 0:  # Only plot if data exists
                    extra_chip_mean = extra_chip_array.mean(axis=0)
                    linestyle = line_styles[i % len(line_styles)]
                    ax2.plot(x, extra_chip_mean, color=extra_chip_color, linestyle=linestyle, linewidth=2)

        ax2.set_ylabel("ChIP-seq Signal intensity", color=chip_color, fontsize=fontsize*1.2)

        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        # Apply scientific formatter to both axes
        ax.yaxis.set_major_formatter(formatter)
        ax2.yaxis.set_major_formatter(formatter)
        # Ensure offset/exponent text (the scientific '10^n' part) uses the same font size
        try:
            ax.yaxis.get_offset_text().set_fontsize(int(fontsize * 1.2))
            ax2.yaxis.get_offset_text().set_fontsize(int(fontsize * 1.2))
        except Exception:
            # some backends/versions may not expose offset_text yet; ignore safely
            pass
        ax.tick_params(axis="y", labelcolor=atac_color, labelsize=int(fontsize * 1.2))
        ax.tick_params(axis="x", labelcolor=atac_color, labelsize=int(fontsize * 1.2))
        ax2.tick_params(axis="y", labelcolor=chip_color, labelsize=int(fontsize * 1.2))

        # Legend if multiple categories and/or an extra signal is present
        has_extra_dnase = extra_dnase_array is not None 
        has_extra_chip = extra_chip_array is not None
        if len(category_names) > 1 or has_extra_dnase or has_extra_chip:
            line_styles_list = ["-", "--", ":", "-."]
            legend_lines = []
            if len(category_names) > 1:
                for i, cat_name in enumerate(category_names):
                    linestyle = line_styles_list[i % len(line_styles_list)]
                    # Create a line that represents this category (showing the linestyle)
                    # Using black color to show it applies to both ATAC (blue) and ChIP (red) signals
                    line = mlines.Line2D([], [], color="black", linestyle=linestyle,
                                        label=cat_name.capitalize(), marker="")
                    legend_lines.append(line)
            if has_extra_dnase:
                legend_lines.append(mlines.Line2D(
                    [], [], color=extra_dnase_color, linestyle="-", linewidth=2,
                    label=extra_dnase_label, marker="",
                ))
            if has_extra_chip:
                legend_lines.append(mlines.Line2D(
                    [], [], color=extra_chip_color, linestyle="-", linewidth=2,
                    label=extra_chip_label, marker="",
                ))
            ax.legend(
                handles=legend_lines,
                loc="upper left",
                bbox_to_anchor=(0.0, 1.0),
                bbox_transform=ax.transAxes,
                fontsize=fontsize,
                frameon=True,
                fancybox=True,
                shadow=False,
                title="FACT Category" if not has_extra_dnase and not has_extra_chip else "Legend",
                title_fontsize=fontsize,
            )

        fig.tight_layout(pad=1.5)

        return fig

    @staticmethod
    def quadrant_signal_profile(
        arrays: List[np.ndarray],
        xAxes: np.ndarray,
        group_labels: List[str],
        colors: Optional[List[str]] = None,
        name: str = "",
        ylabel: str = "Average ChIP-seq signal\n(CPM)",
        fontsize: int = 12,
    ):
        """
        Plot one average signal profile per gene group on a single axis.

        Used to compare a single signal (typically ChIP-seq) across gene
        groups defined by combinations of two other modalities' log2 fold
        changes (e.g. ATAC-high/RNA-high vs. ATAC-high/RNA-low quadrants).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `arrays` | `list[np.ndarray]` | — | One signal matrix per group, each shape `(N, B)`. |
        | `xAxes` | `np.ndarray` | — | X-axis positions (distance from TSS in bp). |
        | `group_labels` | `list[str]` | — | Legend label per group, same length/order as `arrays`. |
        | `colors` | `list[str] \| None` | `None` | Line color per group. Defaults to a 4-color qualitative palette. |
        | `name` | `str` | `""` | Figure title. |
        | `ylabel` | `str` | `"Average ChIP-seq signal\\n(CPM)"` | Y-axis label. |
        | `fontsize` | `int` | `12` | Base font size for labels/legend. |

        ## Returns

        `matplotlib.figure.Figure` — Created matplotlib figure.

        ## Examples

        ```python
        fig = Plotting.quadrant_signal_profile(
            [atac_hi_rna_hi, atac_hi_rna_lo, atac_lo_rna_lo, atac_lo_rna_hi],
            x_axis,
            group_labels=["ATAC high / RNA high", "ATAC high / RNA low",
                          "ATAC low / RNA low", "ATAC low / RNA high"],
            name="ChIP-seq by ATAC/RNA quadrant",
        )
        ```
        """
        x = xAxes
        default_colors = ["#d62728", "#9467bd", "#1f77b4", "#2ca02c"]
        colors = colors if colors is not None else default_colors

        fig = plt.figure(figsize=(9, 7))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.13)

        for arr, label, color in zip(arrays, group_labels, colors):
            if arr.shape[0] == 0:
                continue
            ax.plot(x, arr.mean(axis=0), color=color, linewidth=2, label=f"{label} (n={arr.shape[0]})")

        ax.axvline(0, linestyle=":", color="k", alpha=0.3)
        ax.set_xlabel("Distance from TSS (bp)", fontsize=fontsize * 1.1)
        ax.set_ylabel(ylabel, fontsize=fontsize * 1.1)
        ax.set_title(name, fontsize=fontsize * 1.3)
        ax.tick_params(labelsize=fontsize)
        ax.legend(loc="best", fontsize=fontsize * 0.85, frameon=True)
        fig.tight_layout(pad=1.5)
        return fig

    @staticmethod
    def atac_or_chip_with_rna(
        signal_values: np.ndarray,
        rna: Any,
        xAxis: np.ndarray,
        signalName: str = "",
        divided_by: str = "",
        threshold: float = 1.0,
        padj_col: Optional[str] = None,
        padj_threshold: float = 0,
        subset_by: Optional[np.ndarray] = None,
        subset_order: Optional[Sequence[str]] = None,
        show_heatmap_panels: bool = None,
        line_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        fill_kwargs: Optional[Union[Dict[str, Any], Sequence[Dict[str, Any]]]] = None,
        counts_output_path: Optional[str] = None,
    ):
        """
        Plot a ChIP/ATAC signal heatmap with RNA expression categories overlaid.

        Genes are classified by RNA log2FC: up-regulated (> `threshold`, orange), down-regulated
        (< `-threshold`, brown), or unchanged (gray). Average profiles per category are shown below
        the heatmap. Optionally splits into multiple panels via `subset_by`/`subset_order`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `signal_values` | `np.ndarray` | — | Signal data matrix, shape `(N_genes, N_bins)`. |
        | `rna` | `object` | — | Object exposing `.data` (DataFrame) and `.lfc_column` (str). |
        | `xAxis` | `np.ndarray` | — | X-axis coordinates describing the signal window. |
        | `signalName` | `str` | `""` | Label for the signal axis/title. |
        | `divided_by` | `str` | `""` | Legend title for the RNA category panel. |
        | `threshold` | `float` | `1.0` | log2FC threshold for classifying genes as up/down. |
        | `padj_col` | `str` | `None` | Column in `rna.data` with adjusted p-values. Padj filtering is applied only when this is set and `padj_threshold` is non-zero. |
        | `padj_threshold` | `float` | `0` | When non-zero (and `padj_col` set), a gene is up/down only if `padj < padj_threshold`; all others are unchanged. Must be a float > 0. Empty/non-numeric padj → not significant. |
        | `subset_by` | `np.ndarray` | `None` | Category array to split data into multiple panels. |
        | `subset_order` | `sequence[str]` | `None` | Ordered list of subset category labels. |
        | `line_kwargs` | `dict` or `list[dict]` | `None` | Line plot kwargs forwarded to `imshow_panels`. |
        | `fill_kwargs` | `dict` or `list[dict]` | `None` | Fill kwargs forwarded to `imshow_panels`. |
        | `show_heatmap_panels` | `bool \| None` | `None` | Whether to render the per-status heatmap panels below the profile. |
        | `counts_output_path` | `str` | `None` | If given, write the up/down/unchanged gene counts to this `.txt` file. Only written (and only printed to stdout) when `show_heatmap_panels` is `False` or `None`. |

        ## Returns

        `matplotlib.figure.Figure` — Figure with heatmap and overlaid RNA category profiles.

        ## Raises

        - `ValueError` — If `signal_values` row count does not match RNA row count.

        ## Examples

        ```python
        # Global view: all genes colored by RNA category
        fig = Plotting.atac_or_chip_with_rna(
            signal_matrix, rna_obj, x_axis,
            signalName="ATAC signal",
            threshold=1.0,
        )

        # Split by FACT status
        fig = Plotting.atac_or_chip_with_rna(
            signal_matrix, rna_obj, x_axis,
            signalName="ChIP signal",
            subset_by=fact_labels,
            subset_order=["high", "low"],
            threshold=0.5,
        )
        ```
        """
        
        x = xAxis
        lfc = rna.data[rna.lfc_column].values

        if signal_values.shape[0] != len(lfc):
            raise ValueError(
                f"signal_values has {signal_values.shape[0]} rows but RNA has {len(lfc)} rows. "
                "Both must be aligned before plotting."
            )

        # -----------------------------
        # RNA masks (optionally filtered by adjusted p-value)
        #
        # When both `padj_col` is provided and `padj_threshold` is non-zero, a gene
        # only counts as up/down if its adjusted p-value is below the threshold.
        # Empty/non-numeric padj values are coerced to NaN, and `NaN < threshold`
        # is False, so such genes are treated as not significant (i.e. unchanged).
        # When padj filtering is off the masks are identical to the lfc-only logic.
        # -----------------------------
        apply_padj = bool(padj_col) and padj_threshold != 0
        if apply_padj:
            if isinstance(padj_threshold, bool) or not isinstance(padj_threshold, (int, float)) or padj_threshold <= 0:
                raise ValueError(
                    f"padj_threshold must be a float greater than 0, got {padj_threshold!r}."
                )
            if padj_col not in rna.data.columns:
                raise ValueError(
                    f"padj_col '{padj_col}' not found in expression table columns."
                )
            padj = pd.to_numeric(rna.data[padj_col], errors="coerce").values
            significant = padj < padj_threshold  # NaN/non-numeric -> False
            up_mask = (lfc > threshold) & significant
            down_mask = (lfc < -threshold) & significant
            unchanged_mask = ~(up_mask | down_mask)  # all the remaining
        else:
            up_mask = lfc > threshold
            down_mask = lfc < -threshold
            unchanged_mask = (lfc >= -threshold) & (lfc <= threshold)
        sort_by = signal_values.mean(axis=1)

        # -----------------------------
        # Global matrices aligned with signal_values.
        # Rows that do not belong to the class are filled with NaN.
        # -----------------------------
        up_data_full = np.full(signal_values.shape, np.nan, dtype=float)
        down_data_full = np.full(signal_values.shape, np.nan, dtype=float)
        unchanged_data_full = np.full(signal_values.shape, np.nan, dtype=float)

        up_data_full[up_mask, :] = signal_values[up_mask, :]
        down_data_full[down_mask, :] = signal_values[down_mask, :]
        unchanged_data_full[unchanged_mask, :] = signal_values[unchanged_mask, :]
        extra_bottom_lines = [
            {
                "data": up_data_full,
                "line_kwargs": dict(color="#fe9829", label="up"),
                "fill_kwargs": dict(color="#fe9829", alpha=0.1),
                "conf": 0.5,
            },
            {
                "data": down_data_full,
                "line_kwargs": dict(color="#8e3104", label="down"),
                "fill_kwargs": dict(color="#8e3104", alpha=0.1),
                "conf": 0.5,
            },
            {
                "data": unchanged_data_full,
                "line_kwargs": dict(color=".5", label="unchanged"),
                "fill_kwargs": dict(color=".5", alpha=0.1),
                "conf": 0.5,
            },
        ]


        # -----------------------------
        # Without subsets
        # -----------------------------
        if subset_by is None:
            counts_lines = [
                f"Up-regulated genes: {np.sum(up_mask)}",
                f"Down-regulated genes: {np.sum(down_mask)}",
                f"Unchanged genes: {np.sum(unchanged_mask)}",
            ]

            fig = imshow_panels(
                signal_values,
                x=x,
                vmin=5,
                vmax=99,
                percentile=True,
                conf=0.5,
                subset_by=None,
                subset_order=None,
                sort_by=sort_by,
                line_kwargs=line_kwargs if line_kwargs is not None else dict(color="k", label="All"),
                fill_kwargs=fill_kwargs if fill_kwargs is not None else dict(color="k", alpha=0.1),
                extra_bottom_lines=extra_bottom_lines,
                extra_bottom_ylabel="Average\nenrichment",
                extra_bottom_xlabel="Distance from TSS (bp)",
                extra_bottom_legend_title = divided_by if divided_by else None,
                subplot_params=dict(wspace=0.15, hspace=0.15),
                single_category_title=signalName,
                show_panels=show_heatmap_panels,

            )

        # -----------------------------
        # With subsets
        # -----------------------------
        else:
            subset_by_arr = np.asarray(subset_by, dtype=object)

            if len(subset_by_arr) != signal_values.shape[0]:
                raise ValueError(
                    f"subset_by has length {len(subset_by_arr)} but signal_values has "
                    f"{signal_values.shape[0]} rows. They must match."
                )

            if subset_order is None:
                seen = []
                for c in subset_by_arr.tolist():
                    if c not in seen:
                        seen.append(c)
                groups = seen
            else:
                groups = [g for g in subset_order if g in subset_by_arr]

            n_groups = len(groups)

            line_list = Plotting._ensure_list_kwargs(line_kwargs, n_groups)
            fill_list = Plotting._ensure_list_kwargs(fill_kwargs, n_groups)

            if line_kwargs is None:
                line_list = [dict(color="k", label=g) for g in groups]
            if fill_kwargs is None:
                fill_list = [dict(color="k", alpha=0.1) for _ in groups]

            counts_lines = []
            for g in groups:
                counts_lines.append(g.capitalize() + ":")
                counts_lines.append(f"\t - Up-regulated genes: {np.sum(up_mask & (subset_by_arr == g))}")
                counts_lines.append(f"\t - Down-regulated genes: {np.sum(down_mask & (subset_by_arr == g))}")
                counts_lines.append(f"\t - Unchanged genes: {np.sum(unchanged_mask & (subset_by_arr == g))}")

            fig = imshow_panels(
                signal_values,
                x=x,
                vmin=5,
                vmax=99,
                percentile=True,
                conf=0.5,
                subset_by=subset_by_arr,
                subset_order=list(subset_order) if subset_order is not None else None,
                sort_by=sort_by,
                line_kwargs=line_list,
                fill_kwargs=fill_list,
                extra_bottom_lines=extra_bottom_lines,
                extra_bottom_ylabel="Average\nenrichment",
                extra_bottom_xlabel="Distance from TSS (bp)",
                extra_bottom_legend_kwargs=dict(ncol=3),
                extra_bottom_legend_title=divided_by if divided_by else None,
                show_panels=show_heatmap_panels,
            )

        if fig.ax_hm:
            fig.ax_hm[0].set_ylabel(signalName + " transcripts")
        if fig.ax_ln:
            fig.ax_ln[0].set_ylabel("Average\nenrichment")
        if fig.ax_hm or fig.ax_ln:
            fig.cax.set_ylabel("Enrichment")

        hm_axes = list(getattr(fig, "ax_hm", []))
        ln_axes = list(getattr(fig, "ax_ln", []))
        extra_axes = list(getattr(fig, "ax_extra", []))

        for ax in hm_axes + ln_axes + extra_axes:
            ax.axvline(0, linestyle=":", color="k")

        # Counts are only reported (printed/exported) when heatmap panels are not
        # shown on their own (show_heatmap_panels is False or None).
        if show_heatmap_panels is not True:
            for line in counts_lines:
                print(line)

            if counts_output_path:
                out_dir = os.path.dirname(counts_output_path)
                if out_dir:
                    os.makedirs(out_dir, exist_ok=True)
                with open(counts_output_path, "w") as fh:
                    fh.write("\n".join(counts_lines) + "\n")

        return fig

    @staticmethod
    def all_signal_together_new(
        rna: Any,
        chip_lfc_column: str,
        atac_lfc_column: str,
        title: str = "",
    ):
        """
        Create three 2-D scatter plots comparing RNA, ChIP, and ATAC log2 fold changes.

        Produces a 1×3 figure: RNA vs ChIP, RNA vs ATAC, and ChIP vs ATAC. Returns `None` with
        a warning if any dataset is empty or a required column is missing.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `rna` | `object` | — | Object exposing `.data` (DataFrame) and `.lfc_column` (str). |
        | `chip_lfc_column` | `str` | — | Column name in `rna.data` with ChIP log2 fold change. |
        | `atac_lfc_column` | `str` | — | Column name in `rna.data` with ATAC log2 fold change. |
        | `title` | `str` | `""` | Figure suptitle. Empty string draws no suptitle. |

        ## Returns

        `matplotlib.figure.Figure | None` — Figure with three scatter panels, or `None` if any input is invalid.

        ## Examples

        ```python
        fig = Plotting.all_signal_together_new(rna_obj, "chip_lfc", "atac_lfc")
        if fig is not None:
            fig.savefig("scatter_comparison.png", dpi=150)
        ```
        """
        if rna is None or not hasattr(rna, "data") or rna.data is None or len(rna.data) == 0:
            print("[Pyntegrate WARNING] Empty RNA dataset → skipping all_signal_together_new")
            return None

        if chip_lfc_column not in rna.data.columns:
            print(rna.data.columns)
            print(f"[Pyntegrate WARNING] Column not found: {chip_lfc_column}")
            return None

        if atac_lfc_column not in rna.data.columns:
            print(f"[Pyntegrate WARNING] Column not found: {atac_lfc_column}")
            return None

        rna_lfc = pd.to_numeric(rna.data[rna.lfc_column], errors="coerce").to_numpy()
        chip_lfc = pd.to_numeric(rna.data[chip_lfc_column], errors="coerce").to_numpy()
        atac_lfc = pd.to_numeric(rna.data[atac_lfc_column], errors="coerce").to_numpy()

        valid = (
            np.isfinite(rna_lfc) &
            np.isfinite(chip_lfc) &
            np.isfinite(atac_lfc)
        )

        rna_lfc = rna_lfc[valid]
        chip_lfc = chip_lfc[valid]
        atac_lfc = atac_lfc[valid]

        if len(rna_lfc) == 0:
            print("[Pyntegrate WARNING] No valid points after filtering NaNs")
            return None

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        axes[0].scatter(rna_lfc, chip_lfc, alpha=0.7)
        axes[0].set_xlabel("RNA-seq log2 Fold Change")
        axes[0].set_ylabel("ChIP-seq log2 Fold Change")
        axes[0].set_title("RNA-seq vs ChIP-seq")

        axes[1].scatter(rna_lfc, atac_lfc, alpha=0.7)
        axes[1].set_xlabel("RNA-seq log2 Fold Change")
        axes[1].set_ylabel("ATAC-seq log2 Fold Change")
        axes[1].set_title("RNA-seq vs ATAC-seq")

        axes[2].scatter(chip_lfc, atac_lfc, alpha=0.7)
        axes[2].set_xlabel("ChIP-seq log2 Fold Change")
        axes[2].set_ylabel("ATAC-seq log2 Fold Change")
        axes[2].set_title("ChIP-seq vs ATAC-seq")

        if title:
            fig.suptitle(title)

        plt.tight_layout()
        return fig

    @staticmethod
    def decile_boxplot_rna_scatter(
        chip_data: "pd.DataFrame",
        scatter_data: "pd.DataFrame",
        chip_col: str = "Chip",
        rna_col: str = "RNA",
        group_col: str = "Group",
        chip_ylabel: str = "ChIP (mean signal)",
        chip_range: Optional[Tuple[float, float]] = None,
        rna_ylabel: str = "RNA (log2 fold change)",
        title: str = "Box Plot (ChIP) + Scatter (RNA)",
        xlabel: str = "Deciles (according to ChIP signal)",
        label_fontsize: int = 25,
        figsize: Tuple[int, int] = (12, 8),
    ):
        """
        Dual-axis figure: boxplot of ChIP signal per decile (left y-axis) and RNA stripplot (right y-axis).

        The two axes share the same categorical x-axis (decile groups). ChIP boxes are drawn on top
        of the RNA scatter points so that both signal distributions can be compared simultaneously.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `chip_data` | `pd.DataFrame` | — | Must contain `chip_col` (numeric) and `group_col` (int) columns. Typically produced by `pandas.qcut`. |
        | `scatter_data` | `pd.DataFrame` | — | Must contain `rna_col` (numeric) and `group_col` (int). Rows with NaN RNA values should be dropped before passing. |
        | `chip_col` | `str` | `"Chip"` | Column in `chip_data` with ChIP signal values. |
        | `rna_col` | `str` | `"RNA"` | Column in `scatter_data` with RNA log2 fold-change values. |
        | `group_col` | `str` | `"Group"` | Decile group label column shared by both DataFrames. |
        | `chip_ylabel` | `str` | `"ChIP (mean signal)"` | Left y-axis label. |
        | `chip_range` | `tuple[float, float]` | `None` | Fixed y-axis limits for ChIP signal. |
        | `rna_ylabel` | `str` | `"RNA (log2 fold change)"` | Right y-axis label (displayed in red). |
        | `title` | `str` | `"Box Plot (ChIP) + Scatter (RNA)"` | Figure title. |
        | `xlabel` | `str` | `"Deciles (according to ChIP signal)"` | X-axis label. |
        | `figsize` | `tuple[int, int]` | `(12, 8)` | Figure size in inches `(width, height)`. |
        | `label_fontsize` | `int` | `25` | Base font size for axis labels and title. |

        ## Returns

        `matplotlib.figure.Figure` — Figure with two y-axes sharing the same x-axis.

        ## Examples

        ```python
        fig = Plotting.decile_boxplot_rna_scatter(
            chip_data, scatter_data,
            chip_ylabel="ChIP log2FC", rna_ylabel="RNA log2FC",
        )
        fig.savefig("decile_plot.png", dpi=150)
        ```
        """
        import seaborn as sns


        fig, ax1 = plt.subplots(figsize=figsize)
        ax2 = ax1.twinx()

        # ------ y-limits for RNA (true min/max) ------
        if len(scatter_data) > 0:
            # rna_lo, rna_hi = np.nanpercentile(scatter_data[rna_col], [1, 99])

            rna_lo = float(np.nanmin(scatter_data[rna_col]))
            rna_hi = float(np.nanmax(scatter_data[rna_col]))
        else:
            rna_lo, rna_hi = 0.0, 1.0
        pad_rna = (rna_hi - rna_lo) * 0.10 if rna_hi > rna_lo else 1.0

        # ------ RNA scatter (drawn first, z-order=1) ------
        sns.stripplot(
            x=group_col, y=rna_col, data=scatter_data,
            ax=ax2, color="red", jitter=0.2, alpha=0.35, size=3, zorder=1,
        )
        print("Max RNA signal:", rna_hi)
        print("Min RNA signal:", rna_lo)
        ax2.set_ylabel(rna_ylabel, color="red", fontsize=label_fontsize)
        ax2.set_ylim(rna_lo - pad_rna, rna_hi + pad_rna)
        ax2.ticklabel_format(style="plain", axis="y", useOffset=False)
        ax2.tick_params(axis="y", colors="red", labelsize=label_fontsize * 0.9)

        # ------ ChIP boxplot (drawn on top, z-order=10) ------
        sns.boxplot(
            x=group_col, y=chip_col, data=chip_data,
            ax=ax1, width=0.6, color="black",
            boxprops=dict(alpha=0.6), showfliers=True, zorder=10,
        )
        chip_p1, chip_p99 = np.nanpercentile(chip_data[chip_col], [1, 99])
        chip_min = np.nanmin(chip_data[chip_col])
        chip_max = np.nanmax(chip_data[chip_col])
        chip_lo = min(chip_p1, chip_min)
        chip_hi = max(chip_p99, chip_max)

        if chip_range is not None:
            ax1.set_ylim(chip_range)
        else:
            pad_chip = (chip_hi - chip_lo) * 0.10 if chip_hi > chip_lo else 1.0
            ax1.set_ylim(chip_lo - pad_chip, chip_hi + pad_chip)
        print("Max ChIP signal:", chip_max)
        print("Min ChIP signal:", chip_min)
        ax1.set_ylabel(chip_ylabel, fontsize=label_fontsize)
        ax1.set_xlabel(xlabel, fontsize=label_fontsize)
        ax1.set_title(title + (" (n="+str(len(scatter_data)) + ")") if len(scatter_data) == len(chip_data) else "", fontsize=label_fontsize + 2)
        ax1.ticklabel_format(style="plain", axis="y", useOffset=False)
        ax1.tick_params(axis="y", labelsize=label_fontsize * 0.9)
        # Make x-axis tick labels large and readable
        ax1.tick_params(axis="x", labelsize=int(label_fontsize * 0.95))

        # keep ChIP axis on top visually
        ax1.set_zorder(2)
        ax1.patch.set_visible(False)

        plt.tight_layout()
        return fig

    @staticmethod
    def rna_violin_with_significance(
        rna_values: "np.ndarray",
        significant_mask: "np.ndarray",
        title: str = "",
        ylabel: str = "RNA log2 fold change",
        sig_color: str = "red",
        sig_alpha: float = 0.6,
        sig_size: float = 4,
        figsize: Tuple[int, int] = (4, 8),
        label_fontsize: int = 25,
    ):
        """
        Violin plot of RNA log2 fold-change values with significant genes overlaid as strip points.

        The violin shows the full distribution; `True` entries in `significant_mask` are drawn
        as colored strip points on top. Quartile lines are shown inside the violin.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `rna_values` | `np.ndarray` | — | 1-D array of RNA log2 fold-change values. |
        | `significant_mask` | `np.ndarray` | — | Boolean 1-D array (same length as `rna_values`); `True` entries are drawn as colored points. |
        | `title` | `str` | `""` | Figure title. |
        | `ylabel` | `str` | `"RNA log2 fold change"` | Y-axis label. |
        | `sig_color` | `str` | `"red"` | Color of significant-gene strip points. |
        | `sig_alpha` | `float` | `0.6` | Transparency of strip points (0–1). |
        | `sig_size` | `float` | `4` | Marker size of strip points in points. |
        | `figsize` | `tuple[int, int]` | `(4, 8)` | Figure size in inches `(width, height)`. |
        | `label_fontsize` | `int` | `25` | Base font size for axis labels and title. |

        ## Returns

        `matplotlib.figure.Figure` — Single-axes figure.

        ## Examples

        ```python
        sig_mask = padj_values < 0.1
        fig = Plotting.rna_violin_with_significance(
            rna_lfc, sig_mask,
            title="FACT-high genes – SSRP1 KD",
            sig_color="crimson",
        )
        fig.savefig("violin_significance.png", dpi=300, bbox_inches="tight")
        ```
        """
        import seaborn as sns

        rna_values = np.asarray(rna_values, dtype=float)
        significant_mask = np.asarray(significant_mask, dtype=bool)

        df_all = pd.DataFrame({"lfc": rna_values, "x": 1})
        df_sig = df_all[significant_mask].copy()

        fig, ax = plt.subplots(figsize=figsize)

        sns.violinplot(
            x="x", y="lfc", data=df_all, ax=ax,
            inner="quart", cut=0, linewidth=1, color="lightgray",
        )
        sns.stripplot(
            x="x", y="lfc", data=df_sig, ax=ax,
            color=sig_color, jitter=0.4, alpha=sig_alpha, size=sig_size,
        )

        y_lo, y_hi = np.nanmin(rna_values), np.nanmax(rna_values)
        pad = (y_hi - y_lo) * 0.10 if y_hi > y_lo else 1.0
        ax.set_ylim(y_lo - pad, y_hi + pad)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.set_xlabel("")
        ax.set_xticks([])
        ax.set_title(title, fontsize=label_fontsize)
        ax.ticklabel_format(style="plain", axis="y", useOffset=False)
        ax.tick_params(axis="y", labelsize=label_fontsize)

        plt.tight_layout()
        return fig

    @staticmethod
    def rna_violin_multi_panel(
        combined: pd.DataFrame,
        panel_titles: List[str],
        y_limits: Tuple[float, float],
        ylabel: str = "RNA log2 fold change",
        figure_title: str = "",
        figsize_per_panel: Tuple[float, float] = (4, 9.0),
        label_fontsize: int = 15,
    ):
        """
        Multiple RNA log2FC violins side by side on a single shared axes, without subplot borders.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `combined` | `pd.DataFrame` | — | Long-format frame with columns `"Panel"`, `"RNA"`, `"Significant"`. |
        | `panel_titles` | `list[str]` | — | Panel names in display order; must match values in `combined["Panel"]`. |
        | `y_limits` | `tuple[float, float]` | — | Shared `(ymin, ymax)` for the y-axis. |
        | `ylabel` | `str` | `"RNA log2 fold change"` | Y-axis label. |
        | `figure_title` | `str` | `""` | Figure suptitle. |
        | `figsize_per_panel` | `tuple[float, float]` | `(4, 9.0)` | Approximate `(width, height)` per panel in inches. |
        | `label_fontsize` | `int` | `15` | Font size for axis labels, title, and ticks. |

        ## Returns

        `matplotlib.figure.Figure` — Single-axes figure with one violin per panel.
        """
        import seaborn as sns

        fig, ax = plt.subplots(figsize=(figsize_per_panel[0] * len(panel_titles), figsize_per_panel[1]))
        sns.violinplot(
            data=combined,
            x="Panel",
            y="RNA",
            order=panel_titles,
            ax=ax,
            inner="quart",
            cut=0,
            linewidth=1,
            color="lightgray",
        )
        sns.stripplot(
            data=combined[combined["Significant"]],
            x="Panel",
            y="RNA",
            order=panel_titles,
            ax=ax,
            color="red",
            jitter=0.35,
            alpha=0.6,
            size=4,
        )

        ax.axhline(0, color="black", linewidth=0.8, alpha=0.5)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        ax.set_xlabel("")
        ax.set_ylim(*y_limits)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(False)
        ax.tick_params(axis="x", rotation=0, labelsize=label_fontsize)
        ax.ticklabel_format(style="plain", axis="y", useOffset=False)
        ax.tick_params(axis="y", labelsize=label_fontsize)
        fig.suptitle(figure_title + " (n=" + str(len(combined)) + ")", fontsize=label_fontsize)
        plt.tight_layout()
        return fig

    @staticmethod
    def generate_gene_scale_regions_profile(
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
        skip_zeros: bool = True,
        processes: Optional[int] = None,
        figsize: Tuple[float, float] = (4.7244, 3.1496),
        dpi: int = 300,
        ylabel: str = "ChIP-seq Signal",
        title: str = "",
        legend_location: str = "upper right",
        show: bool = True,
    ) -> Optional[str]:
        """
        Generate a deepTools scale-regions meta-profile figure.

        Thin wrapper delegating to `Pyntegrate.plotutils.generate_gene_scale_regions_profile`,
        which runs `computeMatrix scale-regions` followed by `plotProfile`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `genes_bed_path` | `str` | — | Path to a BED file of gene/region coordinates. May be gzip-compressed. |
        | `signal_paths` | `list[str]` | — | Paths to BigWig signal files, one per sample. |
        | `sample_labels` | `list[str]` | — | Sample labels shown in the plot legend, one per entry in `signal_paths`. |
        | `colors` | `list[str]` | — | Line colors, one per entry in `signal_paths`. |
        | `output_path` | `str` | `None` | Path to save the figure. If `None`, the figure is written to a temporary directory only. |
        | `regions_label` | `str` | `None` | Label for the region set (deepTools `--regionsLabel`). |
        | `before_region_start_length` | `int` | `500` | Distance upstream of the region start to include (bp). |
        | `region_body_length` | `int` | `1000` | Length the gene/region body is scaled to (bp). |
        | `after_region_start_length` | `int` | `500` | Distance downstream of the region end to include (bp). |
        | `bin_size` | `int` | `10` | Bin size for signal aggregation (bp). |
        | `skip_zeros` | `bool` | `True` | Passed to `computeMatrix --skipZeros` to drop all-zero regions. |
        | `processes` | `int` | `None` | Number of processors for `computeMatrix` (defaults to 1). |
        | `figsize` | `tuple[float, float]` | `(4.7244, 3.1496)` | Figure size in inches; converted to cm for deepTools. |
        | `dpi` | `int` | `300` | Output resolution in dots per inch. |
        | `ylabel` | `str` | `"ChIP-seq Signal"` | Y-axis label for the profile plot. |
        | `title` | `str` | `""` | Plot title. |
        | `legend_location` | `str` | `"upper right"` | Legend placement, passed to deepTools `plotProfile`. |
        | `show` | `bool` | `True` | Whether to display the figure interactively. |

        ## Returns

        `str | None` — Absolute path of the saved figure, or `None` when
        `output_path` was not provided.
        """
        return generate_gene_scale_regions_profile(
            genes_bed_path=genes_bed_path,
            signal_paths=signal_paths,
            sample_labels=sample_labels,
            colors=colors,
            output_path=output_path,
            regions_label=regions_label,
            before_region_start_length=before_region_start_length,
            region_body_length=region_body_length,
            after_region_start_length=after_region_start_length,
            bin_size=bin_size,
            skip_zeros=skip_zeros,
            processes=processes,
            figsize=figsize,
            dpi=dpi,
            ylabel=ylabel,
            title=title,
            legend_location=legend_location,
            show=show,
        )

    @staticmethod
    def volcano_with_intensity_alpha(
        df: pd.DataFrame,
        padj_col: str,
        lfc_col: str,
        gene_id_col: str,
        gene_name_col: str | None = None,
        intensity_col: str | None = None,
        log2FoldChange_upper: float = 1.5,
        log2FoldChange_lower: float = -1.5,
        nlog10_upper: float = 2.0,
        text_nlog10: float | None = None,
        text_log2FoldChange: float | None = None,
        title: str | None = None,
        highlight_genes: list[str] | None = None,
        fontsize: int = 12,
        fact_col: str | None = None,
        categories: list[str] | None = None,
    ):
        """
        Volcano plot with optional intensity-based alpha transparency.

        When `fact_col` is `None` or `categories` is empty, delegates to
        `Pyntegrate.plotutils.volcano_plot` once for the whole `df`. Otherwise,
        filters `df` by each requested category and produces one figure per
        category.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `df` | `pd.DataFrame` | — | DataFrame containing expression results. |
        | `padj_col` | `str` | — | Column name for adjusted p-values. |
        | `lfc_col` | `str` | — | Column name for log2 fold-change values. |
        | `gene_id_col` | `str` | — | Column name for gene identifiers. |
        | `gene_name_col` | `str` | `None` | Column name for gene names/labels. |
        | `intensity_col` | `str` | `None` | Column name for intensity values used to scale point alpha. |
        | `log2FoldChange_upper` | `float` | `1.5` | Upper log2FC threshold for the up-regulation guide line. |
        | `log2FoldChange_lower` | `float` | `-1.5` | Lower log2FC threshold for the down-regulation guide line. |
        | `nlog10_upper` | `float` | `2.0` | `-log10(padj)` threshold for the significance guide line. |
        | `text_nlog10` | `float` | `None` | `-log10(padj)` threshold above which genes are auto-labeled. |
        | `text_log2FoldChange` | `float` | `None` | log2FC magnitude threshold above which genes are auto-labeled. |
        | `title` | `str` | `None` | Plot title. Ignored per-category (a generated title is used instead) when `fact_col`/`categories` are set. |
        | `highlight_genes` | `list[str]` | `None` | Gene names to always label, regardless of significance thresholds. |
        | `fontsize` | `int` | `12` | Font size for labels and title. |
        | `fact_col` | `str` | `None` | Column in `df` used to split genes into categories. If `None`, a single figure is produced. |
        | `categories` | `list[str]` | `None` | Category values (matched case-insensitively against `fact_col`) to plot; one figure per category. |

        ## Returns

        `matplotlib.figure.Figure | dict[str, matplotlib.figure.Figure]` —
        A single figure when `fact_col` is `None` or `categories` is empty;
        a dictionary keyed by lower-cased category name otherwise.
        """
        # If no fact_col OR categories empty → single plot (ignore categories)
        if fact_col is None or not categories:
            return volcano_plot(
                df=df,
                padj_col=padj_col,
                lfc_col=lfc_col,
                gene_id_col=gene_id_col,
                gene_name_col=gene_name_col,
                intensity_col=intensity_col,
                log2FoldChange_upper=log2FoldChange_upper,
                log2FoldChange_lower=log2FoldChange_lower,
                nlog10_upper=nlog10_upper,
                text_nlog10=text_nlog10,
                text_log2FoldChange=text_log2FoldChange,
                highlight_genes=highlight_genes,
                fontsize=fontsize,
                title=title,
            )

        # Get all categories and filter if specified
        all_categories = df[fact_col].unique()

        # If categories specified, filter to only those
        if categories:
            categories_to_plot = [c for c in all_categories if str(c).lower() in [s.lower() for s in categories]]
        else:
            categories_to_plot = all_categories

        # Generate plots for each category
        figures = {}

        for category in categories_to_plot:
            df_cat = df[df[fact_col] == category].reset_index(drop=True)
            figures[str(category).lower()] = volcano_plot(
                df=df_cat,
                padj_col=padj_col,
                lfc_col=lfc_col,
                gene_id_col=gene_id_col,
                gene_name_col=gene_name_col,
                intensity_col=intensity_col,
                log2FoldChange_upper=log2FoldChange_upper,
                log2FoldChange_lower=log2FoldChange_lower,
                nlog10_upper=nlog10_upper,
                text_nlog10=text_nlog10,
                text_log2FoldChange=text_log2FoldChange,
                title=f"{padj_col.split()[0].upper()} {str(category).upper()}",
                highlight_genes=highlight_genes,
                fontsize=fontsize,
            )

        return figures

    @staticmethod
    def circplot(
        csv_path: str,
        gene_list_input: List[str],
        chip_col: str,
        rna_col: str,
        atac_col: str,
        chip_basemean_col: Optional[str] = None,
        rna_basemean_col: Optional[str] = None,
        atac_basemean_col: Optional[str] = None,
        save_tsv: bool = False,
        tsv_output_path: Optional[str] = None,
        threshold: float = 1.0,
        allowed_combinations: Optional[List[Tuple[str, str, str]]] = None,
        identifier_col: Optional[str] = None,
        value_transform: Optional[Callable[[float], float]] = None,
        value_label: str = "log2 Fold Change",
    ):
        """
        Generate a circos plot comparing ChIP, ATAC, and RNA log2 fold-change across selected genes.

        The plot has three outer sector groups (ChIP-seq, ATAC-seq, RNA-seq) with one sector per
        gene per group. Inner links connect genes whose signals match user-specified combinations.
        A coolwarm colormap encodes log2FC magnitude.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `csv_path` | `str` | — | Path to CSV/TSV file with gene expression data. |
        | `gene_list_input` | `list[str]` | — | Gene IDs (Ensembl) or symbols to include; minimum 2 required. |
        | `chip_col` | `str` | — | Column name for ChIP log2 fold-change values. |
        | `rna_col` | `str` | — | Column name for RNA log2 fold-change values. |
        | `atac_col` | `str` | — | Column name for ATAC log2 fold-change values. |
        | `chip_basemean_col` | `str` | `None` | Optional column name for ChIP baseMean values exported to TSV. |
        | `rna_basemean_col` | `str` | `None` | Optional column name for RNA baseMean values exported to TSV. |
        | `atac_basemean_col` | `str` | `None` | Optional column name for ATAC baseMean values exported to TSV. |
        | `identifier_col` | `str` | `None` | Column used to match `gene_list_input`. If `None`, the function falls back to automatic matching by Ensembl ID or symbol. |
        | `save_tsv` | `bool` | `False` | If `True`, save gene data to a TSV file. |
        | `tsv_output_path` | `str` | `None` | Path for output TSV file. If `None`, auto-generates name based on condition. |
        | `threshold` | `float` | `1.0` | Threshold value to classify signals as "high" (≥threshold) or "low" (≤-threshold). Given in the same units as the *raw* `chip_col`/`rna_col`/`atac_col` values (e.g. a padj cutoff of `0.05`) — when `value_transform` is set, `threshold` itself is passed through `value_transform` once (raising `ValueError` if that fails or isn't finite) so it stays comparable to the already-transformed values it's checked against. |
        | `allowed_combinations` | `list[tuple[str, str, str]]` | `None` | List of allowed (chip, rna, atac) combinations where each is "high" or "low". If `None`, all 8 possible combinations are allowed (no restriction). |
        | `value_transform` | `Callable[[float], float]` | `None` | Optional function applied to every `chip_col`/`rna_col`/`atac_col` value (and to `threshold`, see above) right after reading it from the CSV (e.g. `lambda v: -np.log10(v)`), before it is used for the rings, connections, and `threshold` comparisons. Does **not** change the `ChIP_{value_label}`/`ATAC_{value_label}`/`RNA_{value_label}` numbers written to the TSV export (`save_tsv=True`) — those stay the raw, as-read values — but the TSV's `ChangeChIP`/`ChangeATAC`/`ChangeRNA` columns are still derived from the transformed values and `threshold`, matching the classification the plot itself used. Does not affect `chip_basemean_col`/`rna_basemean_col`/`atac_basemean_col`. If `None` (default), values are used as read — identical to the previous behavior. A gene/signal whose raw value is missing, non-numeric, or for which `value_transform` raises or returns NaN (e.g. `log10` of 0 or a negative number) is treated as missing (`NaN`), same as any other missing value; `threshold` failing the same way raises instead, since that would silently break every classification. |
        | `value_label` | `str` | `"log2 Fold Change"` | Text describing what's plotted, shown in the figure title as `f"{condition} vs Control - {value_label}"`. Defaults to `"log2 Fold Change"` — identical to the previous hardcoded title. Set this to match `value_transform` when plotting something else (e.g. `"-log10(padj)"` when `chip_col`/`rna_col`/`atac_col` are padj columns transformed with `-log10`). |

        ## Returns

        `plt.Figure` — Circos plot figure when `save_tsv=False`.
        `tuple[plt.Figure, str]` — `(circos_figure, tsv_file_path)` when `save_tsv=True`.

        ## Raises

        - `Exception` — If any gene is not found in the CSV or fewer than two genes are provided.

        ## Examples

        ```python
        genes = ["ENSG00000162494", "ENSG00000186407", "TP53"]
        fig = Plotting.circplot(
            "data.csv", genes,
            "ssrp1_chip_lfc", "ssrp1_rna_lfc", "ssrp1_atac_lfc",
        )
        fig.savefig("circplot_ssrp1.png", dpi=300, bbox_inches="tight")

        fig, tsv_path = Plotting.circplot(
            "data.csv", genes,
            "ctrl_chip_lfc", "ctrl_rna_lfc", "ctrl_atac_lfc",
            save_tsv=True,
        )
        print(f"Data saved to: {tsv_path}")

        # Treat chip_col/rna_col/atac_col as padj columns and plot -log10(padj)
        # instead of the raw values (missing/zero/negative -> NaN, dropped as usual).
        # threshold=0.05 stays in padj units — it's auto-converted to
        # -log10(0.05) ≈ 1.301 internally, so "high" means padj <= 0.05.
        fig = Plotting.circplot(
            "data.csv", genes,
            "ssrp1_chip_padj", "ssrp1_rna_padj", "ssrp1_atac_padj",
            threshold=0.05,
            value_transform=lambda v: -np.log10(v) if v > 0 else float("nan"),
            value_label="-log10(padj)",
        )
        ```
        """
        from collections import OrderedDict
        from pycirclize import Circos
        from matplotlib import cm, colors
        import csv, itertools

        def find_duplicates(items: List[str]) -> List[str]:
            """Find duplicated items."""
            counter = {}
            for item in items:
                if item in counter:
                    counter[item] += 1
                else:
                    counter[item] = 1
            return [item for item, count in counter.items() if count > 1]

        def extract_data(gene_list, chip_name, rna_name, atac_name, file_path, identifier_col=None):
            """Extract data for genes."""
            condition_data, condition_data_to_save, genes, genes_to_save = {}, {}, [], []

            # Use synthetic, guaranteed-unique keys when a *_basemean_col is
            # None (the default), so the three don't all collide on the same
            # `None` dict key and silently drop entries.
            chip_bm_key = chip_basemean_col if chip_basemean_col is not None else "__chip_basemean__"
            rna_bm_key = rna_basemean_col if rna_basemean_col is not None else "__rna_basemean__"
            atac_bm_key = atac_basemean_col if atac_basemean_col is not None else "__atac_basemean__"

            def extract_value(raw_value):
                """Apply `value_transform` (if given) to a raw CSV cell value.

                With no `value_transform`, returns `raw_value` unchanged — identity,
                preserves prior behavior exactly. With a `value_transform`, coerces
                `raw_value` to float first, then calls `value_transform` on it. Any
                failure (missing/non-numeric raw value, or `value_transform` raising
                or returning NaN — e.g. log10 of 0 or a negative number) becomes
                `None`, same as any other missing value elsewhere in this function.
                """
                if value_transform is None:
                    return raw_value
                numeric = pd.to_numeric(raw_value, errors="coerce")
                if pd.isna(numeric):
                    return None
                try:
                    transformed = value_transform(float(numeric))
                except Exception:
                    return None
                if transformed is None or (isinstance(transformed, (int, float)) and not np.isfinite(transformed)):
                    return None
                return transformed

            with open(file_path, "r", encoding="utf-8") as file:
                sample = file.read(4096)
                file.seek(0)
                try:
                    dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
                    delimiter = dialect.delimiter
                except csv.Error:
                    delimiter = "\t" if file_path.lower().endswith(".tsv") else ","

                reader = csv.DictReader(file, delimiter=delimiter)
                for row in reader:
                    identifier_value = (row.get(identifier_col) if identifier_col else (row.get("EnsemblGeneID") or row.get("GeneID") or row.get("ensembl gene ID.4") or row.get("gene name") or row.get("Symbol")))
                    chip_value = extract_value(row.get(chip_name))
                    rna_value = extract_value(row.get(rna_name))
                    atac_value = extract_value(row.get(atac_name))
                    if str(identifier_value).strip() in [str(g).strip() for g in gene_list]:
                        condition_data[identifier_value] = {chip_name: chip_value, rna_name: rna_value, atac_name: atac_value}
                        genes.append(identifier_value)
                        print(f"Added: {identifier_value}")
                    if(chip_value not in (None, "") and rna_value not in (None, "") and atac_value not in (None, "") and
                        (((abs(pd.to_numeric(chip_value, errors="coerce"))> threshold) and
                        (abs(pd.to_numeric(rna_value, errors="coerce")) > threshold)) or
                        ((abs(pd.to_numeric(atac_value, errors="coerce")) > threshold) and
                        (abs(pd.to_numeric(rna_value, errors="coerce")) > threshold)) or
                        ((abs(pd.to_numeric(chip_value, errors="coerce")) > threshold) and
                        (abs(pd.to_numeric(atac_value, errors="coerce")) > threshold))) or
                        str(identifier_value).strip() in [str(g).strip() for g in gene_list]):
                        # "__*_raw__" carry the untransformed CSV values (identical to
                        # row.get(chip_name)/row.get(rna_name)/row.get(atac_name)), kept
                        # alongside the *_transform'd chip_name/rna_name/atac_name entries
                        # so the TSV export can show raw numbers while rings/connections
                        # keep using the transformed ones. See `chip_values`/`rna_values`/
                        # `atac_values` below.
                        condition_data_to_save[identifier_value] = {chip_name: chip_value, chip_bm_key: row.get(chip_basemean_col), rna_name: rna_value, rna_bm_key: row.get(rna_basemean_col), atac_name: atac_value, atac_bm_key: row.get(atac_basemean_col), "__chip_raw__": row.get(chip_name), "__rna_raw__": row.get(rna_name), "__atac_raw__": row.get(atac_name)}
                        genes_to_save.append(identifier_value)

            return genes, genes_to_save, condition_data, condition_data_to_save
           

        def extract_extremes(rna_limits, atac_limits, chip_limits):
            """Extract max/min."""
            maximum = max([chip_limits[1], rna_limits[1], atac_limits[1]])
            minimum = min([chip_limits[0], rna_limits[0], atac_limits[0]])
            return maximum, minimum

        def implement_data(circos, num_genes, rna_values, atac_values, chip_values, cmap, norm, ax=None):
            """Implement data tracks, labels, legend, and colorbar."""
            category_maps = {
                "RNA-seq": rna_values,
                "ATAC-seq": atac_values,
                "ChIP-seq": chip_values,
            }

            show_legend = num_genes > 5
            gene_labels = {}
            label_counter = 1

            for sector in circos.sectors:
                category, gene = sector.name.split(":")
                values = category_maps[category]
                value = values[gene]
                color = cmap(norm(value))

                track = sector.add_track((75, 85), r_pad_ratio=0.05)
                track.axis()
                track.rect(0, 1, facecolor=color, edgecolor="black")

                if show_legend:
                    if gene not in gene_labels:
                        gene_labels[gene] = label_counter
                        label_counter += 1
                    track.text(
                        str(gene_labels[gene]),
                        0.5,
                        size=10,
                        ha="center",
                        va="center",
                        rotation="auto",
                        fontweight="bold",
                    )
                else:
                    track.text(gene, 0.5, size=10, ha="center", va="center", rotation="auto", fontweight="bold",)

            if show_legend and ax is not None:
                legend_elements = [
                    plt.Line2D(
                        [0],
                        [0],
                        marker="o",
                        color="w",
                        label=f"{num}: {gene}",
                        markerfacecolor="gray",
                        markersize=8,
                    )
                    for gene, num in gene_labels.items()
                ]
                legend = ax.legend(
                    handles=legend_elements,
                    title="Genes",
                    loc="upper left",
                    bbox_to_anchor=(-0.25, 1),
                    frameon=True,
                    fontsize=9,
                    title_fontsize=11,
                    labelspacing=0.6,
                    borderpad=0.8,
                )

                legend.get_frame().set_facecolor("whitesmoke")
                legend.get_frame().set_edgecolor("gray")
                legend.get_frame().set_linewidth(0.8)
                legend.get_frame().set_alpha(0.9)
                legend.get_frame().set_boxstyle("round,pad=0.4,rounding_size=2")

            if ax is not None:
                sm = cm.ScalarMappable(cmap=cmap, norm=norm)
                sm.set_array([])
                cbar = plt.colorbar(sm, ax=ax, pad=0.1, fraction=0.04)
                vmin, vmax = norm.vmin, norm.vmax
                sep = (vmax/10)
                mid_ticks = np.arange(np.ceil(vmin / sep) * sep, vmax, sep)
                ticks = np.unique(np.concatenate(([vmin], mid_ticks, [vmax])))
                if len(ticks) >= 2 and (ticks[-1] - ticks[-2]) < sep*0.6:
                    ticks = np.delete(ticks, -2)
                if len(ticks) >= 2 and (ticks[1] - ticks[0]) < sep*0.6:
                    ticks = np.delete(ticks, 1)
                cbar.set_ticks(ticks)
                cbar.set_ticklabels([f"{t:.2f}" for t in ticks])

        def implement_connections(circos, data, cmap, norm, threshold, allowed_combinations):
            """
            Implement connections based on allowed combinations and threshold.

            Logic:
                1. If full combination is allowed (e.g., "high"-"low"-"low"), draw a TRIANGLE
                    connecting all 3 signals with weight 1.5 and color = mean(chip, rna, atac).
                2. Pair connections are evaluated independently from triangle.
                    A pair is linked when it matches any allowed combination for that pair,
                    even if the third signal does not match.

                Example: allowed ("low", "high", "high") and states ("mid", "high", "high"):
                        - RNA-ATAC connects (high-high match)
                        - no triangle (full combination not matched)

            Parameters:
            -----------
            threshold : float
                Threshold for classifying high (≥threshold) or low (≤-threshold)
            allowed_combinations : list of tuple
                List of allowed (chip, rna, atac) combinations
            """
            connections = []
            connected_genes = set()

            for gene, values in data.items():
                chip, rna, atac = (_safe_float(v) for v in values.values())
                def classify(v):
                    if v >= threshold:
                        return "high"
                    elif v <= -threshold:
                        return "low"
                    else:
                        return "mid"

                chip_c = classify(chip)
                rna_c = classify(rna)
                atac_c = classify(atac)
                
                current_combination = (chip_c, rna_c, atac_c)
                
                # ============================================================
                # STEP 1: Draw TRIANGLE (full combination)
                # ============================================================
                if current_combination in allowed_combinations:
                    if (chip > 0 and rna > 0 and atac > 0) or (chip < 0 and rna < 0 and atac < 0):
                        avg_val = (chip + rna + atac) / 3
                        color_triangle = cmap(norm(avg_val))
                    else:
                        color_triangle = "gray"
                    weight_triangle = 1.5
                    
                    triangle_pairs = [
                        ("ChIP-seq", "RNA-seq"),
                        ("ChIP-seq", "ATAC-seq"),
                        ("RNA-seq", "ATAC-seq"),
                    ]
                    
                    for s1, s2 in triangle_pairs:
                        connections.append([
                            f"{s1}:{gene}",
                            f"{s2}:{gene}",
                            weight_triangle,
                            color_triangle
                        ])

                    # Mark this gene as connected (triangle drawn)
                    connected_genes.add(gene)

                    # Important: if triangle is drawn, do not overlay pair links
                    # on the same edges to preserve one uniform color for all sides.
                    continue
                
                # ============================================================
                # STEP 2: Additional pairwise connections
                # (partial match against allowed_combinations)
                # ============================================================

                # Pair 1: ChIP-RNA -> match any (chip_state, rna_state, *) in allowed
                match_chip_rna = any(
                    chip_c == ac_chip and rna_c == ac_rna
                    for ac_chip, ac_rna, _ in allowed_combinations
                )
                if match_chip_rna:
                    if (chip > 0 and rna > 0) or (chip < 0 and rna < 0):
                        avg_pair = (chip + rna) / 2
                        color_pair = cmap(norm(avg_pair))
                    else:
                        color_pair = "gray"
                    connections.append([
                        "ChIP-seq:" + gene,
                        "RNA-seq:" + gene,
                        1.0,
                        color_pair
                    ])
                    connected_genes.add(gene)

                # Pair 2: ChIP-ATAC -> match any (chip_state, *, atac_state) in allowed
                match_chip_atac = any(
                    chip_c == ac_chip and atac_c == ac_atac
                    for ac_chip, _, ac_atac in allowed_combinations
                )
                if match_chip_atac:
                    if (chip > 0 and atac > 0) or (chip < 0 and atac < 0):
                        avg_pair = (chip + atac) / 2
                        color_pair = cmap(norm(avg_pair))
                    else:
                        color_pair = "gray"
                    connections.append([
                        "ChIP-seq:" + gene,
                        "ATAC-seq:" + gene,
                        1.0,
                        color_pair
                    ])
                    connected_genes.add(gene)
                
                # Pair 3: RNA-ATAC -> match any (*, rna_state, atac_state) in allowed
                match_rna_atac = any(
                    rna_c == ac_rna and atac_c == ac_atac
                    for _, ac_rna, ac_atac in allowed_combinations
                )
                if match_rna_atac:
                    if (rna > 0 and atac > 0) or (rna < 0 and atac < 0):
                        avg_pair = (rna + atac) / 2
                        color_pair = cmap(norm(avg_pair))
                    else:
                        color_pair = "gray"
                    connections.append([
                        "RNA-seq:" + gene,
                        "ATAC-seq:" + gene,
                        1.0,
                        color_pair
                    ])
                    connected_genes.add(gene)

            # --- Draw all connections ---
            for source, target, weight, color in connections:
                circos.link(
                    (source, 0.4, 0.6),
                    (target, 0.4, 0.6),
                    color=color,
                    direction=0,
                    lw=weight,
                    alpha=1
                )

            # Return the list of genes that have at least one connection
            return list(connected_genes)

        def get_sector_by_name(circos_obj, name):
            """Get sector by name."""
            for sector in circos_obj.sectors:
                if sector.name == name:
                    return sector
            raise KeyError(f"Sector '{name}' not found.")

        def add_category_labels(circos, ax, sector_groups, radius=100):
            """Add category labels around circos."""
            for category, group in sector_groups.items():
                first_sector = get_sector_by_name(circos, group[0])
                last_sector = get_sector_by_name(circos, group[-1])

                min_deg_lim, max_deg_lim = min(first_sector.deg_lim), max(last_sector.deg_lim)
                mid_angle = (min_deg_lim + max_deg_lim) / 2
                mid_angle_rad = np.radians(mid_angle)

                ax.text(
                    mid_angle_rad,
                    radius,
                    category,
                    ha="center",
                    va="center",
                    fontsize=16,
                    fontweight="bold",
                    color="black",
                    rotation=(360 - mid_angle) if (mid_angle < 90 or mid_angle > 270) else 180 - mid_angle,
                    rotation_mode="anchor",
                )

        def save_table_tsv(
            genes,
            chip_values,
            atac_values,
            rna_values,
            chip_basemean_values=None,
            atac_basemean_values=None,
            rna_basemean_values=None,
            chip_values_classify=None,
            atac_values_classify=None,
            rna_values_classify=None,
            output_path=None,
        ):
            """Save gene values to a TSV file.

            `chip_values`/`atac_values`/`rna_values` are the raw (pre-`value_transform`)
            numbers written to the displayed `ChIP_{value_label}`/`ATAC_{value_label}`/
            `RNA_{value_label}` columns. `chip_values_classify`/`atac_values_classify`/
            `rna_values_classify` (defaulting to the display values themselves when not
            given — i.e. a no-op unless `value_transform` is set) are what
            `ChangeChIP`/`ChangeATAC`/`ChangeRNA` are computed against, together with
            `threshold` — the already-transformed values, matching exactly what
            classified genes as high/low/connected in the circos plot itself.
            """
            chip_values_classify = chip_values_classify if chip_values_classify is not None else chip_values
            atac_values_classify = atac_values_classify if atac_values_classify is not None else atac_values
            rna_values_classify = rna_values_classify if rna_values_classify is not None else rna_values

            # Create DataFrame with gene data (including index)
            data_for_tsv = {
                'Index': list(range(1, len(genes) + 1)),
                'Gene': genes,
                f'ChIP_{value_label}': [float(chip_values_classify[g]) for g in genes],
                'ChIP_baseMean': [float(chip_basemean_values.get(g, np.nan)) for g in genes],
                f'ATAC_{value_label}': [float(atac_values_classify[g]) for g in genes],
                'ATAC_baseMean': [float(atac_basemean_values.get(g, np.nan)) for g in genes],
                f'RNA_{value_label}': [float(rna_values_classify[g]) for g in genes],
                'RNA_baseMean': [float(rna_basemean_values.get(g, np.nan)) for g in genes],
                'ChangeChIP': ['+' if float(chip_values_classify[g]) >= threshold else '-' if float(chip_values_classify[g]) <= -threshold else 'o' for g in genes],
                'ChangeATAC': ['+' if float(atac_values_classify[g]) >= threshold else '-' if float(atac_values_classify[g]) <= -threshold else 'o' for g in genes],
                'ChangeRNA': ['+' if float(rna_values_classify[g]) >= threshold else '-' if float(rna_values_classify[g]) <= -threshold else 'o' for g in genes],
            }


            df = pd.DataFrame(data_for_tsv)

            # Auto-generate filename if not provided
            if output_path is None:
                output_path = f"circplot_data_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.tsv"
            
            # Save to TSV with 3 decimal places for numeric columns
            df.to_csv(
                output_path,
                sep='\t',
                index=False,
                float_format=lambda x: f"{x:.3f}" if pd.notna(x) else "",
            )
            
            return output_path

        
        # If a `value_transform` is given, transform `threshold` the same way
        # so it stays comparable to the already-transformed chip/rna/atac
        # values used everywhere below (row filtering in extract_data,
        # connection classification, and ChangeChIP/ChangeATAC/ChangeRNA in
        # the TSV export). Pass `threshold` in the same units as the raw
        # column values (e.g. a padj cutoff of 0.05) — this converts it into
        # `value_transform`'s output units for you.
        raw_threshold = threshold
        if value_transform is not None:
            _threshold_numeric = pd.to_numeric(threshold, errors="coerce")
            if pd.isna(_threshold_numeric):
                raise ValueError(f"threshold must be numeric, got {threshold!r}.")
            try:
                threshold = float(value_transform(float(_threshold_numeric)))
            except Exception as exc:
                raise ValueError(
                    f"value_transform(threshold) failed for threshold={raw_threshold!r}: {exc}"
                ) from exc
            if not np.isfinite(threshold):
                raise ValueError(
                    f"value_transform({raw_threshold!r}) did not produce a finite number "
                    f"(got {threshold!r}); pick a threshold within value_transform's domain."
                )

        # Validate allowed combinations
        valid_states = ["high", "low"]
        if allowed_combinations is None:
            allowed_combinations = list(itertools.product(valid_states, repeat=3))
        for combo in allowed_combinations:
            if len(combo) != 3 or not all(state in valid_states for state in combo):
                raise ValueError(
                    f"Invalid combination: {combo}. Each combination must be a 3-tuple of 'high' or 'low'"
                )
        # Validate
        dup = find_duplicates(gene_list_input)
        if dup:
            raise Exception(f"Duplicated genes: {dup}")
        if len(gene_list_input) <= 1:
            raise Exception("Insert at least two genes.")

        genes, genes_to_save, condition_data, condition_data_to_save = extract_data(gene_list_input, chip_col, rna_col, atac_col, csv_path, identifier_col=identifier_col)

        num_genes = len(genes)
        categories = ["ATAC-seq", "ChIP-seq", "RNA-seq"]

        sector_dict = OrderedDict()
        for cat in categories:
            for g in genes:
                sector_dict[f"{cat}:{g}"] = 1

        spaces = []
        for i, _ in enumerate(categories):
            for _ in genes:
                spaces.append(0)
            if i <= len(categories) - 1:
                spaces[-1] = 5

        sector_groups = {cat: [f"{cat}:{gene}" for gene in genes] for cat in categories}
        cmap = plt.colormaps["seismic"]

        def _safe_float(value):
            """Convert to float, treating missing/empty basemean values as NaN."""
            return float(value) if value not in (None, "") else float("nan")

        # Raw (pre-`value_transform`) values — these are what gets written to the
        # TSV's ChIP_log2FC/ATAC_log2FC/RNA_log2FC columns; identical to
        # `chip_values`/`rna_values`/`atac_values` when `value_transform` is None.
        chip_values = {k: _safe_float(v.get("__chip_raw__")) for k, v in condition_data_to_save.items()}
        rna_values = {k: _safe_float(v.get("__rna_raw__")) for k, v in condition_data_to_save.items()}
        atac_values = {k: _safe_float(v.get("__atac_raw__")) for k, v in condition_data_to_save.items()}
        # Transformed values — same values used for the rings/connections, kept
        # here so the TSV's ChangeChIP/ChangeATAC/ChangeRNA can still classify
        # against exactly what the plot itself classified against.
        chip_values_transformed = {k: _safe_float(list(v.values())[0]) for k, v in condition_data_to_save.items()}
        rna_values_transformed = {k: _safe_float(list(v.values())[2])  for k, v in condition_data_to_save.items()}
        atac_values_transformed = {k: _safe_float(list(v.values())[4]) for k, v in condition_data_to_save.items()}
        chip_basemean_values = {k: _safe_float(list(v.values())[1]) for k, v in condition_data_to_save.items()}
        rna_basemean_values = {k: _safe_float(list(v.values())[3])  for k, v in condition_data_to_save.items()}
        atac_basemean_values = {k: _safe_float(list(v.values())[5]) for k, v in condition_data_to_save.items()}

        chip_values_selected = {k: _safe_float(list(v.values())[0]) for k, v in condition_data_to_save.items() if k in genes}
        rna_values_selected = {k: _safe_float(list(v.values())[2])  for k, v in condition_data_to_save.items() if k in genes}
        atac_values_selected = {k: _safe_float(list(v.values())[4]) for k, v in condition_data_to_save.items() if k in genes}

        chip_limits = [np.nanmin(list(chip_values_selected.values())), np.nanmax(list(chip_values_selected.values()))]
        rna_limits = [np.nanmin(list(rna_values_selected.values())), np.nanmax(list(rna_values_selected.values()))]
        atac_limits = [np.nanmin(list(atac_values_selected.values())), np.nanmax(list(atac_values_selected.values()))]

        circos = Circos(sector_dict, space=spaces)
        fig, ax = plt.subplots(figsize=(12, 10), subplot_kw=dict(polar=True))

        max_val, min_val = extract_extremes(rna_limits, atac_limits, chip_limits)
        # Single symmetric scale around 0: color intensity is proportional to
        # the actual value magnitude, on the same scale for both signs, so
        # e.g. -0.50 and +0.50 get identical color intensity. 0 stays white.
        # Whichever sign has the larger absolute extent reaches full
        # saturation (blue or red); the other sign only reaches it if its
        # magnitude matches. If one sign is absent entirely, its color never
        # appears since no value crosses the midpoint on that side.
        
        vabs = max(abs(min_val), abs(max_val), 1e-9)
        import matplotlib.colors as mcolors
        print(min_val, max_val, vabs)
        print(1-((vabs-min_val)/(vabs*2)), 1, 256) if vabs == abs(max_val) else print(0, 1-(vabs-abs(max_val))/(vabs*2), 256)
        if vabs == abs(max_val):
            if min_val >= 0:
                min_scale = 0.5
                min_v = 0
            else:
                min_scale = 1-((vabs-min_val)/(vabs*2))
                min_v = min_val
            colors = cmap(np.linspace(min_scale, 1, 256))
        else:
            min_v = min_val
            colors = cmap(np.linspace(0, 1-(vabs-abs(max_val))/(vabs*2), 256))

        cut_cmap = mcolors.LinearSegmentedColormap.from_list('seismic', colors)

        norm = mcolors.Normalize(vmin=min_v, vmax=max_val)

        implement_data(
            circos,
            num_genes,
            rna_values_selected,
            atac_values_selected,
            chip_values_selected,
            cut_cmap,
            norm,
            ax=ax,
        )
        connected_genes = implement_connections(circos, condition_data, cut_cmap, norm, threshold, allowed_combinations)

        circos.plotfig(ax=ax)
        condition_label = atac_col.split()[0].upper() if atac_col else "Condition"
        title = f"{condition_label} vs Control ({value_label})"
        subtitle = (
            f"Threshold = {raw_threshold} (transformed: {threshold:.3g})"
            if value_transform is not None else f"Threshold = {threshold}"
        )
        ax.set_title(f"{title}\n{subtitle}", y=0.93)

        add_category_labels(circos, ax, sector_groups)

        # Save data to TSV if requested — include genes that match the threshold rule
        if save_tsv:
            tsv_path = save_table_tsv(
                genes_to_save, chip_values, atac_values, rna_values,
                chip_basemean_values=chip_basemean_values,
                atac_basemean_values=atac_basemean_values,
                rna_basemean_values=rna_basemean_values,
                chip_values_classify=chip_values_transformed,
                atac_values_classify=atac_values_transformed,
                rna_values_classify=rna_values_transformed,
                output_path=tsv_output_path
            )
            print(f"[Pyntegrate] Data exported to TSV: {tsv_path}")
            return fig, tsv_path
        else:
            return fig

class Homer:
   """
   HOMER CLI wrappers for peak calling and annotation.

   Requires HOMER installed with `makeTagDirectory`, `findPeaks`, and `annotatePeaks.pl` on PATH.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `create_homer_tag_directory()` | Creates a HOMER tagDirectory from a signal file. |
   | `homer_peaks()` | Calls `findPeaks` to identify peaks from a tagDirectory. |
   | `homer_annotate_peaks()` | Annotates peaks with genomic features and nearest genes. |
   """


   @staticmethod
   def create_homer_tag_directory(name: str, file: str) -> None:
       """
       Create a HOMER tagDirectory from a signal file.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `name` | `str` | — | Name/path of the tagDirectory to create. |
       | `file` | `str` | — | Path to the input signal file. |

       ## Returns

       `None`

       ## Examples

       ```python
       Homer.create_homer_tag_directory("myTagDir", "chip.bam")
       ```
       """
       subprocess.run(["makeTagDirectory", name, file])
       return


   @staticmethod
   def homer_peaks(tag_directory: str, style: str = "factor", output: str = "auto") -> None:
       """
       Call HOMER `findPeaks` to identify enriched peaks from a tagDirectory.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `tag_directory` | `str` | — | Path to HOMER tagDirectory. |
       | `style` | `str` | `"factor"` | Peak calling style (e.g., `"factor"`, `"histone"`). |
       | `output` | `str` | `"auto"` | Output filename passed to HOMER via `-o`. |

       ## Returns

       `None`

       ## Examples

       ```python
       Homer.homer_peaks("myTagDir", style="factor", output="peaks.txt")
       ```
       """
       subprocess.run(["findPeaks", tag_directory, "-style", style, "-o", output])
       return


   @staticmethod
   def homer_annotate_peaks(file_directory: str, gene: str, output: str) -> None:
       """
       Annotate peaks with HOMER `annotatePeaks.pl`.

       Writes the annotated peaks table to `output` and HOMER stderr to `logs.log`.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `file_directory` | `str` | — | Path to peaks file (BED or narrowPeak). |
       | `gene` | `str` | — | Genome identifier installed in HOMER (e.g., `"mm10"`, `"hg19"`). |
       | `output` | `str` | — | Output filename for the annotated peaks table. |

       ## Returns

       `None`

       ## Examples

       ```python
       Homer.homer_annotate_peaks("peaks.bed", "mm10", "peaks.annot.txt")
       ```
       """
       with open(output, "w") as out_file, open("logs.log", "w") as logs_file:
           subprocess.run(["annotatePeaks.pl", file_directory, gene], stdout=out_file, stderr=logs_file)

class Filtering:
   """
   Helpers for intersecting and aligning gene sets across modalities (ATAC/ChIP/RNA).

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `filter_shared()` | Returns entries shared across two or more signal arrays, matched by a common key. |
   | `chip_or_atac_genes_not_used_with_rna()` | Keeps signal entries whose gene IDs exist in an RNA results object. |
   """


   @staticmethod
   def filter_shared(by, *signals):

       """
       Return only entries shared across two or more signal arrays, matched by a common key.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `by` | `str` | — | Key used to match entries (e.g., `"gene_name"`, `"id"`). |
       | `*signals` | `list[dict]` | — | Two or more signal arrays of dict entries; each entry must contain the `by` key. |

       ## Returns

       `list[list[dict]]` — One filtered list per input signal, in the same order as `*signals`, aligned to the shared keys.

       ## Examples

       ```python
       atac_good, chip_good = Filtering.filter_shared("gene_name", atac, chip)
       assert len(atac_good) == len(chip_good)
       ```
       """
       dicts = [{entry[by]: entry for entry in signal} for signal in signals]

       common_keys = set.intersection(*(set(d) for d in dicts))

       return [
            [d[key] for key in common_keys]
            for d in dicts
        ]


   @staticmethod
   def chip_or_atac_genes_not_used_with_rna(
       data: Any,
       signal: List[Dict[str, Any]],
       id: str = "gene_name",
       rna_column_name: str = "gene_name",
       delete_no_peaks: bool = False,
   ):
       """
       Keep only signal entries whose gene IDs exist in an RNA results object.

       The RNA object is never mutated; it is returned as-is in the output tuple.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `data` | `Any` | — | ResultsTable-like object exposing `.data` as a pandas DataFrame. |
       | `signal` | `list[dict]` | — | Signal entries (ATAC or ChIP) as dicts containing the key named by `id`. |
       | `id` | `str` | `"gene_name"` | Signal dict key used as gene identifier. When `"id"`, version suffix is stripped at `'.'`. |
       | `rna_column_name` | `str` | `"gene_name"` | RNA DataFrame column to match against. |
       | `delete_no_peaks` | `bool` | `False` | If `True`, first remove entries where `mean(values) == 0`. |

       ## Returns

       `tuple[list[dict], Any]` — `(signal_good, data)` where `data` is the unchanged RNA object.

       ## Examples

       ```python
       signal_good, rna = Filtering.chip_or_atac_genes_not_used_with_rna(rna, chip_array)
       print(len(signal_good))
       ```

       **Output:**
       ```
       3842
       ```
       """
       # 1) keep only valid peaks if requested
       if delete_no_peaks:
           signal_semigood = [x for x in signal if np.mean(x["values"]) != 0]
       else:
           signal_semigood = signal


       # 2) intersect with RNA without altering RNA object
       ids_selected = []
       signal_good = []
       rna_series = data.data[rna_column_name]
       rna_series = rna_series.astype("string")
       rna_ids = set(rna_series.dropna().astype(str).values)


       for x in signal_semigood:
           gid = str(x[id])
           gid = gid.split(".")[0] if id == "id" else gid
           if gid in rna_ids and gid not in ids_selected:
               x[id] = gid
               signal_good.append(x)
               ids_selected.append(gid)

       # No reconstruction of DESeq, no dataframe mutation
       return signal_good, data

class BigWigUtils:
   """
   Helpers for BigWig format conversion.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `bigwigToBed()` | Converts a BigWig file to a BED-like TSV text file. |
   """


   @staticmethod
   def bigwigToBed(bwFile: str, bedNameFile: str) -> None:
       """
       Convert a BigWig file to a BED-like TSV file.

       Output columns: `["chrom", "start", "end", "ID", "value", ""]`.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `bwFile` | `str` | — | Path to the input BigWig file. |
       | `bedNameFile` | `str` | — | Output path for the BED/TSV file. |

       ## Returns

       `None`

       ## Examples

       ```python
       BigWigUtils.bigwigToBed("signal.bigWig", "signal.bed.tsv")
       ```
       """
       import pyBigWig


       bw = pyBigWig.open(bwFile)
       chroms = bw.chroms()


       bed_data = []
       i = 0
       for chrom in chroms:
           intervals = bw.intervals(chrom)
           if intervals is None:
               continue
           for start, end, score in intervals:
               bed_data.append([chrom, start, end, i, score, ""])
               i += 1


       bed_df = pd.DataFrame(bed_data, columns=["chrom", "start", "end", "ID", "value", ""])
       bed_df.to_csv(bedNameFile, sep="\t", header=True, index=False)

class AnnotationAnalysis:
   """
   Integrated analysis combining peak signal data with HOMER and FACT annotations.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `array_gene_name_and_annotation()` | Combines peaks with HOMER annotations and FACT categories for stratified heatmap generation. |
   """


   @staticmethod
   def array_gene_name_and_annotation(
       normalized_subtracted: List[Dict[str, Any]],
       xAxes: np.ndarray,
       fact_table_path: str,
       fact_col: str = "ssrp1 FACT category",
       homer_gene: str = "mm10",
   ) -> None:
       """
       Stratify peak signal by HOMER annotation class and FACT high/low category, then display heatmaps.

       Workflow: (1) build a BED from `normalized_subtracted` and run HOMER `annotatePeaks.pl`;
       (2) map each gene → annotation class; (3) load FACT table and classify each gene as
       high/low/other; (4) for each annotation category, generate three heatmaps (unsorted,
       sorted-by-mean, sorted-by-max-index) split by FACT category.

       Writes `normalized_subtracted.bed`, `normalized_subtracted.bed.txt`, and `logs.log` to
       the current working directory.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `normalized_subtracted` | `list[dict]` | — | List of dicts containing at least `'chr'`, `'start'`, `'end'`, `'gene_name'`, `'values'`. |
       | `xAxes` | `np.ndarray` | — | X-axis positions with shape `(B,)`. |
       | `fact_table_path` | `str` | — | Path to FACT table (CSV/TSV). |
       | `fact_col` | `str` | `"ssrp1 FACT category"` | Column in FACT table for gene classification (`"high"`/`"low"`/other). |
       | `homer_gene` | `str` | `"mm10"` | Genome identifier for HOMER (e.g., `"mm10"`, `"hg19"`). |

       ## Returns

       `None` — Displays figures via `plt.show()`.

       ## Examples

       ```python
       AnnotationAnalysis.array_gene_name_and_annotation(
           normalized_subtracted=subtracted,
           xAxes=x_axis,
           fact_table_path="final_table.tsv",
           fact_col="ssrp1 FACT category",
           homer_gene="mm10",
       )
       ```
       """
       # 1) Build BED for HOMER
       bed_data = []
       i = 0
       for e in normalized_subtracted:
           i += 1
           bed_data.append([e["chr"], e["start"], e["end"], i])


       df_bed = pd.DataFrame(bed_data, columns=["chr", "start", "end", "ID"])
       df_bed.to_csv("normalized_subtracted.bed", sep="\t", header=True, index=False)


       Homer.homer_annotate_peaks(
           file_directory="normalized_subtracted.bed",
           gene=homer_gene,
           output="normalized_subtracted.bed.txt",
       )
       df_homer = pd.read_csv("normalized_subtracted.bed.txt", delimiter="\t")


       # Normalize annotation field (remove parentheses content)
       df_homer["Annotation"] = (
           df_homer["Annotation"].astype(str).str.replace(r"\(.*", "", regex=True)
       )


       # Map gene_name -> annotation (keep first for duplicates)
       gene_to_annotation: Dict[str, str] = {}
       for _, row in df_homer.iterrows():
           gname = str(row.get("Gene Name", "")).strip()
           ann = str(row.get("Annotation", "")).strip()
           if gname and ann and gname not in gene_to_annotation:
               gene_to_annotation[gname] = ann


       unique_annotations = sorted(set(gene_to_annotation.values()))


       # 2) Load FACT table and map gene -> high/low/other
       fact_df = pd.read_csv(fact_table_path, low_memory=False)


       if "Symbol" in fact_df.columns:
           fact_df["Symbol"] = fact_df["Symbol"].astype(str).str.strip()
           fact_df[fact_col] = fact_df[fact_col].astype(str).str.strip().str.lower()
           fact_map = dict(zip(fact_df["Symbol"], fact_df[fact_col]))


           def gene_to_fact(g):
               cls = fact_map.get(str(g).strip(), "other")
               return cls if cls in ("high", "low") else "other"


       else:
           fact_df["EnsemblGeneID"] = fact_df["EnsemblGeneID"].astype(str)
           fact_df[fact_col] = fact_df[fact_col].astype(str).str.strip().str.lower()
           fact_map = dict(zip(fact_df["EnsemblGeneID"], fact_df[fact_col]))


           def gene_to_fact(g):
               cls = fact_map.get(str(g).strip(), "other")
               return cls if cls in ("high", "low") else "other"


       # 3) Group by annotation while keeping values and gene names aligned
       by_ann = {ann: {"values": [], "genes": []} for ann in unique_annotations}


       for entry in normalized_subtracted:
           g = entry.get("gene_name", None)
           if g is None:
               continue
           ann = gene_to_annotation.get(str(g).strip(), None)
           if ann is None:
               continue


           by_ann[ann]["values"].append(entry["values"])
           by_ann[ann]["genes"].append(str(g).strip())


       # 4) For each annotation: filter to high/low and create heatmaps
       x = xAxes
       for ann in by_ann:
           vals = by_ann[ann]["values"]
           genes = by_ann[ann]["genes"]
           if len(vals) == 0:
               continue


           arr = np.asarray(vals)
           subset_by = np.asarray([gene_to_fact(g) for g in genes], dtype=object)


           mask = np.isin(subset_by, ["high", "low"])
           arr_hl = arr[mask, :]
           subset_hl = subset_by[mask]


           if arr_hl.shape[0] == 0:
               continue


           fig1 = Plotting.heatmap_no_sorted(arr_hl, x, subset_by=subset_hl)
           fig2 = Plotting.heatmap_sorted_by_meanValues(arr_hl, x, subset_by=subset_hl)
           fig3 = Plotting.heatmap_sorted_by_maxValueIndex(arr_hl, x, subset_by=subset_hl)


           try:
               fig1.ax_hm_left.set_title(f"{ann} | FACT high (n={np.sum(subset_hl=='high')})")
               fig1.ax_hm_right.set_title(f"{ann} | FACT low (n={np.sum(subset_hl=='low')})")
           except Exception:
               pass

           plt.show()

class FactGeneUtils:
   """
   Utilities for gene ID normalization and FACT category mapping.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `strip_version()` | Removes version suffix from Ensembl gene identifiers. |
   | `get_gene_id_from_feature()` | Extracts gene ID from a GFF/GTF feature. |
   | `build_fact_map()` | Constructs a gene-to-FACT-category mapping from an annotation table. |
   """


   @staticmethod
   def strip_version(x: str) -> str:
       """
       Strip the Ensembl version suffix from a gene identifier.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `x` | `str` | — | Gene identifier, possibly including a version suffix (e.g., `"ENSG000001.12"`). |

       ## Returns

       `str` — Identifier without version suffix. Returns `x` unchanged if it is NaN-like.

       ## Examples

       ```python
       assert FactGeneUtils.strip_version("ENSG000001.12") == "ENSG000001"
       ```
       """
       if pd.isna(x):
           return x
       x = str(x)
       return x.split(".")[0]


   @staticmethod
   def get_gene_id_from_feature(f) -> str | None:
       """
       Attempt to extract a gene ID from a pybedtools GTF/GFF interval.

       Strategy (in order): (1) check `.attrs` for `gene_id`/`geneid`/`gene`; (2) use `.name`
       if it looks like an Ensembl ID; (3) regex search for `gene_id "..."` in the string form.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `f` | `Any` | — | A `pybedtools.Interval`-like object. |

       ## Returns

       `str | None` — Version-stripped gene ID if found, otherwise `None`.

       ## Examples

       ```python
       gid = FactGeneUtils.get_gene_id_from_feature(interval)
       assert gid is None or gid.startswith(("ENSG", "ENSMUSG"))
       ```
       """
       if hasattr(f, "attrs"):
           gid = f.attrs.get("gene_id") or f.attrs.get("geneid") or f.attrs.get("gene")
           if gid:
               return FactGeneUtils.strip_version(gid)


       if hasattr(f, "name") and f.name:
           if str(f.name).startswith(("ENSG", "ENSMUSG")):
               return FactGeneUtils.strip_version(f.name)


       s = str(f)
       m = re.search(r'gene_id "([^"]+)"', s)
       if m:
           return FactGeneUtils.strip_version(m.group(1))


       return None


   @staticmethod
   def build_fact_map(csv_path: str, fact_col: str) -> dict:
       """
       Build a mapping from stripped Ensembl gene ID to FACT category label.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `csv_path` | `str` | — | Path to FACT table (CSV/TSV). Must contain columns `'EnsemblGeneID'` and `fact_col`. |
       | `fact_col` | `str` | — | Column name with FACT category values (e.g., `"high"`, `"low"`). |

       ## Returns

       `dict` — Maps version-stripped Ensembl ID → lowercase category string.

       ## Examples

       ```python
       m = FactGeneUtils.build_fact_map("final_table.tsv", "ssrp1 FACT category")
       assert isinstance(m, dict)
       ```
       """
       df = pd.read_csv(csv_path, sep=None, engine="python")
       df["ENSG_stripped"] = df["EnsemblGeneID"].astype(str).apply(FactGeneUtils.strip_version)
       fact = df[fact_col].astype(str).str.strip().str.lower()
       return dict(zip(df["ENSG_stripped"], fact))


class ChipRnaUtils:
   """
   Static helper utilities for matching, deduplicating, and reshaping
   ChIP-seq signal arrays and RNA-seq differential-expression tables,
   used by `DecilePipeline` and related decile/scatter pipelines.

   ## Methods

   | Method | Description |
   |--------|-------------|
   | `deduplicate_rna_by_id()` | Deduplicate an RNA-seq DataFrame by gene ID, keeping the best row per ID. |
   | `create_log2_fold_change()` | Compute log2(chip1 / chip2) per matching gene name. |
   | `extract_ids()` | Extract ordered IDs from a chip signal list. |
   | `order_rna_by_chip()` | Reorder an RNA DataFrame to match a chip signal list's ID order. |
   | `sort_chip_by_mean_desc()` | Sort chip signal entries by mean(values), descending. |
   | `alter_chip_example()` | Build an altered copy of a chip signal list with block-wise randomized values. |
   | `min_max_chip_mean()` | Compute the min/max of mean(values) across chip signal entries. |
   | `replace_inf_with_clip()` | Replace +inf/-inf in a Series with the series' finite max/min. |
   | `split_into_groups()` | Split a 1-D sequence into `n_groups` contiguous, near-equal groups. |
   """

   @staticmethod
   def deduplicate_rna_by_id(
       rna_df: pd.DataFrame,
       id_col: str,
       strategy: str = "best_padj",
       padj_col: Optional[str] = None,
   ) -> pd.DataFrame:
       """
       Deduplicate an RNA-seq results table by gene ID.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `rna_df` | `pd.DataFrame` | — | RNA-seq results table. |
       | `id_col` | `str` | — | Column holding the gene ID to deduplicate on. |
       | `strategy` | `str` | `"best_padj"` | `"first"` keeps the first row per ID; `"best_padj"` keeps the row with the smallest `padj_col` per ID. |
       | `padj_col` | `str` | `None` | Required when `strategy="best_padj"`: column with adjusted p-values. |

       ## Returns

       `pd.DataFrame` — One row per unique `id_col` value.

       ## Raises

       - `ValueError` — If `strategy` is unknown, or `strategy="best_padj"` without a valid `padj_col`.
       """
       df = rna_df.copy()
       df[id_col] = df[id_col].astype(str)

       if strategy == "first":
           return df.drop_duplicates(subset=[id_col], keep="first")

       if strategy == "best_padj":
           if padj_col is None or padj_col not in df.columns:
               raise ValueError("padj_col is required for strategy='best_padj'")
           df = df.sort_values(by=padj_col, ascending=True, na_position="last")
           return df.drop_duplicates(subset=[id_col], keep="first")

       raise ValueError(f"Unknown strategy: {strategy}")

   @staticmethod
   def create_log2_fold_change(
       chip1: List[Dict[str, Any]],
       chip2: List[Dict[str, Any]],
       *,
       gene_key: str = "gene_name",
       values_key: str = "values",
       nan_value: float = 0.0,
   ) -> List[Dict[str, Any]]:
       """
       Compute log2(chip1 / chip2) per gene_name match.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `chip1` | `list[dict]` | — | Numerator signal entries; each must contain `gene_key` and `values_key`. |
       | `chip2` | `list[dict]` | — | Denominator signal entries, matched by `gene_key`. Mutated: matched entries are popped as they're consumed. |
       | `gene_key` | `str` | `"gene_name"` | Dict key used to match entries between `chip1` and `chip2`. |
       | `values_key` | `str` | `"values"` | Dict key holding the 1-D numeric signal array. |
       | `nan_value` | `float` | `0.0` | Value substituted for any NaN produced by the log2 ratio. |

       ## Returns

       `list[dict]` — One entry per match, with `gene_key` and a `"log2FoldChange"` array.
       Only pairs where both means are `> 0` are included; unmatched/zero-mean entries are silently skipped.
       """
       chip_final: List[Dict[str, Any]] = []

       for e in chip1:
           for i, el in enumerate(chip2):
               if e[gene_key] == el[gene_key]:
                   if np.mean(e[values_key]) > 0 and np.mean(el[values_key]) > 0:
                       real_values = np.log2(e[values_key] / el[values_key])
                       real_values = np.nan_to_num(real_values, nan=nan_value)
                       chip2.pop(i)
                       chip_final.append(
                           {
                               gene_key: e[gene_key],
                               "log2FoldChange": real_values,
                           }
                       )
                   break

       return chip_final

   @staticmethod
   def extract_ids(chip: List[Dict[str, Any]], *, id_key: str = "id") -> List[Any]:
       """Extract the `id_key` value from each chip signal entry, preserving order."""
       return [x[id_key] for x in chip]

   @staticmethod
   def order_rna_by_chip(
       chip: List[Dict[str, Any]],
       rna: pd.DataFrame,
       *,
       chip_id_key: str = "id",
       rna_id_col: str = "GeneID",
   ) -> pd.DataFrame:
       """
       Reorder an RNA-seq DataFrame to match the ID order of a chip signal list.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `chip` | `list[dict]` | — | Chip signal entries defining the desired order. |
       | `rna` | `pd.DataFrame` | — | RNA-seq table to reorder. |
       | `chip_id_key` | `str` | `"id"` | Dict key in `chip` holding the gene ID. |
       | `rna_id_col` | `str` | `"GeneID"` | Column in `rna` holding the matching gene ID. |

       ## Returns

       `pd.DataFrame` — `rna` reindexed to `chip`'s ID order.
       """
       ids = ChipRnaUtils.extract_ids(chip, id_key=chip_id_key)
       return rna.set_index(rna_id_col).loc[ids].reset_index()

   @staticmethod
   def sort_chip_by_mean_desc(
       chip: List[Dict[str, Any]],
       *,
       values_key: str = "values",
       in_place: bool = False,
   ) -> List[Dict[str, Any]]:
       """
       Sort chip signal entries by `mean(values)`, descending.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `chip` | `list[dict]` | — | Chip signal entries; each must contain `values_key`. |
       | `values_key` | `str` | `"values"` | Dict key holding the 1-D numeric signal array. |
       | `in_place` | `bool` | `False` | If `True`, sort the input list itself instead of a copy. Either way, the sorted list is returned. |

       ## Returns

       `list[dict]` — `chip` sorted by descending mean signal.
       """
       src = chip if in_place else list(chip)
       return sorted(src, key=lambda d: float(np.mean(d[values_key])), reverse=True)

   @staticmethod
   def alter_chip_example(
       chip: List[Dict[str, Any]],
       *,
       values_key: str = "values",
       segment_size: int = 10,
       ranges: Optional[Sequence[Tuple[float, float]]] = None,
       seed: Optional[int] = None,
   ) -> List[Dict[str, Any]]:
       """
       Build an altered copy of `chip` with block-wise randomized signal values.

       Entries are grouped into consecutive blocks of `segment_size`; each block
       draws its values from the next `(low, high)` range in `ranges`.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `chip` | `list[dict]` | — | Chip signal entries to alter (not mutated; a deep copy is used). |
       | `values_key` | `str` | `"values"` | Dict key holding the 1-D numeric signal array. |
       | `segment_size` | `int` | `10` | Number of consecutive entries per block. |
       | `ranges` | `list[tuple[float, float]]` | `None` | `(low, high)` uniform-sampling range per block, in order. Defaults to 10 descending ranges from `(1.8, 2.0)` to `(0.0, 0.2)`. |
       | `seed` | `int` | `None` | Seed for reproducible random sampling. |

       ## Returns

       `list[dict]` — New entries with randomized `"values"`, plus the same
       `gene_name`/`chr`/`start`/`end`/`strand`/`id` metadata as the input.
       """
       if seed is not None:
           rng = np.random.default_rng(seed)
           uniform = rng.uniform
       else:
           uniform = np.random.uniform

       if ranges is None:
           ranges = [
               (1.8, 2.0), (1.6, 1.8), (1.4, 1.6), (1.2, 1.4),
               (1.0, 1.2), (0.8, 1.0), (0.6, 0.8), (0.4, 0.6),
               (0.2, 0.4), (0.0, 0.2),
           ]

       chip_copy = copy.deepcopy(chip)
       out: List[Dict[str, Any]] = []

       for i, x in enumerate(chip_copy):
           segment = i // segment_size
           if segment >= len(ranges):
               segment = len(ranges) - 1
           low, high = ranges[segment]

           values_new = uniform(low, high, len(x[values_key]))

           out.append({
               "values": values_new,
               "gene_name": x.get("gene_name"),
               "chr": x.get("chr"),
               "start": x.get("start"),
               "end": x.get("end"),
               "strand": x.get("strand"),
               "id": x.get("id"),
           })

       return out

   @staticmethod
   def min_max_chip_mean(
       chip: List[Dict[str, Any]],
       *,
       values_key: str = "values",
   ) -> Tuple[float, float]:
       """
       Compute the min and max of `mean(values)` across chip signal entries.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `chip` | `list[dict]` | — | Chip signal entries; each must contain `values_key`. |
       | `values_key` | `str` | `"values"` | Dict key holding the 1-D numeric signal array. |

       ## Returns

       `tuple[float, float]` — `(min_mean, max_mean)`. `(0.0, 0.0)` if `chip` is empty.
       """
       means = [float(np.mean(x[values_key])) for x in chip]
       if not means:
           return (0.0, 0.0)
       return (min(means), max(means))

   @staticmethod
   def replace_inf_with_clip(
       series: pd.Series,
       *,
       clip_pos: Optional[float] = None,
       clip_neg: Optional[float] = None,
   ) -> pd.Series:
       """
       Replace +inf/-inf values in a Series with finite values.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `series` | `pd.Series` | — | Numeric series that may contain `inf`/`-inf` (e.g. log2 fold changes from zero counts). |
       | `clip_pos` | `float` | `None` | Value to substitute for `+inf`. Defaults to the series' finite max. |
       | `clip_neg` | `float` | `None` | Value to substitute for `-inf`. Defaults to the series' finite min. |

       ## Returns

       `pd.Series` — Copy of `series` with `+inf`/`-inf` replaced. If the series has no
       finite values at all, both defaults fall back to `0.0`.
       """
       s = series.copy()
       finite = s[np.isfinite(s)]
       if finite.empty:
           pos = 0.0 if clip_pos is None else clip_pos
           neg = 0.0 if clip_neg is None else clip_neg
       else:
           pos = float(finite.max()) if clip_pos is None else float(clip_pos)
           neg = float(finite.min()) if clip_neg is None else float(clip_neg)

       s.loc[np.isposinf(s)] = pos
       s.loc[np.isneginf(s)] = neg
       return s

   @staticmethod
   def split_into_groups(
       values: Sequence[float],
       n_groups: int = 10,
   ) -> List[List[float]]:
       """
       Split a 1-D sequence into `n_groups` contiguous, near-equal-size groups.

       ## Parameters

       | Parameter | Type | Default | Description |
       |-----------|------|---------|-------------|
       | `values` | `Sequence[float]` | — | Values to split, in order. |
       | `n_groups` | `int` | `10` | Number of groups to produce. |

       ## Returns

       `list[list[float]]` — Exactly `n_groups` lists (some may be empty if `len(values) < n_groups`),
       covering `values` in order with sizes differing by at most 1.

       ## Raises

       - `ValueError` — If `n_groups <= 0`.
       """
       values = list(values)
       if n_groups <= 0:
           raise ValueError("n_groups must be > 0")

       if len(values) == 0:
           return [[] for _ in range(n_groups)]

       sizes = [len(values) // n_groups] * n_groups
       for i in range(len(values) % n_groups):
           sizes[i] += 1

       out: List[List[float]] = []
       idx = 0
       for sz in sizes:
           out.append(values[idx:idx + sz])
           idx += sz
       return out


