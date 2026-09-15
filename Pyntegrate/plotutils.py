"""
Plotting and statistics utilities for genomic signal, used throughout `Pyntegrate.pipelines` and `Pyntegrate.SeqSignalAnalysis.Plotting`.

This module is the low-level plotting backend for the higher-level
pipeline classes: it builds heatmap/profile figure layouts, computes
sort orders and confidence intervals, and draws the specific chart
types (volcano plots, hexbin correlations, marginal-histogram
scatters, scaled-region profiles) that those pipelines return.

## Statistics helpers

| Function | Description |
|----------|-------------|
| `ci` | Bootstrap confidence interval for the mean of an array. |
| `nice_log` | `log2`, mapping non-positive values to `0` instead of `-inf`/`NaN`. |
| `tip_zscores` | Z-scores of a 1-D array. |
| `fdrcorrection` | Benjamini-Hochberg FDR correction for p-values. |
| `tip_fdr` | FDR-corrected p-values from TIP z-scores. |
| `prepare_logged` | Log-transform paired `x`/`y` arrays, replacing non-finite values. |

## Heatmap / multi-panel figure layout

| Function | Description |
|----------|-------------|
| `imshow_panels` | Main entry point: lays out one or more heatmaps with optional line plots, FACT category labels, and extra panels. |
| `_get_layout`, `_finalize_x_axes`, `_match_ylim`, `_apply_fact_edge_labels`, `_pick_kwargs` | Internal helpers for `imshow_panels`. |
| `new_shell`, `matrix_and_line_shell` | Build a bare figure/Axes shell for a heatmap (+ optional line/strip plot). |
| `clustered_sortind`, `new_clustered_sortind` | Row sort order for a heatmap via k-means clustering. |
| `input_ip_plots` | IP/control/subtracted heatmap trio with an average-profile line plot. |
| `annotate_peaks` | Convenience wrapper adding peak annotations to an existing plot. |

## Chart types

| Function / class | Description |
|-------------------|-------------|
| `ci_plot` | Line plot of a matrix's mean with a shaded confidence-interval band. |
| `MarginalHistScatter` | Scatter plot with marginal histograms on the top/right axes; used by `ResultsTable.scatter`. |
| `volcano_plot` | Volcano plot (log2 fold-change vs. -log10 p-value) with optional gene highlighting. |
| `hexbin_correlation_scatter` | Hexbin density scatter for correlating two signal/expression columns. |
| `generate_gene_scale_regions_profile` | Meta-gene profile over scaled gene bodies (TSS-to-TES), analogous to deepTools `plotProfile --scaleRegions`. |

## Internal helpers

`_updatecopy` and `_clean` are small dict/array utilities used by
`MarginalHistScatter` and `ResultsTable.scatter`.
"""
import itertools
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import gridspec
from .colormap_adjust import *
from scipy import stats
import pandas as pd
from matplotlib.ticker import FuncFormatter
def _pick_kwargs(kwargs_obj, idx):
    """
    Select the kwargs dict that applies to panel/column `idx`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `kwargs_obj` | `dict \| list[dict] \| tuple[dict] \| None` | — | `None` returns `{}`; a single `dict` is returned as-is (shared by all panels); a `list`/`tuple` of dicts returns the element at `idx`, or the last element if `idx` is out of range. |
    | `idx` | `int` | — | Panel/column index to select from `kwargs_obj` when it is a list/tuple. |

    ## Returns

    `dict` — The kwargs to use for panel `idx`.

    ## Raises

    `TypeError` — If `kwargs_obj` is not `None`, a `dict`, or a `list`/`tuple`.
    """
    if kwargs_obj is None:
        return {}

    if isinstance(kwargs_obj, dict):
        return kwargs_obj

    if isinstance(kwargs_obj, (list, tuple)):
        if len(kwargs_obj) == 0:
            return {}
        if idx < len(kwargs_obj):
            return kwargs_obj[idx]
        return kwargs_obj[-1]

    raise TypeError(
        f"line_kwargs/fill_kwargs must be a dict, list, tuple, or None, not {type(kwargs_obj)}"
    )

def ci_plot(x, arr, conf=0.1, ax=None, line_kwargs=None, fill_kwargs=None):
    """
    Plot the mean of `arr` as a line, with a shaded confidence interval band.

    The mean and confidence bounds are computed with [`ci`][Pyntegrate.plotutils.ci]
    over `arr`'s rows (i.e. across features, for each column/position).

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x` | `array-like` | — | 1-D x values for the plot. |
    | `arr` | `array-like` | — | 2-D array to compute the mean and confidence interval over (rows are averaged). |
    | `conf` | `float` | `0.1` | Confidence level forwarded to [`ci`][Pyntegrate.plotutils.ci]. |
    | `ax` | `matplotlib.axes.Axes \| None` | `None` | Axes to plot on. If `None`, a new figure and axes are created. |
    | `line_kwargs` | `dict \| None` | `None` | Additional kwargs passed to `Axes.plot` for the mean line. |
    | `fill_kwargs` | `dict \| None` | `None` | Additional kwargs passed to `Axes.fill_between` for the confidence band. |

    ## Returns

    `matplotlib.axes.Axes` — The axes the mean line and confidence band were drawn on.
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    line_kwargs = line_kwargs or {}
    fill_kwargs = fill_kwargs or {}

    m, lo, hi = ci(arr, conf)
    ax.plot(x, m, **line_kwargs)
    ax.fill_between(x, lo, hi, **fill_kwargs)
    return ax

def _apply_fact_edge_labels(axes, x):
    """
    Apply shared x-axis tick formatting to a set of axes.

    Sets the x-limits to `(min(x), max(x))` on every axes, then formats
    tick labels so that: the minimum tick shows thousands without "Kb",
    the maximum tick shows thousands with a "Kb" suffix, `0` is labeled
    "TSS", and both extremes plus `0` are guaranteed to appear as ticks.
    Edge tick labels are also left/center-aligned for readability.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `axes` | `list[matplotlib.axes.Axes]` | — | Axes to apply the shared x-limits and tick formatting to. No-op if empty. |
    | `x` | `array-like` | — | x values used to determine the axis extremes (`min`/`max`). |
    """
    if not axes:
        return

    xmin = float(np.min(x))
    xmax = float(np.max(x))

    def _edge_kb_formatter(val, pos=None):
        ival = int(round(val))
        if ival == int(round(xmin)) and abs(ival) >= 1000:
            return f"{ival // 1000}"
        if ival == int(round(xmax)) and abs(ival) >= 1000:
            return f"{ival // 1000}Kb"
        if ival == 0:
            return "TSS"
        if ival % 1000 != 0:
            return ""
        return str(ival)

    formatter = FuncFormatter(_edge_kb_formatter)

    for ax in axes:
        ax.set_xlim(xmin, xmax)

    base_ticks = axes[0].get_xticks()
    extra = [xmin, xmax]
    if xmin <= 0 <= xmax:
        extra.insert(1, 0.0)

    ticks = np.unique(np.r_[base_ticks, extra])
    ticks = ticks[(ticks >= xmin) & (ticks <= xmax)]
    ticks.sort()

    labels = [_edge_kb_formatter(t) for t in ticks]

    for ax in axes:
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(formatter)
        ax.set_xticklabels(labels)

    for ax in axes:
        labs = ax.get_xticklabels()
        if labs:
            labs[0].set_ha("left")
            labs[-1].set_ha("center")


def _finalize_x_axes(hm_axes, ln_axes, extra_axes, x, has_extra_bottom, xlabel=None):
    """
    Apply common x-axis tick formatting and decide which panel row shows labels.

    Calls [`_apply_fact_edge_labels`][Pyntegrate.plotutils._apply_fact_edge_labels]
    on all axes, hides x-tick labels on `hm_axes` (and on `ln_axes` too when
    `has_extra_bottom` is set), and sets `xlabel` on whichever row is the
    bottommost one actually shown.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `hm_axes` | `list[matplotlib.axes.Axes]` | — | Heatmap-row axes. |
    | `ln_axes` | `list[matplotlib.axes.Axes]` | — | Line-plot-row axes. |
    | `extra_axes` | `list[matplotlib.axes.Axes]` | — | Extra bottom-row axes (e.g. FACT tracks). |
    | `x` | `array-like` | — | x values forwarded to `_apply_fact_edge_labels`. |
    | `has_extra_bottom` | `bool` | — | Whether an extra bottom row is present; if so, `ln_axes` tick labels are hidden too and `xlabel` is placed on `extra_axes` instead of `ln_axes`. |
    | `xlabel` | `str \| None` | `None` | Label to set on the bottommost visible row. No label is set if `None`. |
    """
    all_axes = hm_axes + ln_axes + extra_axes
    _apply_fact_edge_labels(all_axes, x)

    for ax in hm_axes:
        plt.setp(ax.get_xticklabels(), visible=False)

    if has_extra_bottom:
        for ax in ln_axes:
            plt.setp(ax.get_xticklabels(), visible=False)

    if xlabel is not None:
        target_axes = extra_axes if has_extra_bottom and extra_axes else ln_axes
        for ax in target_axes:
            ax.set_xlabel(xlabel)


def _match_ylim(axes):
    """
    Set the same y-limits (the min/max across all axes) on every axes in `axes`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `axes` | `list[matplotlib.axes.Axes]` | — | Axes to synchronize. No-op if empty. |
    """
    if not axes:
        return
    ymins, ymaxs = zip(*(ax.get_ylim() for ax in axes))
    ymin, ymax = min(ymins), max(ymaxs)
    for ax in axes:
        ax.set_ylim(ymin, ymax)


def _get_layout(n_main_cols, has_extra_bottom, fig_width, adjust, show_panels=None):
    """
    Compute the GridSpec row count, height ratios, and figure size for `imshow_panels`.

    A single layout function serves both display modes of `imshow_panels`
    (with or without category subsetting).

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `n_main_cols` | `int` | — | Number of main (non-colorbar) columns; used to compute figure width. |
    | `has_extra_bottom` | `bool` | — | Whether an extra bottom row (FACT tracks) is present. |
    | `fig_width` | `float` | — | Per-column width, in inches. |
    | `adjust` | `dict` | — | `subplots_adjust` kwargs, passed straight through in the returned dict. |
    | `show_panels` | `bool \| None` | `None` | `True` shows only the heatmap and line rows; `False` shows only the extra bottom row; `None` shows all rows that apply. |

    ## Returns

    `dict` — Keys `nrows`, `height_ratios`, `figsize`, and `adjust`, ready to
    be used to build a `matplotlib.gridspec.GridSpec` and figure.
    """
    if show_panels is False:
        nrows = 1
        height_ratios = (1,)
        fig_height = 3.5
    elif show_panels is True or not has_extra_bottom:
        nrows = 2
        height_ratios = (4, 1)
        fig_height = 7.0
    else:
        nrows = 3
        height_ratios = (4, 1, 1)
        fig_height = 9.6

    figsize = ((fig_width * n_main_cols)+2.5, fig_height)

    return {
        "nrows": nrows,
        "height_ratios": height_ratios,
        "figsize": figsize,
        "adjust": adjust,
    }


def imshow_panels(
    arr,                          # np.ndarray 2D or dict[str, np.ndarray 2D]
    x=None,
    vmin=5,
    vmax=99,
    percentile=True,
    conf=0.95,
    sort_by=None,                 # if subset_by!=None: length == len(subset_by) (filtered or not, but consistent)
    subset_by=None,               # None => no categories | array => category mode
    subset_order=None,            # None => order of first appearance in subset_by | list/tuple => explicit order
    line_kwargs=None,
    fill_kwargs=None,
    imshow_kwargs=None,
    subplot_params=None,

    # optional
    extra_bottom_lines=None,
    extra_bottom_ylabel=None,
    extra_bottom_xlabel=None,
    extra_bottom_legend_kwargs=None,
    extra_bottom_legend_title="FACT",

    # optional
    single_category_title="",
    figsize=None,

    # panel filter: True => only panels 0 and 1 | False => only panel 2 | None => all
    show_panels=None,

):
    """
    Plot one or more heatmap panels (each with a summary line plot below it).

    Each entry of `arr` (a single 2-D array, or a `dict` of `name -> 2-D
    array` for multiple "FACTs") is drawn as a heatmap row plus, below it,
    a mean +/- confidence-interval line plot. An optional third row can show
    extra summary lines (e.g. FACT tracks) via `extra_bottom_lines`. A 3-D
    array of shape `(1, N, M)` is automatically unwrapped to `(N, M)`.

    Two mutually exclusive modes are supported:

    - **No categories** (`subset_by=None`): one panel per FACT. If given,
      `sort_by` must have length equal to that FACT's row count.
    - **Categories** (`subset_by` given): `subset_by` is a per-row label
      vector for `arr` (same length as `arr`'s rows). One panel is drawn
      per `(FACT, category)` combination. If given, `sort_by` must have the
      same length as `subset_by`. If `subset_order` is given, only those
      categories are plotted (in that order); rows within a category are
      sorted stably by `(sort_by, original index)`.

    `vmin`/`vmax` (absolute or percentile) are shared across all panels. In
    category mode, if `"high"`/`"low"` categories exist, they are used by
    default as the reference for the shared color scale (as in a two-panel
    high/low layout).

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `arr` | `np.ndarray \| dict[str, np.ndarray]` | — | A single 2-D array, or a mapping of FACT name to 2-D array. Each array's rows are features, columns are signal bins. |
    | `x` | `array-like \| None` | `None` | x-axis values (one per column of `arr`). Defaults to `np.arange(n_bins)`. |
    | `vmin` | `float` | `5` | Lower bound for the shared color scale (a percentile if `percentile=True`, else absolute). |
    | `vmax` | `float` | `99` | Upper bound for the shared color scale (a percentile if `percentile=True`, else absolute). |
    | `percentile` | `bool` | `True` | If `True`, interpret `vmin`/`vmax` as percentiles of the pooled data; otherwise as absolute values. |
    | `conf` | `float` | `0.95` | Confidence level forwarded to [`ci_plot`][Pyntegrate.plotutils.ci_plot] for the line-plot rows. |
    | `sort_by` | `array-like \| None` | `None` | Per-row sort key. Length must match the relevant FACT's rows (no-category mode) or `subset_by` (category mode). |
    | `subset_by` | `array-like \| None` | `None` | Per-row category labels for `arr`. `None` selects no-category mode. |
    | `subset_order` | `list \| tuple \| None` | `None` | Explicit category order/subset to plot (category mode only). `None` uses order of first appearance in `subset_by`. |
    | `line_kwargs` | `dict \| list[dict] \| None` | `None` | Kwargs for `Axes.plot` in the line-plot row, per panel (see [`_pick_kwargs`][Pyntegrate.plotutils._pick_kwargs]). |
    | `fill_kwargs` | `dict \| list[dict] \| None` | `None` | Kwargs for `Axes.fill_between` in the line-plot row, per panel. |
    | `imshow_kwargs` | `dict \| None` | `None` | Extra kwargs merged into the heatmap `imshow` call. |
    | `subplot_params` | `dict \| None` | `None` | Kwargs forwarded to `gridspec.GridSpec` (e.g. `wspace`, `hspace`). Defaults to `dict(wspace=0.15, hspace=0.15)`. |
    | `extra_bottom_lines` | `list[dict] \| None` | `None` | Extra summary lines drawn in a third row. Each item is a dict with keys `data` (array), optional `conf`, `line_kwargs`, `fill_kwargs`, and (category mode) `group`. |
    | `extra_bottom_ylabel` | `str \| None` | `None` | y-axis label for the extra bottom row (first column only). |
    | `extra_bottom_xlabel` | `str \| None` | `None` | x-axis label placed on the bottommost visible row. |
    | `extra_bottom_legend_kwargs` | `dict \| None` | `None` | Kwargs overriding the default legend placement for the extra bottom row. |
    | `extra_bottom_legend_title` | `str \| None` | `"FACT"` | Title text drawn to the left of the extra bottom row's legend. |
    | `single_category_title` | `str` | `""` | Title prefix used in no-category mode (and in category mode when there is only one category). |
    | `figsize` | `tuple[float, float] \| None` | `None` | Explicit figure size, in inches. Overrides the size computed by [`_get_layout`][Pyntegrate.plotutils._get_layout]. |
    | `show_panels` | `bool \| None` | `None` | `True` shows only the heatmap and line rows; `False` shows only the extra bottom row (requires `extra_bottom_lines`); `None` shows all applicable rows. |

    ## Returns

    `matplotlib.figure.Figure` — The created figure, with `ax_hm`, `ax_ln`,
    `ax_extra`, `cax`, `gs`, and `facts` attached as attributes (plus
    `groups` in category mode).

    ## Raises

    `ValueError` — If a FACT array is not 2-D, if `show_panels=False` is
    requested without `extra_bottom_lines`, or if `sort_by`/`extra_bottom_lines`
    data lengths are inconsistent with `arr`/`subset_by`.
    """
    fig_width = 2

    if subplot_params is None:
        subplot_params = dict(wspace=0.15, hspace=0.15)

    # --------------------------
    # Helpers
    # --------------------------
    def _as_2d_float(a):
        a = np.asarray(a)
        if a.ndim == 3 and a.shape[0] == 1:
            a = a[0]
        if a.ndim != 2:
            return a
        if a.dtype == object:
            a = a.astype(np.float64)
        else:
            a = a.astype(np.float64, copy=False)
        return a

    def _place_middle_legend(ax):
        handles, labels = ax.get_legend_handles_labels()
        if labels:
            ax.legend(
                loc="lower center",
                bbox_to_anchor=(0.5, 1.0),
                frameon=False,
                fontsize=8,
                borderaxespad=0.2,
                handletextpad=0.3,
                labelspacing=0.3,
            )

    def _collect_bottom_legend_items(extra_axes):
        """
        Collect the unique legend handles and labels from all bottom axes.
        """
        handles_out = []
        labels_out = []
        seen = set()

        for ax in extra_axes:
            handles, labels = ax.get_legend_handles_labels()
            for h, l in zip(handles, labels):
                if l and l not in seen:
                    handles_out.append(h)
                    labels_out.append(l)
                    seen.add(l)

        return handles_out, labels_out


    def _place_bottom_legend(extra_axes):

        if not extra_axes:
            return []

        handles, labels = _collect_bottom_legend_items(extra_axes)
        if not labels:
            return []

        artists = []

        for anchor_ax in extra_axes:

            title_artist = None
            legend_x0 = 0.0
            legend_width = 1.0

            # Left-aligned title text
            if extra_bottom_legend_title is not None:
                title_artist = anchor_ax.text(
                    0.0,
                    1.03,
                    extra_bottom_legend_title,
                    ha="left",
                    va="bottom",
                    fontsize=8 if show_panels is not False else 11,
                    transform=anchor_ax.transAxes,
                )
                legend_x0 = 0.16
                legend_width = 1.0 - legend_x0

            legend_kwargs = dict(
                loc="lower left",
                bbox_to_anchor=(legend_x0, 1.0, legend_width, 0.0),
                mode="expand",
                frameon=False,
                fontsize=8 if show_panels is not False else 11,
                handletextpad=0.2,
                columnspacing=0.75,
                handlelength=0.75,
                ncol=len(labels),
                borderaxespad=0.0,
            )

            if extra_bottom_legend_kwargs is not None:
                legend_kwargs.update(extra_bottom_legend_kwargs)

            leg = anchor_ax.legend(handles, labels, **legend_kwargs)

            artists.append((leg, title_artist))

        return artists
    # --------------------------
    # Normalize FACTs
    # --------------------------
    if isinstance(arr, dict):
        facts = [(str(k), _as_2d_float(v)) for k, v in arr.items()]
    else:
        facts = [("FACT", _as_2d_float(arr))] if subset_by is not None else [(single_category_title, _as_2d_float(arr))]

    for name, a in facts:
        if a.ndim != 2:
            raise ValueError(f"FACT '{name}' is not 2-dimensional (ndim={a.ndim}).")

    n_bins = facts[0][1].shape[1]
    if x is None:
        x = np.arange(n_bins)
    else:
        x = np.asarray(x)

    if line_kwargs is None:
        line_kwargs = dict(color="k")
    if fill_kwargs is None:
        fill_kwargs = dict(color="k", alpha=0.3)

    has_extra_bottom = extra_bottom_lines is not None

    if show_panels is False and not has_extra_bottom:
        raise ValueError(
            "show_panels=False requires extra_bottom_lines to render panel 2."
        )

    # =====================================================
    # NO-CATEGORY MODE
    # =====================================================
    if subset_by is None:
        n_facts = len(facts)
        n_main_cols = n_facts
        ncols = n_main_cols + 1  # + colorbar
        adjust_panels_2 = dict(
            left=0.25,
            right=1.0,
            bottom=0.16,
            top=0.85,
        )

        adjust = dict(
            left=0.29,
            right=0.77,
            bottom=0.08,
            top=0.95,
        )
        if show_panels is False and adjust_panels_2 is not None:
            adjust = adjust_panels_2

        layout = _get_layout(
            n_main_cols=n_main_cols,
            has_extra_bottom=has_extra_bottom,
            fig_width=fig_width,
            adjust=adjust,
            show_panels=show_panels,
        )
        if figsize is not None:
            layout["figsize"] = figsize

        fig = plt.figure(figsize=layout["figsize"])
        gs = gridspec.GridSpec(
            layout["nrows"],
            ncols,
            width_ratios=tuple([1] * n_main_cols + [0.06]),
            height_ratios=layout["height_ratios"],
            **subplot_params
        )

        cax = plt.subplot(gs[:, -1])
        hm_axes, ln_axes = [], []
        extra_axes = []
        mappable_for_cbar = None

        if show_panels is False:
            _row_hm, _row_ln, _row_extra = None, None, 0
        elif show_panels is True:
            _row_hm, _row_ln, _row_extra = 0, 1, None
        else:
            _row_hm, _row_ln = 0, 1
            _row_extra = 2 if has_extra_bottom else None

        flat = np.concatenate([a.ravel() for _, a in facts])
        if percentile:
            vmin_abs = np.percentile(flat, vmin)
            vmax_abs = np.percentile(flat, vmax)
        else:
            vmin_abs, vmax_abs = vmin, vmax

        cmap = smart_colormap(vmin_abs, vmax_abs)
        _imshow_kwargs = dict(
            origin="lower",
            cmap=cmap,
            vmin=vmin_abs,
            vmax=vmax_abs,
            aspect="auto"
        )
        if imshow_kwargs:
            _imshow_kwargs.update(imshow_kwargs)

        for i, (fact_name, a) in enumerate(facts):
            n_rows = a.shape[0]
            idx = np.arange(n_rows)

            if sort_by is not None:
                s = np.asarray(sort_by)
                if s.shape[0] != n_rows:
                    raise ValueError("sort_by must have length n_rows when no categories are given.")
                idx = idx[np.lexsort((idx, s[idx]))]

            arr_sorted = a[idx, :]

            ref_ax = None
            if _row_hm is not None:
                ax_hm = plt.subplot(gs[_row_hm, i])
                ref_ax = ax_hm
            else:
                ax_hm = None

            if _row_ln is not None:
                ax_ln = plt.subplot(gs[_row_ln, i], sharex=ref_ax) if ref_ax else plt.subplot(gs[_row_ln, i])
                if ref_ax is None:
                    ref_ax = ax_ln
            else:
                ax_ln = None

            if _row_extra is not None and has_extra_bottom:
                ax_extra = plt.subplot(gs[_row_extra, i], sharex=ref_ax) if ref_ax else plt.subplot(gs[_row_extra, i])
            else:
                ax_extra = None

            if ax_hm is not None:
                m = ax_hm.imshow(
                    arr_sorted,
                    extent=(x.min(), x.max(), 0, arr_sorted.shape[0]),
                    **_imshow_kwargs
                )
                if mappable_for_cbar is None:
                    mappable_for_cbar = m

            lkw = _pick_kwargs(line_kwargs, i)
            fkw = _pick_kwargs(fill_kwargs, i)
            lkw.setdefault("label", fact_name)

            clean = np.where(np.isfinite(arr_sorted), arr_sorted, np.nan)
            clean = clean[~np.isnan(clean).any(axis=1)]

            if ax_ln is not None and clean.shape[0] > 0:
                ci_plot(x, clean, conf=conf, ax=ax_ln, line_kwargs=lkw, fill_kwargs=fkw)

            if ax_ln is not None:
                _place_middle_legend(ax_ln)
            _title_ax = ax_hm if ax_hm is not None else ax_extra
            if _title_ax is not None:
                _title_pad = 20 if ax_hm is None else 6
                _title_ax.set_title(f"{single_category_title} (n={arr_sorted.shape[0]})", pad=_title_pad)

            if ax_extra is not None:
                for item in extra_bottom_lines:
                    data = np.asarray(item["data"], dtype=float)
                    data = np.where(np.isfinite(data), data, np.nan)
                    data = data[~np.isnan(data).any(axis=1)]
                    if data.shape[0] == 0:
                        continue

                    item_conf = item.get("conf", conf)
                    item_lkw = dict(item.get("line_kwargs", {}))
                    item_fkw = dict(item.get("fill_kwargs", {}))

                    ci_plot(
                        x,
                        data,
                        conf=item_conf,
                        ax=ax_extra,
                        line_kwargs=item_lkw,
                        fill_kwargs=item_fkw,
                    )

            if i == 0:
                title = (single_category_title + " transcripts") if (single_category_title != "") else "Transcripts"
                if ax_hm is not None:
                    ax_hm.set_ylabel(title)
                if ax_ln is not None:
                    ax_ln.set_ylabel("Average enrichment")
                if ax_extra is not None and extra_bottom_ylabel is not None:
                    ax_extra.set_ylabel(extra_bottom_ylabel)
            else:
                if ax_hm is not None:
                    plt.setp(ax_hm.get_yticklabels(), visible=False)
                if ax_ln is not None:
                    ax_ln.set_yticks([])
                if ax_extra is not None:
                    ax_extra.set_yticks([])

            if ax_hm is not None:
                hm_axes.append(ax_hm)
            if ax_ln is not None:
                ln_axes.append(ax_ln)
            if ax_extra is not None:
                extra_axes.append(ax_extra)

        if mappable_for_cbar is not None:
            plt.colorbar(mappable_for_cbar, cax=cax)
        else:
            cax.set_visible(False)

        _match_ylim(ln_axes)
        _match_ylim(extra_axes)

        _finalize_x_axes(
            hm_axes=hm_axes,
            ln_axes=ln_axes,
            extra_axes=extra_axes,
            x=x,
            has_extra_bottom=bool(extra_axes),
            xlabel=extra_bottom_xlabel,
        )

        if extra_axes:
            _place_bottom_legend(extra_axes)

        fig.ax_hm = hm_axes
        fig.ax_ln = ln_axes
        fig.ax_extra = extra_axes
        fig.cax = cax
        fig.gs = gs
        fig.facts = [n for n, _ in facts]
        fig.subplots_adjust(**layout["adjust"])

        return fig

    # =====================================================
    # CATEGORY MODE
    # =====================================================
    subset_by = np.asarray(subset_by, dtype=object)

    if subset_order is None:
        seen = []
        for c in subset_by.tolist():
            if c not in seen:
                seen.append(c)
        groups = seen
    else:
        groups = [g for g in subset_order if g in subset_by]

    idx_by_group = {g: np.nonzero(subset_by == g)[0] for g in groups}

    if sort_by is not None:
        sort_by = np.asarray(sort_by)
        if sort_by.shape[0] != subset_by.shape[0]:
            raise ValueError("sort_by must have the same length as subset_by when categories are given.")

    n_groups = len(groups)
    n_facts = len(facts)
    n_main_cols = n_groups * n_facts
    ncols = n_main_cols + 1
    adjust_panels_2 = dict(
        left=0.16,
        right=1.0,
        bottom=0.16,
        top=0.85,
    )
    
    adjust = dict(
        left=0.18,
        right=0.84,
        bottom=0.08,
        top=0.95,
    )
    if show_panels is False and adjust_panels_2 is not None:
        adjust = adjust_panels_2

    layout = _get_layout(
        n_main_cols=n_main_cols,
        has_extra_bottom=has_extra_bottom,
        fig_width=fig_width,
        adjust=adjust,
        show_panels=show_panels,
    )
    if figsize is not None:
        layout["figsize"] = figsize

    fig = plt.figure(figsize=layout["figsize"])
    gs = gridspec.GridSpec(
        layout["nrows"],
        ncols,
        width_ratios=tuple([1] * n_main_cols + [0.06]),
        height_ratios=layout["height_ratios"],
        **subplot_params
    )

    cax = plt.subplot(gs[:, -1])
    hm_axes, ln_axes = [], []
    extra_axes = []
    mappable_for_cbar = None

    if show_panels is False:
        _row_hm, _row_ln, _row_extra = None, None, 0
    elif show_panels is True:
        _row_hm, _row_ln, _row_extra = 0, 1, None
    else:
        _row_hm, _row_ln = 0, 1
        _row_extra = 2 if has_extra_bottom else None

    ref_groups = [g for g in ("high", "low") if g in idx_by_group]
    if not ref_groups:
        ref_groups = groups

    flat = []
    for _, a in facts:
        for g in ref_groups:
            flat.append(a[idx_by_group[g], :].ravel())
    flat = np.concatenate(flat)

    if percentile:
        vmin_abs = np.percentile(flat, vmin)
        vmax_abs = np.percentile(flat, vmax)
    else:
        vmin_abs, vmax_abs = vmin, vmax

    cmap = smart_colormap(vmin_abs, vmax_abs)
    _imshow_kwargs = dict(
        origin="lower",
        cmap=cmap,
        vmin=vmin_abs,
        vmax=vmax_abs,
        aspect="auto"
    )
    if imshow_kwargs:
        _imshow_kwargs.update(imshow_kwargs)

    only_one_group = (len(groups) == 1)

    for i, (fact_name, a) in enumerate(facts):
        for j, g in enumerate(groups):
            col = i * n_groups + j
            idx = idx_by_group[g]

            if sort_by is not None:
                s = sort_by[idx]
                idx = idx[np.lexsort((idx, s))]

            group_arr = a[idx, :]

            ref_ax = None
            if _row_hm is not None:
                ax_hm = plt.subplot(gs[_row_hm, col])
                ref_ax = ax_hm
            else:
                ax_hm = None

            if _row_ln is not None:
                ax_ln = plt.subplot(gs[_row_ln, col], sharex=ref_ax) if ref_ax else plt.subplot(gs[_row_ln, col])
                if ref_ax is None:
                    ref_ax = ax_ln
            else:
                ax_ln = None

            if _row_extra is not None and has_extra_bottom:
                ax_extra = plt.subplot(gs[_row_extra, col], sharex=ref_ax) if ref_ax else plt.subplot(gs[_row_extra, col])
            else:
                ax_extra = None

            if ax_hm is not None:
                m = ax_hm.imshow(
                    group_arr,
                    extent=(x.min(), x.max(), 0, group_arr.shape[0]),
                    **_imshow_kwargs
                )
                if mappable_for_cbar is None:
                    mappable_for_cbar = m

            lkw = _pick_kwargs(line_kwargs, col)
            fkw = _pick_kwargs(fill_kwargs, col)
            if only_one_group:
                lkw.setdefault("label", f"{fact_name} {g.capitalize()}")

            clean = np.where(np.isfinite(group_arr), group_arr, np.nan)
            clean = clean[~np.isnan(clean).any(axis=1)]

            if ax_ln is not None and clean.shape[0] > 0:
                ci_plot(x, clean, conf=conf, ax=ax_ln, line_kwargs=lkw, fill_kwargs=fkw)

            if ax_ln is not None:
                _place_middle_legend(ax_ln)

            _title_ax = ax_hm if ax_hm is not None else ax_extra
            if _title_ax is not None:
                _title_pad = 20 if ax_hm is None else 6
                if only_one_group and single_category_title != "":
                    _title_ax.set_title(f"{single_category_title} ({group_arr.shape[0]})", pad=_title_pad)
                else:
                    _title_ax.set_title(f"{fact_name} {g.capitalize()} {f'({group_arr.shape[0]})'}", pad=_title_pad)

            if ax_extra is not None:
                for item in extra_bottom_lines:
                    item_group = item.get("group", None)

                    if item_group is not None and item_group != g:
                        continue

                    data = np.asarray(item["data"], dtype=float)

                    if data.ndim == 2 and data.shape[0] == subset_by.shape[0]:
                        data = data[idx_by_group[g], :]
                    elif data.ndim == 2 and data.shape[0] == group_arr.shape[0]:
                        pass
                    elif item_group is not None:
                        continue
                    else:
                        raise ValueError(
                            f"extra_bottom_lines['data'] must have {subset_by.shape[0]} rows "
                            f"(unfiltered, to be filtered by category) or {group_arr.shape[0]} rows "
                            f"(already filtered for group '{g}')."
                        )

                    data = np.where(np.isfinite(data), data, np.nan)
                    data = data[~np.isnan(data).any(axis=1)]

                    if data.shape[0] == 0:
                        continue

                    item_conf = item.get("conf", conf)
                    item_lkw = dict(item.get("line_kwargs", {}))
                    item_fkw = dict(item.get("fill_kwargs", {}))

                    if "label" not in item_lkw:
                        item_lkw["label"] = item_group if item_group is not None else f"extra {g}"

                    ci_plot(
                        x,
                        data,
                        conf=item_conf,
                        ax=ax_extra,
                        line_kwargs=item_lkw,
                        fill_kwargs=item_fkw,
                    )

            if col == 0:
                if ax_hm is not None:
                    ax_hm.set_ylabel("features")
                if ax_ln is not None:
                    ax_ln.set_ylabel("mean")
                if ax_extra is not None and extra_bottom_ylabel is not None:
                    ax_extra.set_ylabel(extra_bottom_ylabel)
            else:
                if ax_hm is not None:
                    plt.setp(ax_hm.get_yticklabels(), visible=False)
                if ax_ln is not None:
                    ax_ln.set_yticks([])
                if ax_extra is not None:
                    ax_extra.set_yticks([])

            if ax_hm is not None:
                hm_axes.append(ax_hm)
            if ax_ln is not None:
                ln_axes.append(ax_ln)
            if ax_extra is not None:
                extra_axes.append(ax_extra)

    if mappable_for_cbar is not None:
        plt.colorbar(mappable_for_cbar, cax=cax)
    else:
        cax.set_visible(False)

    _match_ylim(ln_axes)
    _match_ylim(extra_axes)

    _finalize_x_axes(
        hm_axes=hm_axes,
        ln_axes=ln_axes,
        extra_axes=extra_axes,
        x=x,
        has_extra_bottom=bool(extra_axes),
        xlabel=extra_bottom_xlabel,
    )

    if extra_axes:
        _place_bottom_legend(extra_axes)

    fig.ax_hm = hm_axes
    fig.ax_ln = ln_axes
    fig.ax_extra = extra_axes
    fig.cax = cax
    fig.gs = gs
    fig.facts = [n for n, _ in facts]
    fig.groups = groups
    fig.subplots_adjust(**layout["adjust"])

    return fig



def ci(arr, conf=0.95):
    """
    Compute the column-wise mean and a normal-approximation confidence interval.

    NaNs are ignored (`nanmean`/`nanstd`). The z-score is derived from
    `conf` via the normal quantile function.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `arr` | `array-like` | — | 2-D array; the mean/std/CI are computed per column (axis 0). |
    | `conf` | `float` | `0.95` | Confidence level (e.g. `0.95` for a 95% CI). |

    ## Returns

    `tuple[np.ndarray, np.ndarray, np.ndarray]` — `(mean, lower_bound, upper_bound)`, one value per column of `arr`.
    """
    arr = np.asarray(arr)
    m = np.nanmean(arr, axis=0)
    s = np.nanstd(arr, axis=0, ddof=1)
    n = np.sum(~np.isnan(arr), axis=0)
    n = np.maximum(n, 1)
    # normal approximation
    z = stats.norm.ppf(0.5 + conf / 2)
    se = s / np.sqrt(n)
    lo = m - z * se
    hi = m + z * se
    return m, lo, hi



def nice_log(x):
    """
    Apply a signed log2 scale that also handles negative numbers.

    Computes `log2(|x| + 1)`, restoring the original sign afterwards, so
    values on both sides of zero can be shown on a log-like scale.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x` | `np.ndarray` | — | Array of values, positive and/or negative. |

    ## Returns

    `np.ndarray` — Signed log2-transformed values, same shape as `x`.
    """
    neg = x < 0
    xi = np.log2(np.abs(x) + 1)
    xi[neg] = -xi[neg]
    return xi


def tip_zscores(a):
    """
    Calculate "target identification from profiles" (TIP) z-scores.

    Follows Cheng et al. 2011, Bioinformatics 27(23):3221-3227: each row is
    weighted by the column-wise mean signal, summed into a per-feature
    score, and then converted to a z-score.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `a` | `np.ndarray` | — | 2-D array where each row is the signal for one feature. |

    ## Returns

    `np.ndarray` — One z-score per row (feature) of `a`.
    """
    weighted = a * a.mean(axis=0)
    scores = weighted.sum(axis=1)
    zscores = (scores - scores.mean()) / scores.std()
    return zscores


def _ecdf(x):
    '''
    No-frills empirical CDF used by `fdrcorrection`.
    '''
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)


def fdrcorrection(pvals, alpha=0.05, method='indep'):
    '''
    P-value correction for false discovery rate (FDR).

    Covers Benjamini/Hochberg (for independent or positively correlated
    tests) and Benjamini/Yekutieli (for general or negatively correlated
    tests) — the same methods available in `statsmodels.stats.multitest.multipletests`
    as `method="fdr_bh"` and `method="fdr_by"` respectively.

    NOTE: this function was copied from
    `statsmodels.sandbox.stats.multicomp.fdrcorrection0` (statsmodels
    0.5.0), to avoid requiring all of statsmodels as a dependency just for
    this function.

    If there is prior information on the fraction of true hypotheses, then
    `alpha` should be set to `alpha * m / m_0`, where `m` is the number of
    tests (given by the p-values) and `m_0` is an estimate of the number of
    true hypotheses (see Benjamini, Krieger and Yekutieli).

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `pvals` | `array-like` | — | P-values of the individual tests. |
    | `alpha` | `float` | `0.05` | Error rate (family-wise significance threshold). |
    | `method` | `str` | `'indep'` | `'indep'`/`'poscorr'` (or `'i'`/`'p'`) for Benjamini/Hochberg; `'negcorr'` (or `'n'`) for Benjamini/Yekutieli. |

    ## Returns

    `tuple[np.ndarray, np.ndarray]` — `(rejected, pvals_corrected)`: a
    boolean array (`True` where the null hypothesis is rejected) and the
    array of p-values adjusted for multiple hypothesis testing, both in the
    original order of `pvals`.
    '''
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1./np.arange(1, len(pvals_sorted)+1))  # corrected this
        ecdffactor = _ecdf(pvals_sorted) / cm
#    elif method in ['n', 'negcorr']:
#        cm = np.sum(np.arange(len(pvals)))
#        ecdffactor = ecdf(pvals_sorted)/cm
    else:
        raise ValueError('only indep and necorr implemented')
    reject = pvals_sorted <= ecdffactor*alpha
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1] = 1
    return reject[sortrevind], pvals_corrected[sortrevind]


def tip_fdr(a, alpha=0.05):
    """
    Return FDR-adjusted p-values for the TIP z-scores of `a`.

    Computes [`tip_zscores`][Pyntegrate.plotutils.tip_zscores], converts
    them to p-values via the normal PDF, and FDR-corrects them with
    [`fdrcorrection`][Pyntegrate.plotutils.fdrcorrection].

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `a` | `np.ndarray` | — | 2-D array where each row is the signal for one feature. |
    | `alpha` | `float` | `0.05` | False discovery rate threshold forwarded to `fdrcorrection`. |

    ## Returns

    `np.ndarray` — FDR-adjusted p-values, one per row (feature) of `a`.
    """
    zscores = tip_zscores(a)
    pvals = stats.norm.pdf(zscores)
    rejected, fdrs = fdrcorrection(pvals)
    return fdrs


def prepare_logged(x, y):
    """
    Transform `x` and `y` to a log scale while dealing with zeros.

    This function scales `x` and `y` such that the points that are zero in one
    array are set to the min of the other array.

    When plotting expression data, frequently one sample will have reads in
    a particular feature but the other sample will not.  Expression data also
    tends to look better on a log scale, but log(0) is undefined and therefore
    cannot be shown on a plot.  This function allows these points to be shown,
    piled up along one side of the plot.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x`, `y` | `np.ndarray` | — | Arrays to log2-transform, e.g. expression values from two samples. |

    ## Returns

    `tuple[np.ndarray, np.ndarray]` — `(log2(x), log2(y))`, with
    non-finite entries (from zeros) in either array replaced by the global
    minimum finite value across both arrays.
    """
    xi = np.log2(x)
    yi = np.log2(y)

    xv = np.isfinite(xi)
    yv = np.isfinite(yi)

    global_min = min(xi[xv].min(), yi[yv].min())
    global_max = max(xi[xv].max(), yi[yv].max())

    xi[~xv] = global_min
    yi[~yv] = global_min

    return xi, yi


def new_shell(figsize=(5, 12), strip=False, height_ratios=(4, 1),
              width_ratios=(4, 1), subplot_params=None):
    """
    Build an empty figure with array, line, colorbar, and optional strip axes.

    Similar to [`matrix_and_line_shell`][Pyntegrate.plotutils.matrix_and_line_shell]
    but built from an explicit `matplotlib.gridspec.GridSpec` with
    caller-supplied row/column ratios.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `figsize` | `tuple[float, float]` | `(5, 12)` | Figure size, in inches. |
    | `strip` | `bool` | `False` | If `True`, also create `fig.strip_axes`, sharing the y-axis with `fig.array_axes`. |
    | `height_ratios` | `tuple[float, ...]` | `(4, 1)` | Row height ratios for the `GridSpec`. |
    | `width_ratios` | `tuple[float, ...]` | `(4, 1)` | Column width ratios for the `GridSpec`. |
    | `subplot_params` | `dict \| None` | `None` | Extra kwargs forwarded to `gridspec.GridSpec`. |

    ## Returns

    `matplotlib.figure.Figure` — Figure with `array_axes`, `line_axes`
    (sharing the x-axis with `array_axes`), `cax`, `gs`, and (when
    `strip=True`) `strip_axes` attached as attributes.
    """
    if subplot_params is None:
        subplot_params = {}
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        len(height_ratios),
        len(width_ratios),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        **subplot_params)
    fig.array_axes = plt.subplot(gs[0, 0])
    if strip:
        fig.strip_axes = plt.subplot(gs[0, 1], sharey=fig.array_axes)
    fig.line_axes = plt.subplot(gs[1, 0], sharex=fig.array_axes)
    fig.cax = plt.subplot(gs[1, 1])
    fig.gs = gs
    return fig


def matrix_and_line_shell(figsize=(5, 12), strip=False):
    """
    Helper function to construct an empty figure that has space for a matrix,
    a summary line plot directly below it, a colorbar axis, and an optional
    "strip" axis that parallels the matrix (and shares its y-axis) where data
    can be added to create callbacks.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `figsize` | `tuple[float, float]` | `(5, 12)` | Figure size, `(width, height)`, in inches. |
    | `strip` | `bool` | `False` | If `False`, the returned `strip_ax` is `None` and no strip axes are created. |

    ## Returns

    `tuple` — `(fig, matrix_ax, line_ax, strip_ax, colorbar_ax)`, where
    `line_ax` shares its x-axis with `matrix_ax`, `strip_ax` shares its
    y-axis with `matrix_ax` (or is `None` when `strip=False`), and
    `colorbar_ax` is an empty axes reserved for a colorbar.
    """
    fig = plt.figure(figsize=figsize)

    # Constants to keep track
    if strip:
        STRIP_COLS = 1
    else:
        STRIP_COLS = 0
    ROWS = 4
    COLS = 8 + STRIP_COLS
    MAT_COLS = 7
    MAT_ROWS = 3
    LINE_ROWS = ROWS - MAT_ROWS

    mat_ax = plt.subplot2grid(
        shape=(ROWS, COLS),
        loc=(0, STRIP_COLS),
        rowspan=MAT_ROWS,
        colspan=MAT_COLS,
    )

    line_ax = plt.subplot2grid(
        shape=(ROWS, COLS),
        loc=(MAT_ROWS, STRIP_COLS),
        rowspan=LINE_ROWS,
        colspan=MAT_COLS,
        sharex=mat_ax)

    if strip:
        strip_ax = plt.subplot2grid(
            shape=(ROWS, COLS),
            loc=(0, 0),
            rowspan=MAT_ROWS,
            colspan=STRIP_COLS,
            sharey=mat_ax,
        )
    else:
        strip_ax = None

    cax = plt.subplot2grid(
        shape=(ROWS, COLS),
        loc=(ROWS - MAT_ROWS, MAT_COLS + STRIP_COLS),
        rowspan=1,
        colspan=1,
    )

    fig.subplots_adjust(hspace=0.1, wspace=0.2, right=0.88, left=0.23)
    return fig, mat_ax, line_ax, strip_ax, cax


def clustered_sortind(x,random_state, k=10, scorefunc=None):
    """
    Uses MiniBatch k-means clustering to cluster matrix into groups.

    Each cluster of rows is then sorted by `scorefunc` -- by default, the max
    peak height when all rows in a cluster are averaged, or
    cluster.mean(axis=0).max().

    Returns the index that will sort the rows of `x` and a list of "breaks".
    `breaks` is essentially a cumulative row count for each cluster boundary.
    In other words, after plotting the array you can use axhline on each
    "break" to plot the cluster boundary.

    If `k` is a list or tuple, iteratively try each one and select the best
    with the lowest mean distance from cluster centers.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x` | `array-like` | — | Matrix whose rows are to be clustered. |
    | `random_state` | `int \| np.random.RandomState` | — | Random state forwarded to `MiniBatchKMeans`. |
    | `k` | `int \| list[int] \| tuple[int, ...]` | `10` | Number of clusters, or a list/tuple of candidate cluster counts to try (the one with the lowest mean cluster-transform distance is chosen). |
    | `scorefunc` | `callable \| None` | `None` | Function used to sort rows within a cluster; must accept a single NumPy array argument. Defaults to `cluster.mean(axis=0).max()`. |

    ## Returns

    `tuple[np.ndarray, list[int]]` — `(ind, breaks)`: `ind` is the row
    index that sorts `x` by cluster (clusters ordered by `scorefunc`, then
    rows within a cluster in original order); `breaks` is the cumulative
    row count at each cluster boundary, suitable for `axhline`.

    ## Raises

    `ImportError` — If `scikit-learn` is not installed.
    """
    try:
        from sklearn.cluster import MiniBatchKMeans
    except ImportError:
        raise ImportError('please install scikits.learn for '
                          'clustering.')

    # If integer, do it once and we're done
    if isinstance(k, int):
        best_k = k

    else:
        mean_dists = {}
        for _k in k:
            mbk = MiniBatchKMeans(batch_size=100, compute_labels=True, init='k-means++',
        init_size=None, max_iter=100, max_no_improvement=10, n_clusters=_k,
        n_init=3, random_state=random_state, reassignment_ratio=0.01, tol=0.0,
        verbose=0)
            mbk.fit(x)
            mean_dists[_k] = mbk.transform(x).mean()
        best_k = sorted(mean_dists.items(), key=lambda x: x[1])[0][0]

    ## Refit with the chosen best_k
    mbk = MiniBatchKMeans(batch_size=100, compute_labels=True, init='k-means++',
        init_size=None, max_iter=100, max_no_improvement=10, n_clusters=best_k,
        n_init=3, random_state=random_state, reassignment_ratio=0.01, tol=0.0,
        verbose=0)
    mbk.fit(x)
    k = best_k
    labels = mbk.labels_
    scores = np.zeros(labels.shape, dtype=float)

    if not scorefunc:
        def scorefunc(x):
            return x.mean(axis=0).max()

    for label in range(k):
        ind = labels == label
        score = scorefunc(x[ind, :])
        scores[ind] = score

    pos = 0
    breaks = []
    ind = np.argsort(scores)
    for k, g in itertools.groupby(labels[ind]):
        pos += len(list(g))
        breaks.append(pos)

    return ind, breaks


def new_clustered_sortind(x,random_state, k=10, row_key=None, cluster_key=None):
    """
    Uses MiniBatch k-means clustering to cluster matrix into groups.

    Each cluster of rows is then sorted by `scorefunc` -- by default, the max
    peak height when all rows in a cluster are averaged, or
    cluster.mean(axis=0).max().

    Returns the index that will sort the rows of `x` and a list of "breaks".
    `breaks` is essentially a cumulative row count for each cluster boundary.
    In other words, after plotting the array you can use axhline on each
    "break" to plot the cluster boundary.

    If `k` is a list or tuple, iteratively try each one and select the best
    with the lowest mean distance from cluster centers.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x` | `array-like` | — | Matrix whose rows are to be clustered. |
    | `random_state` | `int \| np.random.RandomState` | — | Random state forwarded to `MiniBatchKMeans`. |
    | `k` | `int \| list[int] \| tuple[int, ...]` | `10` | Number of clusters, or a list/tuple of candidate cluster counts to try (the one with the lowest mean cluster-transform distance is chosen). |
    | `row_key` | `callable \| None` | `None` | Optional sort key for rows within a cluster. Signature `row_key(a)`, where `a` is a 1-D NumPy array (one row of `x`). |
    | `cluster_key` | `callable \| None` | `None` | Optional sort key for ordering clusters. Signature `cluster_key(a)`, where `a` is the 2-D NumPy array of all rows of `x` belonging to that cluster; must return a single comparable value. |

    ## Returns

    `tuple[np.ndarray, np.ndarray]` — `(final_ind, breaks)`: `final_ind`
    is the row index that sorts `x` by cluster (clusters ordered by
    `cluster_key` if given, rows within a cluster ordered by `row_key` if
    given); `breaks` is the cumulative row count at each cluster boundary.

    ## Raises

    `ImportError` — If `scikit-learn` is not installed.
    """
    try:
        from sklearn.cluster import MiniBatchKMeans
    except ImportError:
        raise ImportError('please install scikits.learn for '
                          'clustering.')

    # If integer, do it once and we're done
    if isinstance(k, int):
        best_k = k

    else:
        mean_dists = {}
        for _k in k:
            mbk = MiniBatchKMeans(init='k-means++', n_clusters=_k,random_state=random_state)
            mbk.fit(x)
            mean_dists[_k] = mbk.transform(x).mean()
        best_k = sorted(mean_dists.items(), key=lambda x: x[1])[0][0]

    mbk = MiniBatchKMeans(init='k-means++', n_clusters=best_k,batch_size=100,n_init=3,random_state=random_state)
    mbk.fit(x)
    k = best_k
    labels = mbk.labels_
    scores = np.zeros(labels.shape, dtype=float)

    if cluster_key:
        # It's easier for calling code to provide something that operates on
        # a cluster level, but here it's converted to work on a label level
        # that looks in to the array `x`.
        def _cluster_key(i):
            return cluster_key(x[labels == i, :])
        sorted_labels = sorted(range(k), key=_cluster_key)
    else:
        # Otherwise just use them as-is.
        sorted_labels = range(k)

    if row_key:
        # Again, easier to provide a function to operate on a row.  But here we
        # need it to accept an index
        def _row_key(i):
            return row_key(x[i, :])

    final_ind = []
    breaks = []
    pos = 0
    for label in sorted_labels:
        # which rows in `x` have this label
        label_inds = np.nonzero(labels == label)[0]
        if row_key:
            label_sort_ind = sorted(label_inds, key=_row_key)
        else:
            label_sort_ind = label_inds
        for li in label_sort_ind:
            final_ind.append(li)
        pos += len(label_inds)
        breaks.append(pos)

    return np.array(final_ind), np.array(breaks)


def input_ip_plots(iparr, inputarr, diffed, x, sort_ind,
                   prefix=None, limits1=(None, None), limits2=(None, None),
                   hlines=None, vlines=None):

    """
    All-in-one plotting function to make a 5-panel figure.

    Panels are IP, input, and diffed; plus 2 line plots showing averages.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `iparr`, `inputarr` | `np.ndarray` | — | Arrays constructed by a `genomic_signal` object (IP and input signal, respectively). |
    | `diffed` | `np.ndarray` | — | Difference of `iparr` and `inputarr` (or some other transformation of the two). |
    | `x` | `array-like` | — | x-axis extent, e.g. `np.linspace(-1000, 1000, bins)` for TSS-centered windows, or `np.arange(bins)` for bin IDs. |
    | `sort_ind` | `array-like` | — | Row order applied to all 3 heatmap panels — typically the output of [`clustered_sortind`][Pyntegrate.plotutils.clustered_sortind] or [`tip_zscores`][Pyntegrate.plotutils.tip_zscores]. |
    | `prefix` | `str \| None` | `None` | Prefix used in the IP/input panel titles (e.g. `"{prefix} IP"`). |
    | `limits1` | `tuple` | `(None, None)` | `(vmin, vmax)` passed to `Normalize` for the IP and input panels. Missing bounds are filled from the 1st/99.99...th percentile of the pooled IP+input data. |
    | `limits2` | `tuple` | `(None, None)` | `(vmin, vmax)` passed to `Normalize` for the diffed panel, filled the same way from `diffed`. |
    | `hlines` | `list[tuple] \| None` | `None` | `(position, kwargs)` pairs; a horizontal line is drawn on each of the 3 heatmap panels via `axhline(position, **kwargs)`. Useful for marking cluster boundaries from `clustered_sortind`'s `breaks`. |
    | `vlines` | `list[tuple] \| None` | `None` | `(position, kwargs)` pairs; a vertical line is drawn on each of the 3 heatmap panels via `axvline(position, **kwargs)`. |

    ## Returns

    `matplotlib.figure.Figure` — The 5-panel figure (3 heatmaps, 2 line
    plots, 2 colorbars).
    """

    # global min and max
    gmin = min(iparr.min(), inputarr.min())
    gmax = max(iparr.max(), inputarr.max())

    fig = plt.figure(figsize=(10, 10))

    # 3 arrays, 2 line plots, a gene strip, and 2 colorbars.  Plots share the
    # axes that make sense
    #
    # 3 arrays
    ax1 = plt.subplot2grid(
        (9, 9), (0, 0), colspan=3, rowspan=6)
    ax2 = plt.subplot2grid(
        (9, 9), (0, 3), colspan=3, rowspan=6, sharex=ax1, sharey=ax1)
    ax3 = plt.subplot2grid(
        (9, 9), (0, 6), colspan=3, rowspan=6, sharex=ax1, sharey=ax1)

    # 2 line plots
    ax4 = plt.subplot2grid((9, 9), (6, 3), colspan=3, rowspan=3, sharex=ax1)
    ax5 = plt.subplot2grid((9, 9), (6, 6), colspan=3, rowspan=3, sharex=ax1)

    # 2 colorbars
    cax1 = plt.Axes(fig, rect=(0.05, 0.25, 0.25, 0.025))
    cax2 = plt.Axes(fig, rect=(0.05, 0.15, 0.25, 0.025))

    # For nice imshow axes
    extent = (min(x), max(x), 0, diffed.shape[0])

    cm = matplotlib.cm.gist_gray
    cm.set_bad('k')
    cm.set_over('r')
    cm.set_under('b')

    limits1 = list(limits1)
    limits2 = list(limits2)

    all_base = np.column_stack((iparr.ravel(), inputarr.ravel())).ravel()

    if limits1[0] is None:
        limits1[0] = np.percentile(
            all_base, 1. / all_base.size)
    if limits1[1] is None:
        limits1[1] = np.percentile(
            all_base, 100 - 1. / all_base.size)
    if limits2[0] is None:
        limits2[0] = np.percentile(
            diffed.ravel(), 1. / all_base.size)
    if limits2[1] is None:
        limits2[1] = np.percentile(
            diffed.ravel(), 100 - 1. / all_base.size)

    del all_base

    imshow_kwargs = dict(
        interpolation='nearest',
        aspect='auto',
        cmap=cm,
        norm=matplotlib.colors.Normalize(*limits1),
        extent=extent,
        origin='lower')

    # modify kwargs for diffed (by changing the normalization)
    diffed_kwargs = imshow_kwargs.copy()
    diffed_kwargs['norm'] = matplotlib.colors.Normalize(*limits2)

    # IP
    mappable1 = ax1.imshow(iparr[sort_ind, :], **imshow_kwargs)

    # input
    mappable2 = ax2.imshow(inputarr[sort_ind, :], **imshow_kwargs)

    # diffed
    mappable3 = ax3.imshow((diffed)[sort_ind, :], **diffed_kwargs)

    # IP and input line plot with vertical line
    ax4.plot(x, inputarr.mean(axis=0), color='k', linestyle='--',
             label='input')
    ax4.plot(x, iparr.mean(axis=0), color='k', label='ip')
    ax4.axvline(0, color='k', linestyle=':')

    # Diffed line plot with vertical line
    ax5.plot(x, diffed.mean(axis=0), 'k', label='enrichment')
    ax5.axvline(0, color='k', linestyle=':')

    # Colorbars
    cbar1 = fig.colorbar(mappable1, cax1, orientation='horizontal')
    cbar2 = fig.colorbar(mappable3, cax2, orientation='horizontal')
    fig.add_axes(cax1)
    fig.add_axes(cax2)

    # labeling...
    ax1.set_ylabel('features')
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    ax4.set_xlabel('bp')
    ax4.set_ylabel('mean reads per million mapped reads')
    ax5.set_xlabel('bp')
    cax1.set_xlabel('Reads per million mapped reads')
    cax2.set_xlabel('Enrichment (RPMMR)')

    if prefix is None:
        prefix = ""
    ax1.set_title('%s IP' % prefix)
    ax2.set_title('%s input' % prefix)
    ax3.set_title('Difference')

    # diffed line plot should have y ax on right
    ax5.yaxis.set_ticks_position('right')
    ax5.yaxis.set_label_position('right')
    ax5.set_ylabel('enriched reads per million mapped reads')

    # Legends
    ax4.legend(loc='best', frameon=False)
    ax5.legend(loc='best', frameon=False)

    # Make sure everybody snaps to xmin/xmax
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.axis(xmin=extent[0], xmax=extent[1])

    if not hlines:
        hlines = []
    if not vlines:
        vlines = []

    for ax in [ax1, ax2, ax3]:
        for pos, kwargs in hlines:
            ax.axhline(pos, **kwargs)
        for pos, kwargs in vlines:
            ax.axvline(pos, **kwargs)

    fig.subplots_adjust(bottom=0.05, top=0.95, hspace=0.75, wspace=0.9)

    return fig


def annotate_peaks(pandas_df):
    """
    Plot a pie chart of peak counts by annotation category.

    Strips any parenthetical suffix from the `"Annotation"` column (e.g.
    `"Promoter (<=1kb)"` becomes `"Promoter "`), counts occurrences of each
    resulting category, and plots them as a pie chart with a legend.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `pandas_df` | `pandas.DataFrame` | — | DataFrame with an `"Annotation"` column (as produced by peak annotation tools such as HOMER). Modified in place (the column is overwritten with the stripped values). |

    ## Returns

    `matplotlib.figure.Figure` — The pie chart figure.
    """
    # fig = plt.figure()
    pandas_df["Annotation"] = pandas_df["Annotation"].str.replace(r'\(.*', '', regex=True)
    unique_values = pandas_df['Annotation'].unique()
    counts = pandas_df['Annotation'].value_counts()

    labels = []
    sizes = []
    explodes = []


    for value in unique_values:
        labels.append(str(value))
        size = counts.get(value)
        sizes.append(size)

    labels = ['Not known' if label is None else label for label in labels]
    sizes = [0 if size is None else size for size in sizes]

    fig = plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    plt.pie(sizes,autopct='%1.2f%%',shadow=True, radius=3)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.legend(labels, loc="best", bbox_to_anchor=(1.0, 1.0))
    plt.axis("equal")
    return fig

def _updatecopy(orig, update_with, keys=None, override=False):
    """
    Return a copy of `orig` updated with values from `update_with`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `orig` | `dict` | — | Base dict to copy and update. |
    | `update_with` | `dict` | — | Dict to pull new values from. |
    | `keys` | `iterable \| None` | `None` | If given, only these keys are considered (still only applied if present in `update_with`). Defaults to all of `update_with`'s keys. |
    | `override` | `bool` | `False` | If `False`, existing keys in `orig` are left untouched. If `True`, they are overwritten by `update_with`. |

    ## Returns

    `dict` — A new dict; `orig` is not modified.
    """
    d = orig.copy()
    if keys is None:
        keys = update_with.keys()
    for k in keys:
        if k in update_with:
            if k in d and not override:
                continue
            d[k] = update_with[k]
    return d


def _clean(z):
    """
    Return only the finite values of `z`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `z` | `np.ndarray` | — | Array to filter. |

    ## Returns

    `np.ndarray` — 1-D array of `z`'s finite values (NaN/inf entries removed).
    """
    return z[np.isfinite(z)]


class MarginalHistScatter(object):
    """
    Class to enable incremental appending of scatterplots, each of which
    generate additional marginal histograms.

    Each call to `append` adds a new scatter to `scatter_ax` plus a pair of
    marginal histogram axes (one above, one to the right), carved out of the
    figure via `divider`.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `scatter_ax` | `matplotlib.axes.Axes` | — | The central scatter-plot axes, passed in as the constructor's `ax` argument; new scatters are drawn here via `append`. |
    | `hist_size` | `float` | `0.6` | Size of each new marginal histogram axes, as passed to `mpl_toolkits.axes_grid1.make_axes_locatable`. |
    | `pad` | `float` | `0.05` | Padding between `scatter_ax` and each new marginal histogram axes. |
    | `fig` | `matplotlib.figure.Figure` | — | Figure that `scatter_ax` belongs to. |
    | `divider` | `mpl_toolkits.axes_grid1.axes_divider.AxesDivider` | — | Divider used to carve new marginal histogram axes out of `scatter_ax` on each `append` call. |
    | `top_hists` | `list[matplotlib.axes.Axes]` | — | Marginal histogram axes appended above `scatter_ax`, one per `append` call (when `marginal_histograms=True`). |
    | `right_hists` | `list[matplotlib.axes.Axes]` | — | Marginal histogram axes appended to the right of `scatter_ax`, one per `append` call (when `marginal_histograms=True`). |
    | `xfirst_ax` | `matplotlib.axes.Axes \| None` | — | First top-histogram axes created; reused as the shared y-axis for later top histograms when `hist_share=True`. Stays `None` until the first such call. |
    | `yfirst_ax` | `matplotlib.axes.Axes \| None` | — | First right-histogram axes created; reused as the shared x-axis for later right histograms when `hist_share=True`. Stays `None` until the first such call. |
    | `hxs` | `list[np.ndarray]` | — | Cleaned (NaN-stripped) x-data arrays from each `append` call, in order. |
    | `hys` | `list[np.ndarray]` | — | Cleaned (NaN-stripped) y-data arrays from each `append` call, in order. |
    """
    def __init__(self, ax, hist_size=0.6, pad=0.05):
        """
        Initialize the scatter axes and its (initially empty) histogram bookkeeping.

        See the class docstring's `## Attributes` table for the parameters
        and the attributes they populate.
        """
        self.scatter_ax = ax
        self.fig = ax.figure

        self.divider = make_axes_locatable(self.scatter_ax)

        self.top_hists = []
        self.right_hists = []
        self.hist_size = hist_size
        self.pad = pad
        self.xfirst_ax = None
        self.yfirst_ax = None

        # will hold histogram data
        self.hxs = []
        self.hys = []

    @property
    def xmax(self):
        """
        `float` — Maximum x-value across all data plotted on `scatter_ax`
        so far (`scatter_ax.dataLim.xmax`).
        """
        return self.scatter_ax.dataLim.xmax

    @property
    def ymax(self):
        """
        `float` — Maximum y-value across all data plotted on `scatter_ax`
        so far (`scatter_ax.dataLim.ymax`).
        """
        return self.scatter_ax.dataLim.ymax

    @property
    def xmin(self):
        """
        `float` — Minimum x-value across all data plotted on `scatter_ax`
        so far (`scatter_ax.dataLim.xmin`).
        """
        return self.scatter_ax.dataLim.xmin

    @property
    def ymin(self):
        """
        `float` — Minimum y-value across all data plotted on `scatter_ax`
        so far (`scatter_ax.dataLim.ymin`).
        """
        return self.scatter_ax.dataLim.ymin

    @property
    def limits(self):
        """
        `tuple[float, float, float, float]` — `(xmin, xmax, ymin, ymax)`
        across all data plotted on `scatter_ax` so far.
        """
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def append(self, x, y, scatter_kwargs, hist_kwargs=None, xhist_kwargs=None,
               yhist_kwargs=None, num_ticks=3, labels=None, hist_share=False,
               marginal_histograms=True):
        """
        Add a new scatter to `scatter_ax`, plus marginal histograms for the same data.

        Each call carves a new pair of histogram axes (top and right) out of
        `scatter_ax` via `self.divider`, in addition to drawing the scatter
        itself. NaNs in `x`/`y` are tolerated by `scatter` but stripped
        before histogramming.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `x`, `y` | `array-like` | — | Data to plot. |
        | `scatter_kwargs` | `dict` | — | Keyword arguments passed directly to `Axes.scatter`. |
        | `hist_kwargs` | `dict \| None` | `None` | Keyword arguments passed directly to `Axes.hist`, for both the top and right histograms. |
        | `xhist_kwargs`, `yhist_kwargs` | `dict \| None` | `None` | Margin-specific kwargs for the x or y histogram respectively; used to update a copy of `hist_kwargs` (without overriding keys already set there). |
        | `num_ticks` | `int` | `3` | Number of tick marks on each histogram's frequency axis. |
        | `labels` | `array-like \| None` | `None` | Optional labels set on the scatter `PathCollection` (as `.labels`) so a callback can look them up later. |
        | `hist_share` | `bool` | `False` | If `True`, all top histograms share one y-axis and all right histograms share one x-axis (useful for comparing relative heights without `density=True`). |
        | `marginal_histograms` | `bool` | `True` | If `False`, only the scatter is drawn — no histogram axes are created and the method returns immediately after. |
        """
        scatter_kwargs = scatter_kwargs or {}
        hist_kwargs = hist_kwargs or {}
        xhist_kwargs = xhist_kwargs or {}
        yhist_kwargs = yhist_kwargs or {}
        yhist_kwargs.update(dict(orientation='horizontal'))

        # Plot the scatter
        coll = self.scatter_ax.scatter(x, y, **scatter_kwargs)
        coll.labels = labels

        if not marginal_histograms:
            return

        xhk = _updatecopy(hist_kwargs, xhist_kwargs)
        yhk = _updatecopy(hist_kwargs, yhist_kwargs)

        axhistx = self.divider.append_axes(
            'top', size=self.hist_size,
            pad=self.pad, sharex=self.scatter_ax, sharey=self.xfirst_ax)

        axhisty = self.divider.append_axes(
            'right', size=self.hist_size,
            pad=self.pad, sharey=self.scatter_ax, sharex=self.yfirst_ax)

        axhistx.yaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        axhisty.xaxis.set_major_locator(
            MaxNLocator(nbins=num_ticks, prune='both'))

        if not self.xfirst_ax and hist_share:
            self.xfirst_ax = axhistx

        if not self.yfirst_ax and hist_share:
            self.yfirst_ax = axhisty

        # Keep track of which axes are which, because looking into fig.axes
        # list will get awkward....
        self.top_hists.append(axhistx)
        self.right_hists.append(axhisty)

        # Scatter will deal with NaN, but hist will not.  So clean the data
        # here.
        hx = _clean(x)
        hy = _clean(y)

        self.hxs.append(hx)
        self.hys.append(hy)

        # Only plot hists if there's valid data
        if len(hx) > 0:
            if len(hx) == 1:
                _xhk = _updatecopy(orig=xhk, update_with=dict(bins=[hx[0], hx[0]]), keys=['bins'])
                axhistx.hist(hx, **_xhk)
            else:
                axhistx.hist(hx, **xhk)
        if len(hy) > 0:
            if len(hy) == 1:
                _yhk = _updatecopy(orig=yhk, update_with=dict(bins=[hy[0], hy[0]]), keys=['bins'])
                axhisty.hist(hy, **_yhk)
            else:
                axhisty.hist(hy, **yhk)

        # Turn off unnecessary labels -- for these, use the scatter's axes
        # labels
        for txt in axhisty.get_yticklabels() + axhistx.get_xticklabels():
            txt.set_visible(False)

        for txt in axhisty.get_xticklabels():
            txt.set_rotation(-90)

    def add_legends(self, xhists=True, yhists=False, scatter=True, **kwargs):
        """
        Add a legend to selected sets of axes.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `xhists` | `bool` | `True` | If `True`, include the top (x) marginal histogram axes. |
        | `yhists` | `bool` | `False` | If `True`, include the right (y) marginal histogram axes. |
        | `scatter` | `bool` | `True` | If `True`, include the scatter axes. |
        | `**kwargs` | — | — | Forwarded to `Axes.legend` for each selected axes. |
        """
        axs = []
        if xhists:
            axs.extend(self.top_hists)
        if yhists:
            axs.extend(self.right_hists)
        if scatter:
            axs.append(self.scatter_ax)

        for ax in axs:
            ax.legend(**kwargs)


def volcano_plot(
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
):
    """
    Generate a volcano plot with optional intensity-based alpha transparency.

    Points are colored by up/down/not-significant status based on
    `log2FoldChange_upper`/`_lower` and `nlog10_upper`, with dashed
    threshold lines drawn at each. If `text_nlog10`/`text_log2FoldChange`
    are given (or `highlight_genes` matches), qualifying genes are circled
    and labeled with leader lines to the left/right margins.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `df` | `pd.DataFrame` | — | DataFrame containing expression results. |
    | `padj_col` | `str` | — | Column name for adjusted p-values. |
    | `lfc_col` | `str` | — | Column name for log2 fold change values. |
    | `gene_id_col` | `str` | — | Column name for gene identifiers. |
    | `gene_name_col` | `str \| None` | `None` | Column name for gene names/labels. Falls back to `gene_id_col` if `None` or not present. |
    | `intensity_col` | `str \| None` | `None` | Column name for intensity values; when given, min-max scaled to `[0.08, 1.0]` and used as per-point scatter alpha. |
    | `log2FoldChange_upper` | `float` | `1.5` | Upper log2 fold change threshold (upregulation), drawn as a vertical dashed line. |
    | `log2FoldChange_lower` | `float` | `-1.5` | Lower log2 fold change threshold (downregulation), drawn as a vertical dashed line. |
    | `nlog10_upper` | `float` | `2.0` | `-log10(padj)` significance threshold, drawn as a horizontal dashed line. |
    | `text_nlog10` | `float \| None` | `None` | `-log10(padj)` threshold above which genes are auto-labeled. Requires `text_log2FoldChange` to also be set. |
    | `text_log2FoldChange` | `float \| None` | `None` | `\|log2FoldChange\|` threshold above which genes are auto-labeled. Requires `text_nlog10` to also be set. |
    | `title` | `str \| None` | `None` | Plot title. |
    | `highlight_genes` | `list[str] \| None` | `None` | Gene names (matched case-insensitively against `gene_name_col`) to always label, regardless of the `text_*` thresholds. |
    | `fontsize` | `int` | `12` | Base font size for labels, title, ticks, and gene annotations. |

    ## Returns

    `matplotlib.figure.Figure` — The generated volcano plot figure. If, after
    dropping rows with missing/zero p-values or fold changes, no data
    remains, a figure with a placeholder "no valid data" message is
    returned instead.

    ## Raises

    `ValueError` — If `padj_col`, `lfc_col`, or `gene_id_col` is missing
    from `df`, or if `intensity_col` is given but not present in `df`.
    """
    from sklearn.preprocessing import MinMaxScaler
    from adjustText import adjust_text
    from . import helpers

    df = df.copy()

    required = [padj_col, lfc_col, gene_id_col]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["gene_id"] = df[gene_id_col].astype("string")

    if gene_name_col is not None and gene_name_col in df.columns:
        df["gene_label"] = df[gene_name_col].astype("string")
    else:
        df["gene_label"] = df["gene_id"]

    df["log2FoldChange"] = pd.to_numeric(df[lfc_col], errors="coerce")
    padj = pd.to_numeric(df[padj_col], errors="coerce")
    padj = padj.replace(0, np.nan)
    df["nlog10"] = -np.log10(padj)

    df = df.dropna(subset=["log2FoldChange", "nlog10"]).reset_index(drop=True)

    if df.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        if title is not None:
            ax.set_title(title, fontsize=fontsize, weight="bold")
        ax.text(0.5, 0.5, "No valid data available", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return fig

    df["color_category"] = df[["log2FoldChange", "nlog10"]].apply(
        helpers.map_color,
        axis=1,
        down_log2FoldChange=log2FoldChange_lower,
        upper_log2FoldChange=log2FoldChange_upper,
        upper_nlog10=nlog10_upper,
    )

    if intensity_col is not None:
        if intensity_col not in df.columns:
            raise ValueError(f"Missing intensity column: '{intensity_col}'")

        intensity = pd.to_numeric(df[intensity_col], errors="coerce").fillna(0).to_numpy().reshape(-1, 1)
        alpha = MinMaxScaler().fit_transform(intensity).flatten()
        alpha = np.clip(alpha, 0.08, 1.0)
    else:
        alpha = np.ones(len(df), dtype=float)

    fig, ax = plt.subplots(figsize=(8, 7))

    color_map = {
        "increase": "red",
        "decrease": "blue",
        "not significant": "grey",
    }
    order = ["increase", "decrease", "not significant"]

    for cat in order:
        mask = (df["color_category"] == cat).to_numpy()

        if not mask.any():
            continue

        ax.scatter(
            df.loc[mask, "log2FoldChange"],
            df.loc[mask, "nlog10"],
            color=color_map[cat],
            alpha=alpha[mask],
            s=20,
        )

    ax.axhline(nlog10_upper, zorder=0, c="k", lw=2, ls="--")
    ax.axvline(log2FoldChange_upper, zorder=0, c="k", lw=2, ls="--")
    ax.axvline(log2FoldChange_lower, zorder=0, c="k", lw=2, ls="--")

    _max_label_chars_right = 0
    _max_label_chars_left = 0
    if text_nlog10 is not None and text_log2FoldChange is not None:
        mask_auto = (
            (df["nlog10"] > text_nlog10) &
            (df["log2FoldChange"].abs() > text_log2FoldChange)
        )
    else:
        mask_auto = np.zeros(len(df), dtype=bool)

    if highlight_genes is not None:
        highlight_genes_upper = [g.upper() for g in highlight_genes]
        mask_manual = df["gene_label"].str.upper().isin(highlight_genes_upper).to_numpy()
    else:
        mask_manual = np.zeros(len(df), dtype=bool)

    mask_labels = (mask_auto | mask_manual) & df["gene_label"].notna().to_numpy()
    df_labels = df.loc[mask_labels].copy()

    x_min, x_max = df["log2FoldChange"].min(), df["log2FoldChange"].max()
    y_min, y_max = df["nlog10"].min(), df["nlog10"].max()

    x_span = x_max - x_min if x_max > x_min else 1.0
    y_span = y_max - y_min if y_max > y_min else 1.0

    # Deliberately generous axis margin
    ax.set_xlim(x_min - 0.90 * x_span, x_max + 0.90 * x_span)
    ax.set_ylim(y_min - 0.05 * y_span, y_max + 0.20 * y_span)

    if not df_labels.empty:
        ax.scatter(
            df_labels["log2FoldChange"],
            df_labels["nlog10"],
            s=40,
            facecolors="none",
            edgecolors="black",
            linewidths=1.2,
            zorder=4,
        )

    x_text_right = x_max + 0.70 * x_span
    x_text_left = x_min - 0.70 * x_span

    df_labels_right = df_labels[df_labels["log2FoldChange"] >= 0].sort_values("nlog10")
    df_labels_left = df_labels[df_labels["log2FoldChange"] < 0].sort_values("nlog10")
    if not df_labels_right.empty:
        _max_label_chars_right = int(df_labels_right["gene_label"].str.len().max())
    if not df_labels_left.empty:
        _max_label_chars_left = int(df_labels_left["gene_label"].str.len().max())

    def add_side_annotations(df_side, x_text, ha):
        if df_side.empty:
            return

        n = len(df_side)
        y_positions = np.linspace(
            df_side["nlog10"].min(),
            df_side["nlog10"].max() + 0.18 * y_span,
            n
        )

        for (_, row), y_text in zip(df_side.iterrows(), y_positions):
            x = row["log2FoldChange"]
            y = row["nlog10"]
            label = str(row["gene_label"])

            ann = ax.annotate(
                label,
                xy=(x, y),
                xytext=(x_text, y_text),
                textcoords="data",
                ha=ha,
                va="center",
                fontsize=fontsize*0.8,
                weight="bold",
                zorder=5,
                annotation_clip=False,
                bbox=dict(
                    boxstyle="round,pad=0.25",
                    fc="white",
                    ec="black",
                    alpha=0.95,
                ),
                arrowprops=dict(
                    arrowstyle="-",
                    color="black",
                    lw=1.0,
                    shrinkA=0,
                    shrinkB=0,
                ),
            )
            ann.set_clip_on(False)

    add_side_annotations(df_labels_left, x_text_left, "right")
    add_side_annotations(df_labels_right, x_text_right, "left")

    legend_handles = [
        ax.scatter([], [], color=color_map["increase"], alpha=1.0, s=30, label="increase"),
        ax.scatter([], [], color=color_map["decrease"], alpha=1.0, s=30, label="decrease"),
        ax.scatter([], [], color=color_map["not significant"], alpha=1.0, s=30, label="not significant"),
    ]

    ax.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        bbox_transform=fig.transFigure,
        ncol=3,
        frameon=True,
        edgecolor="black",
        framealpha=1.0,
        fancybox=True,
        prop={"weight": "bold",
              "size": fontsize*0.7},
    )

    _char_w_inch = fontsize * 0.65 / 72
    _right_margin = max(_max_label_chars_right * _char_w_inch + 0.5, 0.5)
    _left_margin = max(_max_label_chars_left * _char_w_inch + 0.5, 0.5)
    _fw = fig.get_figwidth()
    plt.subplots_adjust(
        left=max(_left_margin / _fw, 0.12),
        right=min(1.0 - _right_margin / _fw, 0.90),
        bottom=0.2,
    )

    for spine in ["bottom", "left"]:
        ax.spines[spine].set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=2, labelsize=fontsize*0.8)
    
    ax.set_xlabel(r"$log_{2}$ fold change", size=fontsize)
    ax.set_ylabel(r"-$log_{10}$ FDR", size=fontsize)

    if title is not None:
        ax.set_title(title + " (n=" + str(len(df)) + ")", fontsize=fontsize)

    return fig


def hexbin_correlation_scatter(
    x_vals,
    y_vals,
    bins: int = 60,
    xlabel: str = "",
    ylabel: str = "",
    title: str = "",
    hex_low_color: str = "#D9D9D9",
    hex_high_color: str = "black",
    figsize: tuple[float, float] = (5, 4),
    fontsize: int = 12,
    annotation_fontsize_factor: float = 0.9,
    annotate_offset: float = 0.4,
):
    """
    Hexbin density scatter of two numeric arrays, annotated with Spearman R and n.

    Draws a hexbin density plot (with a colorbar for point counts), overlays
    a faint raw scatter for individual-point context, and annotates the
    Spearman correlation coefficient and sample size in the top-left corner.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x_vals`, `y_vals` | `array-like` | — | Values to plot on the x- and y-axes. |
    | `bins` | `int` | `60` | Hexbin gridsize. |
    | `xlabel`, `ylabel`, `title` | `str` | `""` | Axis labels and plot title. |
    | `hex_low_color`, `hex_high_color` | `str` | `"#D9D9D9"`, `"black"` | Colors for the low- and high-count ends of the hexbin colormap. |
    | `figsize` | `tuple[float, float]` | `(5, 4)` | Figure size, in inches. |
    | `fontsize` | `int` | `12` | Base font size for labels and title; ticks and the annotation scale from it. |
    | `annotation_fontsize_factor` | `float` | `0.9` | Multiplier of `fontsize` used for the R/n annotation text. |
    | `annotate_offset` | `float` | `0.4` | Offset, in data units, from the plot's min-x/max-y corner used to place the R/n annotation. |

    ## Returns

    `matplotlib.figure.Figure` — The hexbin scatter figure.
    """
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.stats import spearmanr

    x_vals = np.asarray(x_vals, dtype=float)
    y_vals = np.asarray(y_vals, dtype=float)
    n = len(x_vals)
    r, _ = spearmanr(x_vals, y_vals)

    cmap = LinearSegmentedColormap.from_list("", [hex_low_color, hex_high_color])

    fig, ax = plt.subplots(figsize=figsize)
    hb = ax.hexbin(x_vals, y_vals, gridsize=bins, cmap=cmap, mincnt=1)
    cbar = plt.colorbar(hb, ax=ax, label="count")
    cbar.ax.tick_params(labelsize=fontsize * 0.8)
    cbar.set_label("count", fontsize=fontsize * 0.9)
    ax.scatter(x_vals, y_vals, alpha=0.05, s=1.6, color="black", linewidths=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.tick_params(axis="both", labelsize=fontsize * 0.8)

    x_annot = float(np.nanmin(x_vals)) + annotate_offset
    y_annot = float(np.nanmax(y_vals)) - annotate_offset
    ax.text(
        x_annot, y_annot, f"R={r:.2f}\nn={n}",
        ha="left", va="top", fontsize=fontsize * annotation_fontsize_factor,
    )

    plt.tight_layout()
    return fig


def generate_gene_scale_regions_profile(
    genes_bed_path: str,
    signal_paths: list[str],
    sample_labels: list[str],
    colors: list[str],
    output_path: str | None = None,
    regions_label: str | None = None,
    before_region_start_length: int = 500,
    region_body_length: int = 1000,
    after_region_start_length: int = 500,
    bin_size: int = 10,
    skip_zeros: bool = True,
    processes: int | None = None,
    figsize: tuple[float, float] = (4.7244, 3.1496),
    dpi: int = 300,
    ylabel: str = "ChIP-seq Signal",
    title: str = "",
    legend_location: str = "upper right",
    show: bool = True,
) -> str | None:
    """
    Generate a deepTools scale-regions meta-profile figure.

    Runs ``computeMatrix scale-regions`` followed by ``plotProfile`` over a
    temporary directory that is cleaned up afterwards. This is the
    function-based counterpart to
    [`GeneScaleRegionsProfilePipeline`][Pyntegrate.pipelines.GeneScaleRegionsProfilePipeline] —
    see that class's `## Attributes` table for full parameter descriptions.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `genes_bed_path` | `str` | — | Path to a BED or BED.GZ file with annotated genes. |
    | `signal_paths` | `list[str]` | — | BigWig files to summarize, one line per sample. |
    | `sample_labels` | `list[str]` | — | Legend labels for each BigWig track. |
    | `colors` | `list[str]` | — | Line colors for each sample. |
    | `output_path` | `str \| None` | `None` | SVG/PNG file to save. If `None`, the figure is written to a temporary file and not returned by path. |
    | `regions_label` | `str \| None` | `None` | Passed as ``--regionsLabel`` to ``plotProfile``; `"{n}"` is replaced with the gene count. Defaults to `"Meta analysis of annotated genes. n={n}"`. |
    | `before_region_start_length` | `int` | `500` | Upstream flank length in bp (``--beforeRegionStartLength``). |
    | `region_body_length` | `int` | `1000` | Scaled gene-body length in bp (``--regionBodyLength``). |
    | `after_region_start_length` | `int` | `500` | Downstream flank length in bp (``--afterRegionStartLength``). |
    | `bin_size` | `int` | `10` | Bin size in bp (``--binSize``). |
    | `skip_zeros` | `bool` | `True` | Pass ``--skipZeros`` to ``computeMatrix``. |
    | `processes` | `int \| None` | `None` | ``--numberOfProcessors`` for ``computeMatrix``; treated as `1` if `None`. |
    | `figsize` | `tuple[float, float]` | `(4.7244, 3.1496)` | Figure size in inches; converted to cm for ``plotProfile``. |
    | `dpi` | `int` | `300` | ``--dpi`` for ``plotProfile``. |
    | `ylabel` | `str` | `"ChIP-seq Signal"` | ``--yAxisLabel`` for ``plotProfile``. |
    | `title` | `str` | `""` | ``--plotTitle`` for ``plotProfile``. |
    | `legend_location` | `str` | `"upper right"` | ``--legendLocation`` for ``plotProfile``; spaces are converted to hyphens. |
    | `show` | `bool` | `True` | Open the saved figure with ``xdg-open`` after saving. |

    ## Returns

    `str | None` — Absolute path of the saved figure, or `None` when `output_path`
    was not provided.

    ## Raises

    `FileNotFoundError` — If the ``computeMatrix`` or ``plotProfile``
    deepTools executables cannot be found on `$PATH` or alongside the
    current Python interpreter.
    """
    import os
    import subprocess
    import shutil
    import tempfile
    import pathlib
    import sys

    def _find_tool(name: str) -> str:
        path = shutil.which(name)
        if path:
            return path
        candidate = pathlib.Path(sys.executable).parent / name
        if candidate.exists():
            return str(candidate)
        raise FileNotFoundError(
            f"{name!r} not found. Install deepTools and make sure it is on "
            "PATH, or install it in the same virtual environment as Pyntegrate."
        )

    compute_matrix = _find_tool("computeMatrix")
    plot_profile = _find_tool("plotProfile")

    tmp_dir = tempfile.mkdtemp(prefix="pyntegrate_gsr_")
    try:
        # computeMatrix does not accept .bed.gz — decompress on the fly
        if genes_bed_path.endswith(".gz"):
            bed_path = os.path.join(tmp_dir, "genes.bed")
            with open(bed_path, "w") as fh:
                subprocess.run(
                    ["gunzip", "-c", genes_bed_path],
                    stdout=fh,
                    check=True,
                )
        else:
            bed_path = genes_bed_path

        with open(bed_path) as _f:
            n_genes = sum(1 for ln in _f if ln.strip() and not ln.startswith("#"))
        # Print number of genes detected in the BED file
        try:
            print(f"computeMatrix: se van a usar {n_genes} genes desde '{genes_bed_path}'")
        except Exception:
            print("computeMatrix: no se pudo determinar el número de genes en el archivo BED")

        matrix_path = os.path.join(tmp_dir, "matrix.gz")

        compute_cmd = [
            compute_matrix, "scale-regions",
            "-S", *signal_paths,
            "-R", bed_path,
            "--beforeRegionStartLength", str(before_region_start_length),
            "--regionBodyLength", str(region_body_length),
            "--afterRegionStartLength", str(after_region_start_length),
            "--binSize", str(bin_size),
            "-o", matrix_path,
            "--numberOfProcessors", str(processes or 1),
        ]
        if skip_zeros:
            compute_cmd.append("--skipZeros")

        subprocess.run(compute_cmd, check=True)

        if output_path:
            out_dir = os.path.dirname(output_path)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            final_output = os.path.abspath(output_path)
        else:
            final_output = os.path.join(tmp_dir, "profile.svg")

        # deepTools plotProfile uses cm; convert from inches
        plot_width_cm = round(figsize[0] * 2.54, 4)
        plot_height_cm = round(figsize[1] * 2.54, 4)

        plot_cmd = [
            plot_profile,
            "-m", matrix_path,
            "--perGroup",
            "--plotType", "lines",
            "--colors", *colors,
            "--samplesLabel", *sample_labels,
            "--plotWidth", str(plot_width_cm),
            "--plotHeight", str(plot_height_cm),
            "--averageType", "mean",
            "--yAxisLabel", ylabel,
            "--legendLocation", legend_location.replace(" ", "-"),
            "--plotTitle", title,
            "--dpi", str(dpi),
            "--outFileName", final_output,
        ]
        final_regions_label = (
            regions_label or "Meta analysis of annotated genes. n={n}"
        ).format(n=f"{n_genes:,}")
        plot_cmd.extend(["--regionsLabel", final_regions_label])

        subprocess.run(plot_cmd, check=True)

        if show:
            subprocess.Popen(["xdg-open", final_output])

        return final_output if output_path else None

    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
