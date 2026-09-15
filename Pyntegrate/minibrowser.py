"""
Spawn mini genome-browser figures for a single feature, using a plugin-style layout.

Point a mini-browser at some data, then call its `plot()` method with a
`pybedtools.Interval`. This spawns a new figure showing that data over
the feature's interval. Each `panels()` override maps a list of created
`Axes` to the methods that draw on them, so a subclass is free to
combine as many data panels (signal tracks, gene models, peaks, ...) as
needed on one figure.

## Classes

| Class | Description |
|-------|-------------|
| `BaseMiniBrowser` | Abstract base: figure/axes bookkeeping, `plot()`/`panels()` contract. |
| `SignalMiniBrowser` | Plots one panel of local coverage per genomic signal object. |
| `ChIPSeqMiniBrowser` | IP vs. control coverage, peaks, and a gene-model track for one region. |
| `GeneModelMiniBrowser` | Adds a gene-model (exon/intron) track below a signal panel. |
| `PeakMiniBrowser` | Adds a peak-interval track below a signal panel. |

## Quick start

```python
from Pyntegrate.minibrowser import SignalMiniBrowser
from Pyntegrate import genomic_signal

ip = genomic_signal("ip.bw", "bigwig")
control = genomic_signal("control.bw", "bigwig")
browser = SignalMiniBrowser([ip, control])
browser.plot(feature)
```

## Writing a new subclass

1. Add one or more methods that accept an `Axes` and a `pybedtools.Interval`
   and return the feature (optionally after modifying it), e.g.:

   ```python
   def my_panel(self, ax, feature):
       gs = self.genomic_signal_objs[0]
       x, y = gs.local_coverage(feature, bins=100)
       ax.plot(x, y)
       ax.axis('tight')
       return feature
   ```

2. Override `panels()` to create `Axes` on `self.fig` (already created by
   `make_fig()`) and return a list of `(ax, method)` tuples mapping each
   `Axes` to the panel method that should draw on it:

   ```python
   def panels(self):
       ax = self.fig.add_subplot(111)
       return [(ax, self.my_panel)]
   ```
"""
from copy import copy
from matplotlib import pyplot as plt
from pybedtools.contrib.plotting import Track
import pybedtools
import gffutils
from gffutils.helpers import asinterval
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
from ._genomic_signal import BamSignal


class BaseMiniBrowser(object):
    """
    Base class for plotting a genomic region.

    This class is designed to be sub-classed, so it really just a shell with
    some example methods filled in. Subclasses typically override `panels()`
    to lay out one or more Axes on `self.fig` and map each one to a method
    that draws on it.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `genomic_signal_objs` | `list` | Genomic signal objects passed at construction; available to panel methods. |
    | `fig` | `matplotlib.figure.Figure` | Figure created by `make_fig()`, immediately before `panels()` is called. |
    | `_all_figures` | `list` | Per-instance list of every figure spawned by `make_fig()`; used by `close_all()`. |
    """
    def __init__(self, genomic_signal_objs):
        self.genomic_signal_objs = genomic_signal_objs
        self._all_figures = []

    def close_all(self):
        """
        Close all figures spawned by this instance.
        """
        for f in self._all_figures:
            plt.close(f)
        self._all_figures = []

    def plot(self, feature):
        """
        Spawn a new figure showing data for `feature`.

        Converts `feature` to a `pybedtools.Interval` if it's a
        `gffutils.Feature`, calls `self.make_fig()`, then calls each
        `(ax, method)` pair returned by `self.panels()` in turn, passing the
        (possibly updated) feature along the chain.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `feature` | `pybedtools.Interval \\| gffutils.Feature` | — | Genomic region to plot. `gffutils.Feature` objects are converted to `pybedtools.Interval` automatically. |

        ## Returns

        `list[matplotlib.axes.Axes]` — The Axes objects created by `panels()`, in the order they were plotted.
        """
        if isinstance(feature, gffutils.Feature):
            feature = asinterval(feature)
        self.make_fig()
        axes = []
        for ax, method in self.panels():
            feature = method(ax, feature)
            axes.append(ax)
        return axes

    def make_fig(self):
        """
        Create a new figure and register it for `close_all()`.

        Called automatically by `plot()`, before `panels()`. Sets `self.fig`
        to a new `matplotlib.figure.Figure` and appends it to this
        instance's `_all_figures` list.
        """
        self.fig = plt.figure(figsize=(8, 4))
        self._all_figures.append(self.fig)

    def panels(self):
        """
        Method to be overridden by subclasses.

        This method should create Axes objects on self.fig and register which
        methods operate on which axes.

        Each method's signature should be func(self, ax, feature), and it will
        take a pybedtools.Interval `feature` and plot it or something derived
        from it onto `ax`. It should return a pybedtools.Interval object,
        usually just the same one that was passed in.

        ## Returns

        `list[tuple[matplotlib.axes.Axes, callable]]` — `(ax, method)` pairs mapping each created Axes to the method that draws on it.
        """
        ax = self.fig.add_subplot(111)
        return [(ax, self.example_panel)]

    def example_panel(self, ax, feature):
        """
        Draw the coordinates of `feature` as text on `ax`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw on. |
        | `feature` | `pybedtools.Interval` | — | Feature whose `chrom:start-stop` text is drawn. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        txt = '%s:%s-%s' % (feature.chrom, feature.start, feature.stop)
        ax.text(0.5, 0.5, txt, transform=ax.transAxes)
        return feature


class ChIPSeqMiniBrowser(BaseMiniBrowser):
    """
    Mini-browser to show IP and control signal, along with optional peaks
    and gene models.

    Builds a stacked-panel figure: IP signal on top, control signal below
    it, and — depending on which optional arguments are given — a peaks
    track and/or a gene-models track underneath. See the module-level
    docstring for details on the mini-browser plugin architecture.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `local_coverage_kwargs` | `dict` | `dict(stranded=False)` | Kwargs forwarded to `ip.local_coverage()` / `control.local_coverage()`. Keep `stranded=False` since the mini-browser always displays the plus strand. |
    | `ip` | genomic signal object | — | IP (treatment) signal; scaled by million mapped reads if it's a `BamSignal`. |
    | `control` | genomic signal object | — | Control (input) signal; scaled the same way as `ip`. |
    | `db` | `gffutils.FeatureDB \\| None` | `None` | Annotation database. Accepts a path string at construction, which is converted to a `FeatureDB`. If set, a gene-models panel is added. |
    | `ip_style` | `dict` | `None` | Kwargs passed to `Axes.fill_between` when drawing the IP panel (e.g. colors, alpha, line styles). |
    | `peaks_style` | `dict` | `None` | Kwargs passed to the `Track`/`PolyCollection` when drawing the peaks panel (e.g. `face_color`, `edgewidth`). |
    | `control_style` | `dict` | `None` | Kwargs passed to `Axes.fill_between` when drawing the control panel. |
    | `peaks` | `pybedtools.BedTool \\| None` | `None` | Called peaks. Accepts a path string or `BedTool` at construction. If set, a peaks panel is added. |
    | `max_bins` | `int` | `500` | Maximum number of bins requested from `local_coverage`, unless overridden by a `bins` key in `local_coverage_kwargs`. |
    | `settings` | `dict` | — | Kwargs forwarded to `gffutils.contrib.plotting.Gene` when drawing gene models (`transcripts`, `cds`, `utrs`, `color`, `featuretype`). Built internally in `__init__`; can be edited between calls to `plot()`. |
    | `genomic_signal_objs` | `list` | — | Set by the base class `__init__` to `[ip, control]`. |
    | `axes` | `dict` | — | Panel name (`'ip'`, `'control'`, `'peaks'`, `'genes'`) mapped to its `Axes`; built by `panels()` based on which of `db`/`peaks` are set. |
    | `gridspec` | `matplotlib.gridspec.GridSpec` | — | Grid layout for the panels, built by `panels()`. Can be used to tweak relative panel sizes. |
    | `_first_ax` | `matplotlib.axes.Axes` | — | The IP panel's Axes; shared x-limits are set on it in `plot()`. |
    | `_last_ax` | `matplotlib.axes.Axes` | — | The bottom-most panel's Axes; gets the x-axis label and tick formatting in `plot()`. |
    | `_current_feature` | `pybedtools.Interval \\| None` | — | The feature most recently passed to `plot()`; used by `coords()`. |

    ## Notes

    After creation, settings can be changed in between calls to the
    `plot()` method.

    For example, the `settings` dictionary, which contains arguments
    passed to the `gffutils.contrib.plotting.Gene` class, can be modified
    to specify which kinds of features should be plotted. The `gridspec`
    attribute contains the `matplotlib.gridspec.GridSpec` object, which
    can be used to modify the relative sizes of the subplots.
    """
    def __init__(self, ip, control, db=None, peaks=None,
                 local_coverage_kwargs=None, ip_style=None,
                 control_style=None, peaks_style=None, max_bins=500):
        super(ChIPSeqMiniBrowser, self).__init__([ip, control])
        self.local_coverage_kwargs = local_coverage_kwargs or dict(stranded=False)
        self.ip = ip
        self.control = control
        if isinstance(db, str):
            db = gffutils.FeatureDB(db)
        self.db = db
        self.ip_style = ip_style or {}
        self.peaks_style = peaks_style or {}
        self.control_style = control_style or {}

        if peaks is not None:
            peaks = pybedtools.BedTool(peaks)
        self.peaks = peaks
        self.max_bins = max_bins

        self.settings = {
            'transcripts': None,
            'cds': ['CDS'],
            'utrs': ['exon'],
            'color': '0.5',
            'featuretype': 'gene',
        }

    def panels(self):
        """
        Build the panel layout and Axes for this browser instance.

        Chooses one of four grid layouts depending on whether `self.db`
        and/or `self.peaks` are set (IP + control, optionally with a peaks
        track and/or a gene-models track), creates the shared-x Axes with
        `matplotlib.gridspec.GridSpec`, and stores them on `self.axes` and
        `self.gridspec`.

        ## Returns

        `list[tuple[matplotlib.axes.Axes, callable]]` — `(ax, method)` pairs in panel order (`ip`, `control`, `peaks`, `genes`, as applicable).
        """
        self.fig.set_facecolor('w')
        method_dispatch = {
            'ip': self.ip_panel,
            'control': self.control_panel,
            'peaks': self.peaks_panel,
            'genes': self.genes_panel}
        if self.db and self.peaks:
            gs = gridspec.GridSpec(4, 1, height_ratios=[1, 1, .15, .5])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            peaks_ax = plt.subplot(gs[2], sharex=ip_ax)
            gene_ax = plt.subplot(gs[3], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'peaks': peaks_ax,
                    'genes': gene_ax,
                    }

        elif self.db and self.peaks is None:
            gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, .5])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            gene_ax = plt.subplot(gs[2], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'genes': gene_ax,
                    }

        elif self.db is None and self.peaks:
            gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, .3])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            peaks_ax = plt.subplot(gs[2], sharex=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    'peaks': peaks_ax,
                    }

        elif self.db is None and self.peaks is None:
            gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1])
            ip_ax = plt.subplot(gs[0])
            control_ax = plt.subplot(gs[1], sharex=ip_ax, sharey=ip_ax)
            axes = {'ip': ip_ax,
                    'control': control_ax,
                    }

        self.axes = axes
        self.gridspec = gs

        # Order shouldn't matter, but just in case...
        mapping = []
        for kind in ['ip', 'control', 'peaks', 'genes']:
            if kind in axes:
                mapping.append((axes[kind], method_dispatch[kind]))

        self._first_ax = mapping[0][0]
        self._last_ax = mapping[-1][0]

        return mapping

    def _bins(self, feature):
        """
        Number of signal bins to request for `feature`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `feature` | `pybedtools.Interval` | — | Feature whose length is compared against `self.max_bins`. |

        ## Returns

        `int` — `min(len(feature), self.max_bins)`.
        """
        return min(len(feature), self.max_bins)

    def _zoomed_feature(self, feature):
        """
        Return a copy of `feature` expanded outward by a zoom fraction.

        Expands `feature` on both ends by
        `self.settings['zoom'] * len(feature) / 2` base pairs.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `feature` | `pybedtools.Interval` | — | Feature to expand. |

        ## Returns

        `pybedtools.Interval` — A copy of `feature` with `start`/`stop` moved outward.
        """
        extra = int(self.settings['zoom'] * len(feature) / 2)
        feature = copy(feature)
        feature.start -= extra
        feature.stop += extra
        return feature

    def ip_panel(self, ax, feature):
        """
        Plot the IP signal for `feature` on `ax`.

        Computes local coverage from `self.ip`, normalizes by million
        mapped reads if `self.ip` is a `BamSignal`, and fills the area
        under the curve on `ax` using `self.ip_style`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the IP panel on. |
        | `feature` | `pybedtools.Interval` | — | Region to compute and plot signal over. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        bins = self._bins(feature)
        x, y = self.ip.local_coverage(
            feature, bins=bins, **self.local_coverage_kwargs)
        if isinstance(self.ip, BamSignal):
            y /= (self.ip.mapped_read_count() / 1e6)
        ax.fill_between(x, y, y2=0, **self.ip_style)
        ax.axis('tight')
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        return feature

    def control_panel(self, ax, feature):
        """
        Plot the control signal for `feature` on `ax`.

        Computes local coverage from `self.control`, normalizes by
        million mapped reads if `self.control` is a `BamSignal`, and
        fills the area under the curve on `ax` using `self.control_style`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the control panel on. |
        | `feature` | `pybedtools.Interval` | — | Region to compute and plot signal over. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        bins = self._bins(feature)
        x, y = self.control.local_coverage(
            feature, bins=bins, **self.local_coverage_kwargs)
        if isinstance(self.control, BamSignal):
            y /= (self.control.mapped_read_count() / 1e6)
        ax.fill_between(x, y, y2=0, **self.control_style)
        ax.axis('tight')
        ax.spines['right'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        return feature

    def peaks_panel(self, ax, feature):
        """
        Plot peaks overlapping `feature` on `ax`.

        Intersects `self.peaks` with `feature`, draws the hits as a
        `pybedtools.contrib.plotting.Track`, and hides the axis
        frame/ticks, leaving only a `"Peaks"` label.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the peaks panel on. |
        | `feature` | `pybedtools.Interval` | — | Region used to select overlapping peaks. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        hits = self.peaks.intersect([feature], u=True)
        track = Track(hits, **self.peaks_style)
        ax.add_collection(track)
        for x in (
            ax.get_yticklines() + ax.get_xticklines() + ax.get_yticklabels()
            + ax.get_xticklabels()
        ):
            x.set_visible(False)
        ax.set_ylabel('Peaks', rotation=0, horizontalalignment='right',
                      verticalalignment='center')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('none')
        ax.set_frame_on(False)

        return feature

    def plot(self, feature):
        """
        Spawn a new figure with IP/control (and optional peaks/genes) panels for `feature`.

        Calls the base class `plot()`, hides the x-tick labels on every
        panel except the bottom-most one, ties the x-limits of the shared
        axes to `feature`, and labels the x-axis with `feature.chrom`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `feature` | `pybedtools.Interval` | — | Genomic region to plot. Stored on `self._current_feature` for later use by `coords()`. |

        ## Returns

        `list[matplotlib.axes.Axes]` — The Axes objects created by `panels()`, in the order they were plotted.
        """
        self._current_feature = feature
        axes = super(ChIPSeqMiniBrowser, self).plot(feature)
        for ax in self.fig.axes:
            if ax is not self._last_ax:
                for txt in ax.get_xticklabels():
                    txt.set_visible(False)
        self.fig.subplots_adjust(wspace=0.05)
        self._first_ax.axis(xmin=feature.start, xmax=feature.stop)
        self._last_ax.set_xlabel(feature.chrom)
        self._last_ax.xaxis.set_major_formatter(FormatStrFormatter("%d"))
        return axes

    def coords(self):
        """
        Current visible genomic coordinates of the browser.

        Reads the x-axis limits of `self._first_ax`, which may have
        changed from `feature`'s original extent if the plot was zoomed
        or panned interactively after `plot()` was called.

        ## Returns

        `str` — `"chrom:start-stop"` for the region currently shown on the IP panel.
        """
        chrom = self._current_feature.chrom
        start, stop, _, _ = self._first_ax.axis()
        return '%s:%s-%s' % (chrom, int(start), int(stop))

    def genes_panel(self, ax, feature):
        """
        Plot gene models near `feature` on `ax`, using `self.db`.

        Queries `self.db` for genes overlapping `feature`, draws each as
        a `gffutils.contrib.plotting.Gene` (highlighting the gene whose
        `ID` matches `feature['ID']`, if any, in a darker color), and
        returns a new interval spanning `feature` plus all nearby genes.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the gene models on. |
        | `feature` | `pybedtools.Interval` | — | Region used to query nearby genes. |

        ## Returns

        `pybedtools.Interval` — A new interval covering `feature` and all nearby genes drawn (strand set to `'.'`).
        """
        from gffutils.contrib.plotting import Gene
        extent = [feature.start, feature.stop]
        nearby_genes = self.db.region(
            (feature.chrom, feature.start, feature.stop), featuretype='gene')
        ybase = 0
        ngenes = 0
        for nearby_gene in nearby_genes:
            # TODO: there should be a better way of identifying which gene is
            # the same as the feature requested.  Might need to expose an "ID"
            # kwarg.
            try:
                if nearby_gene['ID'][0] == feature['ID']:
                    color = '0.2'
                else:
                    color = '0.5'
            except KeyError:
                color = '0.5'
            ngenes += 1
            extent.extend([nearby_gene.start, nearby_gene.stop])
            gene_collection = Gene(
                self.db,
                nearby_gene,
                transcripts=None,
                cds=['CDS'],
                utrs=['exon'],
                ybase=ybase,
                color=color)
            gene_collection.name = nearby_gene.id
            gene_collection.add_to_ax(ax)
            ybase += gene_collection.max_y

        xmin = min(extent)
        xmax = max(extent)

        # 1% padding seems to work well
        padding = (xmax - xmin) * 0.01
        ax.axis('tight')

        # add lines indicating extent of current feature
        vline_kwargs = dict(color='k', linestyle='--')
        ax.axvline(feature.start, **vline_kwargs)
        ax.axvline(feature.stop, **vline_kwargs)

        # Make a new feature to represent the region plus surrounding genes
        interval = pybedtools.create_interval_from_list(feature.fields)
        interval.start = xmin - padding
        interval.stop = xmax + padding
        interval.strand = '.'
        for txt in ax.get_yticklabels():
            txt.set_visible(False)
        for tick in ax.get_yticklines():
            tick.set_visible(False)
        ax.set_ylabel('Genes')
        ax.spines['right'].set_color('None')
        ax.spines['left'].set_color('None')
        ax.spines['top'].set_color('None')
        ax.yaxis.set_ticks_position('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.set_ylabel('Genes', rotation=0, horizontalalignment='right',
                      verticalalignment='center')

        return interval


class SignalMiniBrowser(BaseMiniBrowser):
    """
    Base class for plotting genomic signal.

    Plots genomic signal over a particular area of the genome as a single
    line panel. Designed to be extended by subclasses (see
    `GeneModelMiniBrowser`, `PeakMiniBrowser`), but can stand alone if all
    you want is a simple one-panel browser.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `genomic_signal_objs` | `list` | — | Genomic signal objects (e.g. `Pyntegrate.genomic_signal.BamSignal` instances) to plot; stored by the base class `__init__`. |
    | `local_coverage_kwargs` | `dict` | `None` | Kwargs sent to each signal object's `local_coverage()` method, e.g. `dict(fragment_size=200)`. |
    | `plotting_kwargs` | `list[dict]` | `None` | One kwargs dict per entry in `genomic_signal_objs`, forwarded to `Axes.plot()` in `signal_panel()`, e.g. `[dict(color='r', label='IP'), dict(color='k', label='input')]`. |
    """
    def __init__(self, genomic_signal_objs, local_coverage_kwargs=None,
                 plotting_kwargs=None):
        super(SignalMiniBrowser, self).__init__(genomic_signal_objs)
        self.plotting_kwargs = \
            plotting_kwargs or [{} for i in genomic_signal_objs]
        self.local_coverage_kwargs = local_coverage_kwargs or {}

    def panels(self):
        """
        Add a single panel to the figure.

        ## Returns

        `list[tuple[matplotlib.axes.Axes, callable]]` — A single `(ax, self.signal_panel)` pair.
        """
        ax = self.fig.add_subplot(111)
        return [(ax, self.signal_panel)]

    def signal_panel(self, ax, feature):
        """
        Plot each genomic signal in `self.genomic_signal_objs` as a line on `ax`.

        Each signal's `local_coverage()` is called with
        `self.local_coverage_kwargs`, and plotted with the corresponding
        entry from `self.plotting_kwargs`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the signal lines on. |
        | `feature` | `pybedtools.Interval` | — | Region to compute and plot signal over. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        for gs, kwargs in zip(self.genomic_signal_objs, self.plotting_kwargs):
            x, y = gs.local_coverage(feature, **self.local_coverage_kwargs)
            ax.plot(x, y, **kwargs)
        ax.axis('tight')
        return feature


class GeneModelMiniBrowser(SignalMiniBrowser):
    """
    Mini-browser to show a signal panel on top and gene models on the
    bottom.

    Extends `SignalMiniBrowser` by adding a `gene_panel()` method and
    overriding `panels()` to create two stacked, x-shared Axes: one for
    the signal (top) and one for gene models (bottom).

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `genomic_signal_objs` | `list` | — | Genomic signal objects to plot on the top panel; stored by the base class `__init__`. |
    | `db` | `gffutils.FeatureDB` | — | Annotation database queried by `gene_panel()` for nearby genes. |
    | `local_coverage_kwargs` | `dict` | `None` | Inherited from `SignalMiniBrowser`; kwargs forwarded to each signal object's `local_coverage()` method. Set via `**kwargs`. |
    | `plotting_kwargs` | `list[dict]` | `None` | Inherited from `SignalMiniBrowser`; one dict per entry in `genomic_signal_objs`, forwarded to `Axes.plot()`. Set via `**kwargs`. |
    """
    def __init__(self, genomic_signal_objs, db, **kwargs):
        super(GeneModelMiniBrowser, self).__init__(
            genomic_signal_objs, **kwargs)
        self.db = db

    def panels(self):
        """
        Add 2 panels to the figure: signal on top, gene models on the bottom.

        ## Returns

        `tuple[tuple[matplotlib.axes.Axes, callable], ...]` — `(ax2, self.gene_panel)` and `(ax1, self.signal_panel)`, sharing the x-axis.
        """
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        return (ax2, self.gene_panel), (ax1, self.signal_panel)

    def gene_panel(self, ax, feature):
        """
        Plot gene models near `feature` on `ax`, using `self.db`.

        Queries `self.db` for genes overlapping `feature`, draws each as
        a `gffutils.contrib.plotting.Gene`, adds dashed vertical lines
        marking the extent of `feature`, and returns a new interval
        padded 1% beyond the combined extent of `feature` and the
        nearby genes.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the gene models on. |
        | `feature` | `pybedtools.Interval` | — | Region used to query nearby genes. |

        ## Returns

        `pybedtools.Interval` — A new, padded interval covering `feature` and all nearby genes (strand set to `'.'`).
        """
        from gffutils.contrib.plotting import Gene
        extent = [feature.start, feature.stop]
        nearby_genes = self.db.region(
            (feature.chrom, feature.start, feature.stop), featuretype='gene')
        ybase = 0
        ngenes = 0
        for nearby_gene in nearby_genes:
            ngenes += 1
            extent.extend([nearby_gene.start, nearby_gene.stop])
            gene_collection = Gene(
                self.db,
                nearby_gene,
                transcripts=['mRNA'],
                cds=['CDS'],
                utrs=['exon'],
                ybase=ybase,
                color="0.5", picker=5)
            gene_collection.name = nearby_gene.id
            gene_collection.add_to_ax(ax)
            ybase += gene_collection.max_y

        xmin = min(extent)
        xmax = max(extent)

        # 1% padding seems to work well
        padding = (xmax - xmin) * 0.01
        ax.axis('tight')

        # add lines indicating extent of current feature
        vline_kwargs = dict(color='k', linestyle='--')
        ax.axvline(feature.start, **vline_kwargs)
        ax.axvline(feature.stop, **vline_kwargs)

        # Make a new feature to represent the region plus surrounding genes
        interval = pybedtools.create_interval_from_list(feature.fields)
        interval.start = xmin - padding
        interval.stop = xmax + padding
        interval.strand = '.'
        return interval


class PeakMiniBrowser(SignalMiniBrowser):
    """
    Signal on the top panel, peaks on the bottom panel.

    `bed` is a filename or `pybedtools.BedTool` object (or a list of
    these) used to build a `pybedtools.contrib.plotting.Track` in
    `peak_panel()`.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `genomic_signal_objs` | `list` | — | Genomic signal objects to plot on the top panel; stored by the base class `__init__`. |
    | `bed` | `str \\| pybedtools.BedTool \\| list` | — | Peaks/intervals drawn on the bottom panel by `peak_panel()`. |
    | `local_coverage_kwargs` | `dict` | `None` | Inherited from `SignalMiniBrowser`; kwargs forwarded to each signal object's `local_coverage()` method. Set via `**kwargs`. |
    | `plotting_kwargs` | `list[dict]` | `None` | Inherited from `SignalMiniBrowser`; one dict per entry in `genomic_signal_objs`, forwarded to `Axes.plot()`. Set via `**kwargs`. |
    """
    def __init__(self, genomic_signal_objs, bed, **kwargs):
        super(PeakMiniBrowser, self).__init__(genomic_signal_objs, **kwargs)
        self.bed = bed

    def panels(self):
        """
        Add 2 panels to the figure: signal on top, peaks on the bottom.

        ## Returns

        `tuple[tuple[matplotlib.axes.Axes, callable], ...]` — `(ax1, self.signal_panel)` and `(ax2, self.peak_panel)`, sharing the x-axis.
        """
        ax1 = self.fig.add_subplot(211)
        ax2 = self.fig.add_subplot(212, sharex=ax1)
        return (ax1, self.signal_panel), (ax2, self.peak_panel)

    def peak_panel(self, ax, feature):
        """
        Plot peaks/intervals from `self.bed` overlapping `feature` on `ax`.

        Intersects `self.bed` with `feature` and draws the hits as a
        `pybedtools.contrib.plotting.Track`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ax` | `matplotlib.axes.Axes` | — | Axes to draw the peaks panel on. |
        | `feature` | `pybedtools.Interval` | — | Region used to select overlapping peaks/intervals. |

        ## Returns

        `pybedtools.Interval` — The same `feature` that was passed in.
        """
        bedtool = pybedtools.BedTool(self.bed)
        features = bedtool.intersect([feature], u=True)
        track = Track(features)
        ax.add_collection(track)
        return feature


if __name__ == "__main__":
    import Pyntegrate
    import gffutils
    import pybedtools

    G = gffutils.FeatureDB(
        Pyntegrate.example_filename('Homo_sapiens.GRCh37.66.cleaned.gtf.db'))

    ip = Pyntegrate.genomic_signal(
        Pyntegrate.example_filename('wgEncodeUwTfbsK562CtcfStdAlnRep1.bam'),
        'bam')
    inp = Pyntegrate.genomic_signal(
        Pyntegrate.example_filename('wgEncodeUwTfbsK562InputStdAlnRep1.bam'),
        'bam')
    peaks = pybedtools.BedTool(Pyntegrate.example_filename(
        'wgEncodeUwTfbsK562CtcfStdPkRep1.narrowPeak.gz'))

    plotting_kwargs = [
        dict(color='r', label='IP'),
        dict(color='k', linestyle=':', label='input')]

    local_coverage_kwargs = dict(fragment_size=200)

    b = SignalMiniBrowser(
        [ip, inp],
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    g = GeneModelMiniBrowser(
        [ip, inp],
        G,
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    p = PeakMiniBrowser(
        [ip, inp],
        peaks,
        plotting_kwargs=plotting_kwargs,
        local_coverage_kwargs=local_coverage_kwargs)

    feature = peaks[3]

    b.plot(feature)
    g.plot(feature)
    p.plot(feature)
    plt.show()
