"""
DataFrame-backed wrappers around differential-expression and general tabular results.

`ResultsTable` wraps a `pandas.DataFrame` (typically loaded from a
delimited results file) and adds genomics-aware convenience methods —
annotated scatter plots, gffutils-based feature lookups, boolean
set operations between tables, and RadViz-style visualization — while
still forwarding any undefined attribute access straight to the
underlying DataFrame. `DifferentialExpressionResults` and its
subclasses add conventions for locating the p-value/log-fold-change/
mean-expression columns produced by common DE tools.

## Classes

| Class | Description |
|-------|-------------|
| `ResultsTable` | Base wrapper: DataFrame passthrough, `scatter()`, feature/db integration, set operations. |
| `DifferentialExpressionResults` | Adds `pval_column`/`lfc_column`/`mean_column` conventions, `upregulated()`/`downregulated()`, `volcano()`. |
| `EdgeRResults` | Column-name overrides for edgeR output (`FDR`, `logFC`, `logCPM`). |
| `DESeqResults` | Column-name overrides for DESeq output. |
| `DEseq2Results` | Column-name overrides for DESeq2 output, plus gene-ID version stripping. |
| `LazyDict` | Lazily-instantiated, dict-like cache of `ResultsTable`-like objects keyed by filename. |

## Quick start

```python
from Pyntegrate.results_table import DEseq2Results

res = DEseq2Results("deseq2_results.tsv")
res.scatter(x="baseMean", y="log2FoldChange")
```
"""
from textwrap import dedent
import numpy as np
import pandas
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib

from .plotutils import *
from .plotutils import _updatecopy
from .helpers import *
from matplotlib.transforms import blended_transform_factory
from matplotlib.collections import EventCollection
import gffutils
import pybedtools
from pybedtools import featurefuncs

_base_intro = """%s

The underlying `pandas.DataFrame` is always accessible via the `data`
attribute. Any attribute not explicitly defined by this class is
forwarded to the underlying DataFrame transparently.
"""

_base_doc = _base_intro + """
## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | `str \| pandas.DataFrame` | — | Path to a delimited file, or an existing DataFrame. If a string, it is loaded with `pandas.read_table(data, delimiter=delimiter, **import_kwargs)`. |
| `delimiter` | `str` | `"\\t"` | Field delimiter used by `pandas.read_table()` when `data` is a filename. |
| `db` | `str \| gffutils.FeatureDB \| None` | `None` | Optional gffutils database used by `features()`, `five_prime()`, and other genomic-feature methods. |
| `import_kwargs` | `dict \| None` | `None` | Extra keyword arguments forwarded to `pandas.read_table()` when `data` is a filename. |
"""


class ResultsTable(object):
    __doc__ = _base_intro % dedent(
        """
        Wrapper around a `pandas.DataFrame` with genomic-analysis extras.

        Provides transparent attribute forwarding to the underlying
        DataFrame, along with methods for attaching a gffutils database,
        generating genomic feature intervals, filtering by gene sets,
        and producing annotated scatter plots.

        ## Attributes

        | Attribute | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `data` | `str \| pandas.DataFrame` | — | Path to a delimited file, or an existing DataFrame. If a string, it is loaded with `pandas.read_table(data, delimiter=delimiter, **import_kwargs)` and the result is stored here. Any attribute not defined on this class is forwarded to it via `__getattr__`. |
        | `delimiter` | `str` | `"\\t"` | Field delimiter used by `pandas.read_table()` when `data` is a filename. Constructor-only; not retained as an attribute. |
        | `db` | `str \| gffutils.FeatureDB \| None` | `None` | Path to a `.db` file or a `gffutils.FeatureDB`. Stored as `self.db`, attached via `attach_db()` (called from `__init__`). Enables `features()`, `five_prime()`, `three_prime()`, and `genes_with_peak()`. |
        | `import_kwargs` | `dict \| None` | `None` | Extra keyword arguments forwarded to `pandas.read_table()` when `data` is a filename. Kept, along with `db`, inside `_kwargs`. |
        | `_kwargs` | `dict` | — | Constructor kwargs (`db`, `import_kwargs`) remembered so `update()`, `copy()`, and `align_with()` can reconstruct an instance of the same class. |
        | `_cached_features` | `pybedtools.BedTool \| None` | — | Cache of per-gene features built from `features()`; populated on first call to `genes_with_peak()`. |
        | `marginal` | `MarginalHistScatter` | — | Set by `scatter()`; manages the plot's optional marginal histograms. Only present after `scatter()` has been called. |
        | `_seen` | `numpy.ndarray` | — | Boolean array created by `scatter()` to avoid invoking the click callback more than once for overlapping highlighted genes. |
        """)
    def __init__(self, data, delimiter="\t", db=None, import_kwargs=None):
        if isinstance(data, str):
            import_kwargs = import_kwargs or {}
            data = pandas.read_table(data,delimiter=delimiter,**import_kwargs, low_memory=False)
        if not isinstance(data, pandas.DataFrame):
            raise ValueError("`data` is not a pandas.DataFrame")
        self.data = data

        self._kwargs = dict(db=db, import_kwargs=import_kwargs)
        self.attach_db(db)
        self._cached_features = None

    def __getattr__(self, attr):
        """
        Forward attribute lookups not found on this instance to `self.data`.
        """
        if attr in self.__dict__:
            return getattr(self, attr)
        return getattr(self.data, attr)

    def __getitem__(self, attr):
        """
        Index into the wrapped DataFrame.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `attr` | `str \| Any` | — | Column label, boolean mask, slice, or other pandas indexer. |

        ## Returns

        `pandas.Series | ResultsTable` — If `attr` is a string (a
        column label), returns the underlying `pandas.Series` from
        `self.data`. Otherwise, returns a new instance of the same
        class wrapping `self.data[attr]`, preserving the attached
        database and other constructor kwargs.
        """
        if isinstance(attr, str):
            return self.data.__getitem__(attr)
        else:
            return self.__class__(self.data.__getitem__(attr), **self._kwargs)

    def update(self, dataframe):
        """
        Return a new instance of the same class wrapping `dataframe`.

        Required because pandas fancy indexing (`.iloc`, `.loc`, etc.)
        returns a plain `DataFrame`, not a `ResultsTable`. Use this
        method to preserve the class and its attached database.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `dataframe` | `pandas.DataFrame` | — | Replacement DataFrame to wrap. |

        ## Returns

        A new instance of the same class (or subclass) wrapping `dataframe`.

        ## Examples

        ```python
        subset_df = rt.data[rt.data["baseMean"] > 100]
        rt_subset = rt.update(subset_df)
        print(type(rt_subset).__name__)
        ```

        **Output:**
        ```
        ResultsTable
        ```
        """
        return self.__class__(dataframe, **self._kwargs)

    def copy(self):
        """
        Return a new instance of the same class wrapping a copy of `self.data`.

        For `DEseq2Results` instances, the DESeq2-specific
        configuration (`pval_column`, `lfc_column`, `mean_column`,
        `id_column`, `strip_version_ids`, `id_norm_column`, `verbose`,
        `db`) is read off `self` (via `getattr`, falling back to
        `None`/defaults) and passed explicitly to the constructor. For
        all other classes, the new instance is constructed with
        `db=self.db` and `import_kwargs=self._kwargs`.

        ## Returns

        A new instance of `self.__class__` wrapping a copy of `self.data`.
        """
        data = self.data.copy()

        if self.__class__.__name__ == "DEseq2Results":
            return self.__class__(
                data,
                pval_column=getattr(self, "pval_column", None),
                lfc_column=getattr(self, "lfc_column", None),
                mean_column=getattr(self, "mean_column", None),

                id_column=getattr(self, "id_column", None),
                strip_version_ids=getattr(self, "strip_version_ids", True),
                id_norm_column=getattr(self, "id_norm_column", "GeneID_norm"),
                verbose=getattr(self, "verbose", False),

                db=getattr(self, "db", None),
            )

        import_kwargs = getattr(self, "_kwargs", None)
        return self.__class__(data, db=getattr(self, "db", None), import_kwargs=import_kwargs)


    def __repr__(self):
        """
        Return a string showing the class name and the wrapped DataFrame's repr, indented.
        """
        s = []
        s.append("<%s instance, wrapping the following:"
                 % self.__class__.__name__)
        s.append('')
        s.extend('\t' + i for i in repr(self.data).splitlines(False))
        s.append('>')
        return '\n'.join(s)

    def attach_db(self, db):
        """
        Attach a `gffutils.FeatureDB` to enable genomic-feature methods.

        Can be called after the instance is created to (re-)attach a
        database. Required by `features()`, `five_prime()`,
        `three_prime()`, and `genes_with_peak()`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `db` | `str \| gffutils.FeatureDB \| None` | — | Path to a `.db` file, an existing `FeatureDB`, or `None` to detach. |

        ## Raises

        `ValueError` — If `db` is neither `None`, a string, nor a
        `gffutils.FeatureDB`.

        ## Examples

        ```python
        from Pyntegrate.results_table import ResultsTable

        rt = ResultsTable("deseq2.csv", delimiter=",")
        rt.attach_db("gencode.hg19.db")
        print(type(rt.db).__name__)
        ```

        **Output:**
        ```
        FeatureDB
        ```
        """
        if db is not None:
            if isinstance(db, str):
                db = gffutils.FeatureDB(db)
            if not isinstance(db, gffutils.FeatureDB):
                raise ValueError(
                    "`db` must be a filename or a gffutils.FeatureDB")
        self._kwargs['db'] = db
        self.db = db

    def features(self, ignore_unknown=False):
        """
        Yield a `pybedtools.Interval` for each gene in the DataFrame index.

        Requires a `gffutils.FeatureDB` to be attached (see `attach_db`).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `ignore_unknown` | `bool` | `False` | If `True`, silently skip genes not found in the database. If `False`, raise `FeatureNotFoundError`. |

        ## Yields

        `pybedtools.Interval` — Genomic interval for each gene in the index.

        ## Raises

        `ValueError` — If no database is attached.
        `gffutils.FeatureNotFoundError` — If a gene is not in the database
        and `ignore_unknown=False`.

        ## Examples

        ```python
        from Pyntegrate.results_table import ResultsTable

        rt = ResultsTable("deseq2.csv", delimiter=",", db="gencode.hg19.db")
        for feature in rt.features(ignore_unknown=True):
            print(feature.chrom, feature.start, feature.stop)
            break
        ```

        **Output:**
        ```
        chr1 11868 14409
        ```
        """
        if not self.db:
            raise ValueError("Please attach a gffutils.FeatureDB")
        for i in self.data.index:
            try:
                yield gffutils.helpers.asinterval(self.db[i])
            except gffutils.FeatureNotFoundError:
                if ignore_unknown:
                    continue
                else:
                    raise gffutils.FeatureNotFoundError('%s not found' % i)

    def reindex_to(self, x, attribute="Name"):
        """
        Return a copy whose rows match the feature names/IDs found in `x`, in that order.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `x` | `str \| pybedtools.BedTool` | — | BED/GFF/GTF/VCF file (or `BedTool`) whose features supply the names to reindex to; each feature's `attribute` field (default `"Name"`) becomes one row label. |
        | `attribute` | `str` | `"Name"` | Attribute/field on each feature of `x` containing the name to use as the index. |

        ## Returns

        A new instance of the same class, copied from `self` and then
        reindexed (via `pandas.DataFrame.reindex`) to the names
        extracted from `x`. Rows with no match are filled with NaN;
        rows in `self` not present in `x` are dropped.
        """
        names = [i[attribute] for i in x]
        new = self.copy()
        new.data = new.data.reindex(names)
        return new

    def five_prime(self, upstream=1, downstream=0):
        """
        Return a `pybedtools.BedTool` of the 5′ end of each feature.

        Constructs a BED/GFF file where each entry is the 5′ end of the
        corresponding gene, extended by `upstream` and `downstream` bp.
        Requires an attached `gffutils.FeatureDB`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `upstream` | `int` | `1` | Base pairs upstream of the 5′ end to include. |
        | `downstream` | `int` | `0` | Base pairs downstream of the 5′ end to include. |

        ## Returns

        `pybedtools.BedTool` — Temporary BedTool of 5′ end intervals.
        """
        return pybedtools.BedTool(self.features())\
            .each(featurefuncs.five_prime, upstream, downstream)\
            .saveas()

    def three_prime(self, upstream=0, downstream=1):
        """
        Return a `pybedtools.BedTool` of the 3′ end of each feature.

        Constructs a BED/GFF file where each entry is the 3′ end of the
        corresponding gene, extended by `upstream` and `downstream` bp.
        Requires an attached `gffutils.FeatureDB`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `upstream` | `int` | `0` | Base pairs upstream of the 3′ end to include. |
        | `downstream` | `int` | `1` | Base pairs downstream of the 3′ end to include. |

        ## Returns

        `pybedtools.BedTool` — Temporary BedTool of 3′ end intervals.
        """
        return pybedtools.BedTool(self.features())\
            .each(featurefuncs.three_prime, upstream, downstream)\
            .saveas()

    TSS = five_prime
    TTS = three_prime

    def align_with(self, other):
        """
        Return a copy of this table reindexed to match another object's index.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `other` | `pandas.DataFrame \| pandas.Series \| ResultsTable` | — | Object supplying the index to align to; passed to `pandas.DataFrame.reindex_like()`. |

        ## Returns

        A new instance of the same class wrapping `self.data.reindex_like(other)`.
        """
        return self.__class__(self.data.reindex_like(other), **self._kwargs)

    def genes_in_common(self, other):
        """
        Return the index labels present in both this table and `other`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `other` | `ResultsTable \| pandas.DataFrame` | — | Object whose `.index` is intersected with `self.index`. |

        ## Returns

        `pandas.Index` — Labels found in both `self.index` and `other.index`.
        """
        return self.index & other.index

    def __and__(self, other):
        """Return `self.index & other.index` (labels present in both)."""
        return self.index & other.index

    def __or__(self, other):
        """Return `self.index | other.index` (labels present in either)."""
        return self.index | other.index

    def __sub__(self, other):
        """Return `self.index - other.index` (labels only in `self`)."""
        return self.index - other.index

    def __len__(self):
        """Return the number of rows in `self.data`."""
        return len(self.data)

    def scatter(self, x, y, xfunc=None, yfunc=None,
                xlab=None, ylab=None, genes_to_highlight=None,
                marginal_histograms=False,
                general_kwargs=dict(color="k", alpha=0.2, picker=True),
                general_hist_kwargs=None, offset_kwargs={}, label_kwargs=None,
                ax=None, one_to_one=None, callback=None, xlab_prefix=None,
                ylab_prefix=None, hist_size=0.3, hist_pad=0.0,
                nan_offset=0.015, pos_offset=0.99, linelength=0.01,
                neg_offset=0.005, figure_kwargs=None):
        """
        Do-it-all method for making annotated scatterplots.

        For each of `x` and `y` (column names in `self.data`), values
        are optionally transformed and NaN/+inf/-inf values are drawn
        as "rug" tick marks along the plot edges instead of being
        dropped. Genes can be highlighted in groups with their own
        colors/labels, and optional marginal histograms and a 1:1
        reference line can be added. The resulting collections are
        picker-enabled and connected to `callback` (or
        `self._default_callback`) via `_id_callback`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `x`, `y` | `str` | — | Column names in `self.data` to plot on the x- and y-axes, e.g. `"baseMeanA"`, `"baseMeanB"`. |
        | `xfunc`, `yfunc` | `callable \| None` | `None` | Function applied to the `x`/`y` values before plotting (e.g. `np.log2`). `None` means no transformation. |
        | `xlab`, `ylab` | `str \| None` | `None` | Axis labels. Default is built from `xlab_prefix`/`ylab_prefix`, the column name, and `xfunc.__name__`/`yfunc.__name__` (if given), e.g. `"log2(baseMeanA)"`. |
        | `genes_to_highlight` | `list[tuple] \| None` | `None` | List of `(ind, kwargs)` (or `(ind, kwargs, hist_kwargs)`) tuples, where `ind` selects rows (boolean array or integer index) and `kwargs` updates a copy of `general_kwargs` for those points. `kwargs` may also carry `names` (labels, drawn with `label_kwargs`), `marginal_histograms`, `xhist_kwargs`, and `yhist_kwargs`. |
        | `marginal_histograms` | `bool` | `False` | If `True` (or if a `genes_to_highlight` block sets its own `marginal_histograms=True`), draw marginal histograms for that group of points. |
        | `general_kwargs` | `dict` | `dict(color="k", alpha=0.2, picker=True)` | `matplotlib.pyplot.scatter` kwargs applied to all non-highlighted points. |
        | `general_hist_kwargs` | `dict \| None` | `None` | Kwargs for the marginal histogram of non-highlighted points; `color`/`alpha` are copied from `general_kwargs` if not set. |
        | `offset_kwargs` | `dict` | `{}` | Kwargs for `matplotlib.transforms.offset_copy`, controlling the placement of gene-name labels relative to their points. |
        | `label_kwargs` | `dict \| None` | `None` | `ax.text` kwargs for gene labels; only used for `genes_to_highlight` blocks that supply `names`. Defaults to a right-aligned, italic label with a semi-transparent white box. |
        | `ax` | `matplotlib.axes.Axes \| None` | `None` | Axes to plot onto. If `None`, a new figure (via `figure_kwargs`) and subplot are created. |
        | `one_to_one` | `dict \| None` | `None` | If not `None`, `matplotlib.pyplot.plot` kwargs used to draw a 1:1 reference line spanning the shared data range. |
        | `callback` | `callable \| None` | `None` | Called with a gene ID when a point is clicked. Defaults to `self._default_callback` (prints the row). |
        | `xlab_prefix`, `ylab_prefix` | `str \| None` | `None` | Prepended to the auto-generated `xlab`/`ylab`. |
        | `hist_size` | `float` | `0.3` | Size of the marginal histograms (forwarded to `MarginalHistScatter`). |
        | `hist_pad` | `float` | `0.0` | Spacing between the marginal histograms and the scatter axes. |
        | `nan_offset`, `pos_offset`, `neg_offset` | `float` | `0.015`, `0.99`, `0.005` | Position, in axes-fraction units, of the NaN / +inf / -inf "rug" tick marks. |
        | `linelength` | `float` | `0.01` | Line length of the rug tick marks. |
        | `figure_kwargs` | `dict \| None` | `None` | Kwargs for `plt.figure()` when `ax` is `None`. |

        ## Returns

        `matplotlib.axes.Axes` — The axes the scatterplot was drawn on.

        ## Notes

        Sets `self.marginal` (a `MarginalHistScatter`) and
        `self._seen` (a dedup mask for the click callback) as side
        effects.
        """
        _x = self.data[x]
        _y = self.data[y]

        # Construct defaults---------------------------------------------------
        def identity(x):
            return x.copy()

        # Axis label setup
        if xlab_prefix is None:
            xlab_prefix = ""

        if ylab_prefix is None:
            ylab_prefix = ""

        if xlab is None:
            xlab = x

            if xfunc is not None:
                xlab = xlab_prefix + "%s(%s)" % (xfunc.__name__, str(x))
            else:
                xlab = xlab_prefix + "%s" % (str(x))

        if ylab is None:
            ylab = y
            if yfunc is not None:
                ylab = ylab_prefix + "%s(%s)" % (yfunc.__name__, str(y))
            else:
                ylab = ylab_prefix + "%s" % (str(y))

        if xfunc is None:
            xfunc = identity

        if yfunc is None:
            yfunc = identity

        if general_kwargs is None:
            general_kwargs = {}

        if general_hist_kwargs is None:
            general_hist_kwargs = {}

        if genes_to_highlight is None:
            genes_to_highlight = []

        if ax is None:
            if figure_kwargs is None:
                figure_kwargs = {}
            fig = plt.figure(**figure_kwargs)
            ax = fig.add_subplot(111)

        if label_kwargs is None:
            label_kwargs = dict(
                horizontalalignment='right',
                verticalalignment='center',
                style='italic',
                bbox=dict(facecolor='w', edgecolor='None', alpha=0.5)
            )

        # Clean data ---------------------------------------------------------
        xi = xfunc(_x)
        yi = yfunc(_y)

        # handle inf, -inf, and NaN
        x_is_pos_inf = np.isinf(xi) & (xi > 0)
        x_is_neg_inf = np.isinf(xi) & (xi < 0)
        x_is_nan = np.isnan(xi)
        y_is_pos_inf = np.isinf(yi) & (yi > 0)
        y_is_neg_inf = np.isinf(yi) & (yi < 0)
        y_is_nan = np.isnan(yi)

        # Indexes for valid values
        x_valid = ~(x_is_pos_inf | x_is_neg_inf | x_is_nan)
        y_valid = ~(y_is_pos_inf | y_is_neg_inf | y_is_nan)

        # global min/max
        gmin = max(xi[x_valid].min(), yi[y_valid].min())
        gmax = min(xi[x_valid].max(), yi[y_valid].max())

        # Convert any integer indexes into boolean, and create a new list of
        # genes to highlight.  This handles optional hist kwargs.

        # We'll compile a new list of genes to highlight.
        _genes_to_highlight = []

        for block in genes_to_highlight:
            ind = block[0]

            # Convert to boolean
            if ind.dtype != 'bool':
                new_ind = (np.zeros_like(xi) == 0)
                new_ind[ind] = True
                _genes_to_highlight.append(
                    tuple([new_ind] + list(block[1:]))
                )

            # If it's a DataFrame, we only want the boolean values;
            else:
                if hasattr(ind, 'values'):
                    ind = ind.values
                _genes_to_highlight.append(
                    tuple([ind] + list(block[1:]))
                )

        # Now we remove any genes from in allind (which will be plotted using
        # `general_kwargs`) that will be plotted by genes_to_highlight.  This
        # avoids double-plotting.
        allind = np.zeros_like(xi) == 0
        for block in _genes_to_highlight:
            ind = block[0]
            allind[ind] = False

        # Copy over the color and alpha if they're not specified
        general_hist_kwargs = _updatecopy(
            orig=general_hist_kwargs, update_with=general_kwargs,
            keys=['color', 'alpha'])

        # Put the non-highlighted genes at the beginning of _genes_to_highlight
        # list so we can just iterate over one list
        _genes_to_highlight.insert(
            0,
            (allind, general_kwargs, general_hist_kwargs)
        )

        # Set up the object that will handle the marginal histograms
        self.marginal = MarginalHistScatter(
            ax, hist_size=hist_size, pad=hist_pad)

        # Set up kwargs for x and y rug plots
        rug_x_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transData, ax.transAxes),
        )
        rug_y_kwargs = dict(
            linelength=linelength,
            transform=blended_transform_factory(ax.transAxes, ax.transData),
            orientation='vertical',
        )

        # EventCollection objects need a color as a 3-tuple, so set up
        # a converter here.
        color_converter = matplotlib.colors.ColorConverter().to_rgb

        # Plot the one-to-one line, if kwargs were specified
        if one_to_one:
            ax.plot([gmin, gmax],
                    [gmin, gmax],
                    **one_to_one)

        # Plot 'em all, and label if specified

        # In order to avoid calling the callback function multiple times when
        # we have overlapping genes to highlight (e.g., a gene that is both
        # upregulated AND has a peak), keep track of everything that's been
        # added so far.
        self._seen = np.ones_like(xi) == 0
        for block in _genes_to_highlight:
            ind = block[0]
            kwargs = block[1]

            if len(block) == 3:
                hist_kwargs = block[2]
            else:
                hist_kwargs = {}

            names = kwargs.pop('names', None)
            _marginal_histograms = (
                kwargs.pop('marginal_histograms', False) or
                marginal_histograms)

            updated_kwargs = _updatecopy(
                orig=kwargs, update_with=general_kwargs)

            updated_hist_kwargs = _updatecopy(
                orig=hist_kwargs, update_with=general_hist_kwargs)
            updated_hist_kwargs = _updatecopy(
                orig=updated_hist_kwargs, update_with=kwargs,
                keys=['color', 'alpha'], override=True)

            xhist_kwargs = updated_kwargs.pop('xhist_kwargs', None)
            yhist_kwargs = updated_kwargs.pop('yhist_kwargs', None)

            self.marginal.append(
                xi[ind & x_valid & y_valid],
                yi[ind & x_valid & y_valid],
                scatter_kwargs=dict(**updated_kwargs),
                hist_kwargs=updated_hist_kwargs,
                xhist_kwargs=xhist_kwargs,
                yhist_kwargs=yhist_kwargs,
                marginal_histograms=_marginal_histograms,
            )

            # This is important for callbacks: here we grab the last-created
            # collection,
            coll = self.marginal.scatter_ax.collections[-1]
            coll.df = self.data
            coll.ind = ind & x_valid & y_valid

            color = color_converter(updated_kwargs['color'])
            rug_x_kwargs['color'] = color
            rug_y_kwargs['color'] = color

            # Note: if both x and y are not valid, then they will not be on the
            # plot.
            items = [
                # top rug, y is +inf and x is valid
                (xi, ind & x_valid & y_is_pos_inf, pos_offset, rug_x_kwargs),

                # one of the bottom rugs, where y is NaN
                (xi, ind & x_valid & y_is_nan, nan_offset, rug_x_kwargs),

                # bottom rug, y is -inf
                (xi, ind & x_valid & y_is_neg_inf, neg_offset, rug_x_kwargs),

                # right rug, x is +inf
                (yi, ind & y_valid & x_is_pos_inf, pos_offset, rug_y_kwargs),

                # one of the left rugs; x is NaN
                (yi, ind & y_valid & x_is_nan, nan_offset, rug_y_kwargs),

                # left rug, x is -inf
                (yi, ind & y_valid & x_is_neg_inf, neg_offset, rug_y_kwargs),
            ]
            for values, index, offset, kwargs in items:
                coll = EventCollection(
                    values[index], lineoffset=offset, **kwargs)
                coll.df = self.data
                coll.ind = index
                ax.add_collection(coll)

            if names:
                transOffset = matplotlib.transforms.offset_copy(
                    ax.transData, fig=ax.figure, **offset_kwargs)

                for xii, yii, name in zip(xi[ind], yi[ind], names):
                    ax.text(xii,
                            yii,
                            name,
                            transform=transOffset,
                            **label_kwargs)

        # register callback
        if callback is None:
            callback = self._default_callback

        def wrapped_callback(event):
            for _id in self._id_callback(event):
                callback(_id)

        ax.figure.canvas.mpl_connect('pick_event', wrapped_callback)

        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.axis('tight')

        return ax

    def radviz(self, column_names, transforms=dict(), **kwargs):
        """
        Plot a radviz (radial visualization) projection of multivariate data.

        Useful for exploratory visualization: the variables in
        `column_names` are distributed evenly around the unit circle.
        Each row of `self.data` is treated as attached to each
        variable by a spring whose stiffness is proportional to that
        row's (normalized) value for the variable; a row's plotted
        position is the equilibrium position of all its springs
        combined. Each column is normalized to `[0, 1]` by subtracting
        its minimum and dividing by its range before plotting.

        This is an exploratory plot only — the order of `column_names`
        affects the result, so it's worth trying a few different
        orderings. See [1] for other caveats.

        Adds two columns, `"radviz_x"` and `"radviz_y"`, to
        `self.data`, then calls `self.scatter()` with these as `x`/`y`.
        Additional `kwargs` are forwarded to `self.scatter()` (e.g.
        `genes_to_highlight`), though not every `scatter()` option is
        meaningful for a radviz plot (marginal histograms, for
        instance, are not).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `column_names` | `list[str]` | — | Numeric columns of `self.data` to use; must not contain NaN/inf/-inf values. |
        | `transforms` | `dict` | `{}` | Maps a column name to a transformation applied to it before normalization (e.g. `np.log1p`). Columns not present are used as-is. |
        | `**kwargs` | — | — | Forwarded to `self.scatter()`. |

        ## Returns

        `matplotlib.axes.Axes` — The axes with the radviz plot (a
        circle outline, one anchor point per variable, and the
        projected data points).

        ## Notes

        Adapted from the `pandas.tools.plotting.radviz` function.

        ## References

        1. Hoffman, P.E. et al. (1997) DNA visual and analytic data
           mining. In the Proceedings of the IEEE Visualization.
           Phoenix, AZ, pp. 437-441.
        2. http://www.agocg.ac.uk/reports/visual/casestud/brunsdon/radviz.htm
        3. http://pandas.pydata.org/pandas-docs/stable/visualization.html#radviz
        """
        # make a copy of data
        x = self.data[column_names].copy()

        for k, v in transforms.items():
            x[k] = v(x[k])

        def normalize(series):
            mn = min(series)
            mx = max(series)
            return (series - mn) / (mx - mn)

        df = x.apply(normalize)

        to_plot = []

        n = len(column_names)

        s = np.array([(np.cos(t), np.sin(t))
                      for t in [2.0 * np.pi * (i / float(n))
                                for i in range(n)]])

        for i in range(len(x)):
            row = df.iloc[i].values
            row_ = np.repeat(np.expand_dims(row, axis=1), 2, axis=1)
            to_plot.append((s * row_).sum(axis=0) / row.sum())

        x_, y_ = zip(*to_plot)
        self.data['radviz_x'] = x_
        self.data['radviz_y'] = y_

        ax = self.scatter('radviz_x', 'radviz_y', **kwargs)

        ax.add_patch(patches.Circle((0.0, 0.0), radius=1.0, facecolor='none'))
        for xy, name in zip(s, column_names):
            ax.add_patch(patches.Circle(xy, radius=0.025, facecolor='gray'))
            if xy[0] < 0.0 and xy[1] < 0.0:
                ax.text(xy[0] - 0.025, xy[1] - 0.025, name,
                        ha='right', va='top', size='small')
            elif xy[0] < 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] - 0.025, xy[1] + 0.025, name,
                        ha='right', va='bottom', size='small')
            elif xy[0] >= 0.0 and xy[1] < 0.0:
                ax.text(xy[0] + 0.025, xy[1] - 0.025, name,
                        ha='left', va='top', size='small')
            elif xy[0] >= 0.0 and xy[1] >= 0.0:
                ax.text(xy[0] + 0.025, xy[1] + 0.025, name,
                        ha='left', va='bottom', size='small')

        ax.axis('equal')
        return ax

    def _id_callback(self, event):
        """
        Yield the dataframe index label(s) for a matplotlib pick `event`.

        `event.ind` indexes into the picked artist's x/y data;
        `event.artist.ind` (set in `scatter()`) maps the artist back
        to a boolean/integer index into `self.data`.
        """
        # event.ind is the index into event's x and y data.
        #
        # event.artist.ind is the index of the entire artist into the original
        # dataframe.
        subset_df = event.artist.df.iloc[event.artist.ind]
        for i in event.ind:
            _id = subset_df.index[i]
            yield _id

    def _default_callback(self, i):
        """Print the row of `self.data` at index label `i`."""
        print (self.data.loc[i])

    def strip_unknown_features(self):
        """
        Remove rows whose index is not found in the attached `gffutils.FeatureDB`.

        Typically removes rows like `'ambiguous'`, `'no_feature'`,
        etc., but can also help when the database was built from a
        different annotation than the one used to create this table.

        ## Returns

        A new instance of the same class containing only rows found
        in `self.db`. Returns `self` unchanged if no database is
        attached.
        """
        if not self.db:
            return self
        ind = []
        for i, gene_id in enumerate(self.data.index):
            try:
                self.db[gene_id]
                ind.append(i)
            except gffutils.FeatureNotFoundError:
                pass
        ind = np.array(ind)
        return self.__class__(self.data.iloc[ind], **self._kwargs)

    def genes_with_peak(self, peaks, transform_func=None, split=False,
                        intersect_kwargs=None, id_attribute='ID', *args,
                        **kwargs):
        """
        Return a boolean index of genes that have a peak nearby.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `peaks` | `str \| pybedtools.BedTool` | — | Filename of a BED/GFF/GTF file of intervals, or a `pybedtools.BedTool` directly. |
        | `transform_func` | `callable \| None` | `None` | Applied to each gene feature yielded by `self.features()` (via `.each()`, or `.split()` if `split=True`) before intersecting with `peaks`. E.g. `pybedtools.featurefuncs.TSS` with `upstream=1000, downstream=0` to look for peaks near TSSs. May return an iterable of features (e.g. a gene's exons); each returned feature's `name` field must match a row label in `self.data`. Additional `*args`/`**kwargs` are passed through to `transform_func`. |
        | `split` | `bool` | `False` | If `True` and `transform_func` is given, use `BedTool.split()` instead of `.each()` — needed when `transform_func` yields multiple features per gene. |
        | `intersect_kwargs` | `dict \| None` | `None` | Kwargs passed to `pybedtools.BedTool.intersect()`. |
        | `id_attribute` | `str` | `'ID'` | GFF/GTF attribute holding the gene ID that should match the DataFrame index. Typically `"ID"` for GFF, `"gene_id"` for GTF. |
        | `*args`, `**kwargs` | — | — | Forwarded to `transform_func`. |

        ## Returns

        `numpy.ndarray[bool]` — Boolean mask over `self.data.index`,
        `True` for genes with at least one intersecting peak.

        ## Notes

        Caches the per-gene `pybedtools.BedTool` built from
        `self.features()` in `self._cached_features` on first call.
        """
        def _transform_func(x):
            """
            In order to support transform funcs that return a single feature or
            an iterable of features, we need to wrap it
            """
            result = transform_func(x)
            if isinstance(result, pybedtools.Interval):
                result = [result]
            for i in result:
                if i:
                    yield result

        intersect_kwargs = intersect_kwargs or {}
        if not self._cached_features:
            self._cached_features = pybedtools\
                .BedTool(self.features())\
                .saveas()

        if transform_func:
            if split:
                features = self._cached_features\
                    .split(_transform_func, *args, **kwargs)
            else:
                features = self._cached_features\
                    .each(transform_func, *args, **kwargs)

        else:
            features = self._cached_features

        hits = list(set([i[id_attribute] for i in features.intersect(
            peaks, **intersect_kwargs)]))
        return self.data.index.isin(hits)
    def delete_genes_no_peaks_chip(self,peaks):
        """
        Drop, in place, rows whose gene has an all-zero peak signal.

        For each entry in `peaks`, averages `entry['values']`; if the
        average is exactly `0`, `entry['gene_name']` is treated as
        having no peak and rows in `self.data` with a matching
        `'gene_name'` are dropped. `self.data` is reassigned to the
        filtered DataFrame.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `peaks` | `iterable[dict]` | — | Iterable of mappings, each with a `'gene_name'` key and a `'values'` key (an iterable of numeric peak signal values) for one gene. |

        ## Returns

        `None` — `self.data` is filtered in place.
        """
        new_df = self.data
        gene_no_peaks = []
        for x in peaks:
            cont = 0
            for num in x['values']:
                cont = cont + num
            avg = cont / len(x['values'])

            if avg == 0:
                gene_no_peaks.append(x['gene_name'])
        

        for gene in gene_no_peaks:
            if gene in new_df['gene_name']:
                new_df = new_df[new_df['gene_name'] != gene]
        self.data = new_df


class DifferentialExpressionResults(ResultsTable):

    __doc__ = _base_doc % dedent("""
    A `ResultsTable` subclass for differential expression results.

    Extends `ResultsTable` with significance filters (`upregulated`,
    `downregulated`, `changed`, `unchanged`) and an MA plot. Column
    names for the adjusted p-value, log2 fold-change, and mean
    expression are stored as class attributes and can be overridden by
    subclasses.

    ## Class attributes

    | Attribute | Default | Description |
    |-----------|---------|-------------|
    | `pval_column` | `"padj"` | Column used for significance filtering. |
    | `lfc_column` | `"log2FoldChange"` | Column with log2 fold-change values. |
    | `mean_column` | `"baseMean"` | Column with mean expression values (for MA plot). |

    The constructor also accepts `header_check` (default `True`): if
    `data` is a filename, autodetect and skip any leading comment
    lines (those not starting with `import_kwargs.get('comment', '#')`)
    before reading the file with `pandas.read_table()`.
    """)
    pval_column = 'padj'
    lfc_column = 'log2FoldChange'
    mean_column = 'baseMean'

    def __init__(self, data, db=None, header_check=True, **kwargs):
        import_kwargs = kwargs.pop('import_kwargs', {})
        if header_check and isinstance(data, str):
            comment_char = import_kwargs.get('comment', '#')
            for i, line in enumerate(open(data)):
                if line[0] != comment_char:
                    break
            import_kwargs['skiprows'] = i
        import_kwargs['na_values'] = ['nan']

        import_kwargs['index_col'] = import_kwargs.pop('index_col', 0)
        super(DifferentialExpressionResults, self).__init__(
            data=data, db=db, import_kwargs=import_kwargs, **kwargs)

    def changed(self, thresh=0.05, idx=True):
        """
        Return genes with adjusted p-value ≤ `thresh` (changed in either direction).

        {threshdoc}
        """
        ind = self.data[self.pval_column] <= thresh
        if idx:
            return ind
        return self[ind]

    def unchanged(self, thresh=0.05, idx=True):
        """
        Return genes with adjusted p-value > `thresh` or NaN (not significantly changed).

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] > thresh)
            | np.isnan(self.data[self.pval_column])
        )
        if idx:
            return ind
        return self[ind]

    def enriched(self, thresh=0.05, idx=True):
        """
        Alias for `upregulated`. Return significantly enriched/up-regulated genes.

        {threshdoc}
        """
        return self.upregulated(thresh=thresh, idx=idx)

    def upregulated(self, thresh=0.05, idx=True):
        """
        Return significantly up-regulated genes (padj ≤ `thresh` **and** LFC > 0).

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] <= thresh)
            & (self.data[self.lfc_column] > 0)
        )
        if idx:
            return ind
        return self[ind]

    def downregulated(self, thresh=0.05, idx=True):
        """
        Return significantly down-regulated genes (padj ≤ `thresh` **and** LFC < 0).

        {threshdoc}
        """
        ind = (
            (self.data[self.pval_column] <= thresh)
            & (self.data[self.lfc_column] < 0)
        )
        if idx:
            return ind
        return self[ind]

    def disenriched(self, thresh=0.05, idx=True):
        """
        Alias for `downregulated`. Return significantly disenriched/down-regulated genes.

        {threshdoc}
        """
        return self.downregulated(thresh=thresh, idx=idx)

    def ma_plot(self, thresh, up_kwargs=None, dn_kwargs=None,
                zero_line=None, **kwargs):
        """
        Produce an MA plot (mean expression vs log2 fold-change).

        Plots `mean_column` on the x-axis (log scale) and `lfc_column`
        on the y-axis. Significantly up- and down-regulated genes are
        coloured separately. Additional kwargs are forwarded to
        `self.scatter`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `thresh` | `float` | — | Adjusted p-value threshold; genes ≤ `thresh` are highlighted. |
        | `up_kwargs` | `dict \| None` | `None` | `matplotlib.scatter` kwargs for up-regulated points. Default: red. |
        | `dn_kwargs` | `dict \| None` | `None` | `matplotlib.scatter` kwargs for down-regulated points. Default: blue. |
        | `zero_line` | `dict \| None` | `None` | If not `None`, kwargs for `ax.axhline(0)`. |
        | `**kwargs` | — | — | Forwarded to `self.scatter` (e.g. `genes_to_highlight`). |

        ## Returns

        `matplotlib.axes.Axes` — The axes object containing the plot.

        ## Examples

        ```python
        from Pyntegrate.results_table import DEseq2Results

        res = DEseq2Results(
            "deseq2.csv", delimiter=",",
            id_column="GeneID", pval_column="padj",
            lfc_column="log2FoldChange", mean_column="baseMean",
        )
        ax = res.ma_plot(thresh=0.05, zero_line=dict(color="grey", linestyle="--"))
        ```
        """
        genes_to_highlight = kwargs.pop('genes_to_highlight', [])
        genes_to_highlight.append(
            (self.upregulated(thresh),
             up_kwargs or dict(color='r')))
        genes_to_highlight.append(
            (self.downregulated(thresh),
             dn_kwargs or dict(color='b')))
        if zero_line is None:
            zero_line = {}
        x = self.mean_column
        y = self.lfc_column

        if 'xfunc' not in kwargs:
            kwargs['xfunc'] = np.log
        ax = self.scatter(
            x=x,
            y=y,
            genes_to_highlight=genes_to_highlight,
            **kwargs)
        if zero_line:
            ax.axhline(0, **zero_line)
        return ax

    threshdoc = """
    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `thresh` | `float` | `0.05` | Significance threshold; only features with adjusted p-value ≤ `thresh` are returned. |
    | `idx` | `bool` | `True` | If `True`, return a boolean index array. If `False`, return a new subsetted object of the same class. |

    ## Returns

    `pandas.Series[bool] | DifferentialExpressionResults` — Boolean
    index when `idx=True`; a subsetted results object when `idx=False`.

    ## Examples

    ```python
    # Boolean index of up-regulated genes
    mask = res.upregulated(thresh=0.05, idx=True)
    print(mask.sum(), "up-regulated genes")

    # Subsetted object
    up_res = res.upregulated(thresh=0.05, idx=False)
    print(len(up_res), "up-regulated genes")
    ```
    """
    changed.__doc__ = changed.__doc__.format(threshdoc=threshdoc)
    unchanged.__doc__ = unchanged.__doc__.format(threshdoc=threshdoc)
    enriched.__doc__ = enriched.__doc__.format(threshdoc=threshdoc)
    disenriched.__doc__ = disenriched.__doc__.format(threshdoc=threshdoc)
    upregulated.__doc__ = upregulated.__doc__.format(threshdoc=threshdoc)
    downregulated.__doc__ = downregulated.__doc__.format(threshdoc=threshdoc)


class EdgeRResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from edgeR.

        Just like a DifferentialExpressionResults object, but sets the
        pval_column, lfc_column, and mean_column to the names used in edgeR's
        output.

        ## Class attributes

        | Attribute | Default | Description |
        |-----------|---------|-------------|
        | `pval_column` | `"FDR"` | Column used for significance filtering. |
        | `lfc_column` | `"logFC"` | Column with log2 fold-change values. |
        | `mean_column` | `"logCPM"` | Column with mean expression values (for MA plot). |
        """)
    pval_column = 'FDR'
    lfc_column = 'logFC'
    mean_column = 'logCPM'


class DESeqResults(DifferentialExpressionResults):
    __doc__ = _base_doc % dedent(
        """
        Class for working with results from DESeq.

        Uses the same `pval_column` (`"padj"`), `lfc_column`
        (`"log2FoldChange"`), and `mean_column` (`"baseMean"`) defaults
        inherited from `DifferentialExpressionResults` — DESeq's own
        output already uses these column names, so none are overridden
        here. Adds `colormapped_bedfile()` and `autosql_file()` for
        exporting a padj-colored BED/bigBed track of the results.
        """)
    def colormapped_bedfile(self, genome, cmap=None):
        """
        Create a BED file with padj encoded as color.

        Features are colored by their phred-transformed adjusted
        p-value (`-10 * log10(padj)`); down-regulated features
        (negative `log2FoldChange`) have the score sign flipped, and
        rows with a NaN fold-change get a score of `0`. Requires a
        `gffutils.FeatureDB` to be attached (see `attach_db`).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `genome` | `str` | — | Genome assembly name passed to `pybedtools.BedTool.truncate_to_chrom()`, used to clip features to chromosome boundaries. |
        | `cmap` | `matplotlib.colors.Colormap \| None` | `None` | Colormap used to color-code scores. Defaults to `matplotlib.cm.RdBu_r`, centered at zero. |

        ## Returns

        `pybedtools.BedTool` — Sorted, color-coded, chromosome-truncated
        BED file with extra columns for `padj`, `pval`,
        `log2FoldChange`, and `baseMeanB` (twice).

        ## Raises

        `ValueError` — If no database is attached.
        """
        if self.db is None:
            raise ValueError("FeatureDB required")
        db = self.db

        def scored_feature_generator(d):
            for i in range(len(d)):
                gene_id = (
                    d.data[d.id_column].iloc[i]
                    if getattr(d, 'id_column', None)
                    else d.index[i]
                )
                try:
                    feature = db[gene_id]
                except gffutils.FeatureNotFoundError:
                    raise gffutils.FeatureNotFoundError(gene_id)
                score = -10 * np.log10(d.padj.iloc[i])
                lfc = d.log2FoldChange.iloc[i]
                if np.isnan(lfc):
                    score = 0
                if lfc < 0:
                    score *= -1
                feature.score = str(score)
                feature = featurefuncs.extend_fields(
                    featurefuncs.gff2bed(gffutils.helpers.asinterval(feature)), 9)
                fields = feature.fields[:]
                fields[6] = fields[1]
                fields[7] = fields[2]
                fields.append(str(d.padj.iloc[i]))
                fields.append(str(d.pval.iloc[i]))
                fields.append('%.3f' % d.log2FoldChange.iloc[i])
                fields.append('%.3f' % d.baseMeanB.iloc[i])
                fields.append('%.3f' % d.baseMeanB.iloc[i])
                yield pybedtools.create_interval_from_list(fields)

        x = pybedtools.BedTool(scored_feature_generator(self)).saveas()
        norm = x.colormap_normalize()
        if cmap is None:
            cmap = matplotlib.cm.RdBu_r
        cmap = cmap_center_point_adjust(
            cmap, [norm.vmin, norm.vmax], 0)

        def score_zeroer(f):
            f.score = '0'
            return f
        return x.each(featurefuncs.add_color, cmap=cmap, norm=norm)\
                .sort()\
                .each(score_zeroer)\
                .truncate_to_chrom(genome)\
                .saveas()

    def autosql_file(self):
        """
        Generate an autosql file describing the extra BED fields for bigBed conversion.

        Used together with `colormapped_bedfile()`: the extra columns
        added there (`padj`, `pval`, `log2FoldChange`, `baseMeanA`,
        `baseMeanB`) need an autosql schema to be converted to bigBed,
        so that clicking a feature in a genome browser reports the
        DESeq results.

        ## Returns

        `str` — Path to a temporary file containing the autosql schema.
        """
        fn = pybedtools.BedTool._tmp()

        AUTOSQL = dedent(
            """
            table example
            "output from DESeq"
            (
            string  chrom;  "chromosome"
            uint chromStart; "start coord"
            uint chromEnd; "stop coord"
            string name; "name of feature"
            uint score; "always zero"
            char[1] strand; "+ or - for strand"
            uint    thickStart; "Coding region start"
            uint    thickEnd;  "Coding region end"
            uint reserved; "color according to score"
            string padj; "DESeq adjusted p value"
            string pval; "DESeq raw p value"
            string logfoldchange; "DESeq log2 fold change"
            string basemeana; "DESeq baseMeanA"
            string basemeanb; "DESeq baseMeanB"
        )
        """)

        fout = open(fn, 'w')
        fout.write(AUTOSQL)
        fout.close()
        return fn

class DEseq2Results(DifferentialExpressionResults):
    """
    DESeq2 differential expression results with configurable column names.

    Extends `DifferentialExpressionResults` to support user-defined
    column names for gene IDs, adjusted p-values, log2 fold-change, and
    base-mean expression. Optionally strips Ensembl version suffixes from
    gene IDs (e.g. `ENSG000001.3` → `ENSG000001`) and persists the
    normalised IDs in a new column.

    Configuration is stored in `_kwargs` so that `copy()` and `update()`
    can reconstruct the object correctly.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `data` | `str \| pandas.DataFrame` | — | Path to a CSV/TSV file or an existing DataFrame (see `ResultsTable`). |
    | `pval_column` | `str \| None` | `None` | Column holding adjusted p-values (e.g. `"padj"`); set as an *instance* attribute here, overriding the class-level default from `DifferentialExpressionResults`. |
    | `lfc_column` | `str \| None` | `None` | Column holding log2 fold-change values (e.g. `"log2FoldChange"`); instance override. |
    | `mean_column` | `str \| None` | `None` | Column holding mean expression values (e.g. `"baseMean"`); instance override. |
    | `id_column` | `str \| None` | `None` | Column holding gene identifiers (e.g. `"GeneID"`). |
    | `strip_version_ids` | `bool` | `True` | If `True`, `gene_ids()` strips the Ensembl version suffix (e.g. `.3`) from `id_column` values. |
    | `id_norm_column` | `str` | `"GeneID_norm"` | Name of the column added to `data` (during `__init__`) holding the normalised IDs returned by `gene_ids()`. |
    | `verbose` | `bool` | `False` | If `True`, print configuration details during construction. |

    ## Raises

    `ValueError` — If any of the required columns (`pval_column`,
    `lfc_column`, `mean_column`, `id_column`) are `None` or absent from
    the DataFrame.

    ## Examples

    ```python
    from Pyntegrate.results_table import DEseq2Results

    res = DEseq2Results(
        "deseq2_results.csv",
        delimiter=",",
        id_column="GeneID",
        pval_column="padj",
        lfc_column="log2FoldChange",
        mean_column="baseMean",
        strip_version_ids=True,
        verbose=True,
    )

    print(res.data.shape)
    print(res.upregulated(thresh=0.05, idx=True).sum(), "up-regulated genes")
    ```

    **Output:**
    ```
    (20234, 7)
    342 up-regulated genes
    ```
    """

    def __init__(
        self,
        data,
        pval_column=None,
        lfc_column=None,
        mean_column=None,
        id_column=None,
        strip_version_ids=True,
        id_norm_column="GeneID_norm",
        verbose=False,
        **kwargs,
    ):
        # Allow passing via kwargs (update()/copy() safety)
        pval_column = pval_column or kwargs.pop("pval_column", None)
        lfc_column = lfc_column or kwargs.pop("lfc_column", None)
        mean_column = mean_column or kwargs.pop("mean_column", None)
        id_column = id_column or kwargs.pop("id_column", None)
        strip_version_ids = kwargs.pop("strip_version_ids", strip_version_ids)
        id_norm_column = kwargs.pop("id_norm_column", id_norm_column)
        verbose = kwargs.pop("verbose", verbose)

        # If data already has config (e.g., from update/copy), inherit it
        if hasattr(data, "_kwargs"):
            pval_column = pval_column or data._kwargs.get("pval_column")
            lfc_column = lfc_column or data._kwargs.get("lfc_column")
            mean_column = mean_column or data._kwargs.get("mean_column")
            id_column = id_column or data._kwargs.get("id_column")
            strip_version_ids = data._kwargs.get("strip_version_ids", strip_version_ids)
            id_norm_column = data._kwargs.get("id_norm_column", id_norm_column)
            verbose = data._kwargs.get("verbose", verbose)

        super().__init__(data=data, **kwargs)

        # Save config
        self.pval_column = pval_column
        self.lfc_column = lfc_column
        self.mean_column = mean_column
        self.id_column = id_column
        self.strip_version_ids = strip_version_ids
        self.id_norm_column = id_norm_column
        self.verbose = verbose

        # Persist config for future reconstruction (copy/update)
        self._kwargs.update(
            dict(
                pval_column=pval_column,
                lfc_column=lfc_column,
                mean_column=mean_column,
                id_column=id_column,
                strip_version_ids=strip_version_ids,
                id_norm_column=id_norm_column,
                verbose=verbose,
            )
        )

        # Validate required columns
        missing = []
        for col in [pval_column, lfc_column, mean_column, id_column]:
            if col is None:
                missing.append("<None>")
            elif col not in self.data.columns:
                missing.append(col)

        if missing:
            raise ValueError(
                f"[Pyntegrate ERROR] Missing/invalid columns in DataFrame: {missing}. "
                f"Available columns: {list(self.data.columns)}"
            )

        # Build normalized id column immediately (cheap + prevents mistakes later)
        self.data[self.id_norm_column] = self.gene_ids()

        if self.verbose:
            print("[Pyntegrate] DEseq2Results configuration:")
            print(f"  - id_column        = {self.id_column}")
            print(f"  - id_norm_column   = {self.id_norm_column}")
            print(f"  - strip_version    = {self.strip_version_ids}")
            print(f"  - pval_column      = {self.pval_column}")
            print(f"  - lfc_column       = {self.lfc_column}")
            print(f"  - mean_column      = {self.mean_column}")

    def gene_ids(self):
        """
        Return normalized gene IDs from `id_column`, optionally stripping the version suffix.

        ## Returns

        `pandas.Series[str]` — Values of `self.data[self.id_column]`
        cast to string; if `self.strip_version_ids` is `True`,
        everything from the first `"."` onward is removed (e.g.
        `"ENSG00000001.3"` → `"ENSG00000001"`).
        """
        s = self.data[self.id_column].astype(str)
        if self.strip_version_ids:
            s = s.str.split(".").str[0]
        return s




    def colormapped_bedfile(self, genome, cmap=None):
        """
        Create a BED file with the adjusted p-value encoded as color.

        Features are colored by the phred-transformed value of this
        result's configured `pval_column` (`-10 * log10(...)`);
        down-regulated features (negative `lfc_column` value) have the
        score sign flipped, and rows with a NaN fold-change get a
        score of `0`. Requires a `gffutils.FeatureDB` to be attached
        (see `attach_db`).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `genome` | `str` | — | Genome assembly name passed to `pybedtools.BedTool.truncate_to_chrom()`, used to clip features to chromosome boundaries. |
        | `cmap` | `matplotlib.colors.Colormap \| None` | `None` | Colormap used to color-code scores. Defaults to `matplotlib.cm.RdBu_r`, centered at zero. |

        ## Returns

        `pybedtools.BedTool` — Sorted, color-coded, chromosome-truncated
        BED file with extra columns for this result's configured
        `pval_column`, `lfc_column`, and `mean_column` values.

        ## Raises

        `ValueError` — If no database is attached.
        """
        if self.db is None:
            raise ValueError("FeatureDB required")
        db = self.db

        def scored_feature_generator(d):
            for i in range(len(d)):
                gene_id = (
                    d.data[d.id_column].iloc[i]
                    if getattr(d, 'id_column', None)
                    else d.index[i]
                )
                try:
                    feature = db[gene_id]
                except gffutils.FeatureNotFoundError:
                    raise gffutils.FeatureNotFoundError(gene_id)
                score = -10 * np.log10(d.data[d.pval_column].iloc[i])
                lfc = d.data[d.lfc_column].iloc[i]
                if np.isnan(lfc):
                    score = 0
                if lfc < 0:
                    score *= -1
                feature.score = str(score)
                feature = featurefuncs.extend_fields(
                    featurefuncs.gff2bed(gffutils.helpers.asinterval(feature)), 9)
                fields = feature.fields[:]
                fields[6] = fields[1]
                fields[7] = fields[2]
                fields.append(str(d.data[d.pval_column].iloc[i]))
                fields.append('%.3f' % d.data[d.lfc_column].iloc[i])
                fields.append('%.3f' % d.data[d.mean_column].iloc[i])
                yield pybedtools.create_interval_from_list(fields)

        x = pybedtools.BedTool(scored_feature_generator(self)).saveas()
        norm = x.colormap_normalize()
        if cmap is None:
            cmap = matplotlib.cm.RdBu_r
        cmap = cmap_center_point_adjust(
            cmap, [norm.vmin, norm.vmax], 0)

        def score_zeroer(f):
            f.score = '0'
            return f
        return x.each(featurefuncs.add_color, cmap=cmap, norm=norm)\
                .sort()\
                .each(score_zeroer)\
                .truncate_to_chrom(genome)\
                .saveas()

    def autosql_file(self):
        """
        Generate an autosql file describing the extra BED fields for bigBed conversion.

        Used together with `colormapped_bedfile()`: the extra columns
        added there (this result's configured `pval_column`,
        `lfc_column`, and `mean_column` values) need an autosql schema
        to be converted to bigBed, so that clicking a feature in a
        genome browser reports the DESeq2 results.

        ## Returns

        `str` — Path to a temporary file containing the autosql schema.
        """
        fn = pybedtools.BedTool._tmp()

        AUTOSQL = dedent(
            """
            table example
            "output from DESeq2"
            (
            string  chrom;  "chromosome"
            uint chromStart; "start coord"
            uint chromEnd; "stop coord"
            string name; "name of feature"
            uint score; "always zero"
            char[1] strand; "+ or - for strand"
            uint    thickStart; "Coding region start"
            uint    thickEnd;  "Coding region end"
            uint reserved; "color according to score"
            string pval; "this result's configured pval_column value"
            string log2foldchange; "this result's configured lfc_column value"
            string mean; "this result's configured mean_column value"
        )
        """)

        fout = open(fn, 'w')
        fout.write(AUTOSQL)
        fout.close()
        return fn



class LazyDict(object):
    """
    Dictionary-like object that lazily loads `ResultsTable` objects from files.

    Keys of `fn_dict` become the keys of this object. On first access via
    `__getitem__`, the corresponding file is loaded into a new instance of
    `cls`, reindexed to `index_file`, optionally merged with `extra`, and
    optionally passed through `modifier`. The resulting object is cached so
    subsequent access of the same key returns it directly.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn_dict` | `dict` | — | Keys become the keys of this object; values are filenames loaded into `_cls` objects upon first access. |
    | `dbfn` | `str` | — | Filename of a `gffutils` database, passed as `db=` to each loaded object. |
    | `index` | `list[str]` | — | IDs read from the `index_file` constructor argument, one per line; used to reindex every loaded table to a common index. |
    | `extra` | `pandas.DataFrame \| None` | `None` | Extra columns merged into each loaded table's data (e.g. gene lengths), pre-reindexed to `self.index`. |
    | `modifier` | `callable \| None` | `None` | Applied to each newly-loaded, `extra`-merged object; must return the object to cache and use from then on. |
    | `_dict` | `dict` | — | Cache of already-loaded objects, keyed the same as `fn_dict`. |
    | `_cls` | `ResultsTable` subclass | `DEseq2Results` | Class used to wrap each loaded file (constructor argument name: `cls`). |

    ## Notes

    When a key is accessed for the first time, the workflow is
    `cls(fn, db=dbfn)` → reindex to `index` → merge `extra` → pass to
    `modifier` → cache and return the result. Subsequent access of the
    same key immediately returns the cached, already-processed object.
    """
    def __init__(self, fn_dict, dbfn, index_file, extra=None, cls=DEseq2Results,
                 modifier=None):
        self.fn_dict = fn_dict
        self._dict = {}
        self.dbfn = dbfn
        self.index = [i.strip() for i in open(index_file)]
        if extra is not None:
            self.extra = extra.loc[self.index]
        else:
            self.extra = extra
        self.modifier = modifier
        self._cls = cls

    def __getitem__(self, key):
        """
        Lazily load, process, and cache the object for `key`.

        On first access, loads `self.fn_dict[key]` into a new
        `self._cls` instance (with `db=self.dbfn`), reindexes it to
        `self.index`, merges in `self.extra` if given, and applies
        `self.modifier` if given. The result is cached in `self._dict`
        so subsequent lookups of the same `key` skip all of this and
        return the cached object directly.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | hashable | — | Key into `self.fn_dict` identifying which file to load. |

        ## Returns

        An instance of `self._cls` (or subclass), possibly
        post-processed by `self.modifier`.

        ## Raises

        `KeyError` — If `key` is not present in `self.fn_dict`.
        """
        if key not in self._dict:
            fn = self.fn_dict[key]
            obj = self._cls(fn, db=self.dbfn)
            obj.data = obj.data.loc[self.index]
            if self.extra is not None:
                obj.data = pandas.merge(obj.data, self.extra, left_index=True,
                                        right_index=True)
            if self.modifier:
                obj = self.modifier(obj)

            self._dict[key] = obj

        return self._dict[key]

    def __repr__(self):
        """
        Return a string listing all possible keys (`self.fn_dict`) and the keys already loaded (`self._dict`).
        """
        s = "<%s> with possible keys\n:%s\n" \
            % (self.__class__.__name__, self.fn_dict.keys())
        s += "and existing keys:\n"
        s += repr(self._dict)
        return s

    def keys(self):
        """
        Return the possible keys (`self.fn_dict.keys()`), whether or not they have been loaded yet.
        """
        return self.fn_dict.keys()

    def values(self):
        """
        Return the already-loaded object for each key, in `self.keys()` order.

        ## Returns

        `list` — One entry per key in `self.keys()`, taken directly
        from `self._dict` (does **not** trigger lazy-loading via
        `__getitem__`).

        ## Raises

        `KeyError` — If any key in `self.keys()` has not yet been
        accessed via `__getitem__` (and so is not in `self._dict`).
        """
        return [self._dict[key] for key in self.keys()]

    def items(self):
        """
        Return `(key, value)` pairs for each already-loaded key, in `self.keys()` order.

        ## Returns

        `list[tuple]` — One `(key, value)` pair per key in
        `self.keys()`, taken directly from `self._dict` (does **not**
        trigger lazy-loading via `__getitem__`).

        ## Raises

        `KeyError` — If any key in `self.keys()` has not yet been
        accessed via `__getitem__` (and so is not in `self._dict`).
        """
        return list((key, self._dict[key]) for key in self.keys())
