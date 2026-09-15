"""
General-purpose utilities shared across the Pyntegrate package.

Small, dependency-light helpers for array resampling, chunking
iterables for parallel processing, locating bundled example data, and
parsing/converting genomic interval strings.

## Contents

| Function | Description |
|----------|-------------|
| `rebin` | Linearly interpolate `y` onto a new uniform grid of `nbin` points. |
| `chunker` | Split an iterable into fixed-size chunks (for `multiprocessing.Pool`). |
| `split_feature` | Split a single genomic feature into `n` equal sub-intervals. |
| `tointerval` | Parse a `"chrom:start-stop"`-style string into a `pybedtools.Interval`. |
| `data_dir` | Path to the package's bundled example-data directory. |
| `example_filename` | Resolve the full path of a bundled example file by name. |
| `list_example_files` | List bundled example filenames, optionally filtered by a glob pattern. |
| `map_color` | Map a log2 fold-change value to a color on a red/blue gradient. |
"""
from copy import copy
import os
import glob
import pybedtools
import re
import numpy as np


def rebin(x, y, nbin):
    """
    Linearly interpolate `y` onto a new uniform grid of `nbin` points.

    Constructs a new x-axis of `nbin` evenly-spaced points spanning
    `[x.min(), x.max()]`, then uses linear interpolation to evaluate
    `y` at each new point. Useful for normalising genomic signal
    profiles so that regions of different lengths can be compared or
    averaged.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `x` | `np.ndarray` | — | Original x-axis coordinates, shape `(N,)`. |
    | `y` | `np.ndarray` | — | Original y-axis values, shape `(N,)`. |
    | `nbin` | `int` | — | Number of evenly-spaced points in the output grid. |

    ## Returns

    `tuple[np.ndarray, np.ndarray]` — `(xi, yi)` where `xi` is the new
    x-axis of shape `(nbin,)` and `yi` contains the interpolated y
    values at each point of `xi`, also shape `(nbin,)`.

    ## Examples

    ```python
    import numpy as np
    from Pyntegrate.helpers import rebin

    x = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y = np.array([0.0, 1.0, 0.5, 2.0, 1.5])
    xi, yi = rebin(x, y, nbin=3)
    print(xi)
    print(np.round(yi, 4))
    ```

    **Output:**
    ```
    [0. 2. 4.]
    [0.   0.5  1.5]
    ```
    """
    xi = np.linspace(x.min(), x.max(), nbin)
    return xi, np.interp(xi, x, y)


def chunker(f, n):
    """
    Split an iterable into successive chunks of at most `n` elements.

    The last chunk may contain fewer than `n` elements if the total
    length of `f` is not a multiple of `n`. This is used internally
    to distribute genomic features across parallel worker processes.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `f` | `iterable` | — | Any iterable to be split into chunks. |
    | `n` | `int` | — | Maximum number of elements per chunk. |

    ## Yields

    `tuple` — Successive tuples of at most `n` elements drawn from `f`.
    The last tuple may be shorter.

    ## Examples

    ```python
    from Pyntegrate.helpers import chunker

    result = list(chunker(range(7), 3))
    print(result)
    ```

    **Output:**
    ```
    [(0, 1, 2), (3, 4, 5), (6,)]
    ```
    """
    f = iter(f)
    x = []
    while 1:
        if len(x) < n:
            try:
                x.append(f.__next__())
            except StopIteration:
                if len(x) > 0:
                    yield tuple(x)
                break
        else:
            yield tuple(x)
            x = []


def data_dir():
    """
    Return the absolute path to the bundled example data directory.

    The directory contains small genomic files (BAM, BED, BigWig, etc.)
    that are used by the test suite and the package documentation.

    ## Returns

    `str` — Absolute path to the `test/data/` subdirectory inside the
    installed Pyntegrate package.

    ## Examples

    ```python
    from Pyntegrate.helpers import data_dir

    path = data_dir()
    print(path)
    ```

    **Output:**
    ```
    /path/to/site-packages/Pyntegrate/test/data
    ```
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')


def example_filename(fn):
    """
    Return the absolute path of a bundled example data file.

    Use [`list_example_files`][Pyntegrate.helpers.list_example_files]
    to discover which files are available.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Basename of the requested example file (e.g. `"gdc.bam"`). |

    ## Returns

    `str` — Absolute path to the requested file.

    ## Raises

    `ValueError` — If `fn` does not exist in the example data directory.

    ## Examples

    ```python
    from Pyntegrate.helpers import example_filename

    path = example_filename("gdc.bam")
    print(path)
    ```

    **Output:**
    ```
    /path/to/site-packages/Pyntegrate/test/data/gdc.bam
    ```
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


def list_example_files(pattern=None, full_path=False):
    """
    List example data files bundled with Pyntegrate.

    Scans the internal `test/data/` directory and optionally filters
    results by a glob pattern.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `pattern` | `str \| None` | `None` | Glob pattern to filter filenames (e.g. `"*.bigwig"`). If `None`, all files in the data directory are listed. |
    | `full_path` | `bool` | `False` | If `True`, return absolute paths. If `False`, return bare filenames only. |

    ## Returns

    `list[str]` — List of matching filenames or full paths.

    ## Examples

    ```python
    from Pyntegrate.helpers import list_example_files

    # All example files (basenames only)
    print(list_example_files())
    ```

    **Output:**
    ```
    ['gdc.bam', 'gdc.bam.bai', 'gdc.bed', 'gdc.bigwig', ...]
    ```

    ```python
    # Only BigWig files, with full paths
    print(list_example_files("*.bigwig", full_path=True))
    ```

    **Output:**
    ```
    ['/path/to/Pyntegrate/test/data/gdc.bigwig', ...]
    ```
    """
    if not os.path.isdir(data_dir()):
        raise FileNotFoundError(
            f"Example data directory not found: {data_dir()}. "
            "The bundled metaseq example dataset is not included in this "
            "installation of Pyntegrate."
        )
    if pattern is None:
        files = os.listdir(data_dir())
    else:
        files = glob.glob(os.path.join(data_dir(), pattern))
    if full_path:
        return files
    else:
        return [os.path.basename(i) for i in files]


def split_feature(f, n):
    """
    Split a genomic interval into `n` roughly equal sub-intervals.

    Divides the feature by stepping through it in increments of
    `(stop - start) / n` base pairs. The last sub-interval is
    clipped to the original stop coordinate.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `f` | `pybedtools.Interval` | — | Genomic interval to split. All attributes (chrom, name, strand, etc.) are copied to each piece. |
    | `n` | `int` | — | Number of sub-intervals to create. |

    ## Yields

    `pybedtools.Interval` — Sub-intervals of `f` with adjusted `start`
    and `stop` coordinates and all other attributes preserved.

    ## Raises

    `ValueError` — If `n` is not an integer.

    ## Examples

    ```python
    import pybedtools
    from Pyntegrate.helpers import split_feature

    feature = pybedtools.create_interval_from_list(
        ["chr1", "1000", "2000", "gene1", ".", "+"]
    )
    for piece in split_feature(feature, n=3):
        print(piece.chrom, piece.start, piece.stop)
    ```

    **Output:**
    ```
    chr1 1000 1333
    chr1 1333 1667
    chr1 1667 2000
    ```
    """
    if not isinstance(n, int):
        raise ValueError('n must be an integer')
    orig_feature = copy(f)
    breakpoints = [
        int(round(orig_feature.start + i * (orig_feature.stop - orig_feature.start) / n))
        for i in range(n + 1)
    ]
    for i in range(n):
        piece = copy(orig_feature)
        piece.start = breakpoints[i]
        piece.stop = breakpoints[i + 1]
        yield piece

coord_re = re.compile(
    r"""
    (?P<chrom>.+):
    (?P<start>\d+)-
    (?P<stop>\d+)
    (?:\[(?P<strand>.)\])?""", re.VERBOSE)


def tointerval(s):
    """
    Convert a coordinate string to a `pybedtools.Interval`.

    Pass-through for anything that is already an interval-like object.

    Accepted string formats:

    - `"chrom:start-stop"` — e.g. `"chr1:1000-2000"`
    - `"chrom:start-stop[strand]"` — e.g. `"chr1:1000-2000[+]"`

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `s` | `str \| interval-like` | — | Coordinate string or an existing interval object. Non-string values are returned unchanged. |

    ## Returns

    `pybedtools.Interval` — Parsed interval, or `s` unchanged if it is
    not a string.

    ## Examples

    ```python
    from Pyntegrate.helpers import tointerval

    iv = tointerval("chr1:1000-2000")
    print(iv.chrom, iv.start, iv.stop)

    iv_stranded = tointerval("chrX:500-1500[+]")
    print(iv_stranded.chrom, iv_stranded.start, iv_stranded.stop, iv_stranded.strand)
    ```

    **Output:**
    ```
    chr1 1000 2000
    chrX 500 1500 +
    ```
    """
    if isinstance(s, str):
        m = coord_re.search(s)
        if m.group('strand'):
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
                '.',
                '0',
                m.group('strand')])
        else:
            return pybedtools.create_interval_from_list([
                m.group('chrom'),
                m.group('start'),
                m.group('stop'),
            ])
    return s


def map_color(a, down_log2FoldChange, upper_log2FoldChange, upper_nlog10):
    """
    Classify a gene into an expression-change category for volcano-plot colouring.

    Returns one of three string labels based on whether the gene's
    log2 fold-change and transformed p-value exceed the given thresholds.
    Used internally by volcano-plot functions to assign point colours.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `a` | `tuple[float, float]` | — | A `(log2FoldChange, -log10(padj))` pair for a single gene. |
    | `down_log2FoldChange` | `float` | — | LFC threshold for down-regulation (typically negative, e.g. `-1.5`). |
    | `upper_log2FoldChange` | `float` | — | LFC threshold for up-regulation (e.g. `1.5`). |
    | `upper_nlog10` | `float` | — | `-log10(padj)` threshold for statistical significance. |

    ## Returns

    `str` — One of `"increase"`, `"decrease"`, or `"not significant"`.

    ## Examples

    ```python
    from Pyntegrate.helpers import map_color

    print(map_color((2.0, 5.0), -1.5, 1.5, 2.0))
    print(map_color((-2.5, 4.0), -1.5, 1.5, 2.0))
    print(map_color((0.5, 1.0), -1.5, 1.5, 2.0))
    ```

    **Output:**
    ```
    increase
    decrease
    not significant
    ```
    """
    log2FoldChange, nlog10 = a

    if log2FoldChange > upper_log2FoldChange and nlog10 > upper_nlog10:
        return 'increase'
    elif log2FoldChange < down_log2FoldChange and nlog10 > upper_nlog10:
        return 'decrease'
    return 'not significant'
