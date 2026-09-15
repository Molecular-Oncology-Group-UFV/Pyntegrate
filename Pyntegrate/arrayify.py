"""
Bin whole-genome bigWig signal into fixed-size windows and cache the result to disk.

The single class here, `Binner`, wraps the `bigWigAverageOverBed`
command-line tool (must be available on `PATH`) to convert a bigWig
file into per-chromosome NumPy arrays — bin-midpoint coordinates (`x`)
and a summary metric per bin (`y`) — saved as `.npz` files for fast,
memory-mapped re-use across sessions without recomputing the binning.
"""
import os
import pandas
import pybedtools
import numpy as np
from ._genomic_signal import *


class Binner(object):
    """
    Class to handle converting bigWig files into NumPy arrays.

    `bigWigAverageOverBed` needs to be available on the path.

    The arrays are saved to disk as .npz files, which can then be
    memory-mapped for fast, lightweight re-use. Each .npz file contains
    the coordinates of the bin midpoints (x) and values in each bin (y).

    The class is designed to be set up once, then used many times on
    different bigWig files.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `chromsizes` | `dict` | — | Chromosome sizes for `genome`, as returned by `pybedtools.chromsizes`. |
    | `chroms` | `list[str]` | — | Chromosomes to process — either all keys of `chromsizes` (sorted), or the single `chrom` given to the constructor. |
    | `windowsize` | `int` | — | Bp size of each window. |
    | `window_cache_dir` | `str` | `"."` | Directory where windowed-chromosome BED files are cached. These files are cached to avoid creating them every time, and have the filename pattern `{window_cache_dir}/{chrom}.{windowsize}bp_windows.bed`. |

    The constructor also accepts `genome` (`str`, assembly name such as
    `hg19` or `dm3`) and `chrom` (`str | None`, default `None`), which are
    used to compute `chromsizes` and `chroms` above but are not stored as
    attributes themselves. Two further constructor arguments, `npz_dir`
    (default `"."`) and `metric` (default `"mean0"`), are accepted but
    currently unused.

    ## Examples

    ```python
    b = Binner('mm9', 1000, chrom='chr19')
    b.to_npz('PolII.bigwig')
    pol = np.load('PolII.bigwig.chr19.npz')

    # Assuming matplotlib is installed
    from matplotlib import pyplot as plt
    plt.plot(pol['x'], pol['y'])
    plt.show()
    ```
    """
    def __init__(self, genome, windowsize, chrom=None, window_cache_dir=".",
                 npz_dir='.', metric='mean0'):
        self.chromsizes = pybedtools.chromsizes(genome)
        if chrom is None:
            self.chroms = sorted(self.chromsizes.keys())
        else:
            self.chroms = [chrom]
        self.windowsize = windowsize
        self.window_cache_dir = window_cache_dir

    def make_windows(self, chrom, force=False):
        """
        Create (or reuse) a BED file of fixed-size windows tiling `chrom`.

        The output BED file is cached at
        `{window_cache_dir}/{chrom}.{windowsize}bp_windows.bed`; if it
        already exists and `force` is `False`, the existing path is
        returned without recomputing the windows.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `chrom` | `str` | — | Chromosome name to tile with windows. Must be a key of `self.chromsizes`. |
        | `force` | `bool` | `False` | If `True`, recreate the windows file even if it already exists on disk. |

        ## Returns

        `str` — Path to the (possibly cached) BED file of windows.
        """
        chromsize = self.chromsizes[chrom][1]
        bed = pybedtools.BedTool(
            '{0} 0 {1}'.format(chrom, chromsize),
            from_string=True
        )
        output = os.path.join(
            self.window_cache_dir,
            '%s.%sbp_windows.bed' % (chrom, self.windowsize))
        if force or not os.path.exists(output):
            windows = pybedtools.BedTool()\
                .window_maker(
                    b=bed,
                    i='winnum',
                    w=self.windowsize,
                    output=output,
                )
        return output

    def to_npz(self, bigwig, metric='mean0', outdir=None):
        """
        Bin a bigWig file into fixed-size windows and save the result to disk.

        Uses `bigWigAverageOverBed` to summarize `bigwig` over the cached
        windows for each chromosome in `self.chroms`, then saves one
        `.npz` file per chromosome with the pattern
        `{outdir}/{bigwig}.{chrom}.{windowsize}.{metric}.npz`. Each file
        contains two arrays: `x` (genomic coordinates of the midpoint of
        each window) and `y` (the requested `metric` for each window). It
        can be loaded with `np.load(filename, mmap_mode='r')`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `bigwig` | `str \| BigWigSignal` | — | BigWig data to bin, as a path or an already-open `BigWigSignal` object. |
        | `metric` | `str` | `"mean0"` | Metric to store, as reported by `bigWigAverageOverBed`: `"covered"` (number of bases covered by the bigWig), `"sum"` (sum of values over all covered bases), `"mean0"` (average over all bases, non-covered bases counted as zero), or `"mean"` (average over only the covered bases). |
        | `outdir` | `str \| None` | `None` | Directory to store output files in. If `None`, files are stored in the same directory as `bigwig`. |

        ## Returns

        `list[str]` — Paths to the `.npz` files written, one per
        chromosome in `self.chroms`.
        """
        if isinstance(bigwig, BigWigSignal):
            bigwig = bigwig.fn

        if outdir is None:
            outdir = os.path.dirname(bigwig)

        basename = os.path.basename(bigwig)
        windowsize = self.windowsize

        outfiles = []
        for chrom in self.chroms:
            tmp_output = pybedtools.BedTool._tmp()
            windows = self.make_windows(chrom)

            outfile = os.path.join(
                outdir,
                '{basename}.{chrom}.{windowsize}.{metric}'.format(**locals())
                + '.npz')
            cmds = [
                'bigWigAverageOverBed',
                bigwig,
                windows,
                tmp_output]
            os.system(' '.join(cmds))
            names = ['name', 'size', 'covered', 'sum', 'mean0', 'mean']
            df = pandas.read_table(tmp_output, names=names)
            x = df["size"].cumsum() - df["size"] / 2
            y = df[metric]
            np.savez(outfile, x=x, y=y)
            outfiles.append(outfile)
            del x, y, df
        return outfiles
