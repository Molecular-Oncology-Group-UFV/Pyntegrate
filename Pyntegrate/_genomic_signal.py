"""
Uniform interface for computing genomic signal from multiple file formats.

This module provides classes that wrap BAM, BigWig, BigBed, BED, and
NarrowPeak files behind a single API, so that coverage and count
arrays can be computed with one call regardless of the underlying
file format.

## Supported formats

| Class | Format | Description |
|-------|--------|-------------|
| `BigWigSignal` | BigWig (`.bw`) | Pre-computed continuous coverage tracks |
| `BamSignal` | BAM (`.bam`) | Aligned short reads |
| `BedSignal` | BED / GFF / GTF / VCF | Generic interval-based features |
| `BigBedSignal` | BigBed (`.bb`) | Large indexed binary BED files |
| `NarrowPeakSignal` | NarrowPeak | MACS2/MACS3 peak-call output |

## Quick start

```python
from Pyntegrate import genomic_signal
import pybedtools

sig = genomic_signal("experiment.bw", "bigwig")
genes = pybedtools.BedTool("genes.bed")
array = sig.array(genes, bins=100)   # shape (n_genes, 100)
```

## Extending with new formats

1. Write a new adapter in `Pyntegrate.filetype_adapters`.
2. Subclass `BaseSignal` or `IntervalSignal`, setting `self.adapter`.
3. Register the class in `_registry` at the bottom of this module.
"""

import os
import sys
import subprocess
import pyBigWig

import numpy as np

from .array_helpers import _array, _array_parallel, _local_coverage, \
    _local_count, _count_array, _count_array_parallel
from .filetype_adapters import *
from .helpers import *


def supported_formats():
    """
    Return the file formats supported by Pyntegrate's genomic signal objects.

    ## Returns

    `dict_keys` — Keys of the internal `_registry` dictionary, one entry
    per supported format string (e.g. `"bam"`, `"bigwig"`, `"bed"`, …).

    ## Examples

    ```python
    from Pyntegrate._genomic_signal import supported_formats

    for fmt in supported_formats():
        print(fmt)
    ```

    **Output:**
    ```
    bam
    bed
    gff
    gtf
    vcf
    bigwig
    bigbed
    ```
    """
    return _registry.keys()


def genomic_signal(fn, kind):
    """
    Factory function that instantiates the correct signal class for a file.

    This is the primary entry point for working with genomic signal
    files. It looks up `kind` in the internal registry and returns an
    object of the matching subclass, with the file path stored in `.fn`
    and the format string stored in `.kind`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the genomic signal file. |
    | `kind` | `str` | — | Format string identifying the file type. Case-insensitive. See [`supported_formats`][Pyntegrate._genomic_signal.supported_formats] for valid values. |

    ## Returns

    `BaseSignal` — An instance of the matching signal subclass
    (`BigWigSignal`, `BamSignal`, `BedSignal`, or `BigBedSignal`).

    ## Raises

    `ValueError` — If `kind` is not found in the format registry.

    ## Examples

    ```python
    from Pyntegrate import genomic_signal

    bw = genomic_signal("chip.bw", "bigwig")
    print(type(bw).__name__, bw.kind)

    bam = genomic_signal("aligned.bam", "bam")
    print(type(bam).__name__, bam.kind)
    ```

    **Output:**
    ```
    BigWigSignal bigwig
    BamSignal bam
    ```
    """
    try:
        klass = _registry[kind.lower()]
    except KeyError:
        raise ValueError(
            'No support for %s format, choices are %s'
            % (kind, _registry.keys()))
    m = klass(fn)
    m.kind = kind
    return m


class BaseSignal(object):
    """
    Abstract base class for all genomic signal objects.

    Provides `array` and `local_coverage` methods that work uniformly
    across all supported file formats. Concrete subclasses set
    `self.adapter` to a format-specific reader from
    `Pyntegrate.filetype_adapters`.

    Do not instantiate this class directly — use the
    [`genomic_signal`][Pyntegrate._genomic_signal.genomic_signal]
    factory instead.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the underlying signal file. |
    | `kind` | `str` | Format string set by the factory function. |
    | `adapter` | `BaseAdapter` | File-format-specific reader (set by subclass). |
    """

    def __init__(self, fn):
        """
        Store the path to the signal file.

        See the class docstring for the `fn` attribute description.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the underlying signal file. |
        """
        self.fn = fn

    def array(self, features, processes=None, chunksize=1, ragged=False,
              **kwargs):
        """
        Build an M × N array of genomic signal over a list of features.

        Each row corresponds to one feature in `features`; each column
        corresponds to one signal bin across the feature window. All
        additional keyword arguments are forwarded to `local_coverage`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `features` | `iterable` | — | Interval-like objects (e.g. `pybedtools.Interval`) to compute signal over. |
        | `processes` | `int \| None` | `None` | Number of parallel worker processes. `None` runs sequentially. |
        | `chunksize` | `int` | `1` | Features per worker chunk when `processes` is set. Values of 25–100 work well for most datasets. |
        | `ragged` | `bool` | `False` | If `False` (default), stack rows into a 2-D NumPy array (requires uniform row lengths). If `True`, return a list of 1-D arrays. |
        | `**kwargs` | — | — | Forwarded to `local_coverage` — e.g. `bins`, `fragment_size`, `stranded`. |

        ## Returns

        `np.ndarray | list[np.ndarray]` — Shape `(M, N)` array when
        `ragged=False`; a list of M 1-D arrays when `ragged=True`.

        ## Examples

        ```python
        import pybedtools
        from Pyntegrate import genomic_signal

        sig = genomic_signal("chip.bw", "bigwig")
        genes = pybedtools.BedTool("genes.bed")

        arr = sig.array(genes, bins=100)
        print(arr.shape)
        ```

        **Output:**
        ```
        (5000, 100)
        ```
        """
        if processes is not None:
            arrays = _array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)


        else:
            arrays = _array(self.fn, self.__class__, features, **kwargs)
        if not ragged:
            stacked_arrays = np.row_stack(arrays)
            return stacked_arrays
        else:
            return arrays

    def local_coverage(self, features, *args, **kwargs):
        """
        Compute a binned coverage profile for one or more genomic features.

        Thin wrapper around
        [`_local_coverage`][Pyntegrate.array_helpers._local_coverage]
        using `self.adapter`. When `processes` is set to a positive
        value, a single feature is split into that many chunks, each
        chunk's signal is computed in parallel via `self.array`, and the
        pieces are concatenated back together (optionally rebinned).

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `features` | `interval-like \| str \| list` | — | Feature(s) to compute coverage over. A single feature is required when `processes` is set. |
        | `*args` | — | — | Forwarded to `_local_coverage` (or to `self.array` when `processes` is set). |
        | `processes` | `int \| None` | `None` | Number of parallel worker processes for a single large feature. `None`/`0` runs sequentially through `_local_coverage`. |
        | `bins` | `int \| None` | `None` | Number of bins to rebin the parallel result into. Only used when `processes` is set. |
        | `**kwargs` | — | — | Forwarded to `_local_coverage` (or to `self.array` when `processes` is set). |

        ## Returns

        `tuple` — `(x, profile)` as returned by `_local_coverage` when
        `processes` is not set; `(x, y)` NumPy coordinate/value arrays
        assembled from parallel chunks (rebinned via `rebin` if `bins`
        is given) when `processes` is set.

        ## Raises

        `ValueError` — If `processes` is set and `features` is a list
        or tuple (parallel mode only supports a single feature).
        """
        processes = kwargs.pop('processes', None)
        if not processes:
            return _local_coverage(self.adapter, features, *args, **kwargs)

        if isinstance(features, (list, tuple)):
            raise ValueError(
                "only single features are supported for parallel "
                "local_coverage")

        bins = kwargs.pop('bins', None)

        features = tointerval(features)
        x = np.arange(features.start, features.stop)
        features = list(split_feature(features, processes))
        ys = self.array(
            features, *args, bins=None, processes=processes, ragged=True,
            **kwargs)
        y = np.column_stack(ys).ravel()
        if bins:
            xi, yi = rebin(x, y, bins)
            return xi, yi
        return x, y


class BigWigSignal(BaseSignal):
    """
    Genomic signal interface for BigWig (`.bw`) files.

    BigWig files store pre-computed, continuous coverage values at
    single-base resolution. This class wraps a BigWig adapter and
    inherits the `array` and `local_coverage` methods from `BaseSignal`.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the BigWig file. |
    | `adapter` | `BigWigAdapter` | — | File-format-specific reader wrapping the BigWig file. |
    | `_readcount` | `int \| float \| None` | `None` | Cached mapped-read count. `None` until populated by `mapped_read_count()`. |

    ## Examples

    ```python
    from Pyntegrate import genomic_signal

    sig = genomic_signal("chip.bw", "bigwig")
    print(sig.mapped_read_count())
    ```

    **Output:**
    ```
    15243600
    ```
    """

    def __init__(self, fn):
        """
        Open `fn` as a BigWig file and initialize the read-count cache.

        See the class docstring for attribute descriptions.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BigWig file. |
        """
        super(BigWigSignal, self).__init__(fn)
        self._readcount = None
        self.adapter = BigWigAdapter(fn)

    def mapped_read_count(self, force=False):
        """
        Count the total number of bases covered by the BigWig file.

        Uses pyBigWig's `header()['nBasesCovered']` — the number of
        bases that have *any* interval defined in the file, regardless
        of that interval's value. **This is not the same as the number
        of bases with non-zero signal**: a BigWig that stores explicit
        zero-value intervals (e.g. produced by `bamCoverage` without
        skipping no-coverage regions) counts those bases too, which can
        inflate the result up to the size of the whole genome.

        This value is used as the denominator for CPM normalization
        throughout `pipelines.py` (`value / (mapped_read_count() / 1e6)`)
        — an inflated count here directly shrinks the normalized signal
        the pipelines plot.

        The result is cached in `self._readcount` after the first call.
        If a `<file>.mmr` sidecar file exists, the count is read from
        it directly; otherwise the BigWig is scanned and the result is
        written to `<file>.mmr` for future use.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `force` | `bool` | `False` | If `True`, ignore the cache and re-count from the file. |

        ## Returns

        `int | None` — Total bases covered (`nBasesCovered`), or `None`
        if the file cannot be opened.

        ## Examples

        ```python
        from Pyntegrate import genomic_signal

        sig = genomic_signal("chip.bw", "bigwig")
        n = sig.mapped_read_count()
        print(f"Bases covered: {n:,}")
        ```

        **Output:**
        ```
        Bases covered: 15,243,600
        ```
        """
        if self._readcount and not force:
            return self._readcount

        if os.path.exists(self.fn + '.mmr') and not force:
            for line in open(self.fn + '.mmr'):
                if line.startswith('#'):
                    continue
                self._readcount = float(line.strip())
                return self._readcount

        try:
            bwFile = pyBigWig.open(self.fn)
        except RuntimeError:
            return self._readcount
        mapped_reads = bwFile.header()['nBasesCovered']
        if not os.path.exists(self.fn + '.mmr'):
            fout = open(self.fn + '.mmr', 'w')
            fout.write(str(mapped_reads) + '\n')
            fout.close()

        self._readcount = mapped_reads
        return self._readcount


class IntervalSignal(BaseSignal):
    """
    Abstract base class for interval-based genomic signal files.

    Adds `local_count` and `count_array` methods on top of `BaseSignal`,
    which are meaningful for BED, BAM, and BigBed formats where reads
    or intervals can be counted directly.

    Do not instantiate this class directly — use one of its concrete
    subclasses: `BamSignal`, `BedSignal`, or `BigBedSignal`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the interval-based signal file. |
    """

    def __init__(self, fn):
        """
        Store the path to the interval-based signal file.

        See the class docstring for the `fn` attribute description.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the interval-based signal file. |
        """
        BaseSignal.__init__(self, fn)

    def local_count(self, *args, **kwargs):
        return _local_count(self.adapter, *args, **kwargs)

    local_count.__doc__ = _local_count.__doc__

    def count_array(self, features, processes=None, chunksize=1, **kwargs):
        """
        Count genomic signal occurrences for each feature and return a 1-D array.

        Iterates over `features` and calls `local_count` for each one,
        optionally parallelising the computation. The result is a flat
        NumPy array with one count per feature.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `features` | `iterable` | — | Interval-like objects to count signal in. |
        | `processes` | `int \| None` | `None` | Number of parallel worker processes. `None` runs sequentially. |
        | `chunksize` | `int` | `1` | Features per worker chunk when `processes` is set. |
        | `**kwargs` | — | — | Forwarded to `local_count`. |

        ## Returns

        `np.ndarray` — 1-D integer array of per-feature counts, length
        equal to `len(features)`.

        ## Examples

        ```python
        import pybedtools
        from Pyntegrate import genomic_signal

        bam = genomic_signal("aligned.bam", "bam")
        genes = pybedtools.BedTool("genes.bed")

        counts = bam.count_array(genes)
        print(counts[:5])
        ```

        **Output:**
        ```
        [142  87 310  56 201]
        ```
        """
        if processes is not None:
            arrays = _count_array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)
        else:
            arrays = _count_array(self.fn, self.__class__, features, **kwargs)
        return np.concatenate(arrays)


class NarrowPeakSignal(IntervalSignal):
    """
    Genomic signal interface for NarrowPeak files (MACS2/MACS3 output).

    NarrowPeak is a BED6+4 format produced by peak callers such as
    MACS2. This class provides the same `array`, `local_coverage`,
    `local_count`, and `count_array` interface as other signal classes.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the NarrowPeak file (`.narrowPeak`). |
    | `adapter` | `NarrowPeakAdapter` | — | File-format-specific reader for the NarrowPeak file. |
    """

    def __init__(self, fn):
        """
        Open `fn` as a NarrowPeak file.

        See the class docstring for attribute descriptions.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the NarrowPeak file (`.narrowPeak`). |
        """
        BaseSignal.__init__(self, fn)
        self.adapter = NarrowPeakAdapter(self.fn)


class BamSignal(IntervalSignal):
    """
    Genomic signal interface for BAM (`.bam`) files.

    BAM files contain aligned short reads. This class wraps the
    `BamAdapter` and exposes coverage arrays, read counting, and a
    convenience method to extract the reference genome dictionary from
    the BAM header.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the BAM file. A `.bai` index must be present. |
    | `adapter` | `BamAdapter` | — | File-format-specific reader wrapping the BAM file. |
    | `_readcount` | `int \| float \| None` | `None` | Cached mapped-read count. `None` until populated by `mapped_read_count()`. |

    ## Examples

    ```python
    from Pyntegrate import genomic_signal

    bam = genomic_signal("aligned.bam", "bam")
    print(bam.genome())
    print(bam.mapped_read_count())
    ```

    **Output:**
    ```
    {'chr1': (0, 248956422), 'chr2': (0, 242193529), ...}
    42385761
    ```
    """

    def __init__(self, fn):
        """
        Open `fn` as a BAM file and initialize the read-count cache.

        See the class docstring for attribute descriptions.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BAM file. A `.bai` index must be present. |
        """
        BaseSignal.__init__(self, fn)
        self._readcount = None
        self.adapter = BamAdapter(self.fn)

    def genome(self):
        """
        Build a genome-size dictionary from the BAM file header.

        Returns a `{chrom: (0, length)}` dictionary compatible with
        `pybedtools` genome operations.

        ## Returns

        `dict[str, tuple[int, int]]` — Mapping of chromosome name to
        `(0, length)` tuples, derived from the BAM `@SQ` header lines.

        ## Examples

        ```python
        from Pyntegrate import genomic_signal

        bam = genomic_signal("aligned.bam", "bam")
        genome = bam.genome()
        print(genome["chr1"])
        ```

        **Output:**
        ```
        (0, 248956422)
        ```
        """
        f = self.adapter.fileobj
        d = {}
        for ref, length in zip(f.references, f.lengths):
            d[ref] = (0, length)
        return d

    def mapped_read_count(self, force=False):
        """
        Count the total number of mapped reads in the BAM file.

        Uses `samtools view -c -F 0x4` to count only mapped reads.
        The result is cached in `self._readcount` after the first call.
        If a `<file>.mmr` sidecar file exists the count is read from it
        directly; otherwise samtools is run and the result is written to
        `<file>.mmr` for future use.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `force` | `bool` | `False` | If `True`, ignore the cache and re-run samtools. |

        ## Returns

        `int | None` — Total mapped read count, or `None` if samtools
        returns an error.

        ## Examples

        ```python
        from Pyntegrate import genomic_signal

        bam = genomic_signal("aligned.bam", "bam")
        n = bam.mapped_read_count()
        print(f"Mapped reads: {n:,}")
        ```

        **Output:**
        ```
        Mapped reads: 42,385,761
        ```
        """
        if self._readcount and not force:
            return self._readcount

        if os.path.exists(self.fn + '.mmr') and not force:
            for line in open(self.fn + '.mmr'):
                if line.startswith('#'):
                    continue
                self._readcount = float(line.strip())
                return self._readcount

        cmds = ['samtools', 'view', '-c', '-F', '0x4', self.fn]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if p.returncode != 0:
            sys.stderr.write('samtools says: %s' % stderr.decode(errors='replace'))
            return None
        mapped_reads = int(stdout)

        if not os.path.exists(self.fn + '.mmr'):
            fout = open(self.fn + '.mmr', 'w')
            fout.write(str(mapped_reads) + '\n')
            fout.close()

        self._readcount = mapped_reads
        return self._readcount


class BigBedSignal(IntervalSignal):
    """
    Genomic signal interface for BigBed (`.bb`) files.

    BigBed is an indexed, binary version of BED. This class supports
    the same `array`, `local_coverage`, `local_count`, and
    `count_array` interface as other signal classes.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the BigBed file. |
    | `adapter` | `BigBedAdapter` | — | File-format-specific reader wrapping the BigBed file. |
    | `_readcount` | `int \| float \| None` | `None` | Cached mapped-read count. `None` until populated by `mapped_read_count()`. |
    """

    def __init__(self, fn):
        """
        Open `fn` as a BigBed file and initialize the read-count cache.

        See the class docstring for attribute descriptions.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BigBed file. |
        """
        IntervalSignal.__init__(self, fn)
        self._readcount = None
        self.adapter = BigBedAdapter(fn)

    def mapped_read_count(self, force=False):
        """
        Count the total number of non-zero intervals in the BigBed file.

        Analogous to `BigWigSignal.mapped_read_count`. The result is
        cached and persisted to a `<file>.mmr` sidecar file.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `force` | `bool` | `False` | If `True`, ignore the cache and re-count from the file. |

        ## Returns

        `int | None` — Total non-zero interval count, or `None` if the
        file cannot be opened.

        ## Examples

        ```python
        from Pyntegrate import genomic_signal

        bb = genomic_signal("peaks.bb", "bigbed")
        print(bb.mapped_read_count())
        ```

        **Output:**
        ```
        87432
        ```
        """
        if self._readcount and not force:
            return self._readcount

        if os.path.exists(self.fn + '.mmr') and not force:
            for line in open(self.fn + '.mmr'):
                if line.startswith('#'):
                    continue
                self._readcount = float(line.strip())
                return self._readcount

        mapped_reads = 0
        try:
            bbFile = pyBigWig.open(self.fn)
            bbFile.header()
        except RuntimeError:
            return self._readcount
        for chrom in bbFile.chroms():
            entries = bbFile.entries(chrom, 0, bbFile.chroms(chrom))
            if entries:
                mapped_reads += len(entries)
        if not os.path.exists(self.fn + '.mmr'):
            fout = open(self.fn + '.mmr', 'w')
            fout.write(str(mapped_reads) + '\n')
            fout.close()

        self._readcount = mapped_reads
        return self._readcount


class BedSignal(IntervalSignal):
    """
    Genomic signal interface for BED, GFF, GTF, and VCF files.

    Tabix-indexed interval files are accessed via `BedAdapter`. This
    class supports `array`, `local_coverage`, `local_count`, and
    `count_array` for any format registered under `"bed"`, `"gff"`,
    `"gtf"`, or `"vcf"` in the format registry.

    ## Attributes

    | Attribute | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the BED/GFF/VCF file (must be block-gzip compressed and Tabix-indexed for random access). |
    | `adapter` | `BedAdapter` | — | File-format-specific reader wrapping the BED/GFF/GTF/VCF file. |

    ## Examples

    ```python
    from Pyntegrate import genomic_signal

    bed = genomic_signal("peaks.bed.gz", "bed")
    arr = bed.array(genes, bins=100)
    print(arr.shape)
    ```

    **Output:**
    ```
    (5000, 100)
    ```
    """

    def __init__(self, fn):
        """
        Open `fn` as a BED/GFF/GTF/VCF file.

        See the class docstring for attribute descriptions.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BED/GFF/VCF file (must be block-gzip compressed and Tabix-indexed for random access). |
        """
        IntervalSignal.__init__(self, fn)
        self.adapter = BedAdapter(fn)


_registry = {
    'bam': BamSignal,
    'bed': BedSignal,
    'gff': BedSignal,
    'gtf': BedSignal,
    'vcf': BedSignal,
    'bigwig': BigWigSignal,
    'bigbed': BigBedSignal,
}
