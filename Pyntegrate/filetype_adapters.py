"""
Uniform, low-level file-reading adapters used internally by `Pyntegrate._genomic_signal`.

Each adapter wraps a single genomic file format behind the same minimal
API: the constructor takes a filename (its format is not validated),
`make_fileobj()` builds the underlying format-specific reader, and
`__getitem__` accepts a `pybedtools.Interval` and returns an iterator of
`pybedtools.Interval` objects overlapping it. End users should generally
go through `Pyntegrate._genomic_signal.genomic_signal` rather than these
classes directly.

## Adapters

| Class | Format | Backing library |
|-------|--------|------------------|
| `BaseAdapter` | — | Shared `fn`/`fileobj` machinery; not used directly. |
| `BamAdapter` | BAM (`.bam`) | `pysam` |
| `BedAdapter` | BED / GFF / GTF / VCF | `pybedtools.IntervalFile` |
| `BigBedAdapter` | BigBed (`.bb`) | `bx-python` (`BigBedFile`) |
| `BigWigAdapter` | BigWig (`.bw`) | `bx-python` (`BigWigFile`) |
| `NarrowPeakAdapter` | narrowPeak | `pybedtools` |

## Extending

Subclasses must override `make_fileobj()` (returns the object iterated
over in `__getitem__`) and `__getitem__` (accepts a `pybedtools.Interval`,
returns an iterator of `pybedtools.Interval` objects).
"""
from bx.bbi.bigbed_file import BigBedFile
from bx.bbi.bigwig_file import BigWigFile
from bx.intervals.io import StrandFormatError
import numpy as np
import subprocess
import pysam
import pybedtools
import os
import sys
from textwrap import dedent

strand_lookup = {16: '-', 0: '+'}


class BaseAdapter(object):
    """
    Base class for filetype adapters.

    Provides a uniform, format-agnostic interface for random access to
    genomic interval data. Subclasses must define `__getitem__` (accept a
    `pybedtools.Interval` and return an iterator of `pybedtools.Interval`
    objects) and `make_fileobj` (return the format-specific object to be
    iterated over in `__getitem__`).

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the underlying file. |
    | `fileobj` | `object` | Format-specific file handle/reader returned by `make_fileobj()`. |
    """
    def __init__(self, fn):
        """
        Store the file path and open the format-specific file object.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the underlying file. |
        """
        self.fn = fn
        self.fileobj = None
        self.fileobj = self.make_fileobj()

    def __getitem__(self, key):
        """
        Return an iterator of intervals overlapping `key`.

        Must be implemented by subclasses.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `pybedtools.Interval` | — | Genomic interval to query. |

        ## Returns

        `iterator[pybedtools.Interval]` — Intervals overlapping `key`.
        """
        raise ValueError('Subclasses must define __getitem__')

    def make_fileobj(self):
        """
        Create and return the format-specific file object to read from.

        Must be implemented by subclasses.

        ## Returns

        `object` — Format-specific file handle/reader, later stored as
        `self.fileobj`.
        """
        raise ValueError('Subclasses must define make_fileobj')

class NarrowPeakAdapter(BaseAdapter):
    """
    Adapter that provides random access to narrowPeak files using `pybedtools.BedTool`.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the narrowPeak file (inherited from `BaseAdapter`). |
    | `fileobj` | `pybedtools.BedTool` | Same object as `bedtool`; set via `make_fileobj()` (inherited from `BaseAdapter`). |
    | `bedtool` | `pybedtools.BedTool \| None` | `BedTool` wrapping the narrowPeak file. `None` until lazily created by `make_fileobj()`. |
    | `chroms` | `set[str] \| None` | Chromosome names present in the file. `None` until lazily populated by `make_fileobj()` via `_get_chroms()`. |
    """
    def __init__(self, fn):
        """
        Store the file path and lazily-initialised state, then open the file.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the narrowPeak file. |
        """
        self.bedtool = None
        self.chroms = None
        super(NarrowPeakAdapter, self).__init__(fn)



    def make_fileobj(self):
        """
        Lazily create `self.bedtool` and populate `self.chroms`.

        ## Returns

        `pybedtools.BedTool` — `BedTool` wrapping the narrowPeak file
        (same object as `self.bedtool`).
        """
        if self.bedtool is None:
            self.bedtool = pybedtools.BedTool(self.fn)
        if self.chroms is None:
            self.chroms = self._get_chroms()
        return self.bedtool

    def _get_chroms(self):
        """
        Collect the set of chromosome names present in the narrowPeak file.

        ## Returns

        `set[str]` — Chromosome names found in `self.bedtool`.
        """
        chroms = set()
        for interval in self.bedtool:
            chroms.add(interval.chrom)
        return chroms

    def summarize(self, interval, bins=None, function='mean', zero_inf=True, zero_nan=True):
        """
        Compute binned, aggregated peak values over a genomic interval.

        Intersects `self.bedtool` with `interval`, then aggregates the
        narrowPeak signal value (column 7) of overlapping peaks into
        `bins` evenly-sized bins using `function`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `interval` | `object` | — | Object with `chrom` (`str`), `start` (`int`), and `stop` (`int`) attributes. |
        | `bins` | `int \| None` | `None` | Number of bins. If `None`, uses the length of `interval`. |
        | `function` | `str` | `"mean"` | Aggregation function: `"mean"`, `"sum"`, `"min"`, or `"max"`. |
        | `zero_inf` | `bool` | `True` | If `True`, replace `±inf` values with zero before returning. |
        | `zero_nan` | `bool` | `True` | If `True`, replace `NaN` values with zero before returning. |

        ## Returns

        `np.ndarray` — Array of shape `(bins,)` with the aggregated
        signal value per bin.
        """
        # Set bins to interval length if not specified
        # Ensure the BedTool object is ready
        self.make_fileobj()

        # Rest of the summarize function implementation...
        if bins is None:
            bins = interval.stop - interval.start

        s = np.zeros((bins,))
        seen = np.zeros((bins,), dtype=bool)
        bin_size = (interval.stop - interval.start) / bins
        bin_edges = np.linspace(interval.start, interval.stop, bins + 1)

        query_interval = pybedtools.BedTool(
            [pybedtools.Interval(interval.chrom, interval.start, interval.stop)]
        )
        overlapping_intervals = self.bedtool.intersect(query_interval, wa=True)

        for peak in overlapping_intervals:
            peak_start = max(interval.start, peak.start)
            peak_end = min(interval.stop, peak.end)

            bin_start_idx = int((peak_start - interval.start) // bin_size)
            bin_end_idx = min(int((peak_end - interval.start) // bin_size), bins - 1)

            for i in range(bin_start_idx, bin_end_idx + 1):
                bin_start = bin_edges[i]
                bin_end = bin_edges[i + 1]
                overlap_start = max(peak_start, bin_start)
                overlap_end = min(peak_end, bin_end)
                overlap_length = overlap_end - overlap_start

                if overlap_length > 0:
                    if function == 'mean':
                        s[i] += float(peak[6]) * overlap_length / bin_size
                    elif function == 'sum':
                        s[i] += float(peak[6])
                    elif function == 'min':
                        if not seen[i]:
                            s[i] = float(peak[6])
                            seen[i] = True
                        else:
                            s[i] = min(s[i], float(peak[6]))
                    elif function == 'max':
                        if not seen[i]:
                            s[i] = float(peak[6])
                            seen[i] = True
                        else:
                            s[i] = max(s[i], float(peak[6]))

        if zero_nan:
            s[np.isnan(s)] = 0
        if zero_inf:
            s[np.isinf(s)] = 0

        return s
    
    def __getitem__(self, key):
        """
        Yield intervals in the narrowPeak file that overlap `key`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `object` | — | Object with `chrom` (`str`), `start` (`int`), and `stop` (`int`) attributes. |

        ## Returns

        `iterator[pybedtools.Interval]` — Intervals overlapping `key`,
        each with `file_type` set to `"bed"`. Returns nothing if
        `key.chrom` is not present in `self.chroms`.
        """
        # Create a pybedtools.Interval for the requested region
        if key.chrom not in self.chroms:
            return
        query_interval = pybedtools.BedTool(
            [pybedtools.Interval(key.chrom, key.start, key.stop)]
        )
        overlapping_intervals = self.fileobj.intersect(query_interval, wa=True)

        # Yield each overlapping interval as a pybedtools.Interval
        for interval in overlapping_intervals:
            interval.file_type = 'bed'
            yield interval




class BamAdapter(BaseAdapter):
    """
    Adapter that provides random access to BAM files using Pysam.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the BAM file (inherited from `BaseAdapter`). |
    | `fileobj` | `pysam.Samfile` | Open BAM file handle returned by `make_fileobj()` (inherited from `BaseAdapter`). |
    """
    def __init__(self, fn):
        """
        Store the file path and open the BAM file with Pysam.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BAM file. |
        """
        super(BamAdapter, self).__init__(fn)

    def make_fileobj(self):
        """
        Open the BAM file for reading.

        ## Returns

        `pysam.Samfile` — Open BAM file handle.
        """
        return pysam.Samfile(self.fn, 'rb')

    def __getitem__(self, key):
        """
        Yield aligned-block intervals overlapping `key`.

        Each CIGAR match block (`op == 0`) of every read overlapping
        `key` is yielded as its own interval, with strand derived from
        the SAM flag.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `pybedtools.Interval` | — | Genomic interval to query (`chrom`, `start`, `stop`). |

        ## Returns

        `iterator[pybedtools.Interval]` — One interval per aligned
        match block, each with `file_type` set to `"bed"`.
        """
        iterator = self.fileobj.fetch(
            str(key.chrom),
            key.start,
            key.stop)
        for r in iterator:
            start = r.pos
            curr_end = r.pos
            for op, bp in r.cigar:
                start = curr_end
                curr_end += bp
                if op == 0:
                    interval = pybedtools.Interval(
                        self.fileobj.references[r.rname],
                        start,
                        curr_end,
                        strand=strand_lookup[r.flag & 0x0010])
                    interval.file_type = 'bed'
                    yield interval


class BedAdapter(BaseAdapter):
    """
    Adapter that provides random access to BED files via `pybedtools.IntervalFile`.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the BED file (inherited from `BaseAdapter`). |
    | `fileobj` | `pybedtools.IntervalFile` | Indexed interval file returned by `make_fileobj()` (inherited from `BaseAdapter`). |
    """
    def __init__(self, fn):
        """
        Store the file path and open the BED file as an `IntervalFile`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the BED file. |
        """
        super(BedAdapter, self).__init__(fn)

    def make_fileobj(self):
        """
        Open the BED file as an indexed `pybedtools.IntervalFile`.

        ## Returns

        `pybedtools.IntervalFile` — Indexed interval file.
        """
        return pybedtools.IntervalFile(self.fn)

    def __getitem__(self, key):
        """
        Yield intervals in the BED file overlapping `key`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `pybedtools.Interval` | — | Genomic interval to query. |

        ## Returns

        `iterator[pybedtools.Interval]` — Intervals overlapping `key`.
        """
        for i in self.fileobj.all_hits(key):
            yield i


class BigBedAdapter(BaseAdapter):
    """
    Adapter that provides random access to bigBed files via bx-python.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the bigBed file (inherited from `BaseAdapter`). |
    | `fileobj` | `bx.bbi.bigbed_file.BigBedFile` | Open bigBed file handle returned by `make_fileobj()` (inherited from `BaseAdapter`). |
    """
    def __init__(self, fn):
        """
        Store the file path and open the bigBed file.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the bigBed file. |
        """
        super(BigBedAdapter, self).__init__(fn)

    def make_fileobj(self):
        """
        Open the bigBed file for reading.

        ## Returns

        `bx.bbi.bigbed_file.BigBedFile` — Open bigBed file handle.
        """
        return BigBedFile(open(self.fn, 'rb'))

    def __getitem__(self, key):
        """
        Yield intervals in the bigBed file overlapping `key`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `pybedtools.Interval` | — | Genomic interval to query (`chrom`, `start`, `end`). |

        ## Returns

        `iterator[pybedtools.Interval]` — Intervals overlapping `key`,
        each with `file_type` set to `"bed"`.

        ## Raises

        `NotImplementedError` — If the installed bx-python version does
        not support bigBed files (a workaround via `bigBedToBed` /
        `bedtools bedtobam` is suggested in the error message).
        """
        chrom = key.chrom
        start = key.start
        stop = key.end
        try:
            bx_intervals = self.fileobj.get(chrom, start, stop)
        except StrandFormatError:
            raise NotImplementedError(dedent(
                """
                It appears you have a version of bx-python where bigBed files
                are temporarily unsupported due to recent changes in the
                bx-python dependency. In the meantime, please convert bigBed to
                BAM like this:

                    bigBedToBed {0} tmp.bed
                    bedtools bedtobam -i tmp.bed > {0}.bam

                and create a genomic signal object using this {0}.bam file.
                """.format(self.fn)))
        if bx_intervals is None:
            return
        for i in bx_intervals:
            interval = pybedtools.create_interval_from_list(i.fields)
            interval.file_type = 'bed'
            yield interval


class BigWigAdapter(BaseAdapter):
    """
    Adapter that provides random access to bigWig files via bx-python.

    ## Attributes

    | Attribute | Type | Description |
    |-----------|------|-------------|
    | `fn` | `str` | Path to the bigWig file (inherited from `BaseAdapter`). |
    | `fileobj` | `bx.bbi.bigwig_file.BigWigFile` | Open bigWig file handle returned by `make_fileobj()` (inherited from `BaseAdapter`). |
    """
    def __init__(self, fn):
        """
        Store the file path and open the bigWig file.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `fn` | `str` | — | Path to the bigWig file. |
        """
        super(BigWigAdapter, self).__init__(fn)

    def make_fileobj(self):
        """
        Open the bigWig file for reading.

        ## Returns

        `bx.bbi.bigwig_file.BigWigFile` — Open bigWig file handle.
        """
        return BigWigFile(open(self.fn, 'rb'))

    def __getitem__(self, key):
        """
        Not supported for bigWig files; use `summarize` instead.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `key` | `pybedtools.Interval` | — | Unused. |

        ## Raises

        `NotImplementedError` — Always. BigWig files store pre-computed
        coverage rather than discrete intervals, so random-access
        iteration is not meaningful; use `summarize` to get binned
        coverage values instead.
        """
        raise NotImplementedError(
            "__getitem__ not implemented for %s" % self.__class__.__name__)

    def summarize(self, interval, bins=None, method='summarize',
                  function='mean', zero_inf=True, zero_nan=True):
        """
        Compute binned, aggregated bigWig values over a genomic interval.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `interval` | `object` | — | Object with `chrom` (`str`), `start` (`int`), and `stop` (`int`) attributes. |
        | `bins` | `int \| None` | `None` | Number of bins. If `None`, uses the length of `interval`. |
        | `method` | `str` | `"summarize"` | `"summarize"` and `"get_as_array"` use bx-python; `"ucsc_summarize"` uses `bigWigSummary`. See `Pyntegrate.array_helpers._local_coverage` for further notes. `None` is treated as `"summarize"`. |
        | `function` | `str` | `"mean"` | Aggregation function: `"mean"`, `"sum"`, `"min"`, `"max"`, `"std"`, or `"coverage"`. Ignored when `method="get_as_array"`; `"coverage"` is only valid with `method="ucsc_summarize"`. |
        | `zero_inf` | `bool` | `True` | If `True`, replace `±inf` values with zero before returning. |
        | `zero_nan` | `bool` | `True` | If `True`, replace `NaN` values with zero before returning. |

        ## Returns

        `np.ndarray` — Aggregated signal values, one per bin.

        ## Raises

        `ValueError` — If `function` is not supported by the chosen
        `method`, or if `method` is not one of `"summarize"`,
        `"ucsc_summarize"`, or `"get_as_array"`.
        """

        if method is None:
            method = 'summarize'

        # We may be dividing by zero in some cases, which raises a warning in
        # NumPy based on the IEEE 754 standard (see
        # http://docs.scipy.org/doc/numpy/reference/generated/
        #       numpy.seterr.html)
        #
        # That's OK -- we're expecting that to happen sometimes. So temporarily
        # disable this error reporting for the duration of this method.
        orig = np.geterr()['invalid']
        np.seterr(invalid='ignore')

        if (bins is None) or (method == 'get_as_array'):
            bw = self.fileobj
            chrom = interval.chrom
            if isinstance(chrom, str):
                chrom = chrom.encode()
            s = bw.get_as_array(
                chrom,
                interval.start,
                interval.stop,)
            if s is None:
                s = np.zeros((interval.stop - interval.start,))
            else:
                if zero_nan:
                    s[np.isnan(s)] = 0
                if zero_inf:
                    s[np.isinf(s)] = 0

        elif method == 'ucsc_summarize':
            if function in ['mean', 'min', 'max', 'std', 'coverage']:
                return self.ucsc_summarize(interval, bins, function=function)
            else:
                raise ValueError('function "%s" not supported by UCSC\'s'
                                 'bigWigSummary')

        elif method == 'summarize':
            bw = self.fileobj
            s = bw.summarize(
                interval.chrom,
                interval.start,
                interval.stop, bins)

            if s is None:
                s = np.zeros((bins,))
            else:
                if function == 'sum':
                    s = s.sum_data
                elif function == 'mean':
                    s = s.sum_data / s.valid_count
                    if zero_nan:
                        s[np.isnan(s)] = 0
                elif function == 'min':
                    s = s.min_val
                    if zero_inf:
                        s[np.isinf(s)] = 0
                elif function == 'max':
                    s = s.max_val
                    if zero_inf:
                        s[np.isinf(s)] = 0
                elif function == 'std':
                    s = (s.sum_squares / s.valid_count)
                    if zero_nan:
                        s[np.isnan(s)] = 0
                else:
                    raise ValueError(
                            'function "%s" not supported by bx-python'
                            % function
                    )
            
        else:
            raise ValueError("method '%s' not in [summarize, ucsc_summarize, get_as_array]" % method)

        # Reset NumPy error reporting
        np.seterr(divide=orig)
        return s

    def ucsc_summarize(self, interval, bins=None, function='mean'):
        """
        Summarize a bigWig interval by shelling out to UCSC's `bigWigSummary`.

        ## Parameters

        | Parameter | Type | Default | Description |
        |-----------|------|---------|-------------|
        | `interval` | `object` | — | Object with `chrom` (`str`), `start` (`int`), and `stop` (`int`) attributes. |
        | `bins` | `int \| None` | `None` | Number of bins. If `None`, uses `len(interval)`. |
        | `function` | `str` | `"mean"` | Aggregation type passed to `bigWigSummary -type=`, e.g. `"mean"`, `"min"`, `"max"`, `"std"`, or `"coverage"`. |

        ## Returns

        `np.ndarray` — Array of shape `(bins,)` with the summarized
        value per bin. Bins with no data from `bigWigSummary` remain `0`.

        ## Raises

        `ValueError` — If the `bigWigSummary` subprocess exits with a
        non-zero, non-`None` return code (other than the "no data" case).
        """
        if bins is None:
            bins = len(interval)
        y = np.zeros(bins)
        cmds = [
            'bigWigSummary',
            self.fn,
            interval.chrom,
            str(interval.start),
            str(interval.stop),
            str(bins),
            '-type=%s' % function]
        p = subprocess.Popen(
            cmds,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        def gen():
            try:
                for line in p.stdout:
                    yield line
            finally:
                if p.poll() is None:
                    return
                else:
                    p.wait()
                    err = p.stderr.read().strip()
                    if p.returncode not in (0, None):
                        if err.startswith('no data'):
                            return
                        raise ValueError(
                            "cmds: %s: %s" %
                            (' '.join(cmds), p.stderr.read()))
                    if len(err) != 0:
                        sys.stderr.write(err)

        for line in gen():
            for i, x in enumerate(line.split('\t')):
                try:
                    y[i] = float(x)
                except ValueError:
                    pass
        return np.array(y)
