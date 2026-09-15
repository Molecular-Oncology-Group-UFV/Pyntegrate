"""
Low-level, format-agnostic helpers for extracting genomic signal into NumPy arrays.

These functions power the public API in `Pyntegrate._genomic_signal`
(`BaseSignal.array`, `local_coverage`, `IntervalSignal.local_count`,
`count_array`, etc.) by doing the actual per-feature/per-window
computation — they are not usually called directly by end users.

## Contents

| Function | Description |
|----------|-------------|
| `_local_count` | Count signal items (reads/intervals) within a single genomic interval. |
| `_local_coverage` | Compute a binned coverage vector for one or more genomic intervals. |
| `_array` / `_array_star` / `_array_parallel` | Build per-feature signal arrays, optionally distributed across multiple processes. |
| `_count_array` / `_count_array_star` / `_count_array_parallel` | Build per-feature signal counts, optionally distributed across multiple processes. |
| `remove_duplicates` | Deduplicate `gene_id` entries in a GTF file, overwriting it in place. |

## Exceptions

`ArgumentError` is raised when an argument combination is invalid for the
underlying file format (e.g. passing `use_score` with a BAM file).
"""
import pybedtools
import numpy as np
import multiprocessing
import itertools
from .helpers import *
from .filetype_adapters import *

class ArgumentError(Exception):
    """Raised when an invalid or unsupported combination of arguments is passed for a given file format."""


def _local_count(reader, feature, stranded=False):
    """
    Count genomic signal items (reads or intervals) found within a single
    genomic interval.

    Typically meaningful for BED and BAM files; not intended for BigWig, which
    stores pre-computed coverage values rather than discrete reads.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `reader` | `BaseAdapter` | — | File adapter (e.g. `BamAdapter`) exposing `__getitem__` for random access. |
    | `feature` | `pybedtools.Interval \| str` | — | Genomic interval to query. Strings are parsed by `tointerval`. |
    | `stranded` | `bool` | `False` | If `True`, only count signal items whose strand matches `feature.strand`. |

    ## Returns

    `int` — Number of signal items overlapping `feature`.

    ## Examples

    ```python
    import pybedtools
    from Pyntegrate.filetype_adapters import BamAdapter
    from Pyntegrate.array_helpers import _local_count

    reader = BamAdapter("aligned.bam")
    feature = pybedtools.create_interval_from_list(
        ["chr1", "1000", "5000", ".", ".", "+"]
    )
    count = _local_count(reader, feature, stranded=True)
    print(count)
    ```

    **Output:**
    ```
    142
    ```
    """
    if isinstance(feature, str):
        feature = tointerval(feature)
    if stranded:
        strand = feature.strand
    else:
        strand = '.'

    count = 0
    for al in reader[feature]:
        if stranded and al.strand != strand:
            continue
        count += 1
    return count


def _local_coverage(reader, features, read_strand=None, fragment_size=None,
                    shift_width=0, bins=None, use_score=False, accumulate=True,
                    preserve_total=False, method=None, function="mean",
                    zero_inf=True, zero_nan=True, processes=None,
                    stranded=True, verbose=False):
    """
    Compute a binned coverage vector for one or more genomic intervals.

    Returns `(x, y)` where `x` contains genomic coordinates and `y`
    contains signal values. For a single interval `y` is a 1-D profile;
    for multiple intervals `y` is their concatenation.

    Columns marked **BigWig only** or **Not BigWig** indicate format
    restrictions.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `reader` | `BaseAdapter` | — | File adapter providing random-access to the signal file. |
    | `features` | `str \| interval-like \| list` | — | Single interval, coordinate string (`"chr1:1000-2000"`), or list of intervals. Lists produce a concatenated profile. |
    | `bins` | `int \| list[int] \| None` | `None` | Number of bins for the output. `None` → one value per base pair. For a list of features, provide a list of bin counts of the same length. *(BigWig supported)* |
    | `fragment_size` | `int \| None` | `None` | Extend each read this many bp in the 3′ direction for smoother signal. *(Not BigWig)* |
    | `shift_width` | `int` | `0` | Shift reads this many bp in the 3′ direction before computing coverage. Useful for reconstructing ChIP-seq peaks using the MACS `d` parameter. *(Not BigWig)* |
    | `read_strand` | `str \| None` | `None` | If `"+"` or `"-"`, only count reads on that strand. *(Not BigWig)* |
    | `stranded` | `bool` | `True` | If `True`, reverse the profile for minus-strand features so that all profiles are TSS → TTS. |
    | `use_score` | `bool` | `False` | If `True`, accumulate the BED *score* field instead of read counts. *(Not BigWig, not BAM)* |
    | `accumulate` | `bool` | `True` | If `False`, mark bins as 1/0 rather than summing reads (useful for peak matrices). *(BigWig supported)* |
    | `preserve_total` | `bool` | `False` | If `True`, rescale binned values so that their sum equals the un-binned total. *(BigWig supported, except with `method="ucsc_summarize"`)* |
    | `method` | `str \| None` | `None` | BigWig extraction strategy: `"summarize"` (fast, slight approximation), `"get_as_array"` (exact, slower), or `"ucsc_summarize"` (uses UCSC `bigWigSummary`). *(BigWig only)* |
    | `function` | `str` | `"mean"` | Aggregation function for BigWig bins: `"sum"`, `"mean"`, `"min"`, `"max"`, `"std"`, or `"coverage"` (with `ucsc_summarize`). *(BigWig only)* |
    | `zero_inf` | `bool` | `True` | Replace `±inf` values with zero. *(BigWig only)* |
    | `zero_nan` | `bool` | `True` | Replace `NaN` values with zero. *(BigWig only)* |
    | `processes` | `int \| None` | `None` | Split the feature across multiple processes. |

    ## Returns

    `tuple[np.ndarray, dict]` — `(x, profile)` where `x` is a 1-D array
    of genomic coordinates and `profile` is a dictionary with keys
    `"values"`, `"gene_name"`, `"chr"`, `"start"`, `"end"`, `"strand"`,
    and `"id"`. `"values"` is the concatenation of all features' signal
    when `features` is a list, but the other keys describe only the
    **last** feature in that list — they are only meaningful for
    identifying the feature when `features` is a single interval (the
    normal usage; see all call sites in this codebase). For a list of
    distinct features, call this function once per feature instead of
    relying on the returned metadata to identify which segment of
    `"values"` belongs to which feature.

    ## Notes

    - Minus-strand features are reversed so that index 0 always
      corresponds to the 5′ end of the feature.
    - When `bins` is set and `preserve_total=True`, the returned values
      have units of *total signal per bin*; otherwise units are *density*.

    ## Examples

    ```python
    import pybedtools
    from Pyntegrate.filetype_adapters import BigWigAdapter
    from Pyntegrate.array_helpers import _local_coverage

    adapter = BigWigAdapter("chip.bw")
    feature = pybedtools.create_interval_from_list(
        ["chr1", "1000", "2000", "GENE1", ".", "+", ".", ".", "gene_id \\"GENE1\\""]
    )
    x, profile = _local_coverage(adapter, feature, bins=100)
    print(x.shape, profile["values"].shape)
    ```

    **Output:**
    ```
    (100,) (100,)
    ```
    """
    # bigWig files are handled differently, so we need to know if we're working
    # with one; raise exception if a kwarg was supplied that's not supported.
    is_narrow = False
    if isinstance(reader, BigWigAdapter) :
        is_bigwig = True
        defaults = (
            ('read_strand', read_strand, None),
            ('fragment_size', fragment_size, None),
            ('shift_width', shift_width, 0),
            ('use_score', use_score, False),
            ('preserve_total', preserve_total, False),
        )
        for name, check, default in defaults:
            if (
                ((default is None) and (check is not default))
                or
                (check != default)
            ):
                raise ArgumentError(
                    "Argument '%s' not supported for bigWig" % name)

        if method == 'ucsc_summarize':
            if preserve_total:
                raise ArgumentError(
                    "preserve_total=True not supported when using "
                    "method='ucsc_summarize'")
    else:
        is_bigwig = False

    if isinstance(reader, BamAdapter):
        if use_score:
            raise ArgumentError("Argument 'use_score' not supported for "
                                "bam")
    if isinstance(reader,NarrowPeakAdapter ):
        is_narrow = True
    # e.g., features = "chr1:1-1000"
    if isinstance(features, str):
        features = tointerval(features)

    if not ((isinstance(features, list) or isinstance(features, tuple))):
        if bins is not None:
            if not isinstance(bins, int):
                raise ArgumentError(
                    "bins must be an int, got %s" % type(bins))
        features = [features]
        bins = [bins]
    else:
        if bins is None:
            bins = [None for i in features]
        if not len(bins) == len(features):
            raise ArgumentError(
                "bins must have same length as feature list")

    # nomenclature:
    #   "window" is region we're getting data for
    #   "alignment" is one item in that region
    #
    profiles = []
    xs = []
    for window, nbin in zip(features, bins):
        window = tointerval(window)
        chrom = window.chrom
        start = window.start
        stop = window.stop
        strand = window.strand


        if not is_bigwig:
            # Extend the window to catch reads that would extend into the
            # requested window
            _fs = fragment_size or 0
            padded_window = pybedtools.Interval(
                chrom,
                max(start - _fs - shift_width, 0),
                stop + _fs + shift_width,
            )
            window_size = stop - start

            # start off with an array of zeros to represent the window
            profile = np.zeros(window_size, dtype=float)
            for interval in reader[padded_window]:
                if read_strand:
                    if interval.strand != read_strand:
                        continue

                # Shift interval by modeled distance, if specified.
                if shift_width:
                    if interval.strand == '-':
                        interval.start -= shift_width
                        interval.stop -= shift_width
                    else:
                        interval.start += shift_width
                        interval.stop += shift_width

                # Extend fragment size from 3'
                if fragment_size:
                    if interval.strand == '-':
                        interval.start = interval.stop - fragment_size
                    else:
                        interval.stop = interval.start + fragment_size

                # Convert to 0-based coords that can be used as indices into
                # array
                start_ind = interval.start - start

                # If the feature goes out of the window, then only include the
                # part that's inside the window
                start_ind = max(start_ind, 0)

                # Same thing for stop
                stop_ind = interval.stop - start
                stop_ind = min(stop_ind, window_size)


                # Skip if the feature is shifted outside the window. This can
                # happen with large values of `shift_width`.
                if start_ind >= window_size or stop_ind < 0:
                    continue

                # Finally, increment profile
                if use_score:
                    score = float(interval.score)
                else:
                    score = 1

                if accumulate:
                    if preserve_total:
                        profile[start_ind:stop_ind] += (
                            score / float((stop_ind - start_ind)))
                    else:
                        profile[start_ind:stop_ind] += score

                else:
                    profile[start_ind:stop_ind] = score
        else:  # it's a bigWig
            gene_name = window.name
            profile = reader.summarize(
                window,
                method=method,
                function=function,
                bins=(nbin or len(window)),
                zero_inf=zero_inf,
                zero_nan=zero_nan,
                )

        # If no bins, return genomic coords
        if (nbin is None):
            x = np.arange(start, stop)

        # Otherwise do the downsampling; resulting x is stll in genomic
        # coords
        else:
            if preserve_total:
                total = float(profile.sum())
            if not is_bigwig or method == 'get_as_array':
                ## rebin() function, defined in helpers.py
                xi, profile = rebin(
                    x=np.arange(start, stop), y=profile, nbin=nbin)
                if not accumulate:
                    profile[profile != 0] = 1
                x = xi

            else:
                x = np.linspace(start, stop - 1, nbin)

        # Minus-strand profiles should be flipped left-to-right.
        if stranded and strand == '-':
            profile = profile[::-1]
        xs.append(x)
        if preserve_total:
            scale = profile.sum() / total
            profile /= scale
        profiles.append(profile)


    stacked_xs = np.hstack(xs)
    stacked_profiles = np.hstack(profiles)
    if(window.name):
        name = window.name
    else:
        name = "."

    gene_id = window.attrs.get("gene_id")
    id = gene_id if gene_id else "."
    stacked_profiles = {"values": stacked_profiles, 
                        "gene_name": name, 
                        "chr": window.chrom, 
                        "start": window.start, 
                        "end": window.end, 
                        "strand": window.strand,
                        "id": id}
    del xs
    del profiles
    return stacked_xs, stacked_profiles


def _array_parallel(fn, cls, genelist, chunksize=250, processes=1, **kwargs):
    """
    Build per-gene signal arrays in parallel using a multiprocessing pool.

    Splits `genelist` into chunks of size `chunksize` and distributes
    them across `processes` workers. Each worker calls `_array` on its
    chunk and the results are collected but **not** concatenated here —
    the caller is responsible for stacking.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the genomic signal file. |
    | `cls` | `type` | — | Signal class (e.g. `BigWigSignal`) used to open `fn` inside each worker. |
    | `genelist` | `iterable` | — | Iterable of `pybedtools.Interval` objects to process. |
    | `chunksize` | `int` | `250` | Features per worker chunk. Values of 25–100 work well on 8 cores. |
    | `processes` | `int` | `1` | Number of parallel worker processes. |
    | `**kwargs` | — | — | Forwarded to `_array` (e.g. `bins`, `fragment_size`). |

    ## Returns

    `list[list[dict]]` — Nested list: one sub-list per chunk, each
    containing the `_local_coverage` result dictionaries for that chunk.
    """
    # genelist = remove_duplicates(genelist)
    pool = multiprocessing.Pool(processes)
    chunks = list(chunker(genelist, chunksize))
    # pool.map can only pass a single argument to the mapped function, so you
    # need this trick for passing multiple arguments; idea from
    # http://stackoverflow.com/questions/5442910/
    #               python-multiprocessing-pool-map-for-multiple-arguments

    results = pool.map(
        func=_array_star,
        iterable=zip(
            itertools.repeat(fn),
            itertools.repeat(cls),
            chunks,
            itertools.repeat(kwargs)))
    
    pool.close()
    pool.join()
    
    return results

def _count_array_parallel(fn, cls, genelist, chunksize=250, processes=1, **kwargs):
    """
    Parallel version of `_count_array` using a multiprocessing pool.

    Splits `genelist` into chunks of size `chunksize` and distributes
    them across `processes` workers, each calling `_count_array` on its
    chunk.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the genomic signal file. |
    | `cls` | `type` | — | Adapter class (e.g. `BamAdapter`) used to read `fn`. |
    | `genelist` | `iterable` | — | Iterable of genomic features (intervals) to count signal in. |
    | `chunksize` | `int` | `250` | Number of features per worker chunk. |
    | `processes` | `int` | `1` | Number of parallel worker processes. |
    | `**kwargs` | — | — | Forwarded to `_count_array`. |

    ## Returns

    `list[list]` — Nested list of per-feature counts, one sub-list per chunk.
    """
    pool = multiprocessing.Pool(processes)
    chunks = list(chunker(genelist, chunksize))
    results = pool.map(
        func=_count_array_star,
        iterable=zip(
            itertools.repeat(fn),
            itertools.repeat(cls),
            chunks,
            itertools.repeat(kwargs)))
    pool.close()
    pool.join()
    return results


def _count_array_star(args):
    """
    Unpack `args` and delegate to `_count_array`.

    Required by `multiprocessing.Pool.map`, which only accepts
    single-argument callables.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `args` | `tuple` | — | `(fn, cls, genelist, kwargs)` as produced by `_count_array_parallel`. |

    ## Returns

    `list` — Per-feature count list returned by `_count_array`.
    """
    fn, cls, genelist, kwargs = args
    return _count_array(fn, cls, genelist, **kwargs)

def _count_array(fn, cls, genelist, **kwargs):
    """
    Count genomic signal occurrences for each feature in `genelist`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the genomic signal file. |
    | `cls` | `type` | — | Adapter class used to open `fn` (must expose a `local_count` method). |
    | `genelist` | `iterable` | — | Iterable of genomic features (intervals) to count signal in. |
    | `**kwargs` | — | — | Forwarded to `cls.local_count`. |

    ## Returns

    `list` — List of count values, one per feature in `genelist`.
    """
    reader = cls(fn)
    _local_count_func = cls.local_count
    biglist = []
    for gene in genelist:
        c = _local_count_func(
            reader, gene, **kwargs)
        biglist.append(c)
    return biglist

def _array_star(args):
    """
    Unpack `args` and delegate to `_array`.

    Required by `multiprocessing.Pool.map`, which only accepts
    single-argument callables.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `args` | `tuple` | — | `(fn, cls, genelist, kwargs)` as produced by `_array_parallel`. |

    ## Returns

    `list` — Result returned by `_array`.
    """
    fn, cls, genelist, kwargs = args
    return _array(fn, cls, genelist, **kwargs)


def _array(fn, cls, genelist, **kwargs):
    """
    Build a "meta-feature" array of per-bin signal values for each feature.

    Calls `cls.local_coverage` on each feature in `genelist` and
    collects the resulting per-bin values.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fn` | `str` | — | Path to the genomic signal file. |
    | `cls` | `type` | — | Signal class (e.g. `BigWigSignal`) used to open `fn`. |
    | `genelist` | `iterable` | — | Iterable of genomic features (intervals) to compute signal over. |
    | `**kwargs` | — | — | Forwarded to `cls.local_coverage` (e.g. `bins`, `fragment_size`). |

    ## Returns

    `list[np.ndarray]` — One 1-D array of per-bin signal values per
    feature in `genelist`.
    """
    reader = cls(fn)

    _local_coverage_func = cls.local_coverage
    biglist = []
    if 'bins' in kwargs:
        if isinstance(kwargs['bins'], int):
            kwargs['bins'] = [kwargs['bins']]

    for gene in genelist:
        if not isinstance(gene, (list, tuple)):
            gene = [gene]
        coverage_x, coverage_y = _local_coverage_func(
            reader, gene, **kwargs)
        biglist.append(coverage_y)
    return biglist


def remove_duplicates(fileName):
    """
    Remove duplicate gene entries from a GTF file, overwriting it in-place.

    Reads the file line-by-line, keeps the first occurrence of each
    unique `gene_id` (Ensembl version suffix stripped), and writes the
    result back to the same path. Lines whose `gene_id` cannot be
    parsed are silently dropped.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `fileName` | `str` | — | Path to the GTF file. **The file is overwritten.** |

    ## Returns

    `None`

    ## Notes

    Deduplication key is the `gene_id` attribute after stripping the
    Ensembl version suffix (e.g. `ENSG000001.3` → `ENSG000001`).

    ## Examples

    ```python
    from Pyntegrate.array_helpers import remove_duplicates

    # Deduplicate a GTF before building a gffutils database
    remove_duplicates("gencode.v44.annotation.gtf")
    ```
    """
    seen_gene_ids = set()
    output_lines = []

    with open(fileName, 'r') as f_in:
        for line in f_in:
            fields = line.strip().split('\t')

            attributes = fields[-1]
            gene_id = None

            if 'gene_id "' in attributes:
                try:
                    gene_id = attributes.split('gene_id "')[1].split('"; ')[0]
                    gene_id = gene_id.split('.')[0]
                except IndexError:
                    gene_id = None

            if gene_id:
                if gene_id not in seen_gene_ids:
                    seen_gene_ids.add(gene_id)
                    output_lines.append(line)

    with open(fileName, 'w') as f_out:
        f_out.writelines(output_lines)

        