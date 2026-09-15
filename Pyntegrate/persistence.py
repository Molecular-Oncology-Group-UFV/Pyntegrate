"""
Tools for persisting genomic features and signal arrays across sessions.

Provides two complementary functions —
[`save_features_and_arrays`][Pyntegrate.persistence.save_features_and_arrays]
and
[`load_features_and_arrays`][Pyntegrate.persistence.load_features_and_arrays]
— that write and read a pair of files (a BED/GFF feature file and a
NumPy `.npz` archive) so that expensive array computations do not need
to be repeated in every analysis session.
"""

import os
import pybedtools
import numpy as np


def load_features_and_arrays(prefix, mmap_mode='r'):
    """
    Load genomic features and signal arrays saved by
    [`save_features_and_arrays`][Pyntegrate.persistence.save_features_and_arrays].

    Reads two files produced by `save_features_and_arrays`:

    - **`prefix.features`** — BED/GFF features as a `pybedtools.BedTool`.
    - **`prefix.npz`** — NumPy archive containing named signal arrays.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `prefix` | `str` | — | File path prefix used when the data were originally saved. The suffixes `.features` and `.npz` are appended automatically. |
    | `mmap_mode` | `str \| None` | `"r"` | Memory-mapping mode forwarded to `numpy.load`. One of `None`, `"r+"`, `"r"`, `"w+"`, or `"c"`. Use `"r"` (read-only) for safe concurrent access. |

    ## Returns

    `tuple[pybedtools.BedTool, numpy.lib.npyio.NpzFile]` — A 2-tuple
    `(features, arrays)` where `features` is a `pybedtools.BedTool`
    and `arrays` is a NumPy `.npz` archive whose keys are signal names.

    ## Examples

    ```python
    from Pyntegrate.persistence import load_features_and_arrays

    features, arrays = load_features_and_arrays("my_experiment")

    # Access a named signal array
    ip_array = arrays["ip"]
    print(ip_array.shape)
    ```

    **Output:**
    ```
    (5000, 100)
    ```
    """
    features = pybedtools.BedTool(prefix + '.features')
    arrays = np.load(prefix + '.npz', mmap_mode=mmap_mode, allow_pickle=True)
    return features, arrays


def save_features_and_arrays(features, arrays, prefix, compressed=False,
                             link_features=False, overwrite=False):
    """
    Persist genomic features and signal arrays to disk for later reuse.

    Writes two files starting with `prefix`:

    - **`prefix.features`** — The genomic features in BED/GFF format
      (or a symlink when `link_features=True`).
    - **`prefix.npz`** — A NumPy archive containing the named signal
      arrays (compressed if `compressed=True`).

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `features` | `iterable \| pybedtools.BedTool \| str` | — | Genomic features whose rows correspond to rows in each array. Usually the same features used to create the arrays. |
    | `arrays` | `dict[str, np.ndarray]` | — | Dictionary mapping signal names to NumPy arrays. Keys become the array names in the `.npz` file. |
    | `prefix` | `str` | — | Output path prefix. The suffixes `.features` and `.npz` are appended automatically. |
    | `compressed` | `bool` | `False` | If `True`, use `numpy.savez_compressed` instead of `numpy.savez`. Saves disk space at the cost of slower subsequent reads. |
    | `link_features` | `bool` | `False` | If `True`, create a symlink to the original features file instead of copying it. `features` must be a file-backed `pybedtools.BedTool` or a path string. |
    | `overwrite` | `bool` | `False` | If `True` and `link_features=True`, pass `-f` to `ln -s` to overwrite an existing symlink. |

    ## Examples

    ```python
    import numpy as np
    import pybedtools
    from Pyntegrate.persistence import (
        save_features_and_arrays,
        load_features_and_arrays,
    )

    features = pybedtools.BedTool("genes.bed")
    arrays = {
        "ip":    np.random.rand(100, 200),
        "input": np.random.rand(100, 200),
    }

    save_features_and_arrays(features, arrays, prefix="experiment1")

    # Reload in the next session
    feat, arr = load_features_and_arrays("experiment1")
    print(arr["ip"].shape)
    ```

    **Output:**
    ```
    (100, 200)
    ```
    """
    if link_features:
        if isinstance(features, pybedtools.BedTool):
            assert isinstance(features.fn, str)
            features_filename = features.fn
        else:
            assert isinstance(features, str)
            features_filename = features

        if overwrite:
            force_flag = '-f'
        else:
            force_flag = ''

        cmds = [
            'ln', '-s', force_flag, os.path.abspath(features_filename), prefix + '.features']
        os.system(' '.join(cmds))
    else:
        pybedtools.BedTool(features).saveas(prefix + '.features')

    if compressed:
        np.savez_compressed(
            prefix,
            **arrays)
    else:
        np.savez(prefix, **arrays)
