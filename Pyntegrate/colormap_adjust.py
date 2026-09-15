"""
Custom colormap utilities for genomic signal visualisation.

Provides helpers for creating diverging colormaps centred on zero,
power-law adjustments, and centre-point shifting. The main utility is
[`smart_colormap`][Pyntegrate.colormap_adjust.smart_colormap], which
automatically matches colour saturation between the positive and
negative halves of asymmetric data ranges.

`cmap_powerlaw_adjust`, `cmap_center_adjust`, and
`cmap_center_point_adjust` are adapted from
<https://sites.google.com/site/theodoregoetz/notes/matplotlib_colormapadjust>.
`color_test` is a small debugging helper that renders a solid-colour
figure for a given colour specifier.
"""

import math
import copy
import numpy
import numpy as np
from matplotlib import pyplot, colors, cm
import matplotlib
import colorsys


def color_test(color):
    """
    Display a solid-colour figure for quick visual debugging.

    Opens a 2 × 2 inch matplotlib figure filled uniformly with `color`.
    Useful when exploring hex codes, RGB tuples, or NumPy colour arrays.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `color` | `str \| tuple \| np.ndarray` | — | Any matplotlib-compatible colour specifier. NumPy arrays are ravelled before use. |

    ## Examples

    ```python
    import matplotlib.pyplot as plt
    from Pyntegrate.colormap_adjust import color_test

    color_test("#1a73e8")   # solid Google-blue figure
    plt.show()
    ```
    """
    if isinstance(color, np.ndarray):
        color = color.ravel()
    fig = pyplot.figure(figsize=(2, 2), facecolor=color)


def smart_colormap(vmin, vmax, color_high='#b11902', hue_low=0.6):
    """
    Build a diverging colormap centred on zero for asymmetric data.

    Creates a white-centred, three-stop colormap where the saturation
    and value of the "low" colour mirror those of `color_high` at the
    symmetric counterpart, then the hue is shifted to `hue_low`. This
    ensures that colours at equal distance from zero appear equally
    vivid even when `|vmin| ≠ |vmax|`.

    If both extremes are effectively zero the function returns the
    string `"Greys"` as a safe fallback.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `vmin` | `float` | — | Minimum value in the data to be plotted. |
    | `vmax` | `float` | — | Maximum value in the data to be plotted. |
    | `color_high` | `str` | `"#b11902"` | Matplotlib colour for the high (positive) end of the scale. |
    | `hue_low` | `float` | `0.6` | HSV hue `[0, 1]` for the low (negative) end. `0.6` gives a muted blue. |

    ## Returns

    `matplotlib.colors.LinearSegmentedColormap | str` — A 2048-colour
    diverging colormap, or `"Greys"` when both extremes are zero.

    ## Examples

    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from Pyntegrate.colormap_adjust import smart_colormap

    cmap = smart_colormap(vmin=-3, vmax=9)
    data = np.random.uniform(-3, 9, (50, 50))

    fig, ax = plt.subplots()
    im = ax.imshow(data, cmap=cmap, aspect="auto", vmin=-3, vmax=9)
    fig.colorbar(im, ax=ax)
    plt.title("Asymmetric diverging colormap")
    plt.show()
    ```
    """
    orig_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'test', ['#FFFFFF', color_high], N=2048)

    vmin = float(vmin)
    vmax = float(vmax)
    mx = max([vmin, vmax])
    mn = min([vmin, vmax])
    if mx == 0 or mn == 0 or np.isclose(mx, 0):
        return "Greys"

    frac = abs(mn / mx)

    rgb = orig_cmap(frac)[:-1]

    hsv = list(colorsys.rgb_to_hsv(*rgb))
    hsv[0] = hue_low
    new_rgb = colorsys.hsv_to_rgb(*hsv)

    zeropoint = abs(-vmin / (vmax - vmin))

    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'test', [(0, new_rgb), (zeropoint, '#FFFFFF'), (1, color_high)],
        N=2048)

    return new_cmap


def cmap_powerlaw_adjust(cmap, a):
    """
    Rescale colormap indices via a power-law transformation.

    Applies `new_index = old_index ** a` to the colour lookup table,
    compressing or expanding whichever end of the colormap has the most
    visual resolution. Values `a < 0` return the input unchanged.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `cmap` | `matplotlib.colors.Colormap` | — | Input colormap instance (e.g. `cm.jet`). |
    | `a` | `float` | — | Power-law exponent. `a > 1` compresses the low end; `a < 1` compresses the high end. `a < 0` returns the input unchanged. |

    ## Returns

    `matplotlib.colors.LinearSegmentedColormap` — A new 1024-colour
    colormap with power-law adjusted index positions.

    ## Examples

    ```python
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    import numpy as np
    from Pyntegrate.colormap_adjust import cmap_powerlaw_adjust

    data = np.random.exponential(scale=1.0, size=(40, 40))
    adjusted = cmap_powerlaw_adjust(cm.viridis, a=0.5)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    axes[0].imshow(data, cmap=cm.viridis)
    axes[0].set_title("Original")
    axes[1].imshow(data, cmap=adjusted)
    axes[1].set_title("Power-law adjusted (a=0.5)")
    plt.show()
    ```
    """
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x: (x[0] ** a, x[1], x[2])
    for key in ('red', 'green', 'blue'):
        cdict[key] = list(map(fn, cdict[key]))
        cdict[key].sort()
        assert (cdict[key][0][0] >= 0 and cdict[key][-1][0] <= 1), \
            "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('colormap', cdict, 1024)


def cmap_center_adjust(cmap, center_ratio):
    """
    Shift the perceptual midpoint of a colormap to a new relative position.

    Converts `center_ratio` to a power-law exponent via a log
    transform and delegates to
    [`cmap_powerlaw_adjust`][Pyntegrate.colormap_adjust.cmap_powerlaw_adjust].
    The colour that was previously at index 0.5 will now appear at
    `center_ratio`. Returns the input unchanged when `center_ratio`
    is outside `(0, 1)`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `cmap` | `matplotlib.colors.Colormap` | — | Input colormap instance. |
    | `center_ratio` | `float` | — | New position of the old midpoint colour, in `(0, 1)`. `> 0.5` pulls the centre toward the high end; `< 0.5` toward the low end. |

    ## Returns

    `matplotlib.colors.LinearSegmentedColormap` — Adjusted colormap, or
    the input unchanged if `center_ratio` is not in `(0, 1)`.

    ## Examples

    ```python
    import matplotlib.cm as cm
    from Pyntegrate.colormap_adjust import cmap_center_adjust

    # Pull the midpoint 25 % toward the high end
    adjusted = cmap_center_adjust(cm.seismic, center_ratio=0.75)
    ```
    """
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    return cmap_powerlaw_adjust(cmap, a)


def cmap_center_point_adjust(cmap, range, center):
    """
    Shift a colormap so that a specific data value maps to the midpoint colour.

    Converts the absolute `center` value to a fractional position
    within `range` and calls
    [`cmap_center_adjust`][Pyntegrate.colormap_adjust.cmap_center_adjust].
    Typical use: force zero (or another anchor) to map to white in a
    diverging colormap. Returns the input unchanged when `center` lies
    outside `range`.

    ## Parameters

    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | `cmap` | `matplotlib.colors.Colormap` | — | Input colormap instance. |
    | `range` | `tuple[float, float]` | — | `(min, max)` of the data being plotted. |
    | `center` | `float` | — | Data value that should correspond to the colormap midpoint. Must satisfy `range[0] < center < range[1]`. |

    ## Returns

    `matplotlib.colors.LinearSegmentedColormap` — Adjusted colormap, or
    the input unchanged when `center` is outside `range`.

    ## Examples

    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from Pyntegrate.colormap_adjust import cmap_center_point_adjust

    # Force 0 to be white in a dataset ranging from -2 to +8
    cmap = cmap_center_point_adjust(cm.seismic, range=(-2, 8), center=0)

    data = np.random.uniform(-2, 8, (40, 40))
    plt.imshow(data, cmap=cmap, vmin=-2, vmax=8, aspect="auto")
    plt.colorbar()
    plt.title("Zero centred on white")
    plt.show()
    ```
    """
    if not ((range[0] < center) and (center < range[1])):
        return cmap
    return cmap_center_adjust(
        cmap,
        abs(center - range[0]) / abs(range[1] - range[0]))


if __name__ == '__main__':
    def func3(x, y):
        return (1 - x / 2 + x ** 5 + y ** 3) * numpy.exp(-x ** 2 - y ** 2)
    x = numpy.linspace(-3.0, 3.0, 60)
    y = numpy.linspace(-3.0, 3.0, 60)
    X, Y = numpy.meshgrid(x, y)
    Z = func3(X, Y)
    extent = [x[0], x[-1], y[0], y[-1]]

    plotkwargs = {
        'extent': extent,
        'origin': 'lower',
        'interpolation': 'nearest',
        'aspect': 'auto'}

    fig = pyplot.figure(figsize=(8, 3))
    fig.subplots_adjust(left=.05, bottom=.11, right=.94, top=.83, wspace=.35)
    ax = [fig.add_subplot(1, 3, i) for i in range(1, 4, 1)]

    cmap = cm.seismic

    plt = ax[0].imshow(Z, cmap=cmap, **plotkwargs)
    cb = ax[0].figure.colorbar(plt, ax=ax[0])
    ax[0].set_title('cmap: seismic')

    plt = ax[1].imshow(Z, cmap=cmap_center_adjust(cmap, 0.75), **plotkwargs)
    cb = ax[1].figure.colorbar(plt, ax=ax[1])
    ax[1].set_title('center raised by 25%')

    plt = ax[2].imshow(
        Z,
        cmap=cmap_center_point_adjust(
            cmap, [numpy.min(Z), numpy.max(Z)], 0),
        **plotkwargs)
    cb = ax[2].figure.colorbar(plt, ax=ax[2])
    ax[2].set_title('center set to 0')

    pyplot.show()
