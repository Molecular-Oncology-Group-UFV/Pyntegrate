"""
Pyntegrate — integrate and visualize ChIP-seq, ATAC-seq, and RNA-seq genomic signal.

Provides a uniform interface for reading genomic signal from BAM,
BigWig, BigBed, BED, and NarrowPeak files (`genomic_signal`), tabular
differential-expression results (`results_table`), and a set of
high-level `pipelines` classes that combine both into ready-to-display
`matplotlib` figures (heatmaps, profiles, volcano plots, deciles, and
more) for ChIP-seq/ATAC-seq/RNA-seq integration analyses.

See `Pyntegrate.pipelines` for the highest-level entry points.
"""
import os
import sys
import time
from .helpers import *
from .helpers import data_dir, example_filename
from ._genomic_signal import genomic_signal
from .plotutils import *
from .colormap_adjust import *
from .results_table import *
from .version import __version__
from .persistence import *
from .SeqSignalAnalysis import *
from . import pipelines
