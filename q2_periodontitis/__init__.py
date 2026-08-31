"""
q2-periodontitis
================
QIIME 2 plugin for subgingival biofilm ecological characterization,
Socransky complexes collapse, and dysbiosis indices.
"""

from ._version import __version__
from ._methods import collapse_complexes, delta_indices
from ._visualizer import summarize

__all__ = [
    "__version__",
    "collapse_complexes",
    "delta_indices",
    "summarize",
]
