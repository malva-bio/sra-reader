"""
sra_reader: Direct SRA file reader for single-cell experiments.

Bypasses fasterq-dump by streaming reads directly from .sra files
using vdb-dump, with a Cython-optimized hot loop for line parsing.
"""

from .reader import SRAReader, SRASpot, SRASegment
from .layout import ReadLayout, detect_layout
from .utils import sra_info, check_dependencies

__all__ = [
    "SRAReader",
    "SRASpot",
    "SRASegment",
    "ReadLayout",
    "detect_layout",
    "sra_info",
    "check_dependencies",
]
