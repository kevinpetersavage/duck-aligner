"""
Duck Aligner package.
"""

from duck_aligner.db import connect, register
from duck_aligner.aligner import align, compute_align_score
from duck_aligner.fastq import parse_fastq_records

__all__ = [
    "connect",
    "register",
    "align",
    "compute_align_score",
    "parse_fastq_records",
]
