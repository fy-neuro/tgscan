"""tgscan — transgene hitchhiker gene screener.

Public API:
    from tgscan import verify, parse_matrix, cis_test, GtfIndex
"""
from .models import AnalysisResult, CisResult, VerifyResult, GeneLocation
from .gtf import GtfIndex, get_gtf_index
from .parsers import parse_matrix, numeric_sample_cols, detect_format
from .analysis import analyze
from .cis import cis_test
from .runner import verify, verify_batch
from . import catalog

__version__ = "0.2.0"

__all__ = [
    'verify', 'verify_batch', 'parse_matrix', 'numeric_sample_cols', 'detect_format',
    'analyze', 'cis_test',
    'GtfIndex', 'get_gtf_index',
    'AnalysisResult', 'CisResult', 'VerifyResult', 'GeneLocation',
    'catalog', '__version__',
]
