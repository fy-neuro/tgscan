"""Test parsers module."""
import os
import pytest
from tgscan.parsers import detect_format, parse_matrix, numeric_sample_cols


def test_detect_format():
    assert detect_format('/tmp/x.xlsx') == 'xlsx'
    assert detect_format('/tmp/x.xlsx.gz') == 'xlsx_gz'
    assert detect_format('/tmp/x.xls.gz') == 'xls_gz'
    assert detect_format('/tmp/x.csv.gz') == 'csv'
    assert detect_format('/tmp/x.tsv') == 'tsv'
    assert detect_format('/tmp/x.txt.gz') == 'txt'
    assert detect_format('/tmp/x.h5ad') == 'h5ad'
    assert detect_format('/tmp/x_RAW.tar') == 'raw_tar'
    assert detect_format('/tmp/x_unknown') == 'unknown'


def test_parse_matrix_tsv(nfil3_matrix):
    df, fmt = parse_matrix(nfil3_matrix)
    assert fmt in ('csv', 'tsv')  # tsv.gz detected as tsv by extension
    assert len(df) > 1000  # background + key genes
    # check key columns present
    assert 'EnsemblID' in df.columns
    assert 'Symbol' in df.columns
    # sample columns
    sample_cols = [c for c in df.columns if c.startswith('Sample_')]
    assert len(sample_cols) == 56


def test_numeric_sample_cols(nfil3_matrix):
    df, _ = parse_matrix(nfil3_matrix)
    sample_cols = numeric_sample_cols(df, df.columns[0])
    # Symbol column should not be in samples
    assert 'Symbol' not in sample_cols
    # EnsemblID column should not be in samples (id_col)
    assert 'EnsemblID' not in sample_cols
    # sample columns should be there
    assert len(sample_cols) == 56


def test_unknown_format_raises():
    with pytest.raises(Exception):
        parse_matrix('/tmp/nonexistent_file.xyz')
