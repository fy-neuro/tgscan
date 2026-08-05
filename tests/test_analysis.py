"""Test analysis + cis modules."""
import pytest
import numpy as np
import pandas as pd
from tgscan.parsers import parse_matrix, numeric_sample_cols
from tgscan.analysis import analyze
from tgscan.cis import cis_test
from tgscan.gtf import get_gtf_index


def test_analyze_basic(nfil3_matrix, mini_gtf):
    """Test analyze() returns sensible result."""
    df, _ = parse_matrix(nfil3_matrix)
    gtf = get_gtf_index(mini_gtf)
    id_col = df.columns[0]
    sample_cols = numeric_sample_cols(df, id_col)

    drv_row, _ = gtf.find_in_matrix(df, 'Nfil3', 'ENSMUSG00000056749')
    cand_row, _ = gtf.find_in_matrix(df, 'Auh', 'ENSMUSG00000021460')
    assert drv_row is not None
    assert cand_row is not None

    result = analyze(df, drv_row, cand_row, sample_cols)
    assert result.n_samples == 56
    assert result.n_genes > 100
    assert -1 <= result.r <= 1
    assert 0 <= result.percentile <= 100


def test_cis_finds_driver_in_top(nfil3_matrix, mini_gtf):
    """cis_test should rank driver #1."""
    df, _ = parse_matrix(nfil3_matrix)
    gtf = get_gtf_index(mini_gtf)
    id_col = df.columns[0]
    sample_cols = numeric_sample_cols(df, id_col)
    result = cis_test(df, 'ENSMUSG00000056749', sample_cols, gtf, id_col=id_col)
    assert result.n_genes_genome_wide > 100
    assert 'Nfil3' in result.top_cis_genes
