"""Regression tests from field bugs found 2026-08-15..17 (Task25 pipeline)."""
import numpy as np
import pandas as pd
import pytest

from tgscan.parsers import _merge_dup_value_cols, numeric_sample_cols
from tgscan.design import detect_design_issues
from tgscan.catalog import stats, construct_gate


def test_numeric_leading_mgi_symbols():
    """NO_VALID_R root cause: numeric-leading MGI symbols (0610005C13Rik) must
    still be usable as the id column (symbol lookup must not depend on ^[A-Z])."""
    df = pd.DataFrame({
        "gene": ["0610005C13Rik", "0610007P14Rik", "Sox17", "Tgdrv1"],
        "s1": [1.0, 200.0, 5.0, 900.0],
        "s2": [2.0, 180.0, 3.0, 800.0],
    })
    cols = numeric_sample_cols(df, "gene")
    assert cols == ["s1", "s2"]


def test_cufflinks_duplicate_fpkm_cols_merged():
    """GSE71463-style Cufflinks files repeat the FPKM header (FPKM, FPKM, FPKM)
    which pandas dedups to FPKM/FPKM.1/FPKM.2 — must collapse to one column."""
    df = pd.DataFrame({
        "tracking_id": ["g1", "g2"],
        "FPKM": [1.0, 2.0],
        "FPKM.1": [3.0, 4.0],
        "FPKM.2": [6.0, 8.0],
    })
    out = _merge_dup_value_cols(df)
    assert "FPKM" in out.columns and "FPKM.1" not in out.columns
    assert out["FPKM"].tolist() == pytest.approx([10.0 / 3, 14.0 / 3])


def test_dup_merge_leaves_sample_cols_alone():
    """Legit sample columns sharing a prefix (tg.1, tg.2) must NOT be merged."""
    df = pd.DataFrame({
        "gene": ["g1", "g2"],
        "tg.1": [1.0, 1.0],
        "tg.2": [2.0, 2.0],
    })
    out = _merge_dup_value_cols(df)
    assert "tg.1" in out.columns and "tg.2" in out.columns


def test_zscore_matrix_flagged():
    """z-score normalized matrices (many negative values) must be flagged."""
    rng = np.random.default_rng(0)
    z = rng.normal(0, 1, (100, 4))
    df = pd.DataFrame(z, columns=[f"s{i}" for i in range(4)])
    df.insert(0, "gene", [f"g{i}" for i in range(100)])
    issues = detect_design_issues(df, [f"s{i}" for i in range(4)])
    assert any(i.startswith("zscore_normalized") for i in issues)


def test_counts_matrix_not_zscore_flagged():
    rng = np.random.default_rng(0)
    counts = rng.poisson(50, (100, 4)).astype(float)
    df = pd.DataFrame(counts, columns=[f"s{i}" for i in range(4)])
    df.insert(0, "gene", [f"g{i}" for i in range(100)])
    issues = detect_design_issues(df, [f"s{i}" for i in range(4)])
    assert not any(i.startswith("zscore_normalized") for i in issues)


def test_catalog_bundled_resource():
    """After pip install, the catalog must load from package resources
    (regression: v0.1.0 pointed outside site-packages and crashed)."""
    s = stats()
    assert s["total"] >= 26
    assert s["confirmed"] >= 21


def test_construct_gate():
    """Empirical blacklist gate: exact hits, literal-asterisk MGI symbols,
    WARN vs EXCLUDE separation."""
    from tgscan.catalog import construct_warning
    assert construct_gate("Tg(Adipoq-cre)1Evdr") is not None          # exact EXCLUDE
    assert construct_gate("Tg(Trpm5-EGFP)#Sdmk") is not None          # '#' inside symbol
    assert construct_gate("Tg(Dkk3-cre)D9Tfur") is None               # not blacklisted
    # '*' is a LITERAL MGI character (mutant marker), never a wildcard:
    assert construct_gate("Tg(Mpz*)CDah") is None                     # Mpz is WARN, not gate
    assert construct_warning("Tg(Mpz*)CDah") is not None              # ...but it warns
    assert construct_warning("Tg(Itgax-cre)1-1Reiz") is not None      # exact WARN
