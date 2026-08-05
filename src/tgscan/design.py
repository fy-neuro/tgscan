"""Design issue pre-screening — flag known unusable GEO designs."""
from __future__ import annotations
from typing import List
import pandas as pd


def detect_design_issues(df: pd.DataFrame, sample_cols) -> List[str]:
    """Detect known design problems that make a dataset unusable.

    Args:
        df: parsed matrix
        sample_cols: sample column names

    Returns:
        List of issue strings (empty if no issues).
    """
    issues = []
    if len(df.columns) == 0:
        return ['empty_matrix']
    id_col = df.columns[0]
    # miRNA-only matrix (most IDs start with 'mmu-')
    sample_ids = df[id_col].astype(str).head(100)
    mirna_frac = sample_ids.str.startswith('mmu-').mean()
    if mirna_frac > 0.5:
        issues.append('mirna_only')
    # too few genes (likely strict filtering, ncRNA/pseudogenes excluded)
    if len(df) < 15000:
        issues.append(f'filter_too_strict({len(df)})')
    # too few samples for r
    if len(sample_cols) < 3:
        issues.append(f'too_few_samples({len(sample_cols)})')
    return issues
