"""Stage 1: Pearson correlation + percentile analysis."""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats
from .models import AnalysisResult


def analyze(df: pd.DataFrame, driver_row: pd.Series, candidate_row: pd.Series,
            sample_cols) -> AnalysisResult:
    """Compute Pearson r between driver and candidate, plus percentile in genome-wide distribution.

    Args:
        df: full gene x sample matrix
        driver_row: row of the driver gene
        candidate_row: row of the candidate gene
        sample_cols: list of sample column names

    Returns:
        AnalysisResult with r, percentile, mean_r, status.
    """
    drv = np.array([float(driver_row[c]) for c in sample_cols], dtype=float)
    cand = np.array([float(candidate_row[c]) for c in sample_cols], dtype=float)

    if np.std(drv) == 0:
        return AnalysisResult(r=float('nan'), p_value=float('nan'), percentile=float('nan'),
                              mean_r=float('nan'), n_samples=len(sample_cols), n_genes=0,
                              status='ZERO_VARIANCE_DRIVER')
    if np.std(cand) == 0:
        return AnalysisResult(r=float('nan'), p_value=float('nan'), percentile=float('nan'),
                              mean_r=float('nan'), n_samples=len(sample_cols), n_genes=0,
                              status='ZERO_VARIANCE_CANDIDATE')

    r_val, p_val = stats.pearsonr(drv, cand)

    # genome-wide r distribution
    all_r = []
    for _, row in df.iterrows():
        vals = np.array([float(row[c]) for c in sample_cols], dtype=float)
        if np.std(vals) == 0 or np.nansum(vals) < 10:
            continue
        try:
            rr, _ = stats.pearsonr(drv, vals)
            if not np.isnan(rr):
                all_r.append(rr)
        except Exception:
            continue
    arr = np.array(all_r)
    if len(arr) == 0:
        return AnalysisResult(r=float(r_val), p_value=float(p_val), percentile=float('nan'),
                              mean_r=float('nan'), n_samples=len(sample_cols), n_genes=0,
                              status='NO_BG')

    pct = float((arr < r_val).mean() * 100)
    mean_r = float(arr.mean())

    if mean_r > 0.3:
        verdict = 'BACKGROUND_TOO_HIGH'
    elif pct > 95:
        verdict = 'HIGH_CONFIDENCE'
    elif pct > 90:
        verdict = 'MODERATE'
    else:
        verdict = 'NO_SIGNAL'

    return AnalysisResult(
        r=float(r_val), p_value=float(p_val), percentile=pct,
        mean_r=mean_r, n_samples=len(sample_cols), n_genes=len(arr),
        status=verdict,
    )
