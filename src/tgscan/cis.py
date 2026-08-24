"""Stage 2: cis-enrichment hypergeometric test.

Real hitchhikers must be physically within the BAC construct (cis to driver).
If top-correlated genes are scattered across the genome (trans), the high r is
likely a cell-type composition artifact, not a hitchhiker.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom
from typing import Optional
from .models import CisResult
from .gtf import GtfIndex


def _find_eid_in_row(row, id_col, sym_to_eid):
    """Extract Ensembl ID from a row by checking id_col then Symbol column."""
    v = str(row[id_col])
    if v.startswith('ENSMUSG') or v.startswith('ENSG'):
        return v.split('.')[0]
    # try Symbol-like → Ensembl via GtfIndex
    eid = sym_to_eid(v) if callable(sym_to_eid) else None
    return eid


def cis_test(df: pd.DataFrame, driver_eid: str, sample_cols, gtf: GtfIndex,
             cis_window_kb: int = 1000, top_k=(10, 50, 100),
             id_col: Optional[str] = None) -> CisResult:
    """Run cis-enrichment test for a driver gene.

    Args:
        df: gene x sample matrix
        driver_eid: Ensembl ID of the driver gene (must be in GTF)
        sample_cols: list of sample column names
        gtf: GtfIndex for gene locations
        cis_window_kb: window around driver to consider as cis (default ±1 Mb)
        top_k: tuple of Top-K values to test (default 10, 50, 100)
        id_col: name of the gene ID column in df (auto-detected if None)

    Returns:
        CisResult with hypergeometric p-values for each Top-K and a verdict.
    """
    drv_loc = gtf.get_location(driver_eid)
    if drv_loc is None:
        raise ValueError(f"Driver {driver_eid} not in GTF")

    if id_col is None:
        id_col = df.columns[0]

    # find driver row in matrix (to compute genome-wide r)
    drv_row, _ = gtf.find_in_matrix(df, drv_loc.symbol, driver_eid)
    if drv_row is None:
        raise ValueError(f"Driver {drv_loc.symbol} ({driver_eid}) not in matrix")
    drv_vals = np.array([float(drv_row[c]) for c in sample_cols], dtype=float)
    if np.std(drv_vals) == 0:
        raise ValueError("Driver has zero variance")

    # identify which columns contain Ensembl IDs vs symbols
    sample_first = df[id_col].astype(str).head(20)
    eid_col = sym_col = None
    if sample_first.str.startswith(('ENSMUSG', 'ENSG')).any():
        eid_col = id_col
    else:
        sym_col = id_col
    for c in df.columns:
        if c == id_col:
            continue
        cs = str(c).lower()
        if cs in ('ensembl', 'ensembl_id', 'gene_id', 'ensembl_gene_id'):
            eid_col = c
        elif cs in ('symbol', 'gene_symbol', 'gene_name', 'genename', 'genes', 'gene'):
            sym_col = c

    # compute r genome-wide + locate each gene
    res = []
    for _, row in df.iterrows():
        # get eid
        eid = None
        if eid_col is not None:
            v = str(row[eid_col])
            if v.startswith(('ENSMUSG', 'ENSG')):
                eid = v.split('.')[0]
        if eid is None and sym_col is not None:
            eid = gtf.symbol_to_ensembl(str(row[sym_col]))
        if eid is None:
            # fallback: scan all columns
            for col in df.columns:
                if col == id_col:
                    continue
                v = str(row[col])
                if v.startswith(('ENSMUSG', 'ENSG')):
                    eid = v.split('.')[0]
                    break
        if eid is None:
            continue
        loc = gtf.get_location(eid)
        if loc is None:
            continue
        vals = np.array([float(row[c]) for c in sample_cols], dtype=float)
        if np.std(vals) == 0 or vals.sum() < 10:
            continue
        try:
            r_val, _ = stats.pearsonr(drv_vals, vals)
        except Exception:
            continue
        if np.isnan(r_val):
            continue
        dist = ((loc.start + loc.end) // 2 - (drv_loc.start + drv_loc.end) // 2) / 1000 \
            if loc.chrom == drv_loc.chrom else None
        res.append({'sym': loc.symbol, 'eid': eid, 'chrom': loc.chrom,
                    'r': r_val, 'dist_kb': dist})

    if not res:
        raise ValueError("No genes with valid r — check that matrix has Ensembl IDs or symbols matching the GTF")

    rdf = pd.DataFrame(res).sort_values('r', ascending=False).reset_index(drop=True)
    rdf['rank'] = rdf.index + 1

    # cis window
    near = rdf[(rdf['chrom'] == drv_loc.chrom) & (rdf['dist_kb'].abs() <= cis_window_kb)]
    N = len(rdf)
    K = len(near)

    # hypergeom for each Top-K
    cis_results = {}
    for k in top_k:
        n_cis = int(((rdf['rank'] <= k) & (rdf['chrom'] == drv_loc.chrom) &
                     (rdf['dist_kb'].abs() <= cis_window_kb)).sum())
        expected = k * K / N if N > 0 else 0.0
        pval = float(1 - hypergeom.cdf(n_cis - 1, N, K, k)) if n_cis > 0 else 1.0
        fold = n_cis / expected if expected > 0 else float('inf')
        cis_results[k] = {'n_cis': n_cis, 'expected': expected,
                          'fold': fold, 'p': pval}

    # best p
    best_p = min(cis_results[k]['p'] for k in top_k)

    # verdict
    if best_p < 1e-3:
        verdict = 'CONFIRMED'
    elif best_p < 1e-2:
        verdict = 'CANDIDATE'
    else:
        verdict = 'WEAK_OR_FALSE_POSITIVE'

    # B-side (08-24): exclude the driver's self row (r=1.0, rank 1, dist 0 —
    # a guaranteed cis slot any gene would get) and re-rank. Disclosed alongside
    # the A-side per PI 2026-08-24 option-2 ruling; verdict thresholds unchanged.
    rdf_b = rdf[rdf['eid'] != driver_eid].reset_index(drop=True)
    rdf_b['rank'] = rdf_b.index + 1
    near_b = rdf_b[(rdf_b['chrom'] == drv_loc.chrom) & (rdf_b['dist_kb'].abs() <= cis_window_kb)]
    n_b, k_b = len(rdf_b), len(near_b)
    best_p_b = 1.0
    if n_b > 0:
        for k in top_k:
            n_cis_b = int(((rdf_b['rank'] <= k) & (rdf_b['chrom'] == drv_loc.chrom) &
                           (rdf_b['dist_kb'].abs() <= cis_window_kb)).sum())
            p_b = float(1 - hypergeom.cdf(n_cis_b - 1, n_b, k_b, k)) if n_cis_b > 0 else 1.0
            best_p_b = min(best_p_b, p_b)
    verdict_b = ('CONFIRMED' if best_p_b < 1e-3 else
                 'CANDIDATE' if best_p_b < 1e-2 else 'WEAK_OR_FALSE_POSITIVE')

    # top cis genes (limit 5)
    top_cis = near.sort_values('r', ascending=False).head(5)['sym'].astype(str).tolist()

    return CisResult(
        cis_top10_p=cis_results.get(10, {}).get('p', float('nan')),
        cis_top50_p=cis_results.get(50, {}).get('p', float('nan')),
        cis_top100_p=cis_results.get(100, {}).get('p', float('nan')),
        best_p=best_p,
        verdict=verdict,
        top_cis_genes=top_cis,
        n_cis_1mb=K,
        n_genes_genome_wide=N,
        fold_top10=cis_results.get(10, {}).get('fold', 0.0),
        fold_top50=cis_results.get(50, {}).get('fold', 0.0),
        fold_top100=cis_results.get(100, {}).get('fold', 0.0),
        best_p_excl_driver=best_p_b,
        verdict_excl_driver=verdict_b,
    )
