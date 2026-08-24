"""Design gate — screen unusable GEO designs BEFORE correlation analysis.

v0.3 hard requirement (ttyh3_cerkl_review.md 2026-08-20, PI-endorsed):
FACS-sorted / single-cell / TRAP-IP designs manufacture spurious correlations
(sign flips, fake high r). The old guard only checked background dirtiness;
design type is now a first-class gate.

Two complementary layers:
  A. matrix-intrinsic detection (this module; always-on at verify time)
  B. curated registry data/known_designs.tsv (GSE -> verdict, from GEO
     metadata review; authoritative for known datasets)

Blocking designs:
  facs_sorted       sorted fractions (GFP+/GFP-) as sample groups; r measures
                    sorting structure, not co-regulation. Ground truth:
                    GSE83356 cross-group r=-0.77 while within-GFP+ r=+0.91.
  single_cell_matrix cells-as-samples; r measures cell-type composition
                    (GSE115934: P1/P12/P100 hair cells as "samples").
  trap_rna_ip       TRAP/RiboTag/IP-enriched RNA is not total RNA
                    (DECISIONS 2026-07-25; GSE127845 was also an allele
                    mismatch: Fmr1 KO, no Neurod1-cre at all).
  zscore_normalized per-row z-scaled export; Pearson dominated by scaling.
  mirna_only / filter_too_strict / too_few_samples  (pre-existing checks).
"""
from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd

DATA_DIR = Path(__file__).parent / "data"
KNOWN_DESIGNS = DATA_DIR / "known_designs.tsv"

# issues that make correlation evidence unusable -> BLOCK
BLOCKING_ISSUES = {
    'mirna_only', 'facs_sorted', 'single_cell_matrix', 'trap_rna_ip',
    'zscore_normalized', 'too_few_samples', 'filter_too_strict',
}

# ---------------------------------------------------------------- layer A

# fluorescence markers used as sort reporters
_FLUOR = r'(?:GFP|EGFP|YFP|RFP|tdTomato|mCherry|Venus|DsRed|Tomato|ZsGreen)'
# Fraction markers: "GFP+_1" / "GFP-2" / "EGFP pos" / "GFP-negative" /
# "Venus_plus_22" / "Venus_minus_1". The '-'/'+' attach to the marker
# (optionally via space/underscore) so plain line ids don't misread.
# Pair rule: facs_sorted fires only when BOTH a positive and a negative
# fraction are present (GSE83356 structure). A SINGLE-arm "Venus_plus" may
# just be genotyping (GSE211929: K_WT_Venus_plus, no minus arm, cis passed)
# — single-arm labels do not create cross-fraction sign-flip structure and
# are NOT blocked (flagged for metadata review instead).
_SEP = r'[ _\-]*'
_POS_PAT = re.compile(rf'{_FLUOR}{_SEP}(?:\+|pos(itive)?|plus)(?![a-z])', re.I)
_NEG_PAT = re.compile(rf'{_FLUOR}{_SEP}(?:-(?![a-z])|neg(ative)?|minus)(?![a-z])', re.I)
_FACS_PAT = re.compile(r'\b(?:facs|flow[- ]?sort|sorted)\b', re.I)
# Sort-fraction bin naming: "E18/383_SF10-1", "P5/393_SF17-16" — numbered bins
# of a fluorescence sort (GSE90860, Zbtb16 case 08-24). Require >=3 distinct
# bins so a stray SF token cannot fire the gate.
# TRAP / ribosome IP vocabulary (not total RNA)
_TRAP_PAT = re.compile(r'(?<![A-Za-z0-9])(TRAP|RiboTag|Ribo-?tag|L10a|IgG|IP|input)(?![A-Za-z0-9])', re.I)
_SC_PAT = re.compile(r'(?:^|[^a-z])(?:scRNA|single[- ]?cell|10x)(?:$|[^a-z])', re.I)
_SFBIN_PAT = re.compile(r'(?:^|[/_ \-])(SF\d+)(?:[-_ ]\d+)?(?=$|[/_ \-|])', re.I)


def _detect_facs(sample_cols) -> bool:
    """Sorted-fraction design: paired +/- marker groups, or explicit FACS."""
    names = [str(c) for c in sample_cols]
    joined = ' | '.join(names)
    has_pos = any(_POS_PAT.search(n) for n in names)
    has_neg = any(_NEG_PAT.search(n) for n in names)
    if has_pos and has_neg:
        return True
    if _SFBIN_PAT.search(joined):
        bins = set(_SFBIN_PAT.findall(joined))
        if len(bins) >= 3:
            return True
    return bool(_FACS_PAT.search(joined))


def _detect_trap(sample_cols) -> bool:
    names = [str(c) for c in sample_cols]
    return any(_TRAP_PAT.search(n) for n in names)


def _detect_single_cell(df: pd.DataFrame, sample_cols) -> bool:
    """Cells-as-samples signature: very sparse matrix with many columns.

    Bulk count matrices rarely exceed ~60% zeros genome-wide, while a
    cell-as-sample matrix (or un-pseudobulked 10x export) has >80% zeros.
    Large bulk cohorts also almost never ship unfiltered; n>=24 columns
    plus >80% zeros is a reliable single-cell signature.
    """
    names = [str(c) for c in sample_cols]
    if any(_SC_PAT.search(n) for n in names):
        return True
    if len(sample_cols) < 24:
        return False
    try:
        probe_cols = list(sample_cols[:12])
        sub = df[probe_cols].head(3000).apply(pd.to_numeric, errors='coerce')
        zero_frac = float((sub == 0).sum().sum() / sub.size)
        return zero_frac > 0.80
    except Exception:
        return False


def detect_design_issues(df: pd.DataFrame, sample_cols) -> List[str]:
    """Detect design problems that make a dataset unusable.

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
    # z-score normalized matrix: raw counts are never negative; >30% negative
    # values across the first sample columns implies per-row z-scaling (scRNA
    # exports) — Pearson on such matrices is dominated by normalization
    if sample_cols:
        try:
            probe = df[sample_cols[:3]].apply(pd.to_numeric, errors='coerce')
            neg_frac = float((probe < 0).sum().sum() / probe.size)
            if neg_frac > 0.30:
                issues.append(f'zscore_normalized({neg_frac:.0%}_negative)')
        except Exception:
            pass
    # v0.3 design gate
    if sample_cols:
        if _detect_facs(sample_cols):
            issues.append('facs_sorted')
        if _detect_trap(sample_cols):
            issues.append('trap_rna_ip')
        if _detect_single_cell(df, sample_cols):
            issues.append('single_cell_matrix')
    return issues


# ---------------------------------------------------------------- layer B

@lru_cache(maxsize=1)
def _load_registry() -> dict:
    """GSE -> {design, verdict, evidence, source}; empty dict if file missing."""
    reg = {}
    if KNOWN_DESIGNS.exists():
        with open(KNOWN_DESIGNS) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                gse = (r.get('gse') or '').strip()
                if gse:
                    reg[gse] = r
    return reg


def lookup_registry(gse: Optional[str]) -> Optional[dict]:
    """Curated design record for a GSE (None if not audited)."""
    if not gse:
        return None
    return _load_registry().get(gse.strip().upper())


def registry_verdict(gse: Optional[str]) -> Optional[str]:
    """'BLOCKED' | 'WEAK' | 'OK' | None (not audited)."""
    rec = lookup_registry(gse)
    return (rec.get('verdict') or '').strip().upper() if rec else None


# ---------------------------------------------------------------- gate

def screen_design(df: pd.DataFrame, sample_cols, gse: Optional[str] = None) -> dict:
    """Combine matrix-intrinsic detection with the curated registry.

    Returns dict(verdict='BLOCK'|'WEAK'|'PASS',
                 issues=[...], registry=<record or None>, source='matrix'|'registry')
    Registry 'OK' overrides matrix false-positives (escape hatch); registry
    BLOCKED/WEAK wins over matrix PASS (GEO metadata sees what the matrix
    cannot).
    """
    issues = detect_design_issues(df, sample_cols)
    reg = lookup_registry(gse)
    reg_verdict = (reg.get('verdict') or '').strip().upper() if reg else None

    if reg_verdict == 'OK':
        issues = [i for i in issues
                  if not i.startswith(('facs_sorted', 'trap_rna_ip', 'single_cell_matrix'))]
    if reg_verdict == 'BLOCKED':
        verdict = 'BLOCK'
    elif any(i.split('(')[0] in BLOCKING_ISSUES for i in issues):
        verdict = 'BLOCK'
    elif reg_verdict == 'WEAK':
        verdict = 'WEAK'
    else:
        verdict = 'PASS'
    return {'verdict': verdict, 'issues': issues,
            'registry': reg,
            'source': 'registry' if reg_verdict in ('BLOCKED', 'WEAK') else 'matrix'}


def design_issue_is_blocking(issue: str) -> bool:
    """For evidence-store rows: does this design_issues string block usability?

    Store cells may carry e.g. 'facs_sorted' or 'single_cell; n2' — any
    blocking keyword inside blocks.
    """
    return any(k in str(issue) for k in (
        'facs', 'single_cell', 'trap', 'zscore', 'mirna', 'blocked',
        'too_few', 'filter_too_strict'))
