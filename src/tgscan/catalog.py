"""Known hitchhiker catalog — bundled inside the package (importlib.resources)."""
from __future__ import annotations
from functools import lru_cache
from typing import Optional
from importlib import resources
import pandas as pd


@lru_cache(maxsize=1)
def _load_catalog() -> pd.DataFrame:
    ref = resources.files('tgscan') / 'data' / 'known_hitchhikers.tsv'
    with resources.as_file(ref) as path:
        return pd.read_csv(path, sep='\t', keep_default_na=False, na_values=[''])


@lru_cache(maxsize=1)
def _load_excluded() -> pd.DataFrame:
    """Empirically excluded constructs (promoter cassettes, single-copy, etc.).
    Sources: PI decisions in DECISIONS.md + Task25 audits. Format: allele, action, reason.
    NOTE: '#' is legal inside MGI allele symbols (Tg(...)#Litt), so comment
    lines are filtered manually instead of via read_csv(comment='#')."""
    ref = resources.files('tgscan') / 'data' / 'excluded_constructs.tsv'
    with resources.as_file(ref) as path:
        lines = [ln for ln in path.read_text().splitlines()
                 if ln.strip() and not ln.lstrip().startswith('#')]
        from io import StringIO
        return pd.read_csv(StringIO('\n'.join(lines)), sep='\t', header=None,
                           names=['allele', 'action', 'reason'],
                           keep_default_na=False, na_values=[])


def _blacklist_hit(allele: str, action_want: str) -> Optional[str]:
    """Exact allele-symbol match only. NOTE: '*' is a literal character in MGI
    nomenclature (e.g. Tg(Prnp*)CDah marks a mutant transgene) — it is NOT a
    wildcard here. That collision bit us once; don't reintroduce it."""
    df = _load_excluded()
    hit = df[(df['allele'] == allele) & (df['action'] == action_want)]
    if len(hit):
        return str(hit['reason'].iloc[0])
    return None


def construct_gate(allele: str) -> Optional[str]:
    """Return exclusion reason if allele is EXCLUDE-blacklisted (exact or
    wildcard like 'Tg(Mpz*)'), else None. WARN entries do not gate; use
    construct_warning() for those."""
    return _blacklist_hit(allele, 'EXCLUDE')


def construct_warning(allele: str) -> Optional[str]:
    """Return WARN reason if allele carries a verify-before-use warning."""
    return _blacklist_hit(allele, 'WARN')


def list_confirmed() -> pd.DataFrame:
    """Return all confirmed hitchhikers."""
    df = _load_catalog()
    return df[df['status'] == 'confirmed'].reset_index(drop=True)


def list_candidates() -> pd.DataFrame:
    """Return all candidate hitchhikers (awaiting independent validation)."""
    df = _load_catalog()
    return df[df['status'] == 'candidate'].reset_index(drop=True)


def list_all() -> pd.DataFrame:
    """Return the full catalog."""
    return _load_catalog()


def search(driver: Optional[str] = None, gene: Optional[str] = None) -> pd.DataFrame:
    """Search catalog by driver gene and/or hitchhiker gene."""
    df = _load_catalog()
    mask = pd.Series([True] * len(df))
    if driver:
        mask &= df['driver'].str.lower() == driver.lower()
    if gene:
        mask &= df['gene'].str.lower() == gene.lower()
    return df[mask].reset_index(drop=True)


def stats() -> dict:
    """Summary statistics of the catalog."""
    df = _load_catalog()
    return {
        'total': len(df),
        'confirmed': (df['status'] == 'confirmed').sum(),
        'candidate': (df['status'] == 'candidate').sum(),
        'unique_drivers': df['driver'].nunique(),
        'unique_hitchhikers': df['gene'].nunique(),
    }
