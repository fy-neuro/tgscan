"""Known hitchhiker catalog — bundled with the package."""
from __future__ import annotations
import os
from functools import lru_cache
from typing import Optional
import pandas as pd


CATALOG_TSV = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
                            'data', 'known_hitchhikers.tsv')


@lru_cache(maxsize=1)
def _load_catalog() -> pd.DataFrame:
    if not os.path.exists(CATALOG_TSV):
        raise FileNotFoundError(f"Catalog TSV not found: {CATALOG_TSV}")
    return pd.read_csv(CATALOG_TSV, sep='\t', keep_default_na=False, na_values=[''])


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
