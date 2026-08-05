"""GTF loading and gene location lookup."""
from __future__ import annotations
import gzip, re
from functools import lru_cache
from typing import Optional, Dict, Tuple
from .models import GeneLocation


class GtfIndex:
    """GTF gene index with symbol/Ensembl lookup."""

    def __init__(self, gtf_path: str):
        self.gtf_path = gtf_path
        self._sym_to_eid: Dict[str, str] = {}
        self._gene_locs: Dict[str, GeneLocation] = {}
        self._loaded = False

    def _load(self):
        if self._loaded:
            return
        opener = gzip.open if self.gtf_path.endswith('.gz') else open
        with opener(self.gtf_path, 'rt') as f:
            for line in f:
                if '\tgene\t' not in line:
                    continue
                parts = line.split('\t')
                try:
                    start, end = int(parts[3]), int(parts[4])
                except (ValueError, IndexError):
                    continue
                eid = sym = biotype = None
                for kv in parts[-1].split(';'):
                    kv = kv.strip()
                    if kv.startswith('gene_id'):
                        m = re.search(r'"(ENSMUSG\S+|ENSG\S+)"', kv)
                        eid = m.group(1).split('.')[0] if m else None
                    elif kv.startswith('gene_name'):
                        m = re.search(r'"([^"]+)"', kv)
                        sym = m.group(1) if m else None
                    elif kv.startswith('gene_biotype'):
                        m = re.search(r'"([^"]+)"', kv)
                        biotype = m.group(1) if m else None
                if not (eid and sym):
                    continue
                self._sym_to_eid.setdefault(sym, eid)
                self._gene_locs[eid] = GeneLocation(
                    chrom=parts[0], start=start, end=end, symbol=sym,
                    ensembl_id=eid, biotype=biotype or ''
                )
        self._loaded = True
        return

    def symbol_to_ensembl(self, symbol: str) -> Optional[str]:
        self._load()
        return self._sym_to_eid.get(symbol)

    def get_location(self, ensembl_id: str) -> Optional[GeneLocation]:
        self._load()
        eid = ensembl_id.split('.')[0]
        return self._gene_locs.get(eid)

    def find_in_matrix(self, df, name: str, ensembl_id: Optional[str] = None,
                        entrez_id: Optional[str] = None):
        """Find gene row in a matrix DataFrame.

        Tries: (1) exact symbol match in any column,
               (2) Ensembl ID prefix match,
               (3) Entrez ID exact match.

        Returns (row, col_name) or (None, None).
        """
        # 1. exact symbol match
        for col in df.columns:
            s = df[col].astype(str).str.strip()
            mask = (s == name)
            if mask.any():
                return df[mask].iloc[0], col
        # 2. Ensembl prefix
        if ensembl_id:
            prefix = ensembl_id.split('.')[0][:18]
            for col in df.columns:
                s = df[col].astype(str).str.strip()
                mask = s.str.startswith(prefix)
                if mask.any():
                    return df[mask].iloc[0], col
        # 3. Entrez
        if entrez_id:
            for col in df.columns:
                s = df[col].astype(str).str.strip()
                mask = (s == str(entrez_id))
                if mask.any():
                    return df[mask].iloc[0], col
        return None, None


# Module-level cache so users can re-use the same GTF across calls
@lru_cache(maxsize=4)
def get_gtf_index(gtf_path: str) -> GtfIndex:
    return GtfIndex(gtf_path)
