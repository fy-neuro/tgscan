"""Dataclass result types for tgscan."""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import Optional, List


@dataclass
class AnalysisResult:
    """Stage 1: Pearson correlation + percentile analysis."""
    r: float
    p_value: float
    percentile: float
    mean_r: float  # genome-wide background mean
    n_samples: int
    n_genes: int  # number of genes used in genome-wide distribution
    status: str  # HIGH_CONFIDENCE / MODERATE / NO_SIGNAL / BACKGROUND_TOO_HIGH
    bg_sd: float = float('nan')   # background sd (design C)
    z_abs: float = float('nan')   # absolute background z-score (design C)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class CisResult:
    """Stage 2: cis-enrichment hypergeometric test."""
    cis_top10_p: float
    cis_top50_p: float
    cis_top100_p: float
    best_p: float
    verdict: str  # CONFIRMED / CANDIDATE / WEAK_OR_FALSE_POSITIVE
    top_cis_genes: List[str] = field(default_factory=list)
    n_cis_1mb: int = 0
    n_genes_genome_wide: int = 0
    fold_top10: float = 0.0
    fold_top50: float = 0.0
    fold_top100: float = 0.0
    # B-side convention (08-24, PI option 2 disclosure): driver self-correlation
    # (r=1.0, rank 1) excluded from ranking — guards the guaranteed cis slot.
    best_p_excl_driver: float = float('nan')
    verdict_excl_driver: str = ''

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class VerifyResult:
    """Combined result of verify() — Stage 1 + Stage 2 (if available)."""
    gene: str
    driver: str
    matrix_format: str
    analysis: Optional[AnalysisResult] = None
    cis: Optional[CisResult] = None
    error: Optional[str] = None
    design_issues: List[str] = field(default_factory=list)
    design_verdict: Optional[str] = None  # BLOCK / WEAK / PASS (design gate)

    @property
    def verdict(self) -> str:
        """Final verdict — design gate, then cis, then stage1."""
        if self.design_verdict == 'BLOCK':
            return "BLOCKED_DESIGN"
        if self.cis is not None:
            return self.cis.verdict
        if self.analysis is not None:
            return self.analysis.status
        return "ERROR"

    @property
    def r(self) -> Optional[float]:
        return self.analysis.r if self.analysis else None

    @property
    def percentile(self) -> Optional[float]:
        return self.analysis.percentile if self.analysis else None

    @property
    def cis_p(self) -> Optional[float]:
        return self.cis.best_p if self.cis else None

    def to_dict(self) -> dict:
        d = {
            'gene': self.gene, 'driver': self.driver,
            'matrix_format': self.matrix_format,
            'verdict': self.verdict, 'error': self.error,
            'design_verdict': self.design_verdict,
            'design_issues': self.design_issues,
        }
        if self.analysis:
            for k, v in self.analysis.to_dict().items():
                d[f'analysis_{k}'] = v
        if self.cis:
            for k, v in self.cis.to_dict().items():
                d[f'cis_{k}'] = v
        return d


@dataclass
class GeneLocation:
    """Gene location from GTF."""
    chrom: str
    start: int
    end: int
    symbol: str
    ensembl_id: str
    biotype: str = ''
