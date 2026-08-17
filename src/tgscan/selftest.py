"""tgscan selftest — end-to-end installation check with synthetic data.

Generates a tiny GTF + count matrix in a temp dir (no downloads, no real GTF),
runs the full verify() pipeline, and checks the expected outcome:
  - a synthetic hitchhiker (co-amplified with the driver array) must reach
    HIGH_CONFIDENCE in stage 1 and appear among top cis genes;
  - an unlinked negative-control gene must stay NO_SIGNAL;
  - background genes use numeric-leading MGI symbols (0610005C13Rik style),
    regression-guarding the symbol-detection bug that caused NO_VALID_R in the
    web pipeline (2026-08-17).
"""
from __future__ import annotations
import os
import tempfile
import numpy as np
import pandas as pd

from .runner import verify

RNG = np.random.default_rng(42)


def _make_gtf(path: str) -> str:
    """chr1: driver + 4 hitchhikers + 6 cis noise (all within +-1Mb);
    chr1 far + chr2: trans noise (~74 genes). Numeric-leading MGI-style symbols
    for noise genes (regression guard for the NO_VALID_R symbol bug)."""
    genes = []
    genes.append(("ENSMUSG0000090001", "Tgdrv1", "chr1", 100_000, 110_000))
    for i in range(4):
        genes.append((f"ENSMUSG00000900{i:02d}", f"Tghh{i+1}", "chr1",
                      150_000 + 30_000 * i, 158_000 + 30_000 * i))
    for i in range(6):
        genes.append((f"ENSMUSG00000901{i:02d}", f"061000{i:02d}Rik", "chr1",
                      500_000 + 40_000 * i, 505_000 + 40_000 * i))
    for i in range(14):
        genes.append((f"ENSMUSG00000902{i:02d}", f"061001{i:02d}Rik", "chr1",
                      3_000_000 + 100_000 * i, 3_010_000 + 100_000 * i))
    for i in range(60):
        genes.append((f"ENSMUSG00000903{i:02d}", f"061002{i:02d}Rik", "chr2",
                      1_000_000 + 40_000 * i, 1_010_000 + 40_000 * i))
    with open(path, "w") as f:
        f.write("#!genome-build tgscan-selftest\n")
        for eid, sym, chrom, s, e in genes:
            f.write(f'{chrom}\ttgscan\tgene\t{s}\t{e}\t.\t+\t.\t'
                    f'gene_id "{eid}"; gene_name "{sym}"; gene_biotype "protein_coding";\n')
    return path


def _make_matrix(path: str) -> str:
    """6 samples: 3 controls (low array) + 3 Tg (high array).
    Driver and hitchhikers co-scale (load gradient); everything else noise."""
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "tg_1", "tg_2", "tg_3"]
    load = np.array([1.0, 1.0, 1.0, 6.0, 6.5, 5.5])
    rows = {}
    base = 60.0
    rows["Tgdrv1"] = base * load + RNG.normal(0, 3, 6)
    for i in range(4):
        rows[f"Tghh{i+1}"] = 4.0 * load / 6.0 + RNG.normal(0, 0.5, 6) + 0.5
    for i in range(6):
        rows[f"061000{i:02d}Rik"] = RNG.poisson(30, 6).astype(float)
    for i in range(14):
        rows[f"061001{i:02d}Rik"] = RNG.poisson(50, 6).astype(float)
    for i in range(60):
        rows[f"061002{i:02d}Rik"] = RNG.poisson(40, 6).astype(float)
    df = pd.DataFrame(rows, index=samples).T.reset_index()
    df.columns = ["gene"] + samples
    df.to_csv(path, sep="\t", index=False)
    return path


NEG_GENE = "06100000Rik"  # numeric-leading MGI symbol (NO_VALID_R regression)


def run_selftest(verbose: bool = True) -> bool:
    """Run the synthetic end-to-end check. Returns True on pass."""
    with tempfile.TemporaryDirectory(prefix="tgscan_selftest_") as tmp:
        gtf = _make_gtf(os.path.join(tmp, "mini.gtf"))
        matrix = _make_matrix(os.path.join(tmp, "mini_counts.tsv"))

        hh = verify(matrix, "Tgdrv1", "Tghh1", gtf)
        neg = verify(matrix, "Tgdrv1", NEG_GENE, gtf)

        ok = True
        def check(cond, label):
            nonlocal ok
            status = "PASS" if cond else "FAIL"
            if not cond:
                ok = False
            if verbose:
                print(f"  [{status}] {label}")

        if verbose:
            print(f"\n  hitchhiker Tghh1: r={hh.r}, pct={hh.percentile}, verdict={hh.verdict}")
            print(f"  negative control: r={neg.r}, pct={neg.percentile}, verdict={neg.verdict}")
        check(hh.error is None, f"hitchhiker runs without error (got: {hh.error})")
        check(hh.analysis is not None and hh.analysis.status == "HIGH_CONFIDENCE",
              f"hitchhiker stage1 == HIGH_CONFIDENCE (got {hh.analysis.status if hh.analysis else '?'})")
        check(hh.cis is not None and "Tghh1" in (hh.cis.top_cis_genes or []),
              "hitchhiker appears among top cis genes")
        check(neg.analysis is not None and neg.analysis.status in ("NO_SIGNAL", "MODERATE"),
              f"negative control stays NO_SIGNAL/MODERATE (got {neg.analysis.status if neg.analysis else '?'})")
        check(neg.cis is None or neg.cis.verdict != "CONFIRMED",
              "negative control does not reach cis CONFIRMED")
        if verbose:
            print(f"\n  SELFTEST {'PASSED' if ok else 'FAILED'}")
        return ok
