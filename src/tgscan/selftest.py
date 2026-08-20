"""tgscan selftest — end-to-end installation check with synthetic data.

Generates a tiny GTF + count matrix in a temp dir (no downloads, no real GTF),
runs the full verify() pipeline, and checks the expected outcome:
  - a synthetic hitchhiker (co-amplified with the driver array) must reach
    HIGH_CONFIDENCE in stage 1 and appear among top cis genes;
  - an unlinked negative-control gene must stay NO_SIGNAL;
  - background genes use numeric-leading MGI symbols (0610005C13Rik style),
    regression-guarding the symbol-detection bug that caused NO_VALID_R in the
    web pipeline (2026-08-17);
  - v0.3 design gate: FACS-fraction / TRAP / single-cell designs are BLOCKED
    before correlation (ground truth: GSE83356/115934/127845); a clean bulk
    design passes; the curated registry blocks GSE83356 by accession.
"""
from __future__ import annotations
import os
import tempfile
import numpy as np
import pandas as pd

from .runner import verify
from .design import detect_design_issues, screen_design, lookup_registry

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


def _design_gate_matrix(path: str, sample_names, n_bg: int = 20) -> str:
    """Dense small matrix with caller-chosen sample names (design signal only).

    n_bg default 20 keeps files tiny; pass n_bg >= 15000 to build a matrix
    that also clears filter_too_strict (needed for pure registry tests).
    """
    n = len(sample_names)
    rows = {"Tgdrv1": RNG.poisson(100, n).astype(float),
            "Tghh1": RNG.poisson(40, n).astype(float)}
    for i in range(n_bg):
        rows[f"0610{i:05d}Rik"] = RNG.poisson(30, n).astype(float)
    df = pd.DataFrame(rows, index=sample_names).T.reset_index()
    df.columns = ["gene"] + list(sample_names)
    df.to_csv(path, sep="\t", index=False)
    return path


def _single_cell_matrix(path: str) -> str:
    """30 'cells' as samples, 90% zeros — cells-as-sample signature."""
    n = 30
    names = [f"cell_{i}" for i in range(n)]
    rows = {}
    for i in range(80):
        v = RNG.poisson(0.08, n).astype(float)  # ~92% zeros
        rows[f"061004{i:02d}Rik"] = v
    rows["Tgdrv1"] = RNG.poisson(0.4, n).astype(float)
    rows["Tghh1"] = RNG.poisson(0.4, n).astype(float)
    df = pd.DataFrame(rows, index=names).T.reset_index()
    df.columns = ["gene"] + names
    df.to_csv(path, sep="\t", index=False)
    return path


NEG_GENE = "06100000Rik"  # numeric-leading MGI symbol (NO_VALID_R regression)


def run_selftest(verbose: bool = True) -> bool:
    """Run the synthetic end-to-end check. Returns True on pass."""
    with tempfile.TemporaryDirectory(prefix="tgscan_selftest_") as tmp:
        gtf = _make_gtf(os.path.join(tmp, "mini.gtf"))
        matrix = _make_matrix(os.path.join(tmp, "mini_counts.tsv"))

        # legacy pipeline checks run with the design gate off: the synthetic
        # matrix is deliberately tiny (filter_too_strict) — gate tests below
        # use dedicated matrices.
        hh = verify(matrix, "Tgdrv1", "Tghh1", gtf, gate=False)
        neg = verify(matrix, "Tgdrv1", NEG_GENE, gtf, gate=False)

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

        # ---- v0.3 design gate (positive/negative controls) ----
        facs_m = _design_gate_matrix(os.path.join(tmp, "facs.tsv"),
                                     ["GFP+1", "GFP+2", "GFP-1", "GFP-2", "unsorted_1", "unsorted_2"])
        trap_m = _design_gate_matrix(os.path.join(tmp, "trap.tsv"),
                                     ["CA1_TRAP_1", "CA1_TRAP_2", "CA1_input_1", "CA1_input_2"])
        clean_m = _design_gate_matrix(os.path.join(tmp, "clean.tsv"),
                                      ["AI7685", "AI7686", "AI7687", "AI7688", "AI7689", "AI7690"],
                                      n_bg=15000)
        sc_m = _single_cell_matrix(os.path.join(tmp, "sc.tsv"))

        from .parsers import parse_matrix, numeric_sample_cols

        def _issues(path):
            df, _ = parse_matrix(path)
            sc = numeric_sample_cols(df, df.columns[0])
            return detect_design_issues(df, sc), df, sc

        facs_i, _, _ = _issues(facs_m)
        check('facs_sorted' in facs_i, f"FACS +/- pair detected (got {facs_i})")

        # paired plus/minus words are equally sorting; single-arm is NOT
        pm_m = _design_gate_matrix(os.path.join(tmp, "pm.tsv"),
                                   ["Venus_plus_22", "Venus_plus_24", "Venus_minus_1", "Venus_minus_2"])
        pm_i, _, _ = _issues(pm_m)
        check('facs_sorted' in pm_i, f"paired Venus_plus/minus detected (got {pm_i})")
        single_m = _design_gate_matrix(os.path.join(tmp, "single.tsv"),
                                       ["K_WT_Venus_plus_22", "K_WT_Venus_plus_24", "F_Ko16", "F_Ko29"])
        single_i, _, _ = _issues(single_m)
        check('facs_sorted' not in single_i,
              f"single-arm Venus_plus (genotype-style label, GSE211929) NOT blocked (got {single_i})")
        trap_i, _, _ = _issues(trap_m)
        check('trap_rna_ip' in trap_i, f"TRAP/input design detected (got {trap_i})")
        sc_i, _, _ = _issues(sc_m)
        check('single_cell_matrix' in sc_i, f"sparse cells-as-samples detected (got {sc_i})")
        clean_i, _, _ = _issues(clean_m)
        check(not any(k in ' '.join(clean_i) for k in ('facs', 'trap', 'single_cell')),
              f"clean bulk design passes design detectors (got {clean_i})")

        blocked = verify(facs_m, "Tgdrv1", "Tghh1", gtf)
        check(blocked.design_verdict == 'BLOCK' and 'facs_sorted' in blocked.design_issues
              and blocked.verdict == "BLOCKED_DESIGN",
              f"verify() blocks FACS design end-to-end (got {blocked.design_verdict}: {blocked.design_issues})")

        df_clean, _ = parse_matrix(clean_m)
        sc_clean = numeric_sample_cols(df_clean, df_clean.columns[0])
        clean_issues = detect_design_issues(df_clean, sc_clean)
        check(clean_issues == [],
              f"clean 15k-gene bulk matrix has zero design issues (got {clean_issues})")
        reg_ok = screen_design(df_clean, sc_clean, gse="GSE83356")
        check(reg_ok['verdict'] == 'BLOCK' and reg_ok['source'] == 'registry',
              "registry alone blocks GSE83356 on an otherwise clean matrix")
        check(lookup_registry("GSE94145") is None,
              "unaudited GSE has no registry opinion")
        check(lookup_registry("GSE146304")['verdict'].upper() == 'WEAK',
              "registry carries WEAK verdicts")
        if verbose:
            print(f"\n  SELFTEST {'PASSED' if ok else 'FAILED'}")
        return ok
