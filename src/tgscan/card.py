"""Evidence card generator — per-(driver, gene) evidence dossier for tgscan v0.3.

Replaces the one-word verdict with a transparent card: every dataset, pooled
statistics, heterogeneity, structural capture, and automatic flags.
Aggregation only (no GTF / no matrix parsing required).
"""
from __future__ import annotations

import csv
import json
import math
from typing import Optional

from .design import design_issue_is_blocking


# ---------------------------------------------------------------- math (ported from Task25/p2_meta_analysis.py, assertion-tested there)

def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _chi2_sf(x: float, k: int) -> float:
    """Regularized upper incomplete gamma Q(k/2, x/2) — series / continued fraction."""
    a, xx = k / 2, x / 2
    if xx < a + 1:
        term = 1.0 / a
        s, n_ = term, a
        for _ in range(500):
            n_ += 1
            term *= xx / n_
            s += term
            if abs(term) < abs(s) * 1e-14:
                break
        P = s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
        return max(0.0, min(1.0, 1.0 - P))
    b = xx + 1 - a
    c = 1e300
    d = 1 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-12:
            break
    Q = math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h
    return max(0.0, min(1.0, Q))


def fisher_z_pool(rs, ns):
    """Fixed-effect pooling; returns (pooled_r, lo, hi, I2_pct) or (None,)*4."""
    if not rs:
        return None, None, None, None
    zs = [math.atanh(min(max(r, -0.9995), 0.9995)) for r in rs]
    ws = [n - 3 for n in ns]
    zbar = sum(w * z for w, z in zip(ws, zs)) / sum(ws)
    se = 1 / math.sqrt(sum(ws))
    pooled = math.tanh(zbar)
    lo, hi = math.tanh(zbar - 1.96 * se), math.tanh(zbar + 1.96 * se)
    assert min(rs) - 1e-9 <= pooled <= max(rs) + 1e-9, "pooled r out of input range"
    i2 = None
    if len(rs) >= 2:
        Q = sum(w * (z - zbar) ** 2 for w, z in zip(ws, zs))
        dfree = len(rs) - 1
        i2 = max(0.0, (Q - dfree) / Q * 100) if Q > 0 else 0.0
    return pooled, lo, hi, i2


def fisher_combine_p(ps):
    if not ps:
        return None
    X = sum(-2 * math.log(max(p, 1e-300)) for p in ps)
    p = _chi2_sf(X, 2 * len(ps))
    assert p <= min(ps) + 1e-9, "Fisher combined p must be <= min input"
    return p


# ---------------------------------------------------------------- card

def load_store(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def make_card(driver: str, gene: str, store_rows, ssot_row: Optional[dict] = None) -> dict:
    rows = [r for r in store_rows
            if r.get("driver") == driver and r.get("candidate") == gene]
    ds = []
    usable_rs, usable_ns, observed_rs = [], [], []
    cis_by_gse = {}
    for r in sorted(rows, key=lambda x: (x.get("gse") or "")):
        rv, n = _f(r.get("r")), _f(r.get("n_samples"))
        status = (r.get("status") or "").strip()
        di = (r.get("design_issues") or "").strip()
        design_bad = design_issue_is_blocking(di)
        clean = status in ("HIGH_CONFIDENCE", "SSOT_AGGREGATE") and not design_bad
        usable = clean and rv is not None and (n or 0) >= 4 and abs(rv) < 0.9995
        # 状态不洁但 r 有效的行: 展示+参与符号翻转观察, 不进合并池
        if rv is not None:
            observed_rs.append(rv)
        note = ""
        if r.get("source_file") == "ssot_aggregate":
            note = "ssot_aggregate"
        elif design_bad:
            note = f"design:{di}"
        elif not clean and status:
            note = status
        elif rv is None:
            note = "r_unavailable"
        elif (n or 0) < 4:
            note = "n_lt4"
        elif abs(rv) >= 0.9995:
            note = "r_degenerate"
        if usable:
            usable_rs.append(rv)
            usable_ns.append(int(n))
        cp = _f(r.get("cis_best_p"))
        # cis 只从设计干净的行合并——设计死行的 cis p 同样是伪相关产物
        # (Cerkl 案例 08-20: usable k=0 但 headline cis 1.28e-05 全来自 scRNA 行)
        if cp is not None and not design_bad:
            cis_by_gse[r.get("gse")] = cp
        ds.append({"gse": r.get("gse"), "allele": r.get("allele"),
                   "n": int(n) if n is not None else None,
                   "r": rv, "pct": _f(r.get("pct")),
                   "mean_r": _f(r.get("mean_r")),
                   "cis_best_p": cp, "usable": usable, "note": note})

    pooled_r, lo, hi, i2 = fisher_z_pool(usable_rs, usable_ns)
    sign_flip = len(observed_rs) >= 2 and min(observed_rs) < 0 < max(observed_rs)
    spread = (max(usable_rs) - min(usable_rs)) if len(usable_rs) >= 2 else None
    cis_combined = fisher_combine_p(sorted(set(cis_by_gse.values()))) if cis_by_gse else None

    flags, labels = [], []
    k = len(usable_rs)
    if k == 1:
        flags.append("single-dataset")
    if k and any(n < 6 for n in usable_ns):
        flags.append("small-n")
    if spread is not None and spread > 0.4:
        flags.append(f"r-spread {spread:.2f}")
    if any((d["r"] is not None and not d["usable"] and d["note"]
            and d["note"] not in ("ssot_aggregate",)) for d in ds):
        pass  # 状态不洁行已在数据集行显示 note
    if sign_flip:
        flags.append("sign-flip")
    if i2 is not None and i2 > 75:
        flags.append(f"I2 {i2:.0f}%")
    if rows and all(r.get("source_file") == "ssot_aggregate" for r in rows):
        flags.append("lineage-gap")
    if ssot_row is not None and _f(ssot_row.get("gene_in_bac_pct")) == 0 \
            and ssot_row.get("evidence_level") != "L1":
        flags.append("capture-unverified")
    # labels
    if not rows and ssot_row is not None and ssot_row.get("evidence_level") == "L1":
        labels.append("L1-channel: full-depth re-alignment evidence (no matrix rows by design)")
    if k >= 2 and not sign_flip and pooled_r is not None and pooled_r >= 0.5 \
            and (cis_combined is not None and cis_combined < 1e-3):
        labels.append("replicated")
    if sign_flip or (i2 is not None and i2 > 75):
        labels.append("needs-review")

    structure = None
    if ssot_row is not None:
        structure = {"in_bac_pct": _f(ssot_row.get("gene_in_bac_pct")),
                     "bac_clone": ssot_row.get("bac_clone"),
                     "dist_kb": _f(ssot_row.get("gene_distance_kb"))}
    return {"gene": gene, "driver": driver,
            "ssot": {k: ssot_row.get(k) for k in
                     ("status", "evidence_level", "dataset")} if ssot_row else None,
            "datasets": ds, "n_usable": k,
            "pooled_r": pooled_r, "r_ci95": [lo, hi], "I2_pct": i2,
            "r_spread": spread, "sign_flip": sign_flip,
            "cis_combined_p": cis_combined,
            "flags": flags, "labels": labels, "structure": structure}


def _fmt(x, nd=3):
    if x is None:
        return "–"
    return f"{x:.{nd}f}"


def render_text(card: dict) -> str:
    L = []
    head = f"{card['driver']} → {card['gene']}"
    if card.get("ssot"):
        head += f"   [{card['ssot'].get('status')}/{card['ssot'].get('evidence_level')}]"
    L.append(head)
    L.append("─" * max(46, len(head)))
    if not card["datasets"]:
        L.append("  (no per-dataset matrix evidence rows)")
    for d in card["datasets"]:
        bits = [f"{d['gse']}", f"n={d['n']}", f"r={_fmt(d['r'])}",
                f"pct={_fmt(d['pct'], 1)}", f"bg={_fmt(d['mean_r'], 2)}"]
        if d["cis_best_p"] is not None:
            bits.append(f"cis={d['cis_best_p']:.1e}")
        if d["note"]:
            bits.append(f"[{d['note']}]")
        L.append("  " + "  ".join(bits))
    L.append("")
    L.append(f"  pooled r = {_fmt(card['pooled_r'])} "
             f"[{_fmt(card['r_ci95'][0])}, {_fmt(card['r_ci95'][1])}]"
             f"   I² = {_fmt(card['I2_pct'], 0) if card['I2_pct'] is not None else '–'}%"
             f"   usable k = {card['n_usable']}")
    if card["cis_combined_p"] is not None:
        L.append(f"  cis combined (Fisher) = {card['cis_combined_p']:.2e}")
    if card["structure"]:
        s = card["structure"]
        L.append(f"  structure: in_bac={_fmt(s['in_bac_pct'], 0)}%  "
                 f"clone={s['bac_clone']}  dist={_fmt(s['dist_kb'], 1)} kb")
    if card["flags"]:
        L.append("  ⚠ flags: " + ", ".join(card["flags"]))
    if card["labels"]:
        L.append("  ▸ " + " | ".join(card["labels"]))
    return "\n".join(L)


def render_json(card: dict) -> str:
    return json.dumps(card, indent=1, ensure_ascii=False)


CARD_TSV_FIELDS = ["gene", "driver", "ssot_status", "evidence_level", "n_usable",
                   "pooled_r", "r_ci_lo", "r_ci_hi", "I2_pct", "r_spread",
                   "sign_flip", "cis_combined_p", "flags", "labels"]


def render_tsv_row(card: dict) -> str:
    vals = [card["gene"], card["driver"],
            (card.get("ssot") or {}).get("status", ""),
            (card.get("ssot") or {}).get("evidence_level", ""),
            card["n_usable"], _fmt(card["pooled_r"]), _fmt(card["r_ci95"][0]),
            _fmt(card["r_ci95"][1]), _fmt(card["I2_pct"], 0) if card["I2_pct"] is not None else "",
            _fmt(card["r_spread"], 2) if card["r_spread"] is not None else "",
            int(card["sign_flip"]),
            f"{card['cis_combined_p']:.2e}" if card["cis_combined_p"] is not None else "",
            ";".join(card["flags"]), ";".join(card["labels"])]
    return "\t".join(str(v) for v in vals)
