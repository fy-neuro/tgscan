"""score.py — v0.3 provisional probability scorer (transparent logistic).

设计声明（防"自创框架"）: 本模块是工程工具, 不是新判定教条——
权重系数全部公开, 校准于项目裁决集(SSOT vs store negatives),
LOO-CV 报告, PI 定调前仅作排序参考。输出含组件分解, 逐项可审计。

特征(v1, 裸统计量——2026-08-20 三变体对照后的诚实选择):
  x1 = pooled_r(usable; 无则单源 r)
  x2 = min(-log10(cis_combined), 8)/8
  x3 = min(n_usable, 3)/3
  x4 = in_bac(1/0/0.5-unknown)
  x5 = sign_flip(1/0)  [负向]
  x6 = min(I2,100)/100 [负向]

对照实验(BAYESIAN_DESIGN.md §七, 同一裁决集 141 对):
  v1 裸特征      LOO AUC 0.977   ← 保留
  v2 P_theta     LOO AUC 0.882   (后验饱和: 0.99/0.998/1.0 挤压排序粒度)
  v3 后验 z 值   LOO AUC 0.912   (不饱和但设计闸门清场后无增益)
结论: 设计闸门已移除"高 r 但假"的案例, 裸 r 在干净负例集上判别力最强;
贝叶斯通道(P_theta/τ²/signP)的价值在 card 报告层(连续谱/异质性/sign
怀疑度), 不在本校准器。两者互补, 不是替代。
"""
from __future__ import annotations

import csv
import math
from typing import Optional

from . import card as card_mod

FEATS = ["pooled_r", "cis_strength", "n_usable_n", "in_bac", "sign_flip", "i2_n"]


def features_from_card(c: dict) -> Optional[list]:
    if c["pooled_r"] is None and not c["datasets"]:
        return None
    r = c["pooled_r"]
    if r is None:  # 无可用合并时取任一有效 r 的最大绝对值方向保守值
        vals = [d["r"] for d in c["datasets"] if d["r"] is not None]
        r = max(vals) if vals else 0.0
    cis = c["cis_combined_p"]
    cis_strength = min(-math.log10(max(cis, 1e-8)), 8) / 8 if cis else 0.0
    st = c.get("structure") or {}
    in_bac = 0.5 if st.get("in_bac_pct") is None else (1.0 if st["in_bac_pct"] > 0 else 0.0)
    return [r, cis_strength, min(c["n_usable"], 3) / 3, in_bac,
            1.0 if c["sign_flip"] else 0.0,
            min(c["I2_pct"] or 0, 100) / 100 if c["I2_pct"] is not None else 0.0]


def fit_logistic(X, y, lam=1.0, iters=8000, lr=0.15):
    n, d = len(X), len(X[0])
    mu = [sum(row[j] for row in X) / n for j in range(d)]
    sd = [max(1e-9, (sum((row[j] - mu[j]) ** 2 for row in X) / n) ** .5) for j in range(d)]
    Xz = [[(row[j] - mu[j]) / sd[j] for j in range(d)] for row in X]
    w, b = [0.0] * d, 0.0
    for _ in range(iters):
        pr = [1 / (1 + math.exp(-(sum(a * c_ for a, c_ in zip(Xz[i], w)) + b))) for i in range(n)]
        gw = [sum(Xz[i][j] * (pr[i] - y[i]) for i in range(n)) / n + lam * w[j] / n for j in range(d)]
        gb = sum(pr[i] - y[i] for i in range(n)) / n
        w = [wj - lr * gj for wj, gj in zip(w, gw)]
        b -= lr * gb
    return {"w": w, "b": b, "mu": mu, "sd": sd}


def apply_model(m, x):
    z = m["b"] + sum(m["w"][j] * (x[j] - m["mu"][j]) / m["sd"][j] for j in range(len(x)))
    return 1 / (1 + math.exp(-z))


def auc(pairs):
    lab = [p[0] for p in pairs]
    sc = [p[1] for p in pairs]
    order = sorted(range(len(sc)), key=lambda i: sc[i])
    ranks, i = [0.0] * len(sc), 0
    while i < len(sc):
        j = i
        while j + 1 < len(sc) and sc[order[j + 1]] == sc[order[i]]:
            j += 1
        rk = (i + j) / 2 + 1
        for t in range(i, j + 1):
            ranks[order[t]] = rk
        i = j + 1
    n1, n0 = sum(lab), len(lab) - sum(lab)
    if n1 == 0 or n0 == 0:
        return None
    return (sum(r for r, l in zip(ranks, lab) if l) - n1 * (n1 + 1) / 2) / (n1 * n0)


def build_training(store_rows, ssot_rows):
    """positives = SSOT confirmed+candidate; negatives = store 中状态 NO_SIGNAL 的配对。"""
    pos_keys = {(r["driver"], r["gene"]) for r in ssot_rows}
    X, y, meta = [], [], []
    neg_seen = set()
    for r in store_rows:
        key = (r.get("driver"), r.get("candidate"))
        if key in pos_keys:
            continue
        st = (r.get("status") or "")
        if st != "NO_SIGNAL" or key in neg_seen:
            continue
        neg_seen.add(key)
        c = card_mod.make_card(key[0], key[1], store_rows)
        f = features_from_card(c)
        if f is None:
            continue
        X.append(f); y.append(0.0); meta.append((key[0], key[1], "negative"))
    for r in ssot_rows:
        c = card_mod.make_card(r["driver"], r["gene"], store_rows, r)
        f = features_from_card(c)
        if f is None:
            continue
        X.append(f); y.append(1.0); meta.append((r["driver"], r["gene"], r["status"]))
    return X, y, meta


def loo_report(X, y):
    out = []
    for i in range(len(X)):
        tr = [k for k in range(len(X)) if k != i]
        m = fit_logistic([X[k] for k in tr], [y[k] for k in tr])
        out.append((y[i], apply_model(m, X[i])))
    return auc(out)


def loo_predictions(X, y):
    """LOO 的 (y_true, p_pred) 序列——校准曲线的输入。"""
    out = []
    for i in range(len(X)):
        tr = [k for k in range(len(X)) if k != i]
        m = fit_logistic([X[k] for k in tr], [y[k] for k in tr])
        out.append((y[i], apply_model(m, X[i])))
    return out


def calibration_table(preds, n_bins: int = 5):
    """可靠性表：[(bin_lo, bin_hi, n, mean_pred, observed_frac), ...]。

    诚实口径：正例仅 ~29 个，bin 内样本小，observed_frac 噪声大——
    论文中只作可靠性图参考，不作判定依据。
    """
    edges = [i / n_bins for i in range(n_bins + 1)]
    tab = []
    for i in range(n_bins):
        grp = [(y, p) for y, p in preds if edges[i] <= p < edges[i + 1] or
               (i == n_bins - 1 and p == 1.0)]
        if not grp:
            tab.append((edges[i], edges[i + 1], 0, None, None))
            continue
        ys = [g[0] for g in grp]
        tab.append((edges[i], edges[i + 1], len(grp),
                    sum(p for _, p in grp) / len(grp), sum(ys) / len(ys)))
    return tab


def score_pair(driver, gene, store_rows, ssot_row, model):
    c = card_mod.make_card(driver, gene, store_rows, ssot_row)
    f = features_from_card(c)
    if f is None:
        return c, None, None
    p = apply_model(model, f)
    contrib = {FEATS[j]: round(model["w"][j] * (f[j] - model["mu"][j]) / model["sd"][j], 3)
               for j in range(len(f))}
    return c, p, contrib
