"""bayes.py — Bayesian evidence channels for tgscan (BAYESIAN_DESIGN.md).

设计声明（防"自创框架"）：本模块是透明证据通道的实现，不是新判定教条。
所有先验显式披露（PRIOR_*），先验敏感性可由 prior_sensitivity() 复算；
θ 后验用正态共轭闭式解（k=2-4 数据集规模下无需 MCMC），系数/公式全部公开。

设计 A — 分层随机效应合并（DerSimonian-Laird τ² + 正态后验）:
    z_i = atanh(r_i) ~ N(θ, 1/(n_i-3) + τ²);  θ ~ N(PRIOR_MU, PRIOR_SD²)
    输出 pooled r / 95% CI / P(θ>0.5)（真相关超过 0.5 的后验概率）
设计 D — sign 一致性后验:
    数据集级 P(符号为正) = π ~ Beta(1,1)；P(π>0.5 | 观测符号) 为连续怀疑分
    （全正 k=4 → 0.969；一半反向 → 0.500——替代二元 sign-flip flag 的刻度）
"""
from __future__ import annotations

import math
from typing import List, Optional, Sequence, Tuple

# ---- 披露的先验（改动必须同步 BAYESIAN_DESIGN.md 先验敏感性表）----
PRIOR_MU = 0.0      # θ 先验均值（z 尺度）
PRIOR_SD = 1.5      # θ 先验标准差（弱信息：r 尺度约 [-0.91, 0.91] 之外才受抑制）
THETA_THRESHOLD = 0.5   # P(θ > 0.5)：报告"真相关 >0.5"的后验概率


def _clamp_r(r: float) -> float:
    return min(max(r, -0.9995), 0.9995)


def _norm_sf(x: float) -> float:
    return 0.5 * math.erfc(x / math.sqrt(2))


def dl_random_effects(rs: Sequence[float], ns: Sequence[int]
                      ) -> Tuple[float, float, float]:
    """DerSimonian-Laird 随机效应合并。

    Returns:
        (tau2, z_hat, se) — 异质性方差、合并 z 后验均值、标准误
    """
    zs = [math.atanh(_clamp_r(r)) for r in rs]
    w = [n - 3 for n in ns]
    sw = sum(w)
    zf = sum(wi * zi for wi, zi in zip(w, zs)) / sw
    Q = sum(wi * (zi - zf) ** 2 for wi, zi in zip(w, zs))
    k = len(rs)
    if k > 1:
        c = sw - sum(wi * wi for wi in w) / sw
        tau2 = max(0.0, (Q - (k - 1)) / c) if c > 0 else 0.0
    else:
        tau2 = 0.0
    w_re = [1.0 / (1.0 / wi + tau2) for wi in w]
    z = sum(wi * zi for wi, zi in zip(w_re, zs)) / sum(w_re)
    se = 1.0 / math.sqrt(sum(w_re))
    return tau2, z, se


def posterior_theta(rs: Sequence[float], ns: Sequence[int],
                    prior_sd: float = PRIOR_SD,
                    threshold: float = THETA_THRESHOLD
                    ) -> dict:
    """设计 A 主入口：θ 后验（正态共轭）。

    Returns:
        dict(pooled_r, ci_lo, ci_hi, tau2, se, P_theta)
    """
    tau2, z_hat, se = dl_random_effects(rs, ns)
    z0 = math.atanh(threshold)
    prec = 1.0 / prior_sd**2 + 1.0 / se**2
    mu = (PRIOR_MU / prior_sd**2 + z_hat / se**2) / prec
    sd = math.sqrt(1.0 / prec)
    p_theta = _norm_sf((z0 - mu) / sd)
    return {
        'pooled_r': math.tanh(mu),
        'ci_lo': math.tanh(mu - 1.96 * sd),
        'ci_hi': math.tanh(mu + 1.96 * sd),
        'tau2': tau2, 'se': se, 'P_theta': p_theta,
    }


def prior_sensitivity(rs: Sequence[float], ns: Sequence[int],
                      sds: Optional[List[float]] = None) -> List[Tuple[float, float]]:
    """先验敏感性表：(prior_sd, P_theta) 序列——报告用，防"先验偷跑"。"""
    return [(sd, posterior_theta(rs, ns, prior_sd=sd)['P_theta'])
            for sd in (sds or [0.5, 1.0, 1.5, 3.0, 10.0])]


def sign_consistency(signs: Sequence[float]) -> Optional[float]:
    """设计 D：P(π>0.5 | 观测符号)，π ~ Beta(1,1) 先验。

    数值积分（k 为数据集个数，最多十余项，梯形即可稳定）。
    None 当观测 <2（无信息）。
    """
    s = sum(1 for x in signs if x > 0)
    f = len(signs) - s
    if len(signs) < 2:
        return None
    a, b = 1 + s, 1 + f
    logB = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    n = 4000
    h = 0.5 / n
    dens = lambda t: t ** (a - 1) * (1 - t) ** (b - 1)
    cdf05 = sum(dens(h * (i + 0.5)) for i in range(n)) * h
    return 1.0 - cdf05 / math.exp(logB)
