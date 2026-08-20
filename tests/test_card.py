"""test_card.py — evidence card generator unit tests (synthetic fixture only).

数学性质断言沿用 Task25/p2_meta_analysis 的验证口径:
- pooled r 必落在输入 min/max 内
- Fisher 合并 p ≤ 最小单源 p
- 状态不洁(HIGH_CONFIDENCE 之外)的行不入池但参与 sign-flip 观察
"""
import math
from tgscan import card


def _store(tmp_path, rows):
    import csv
    p = tmp_path / "store.tsv"
    fields = ["gse", "allele", "driver", "candidate", "dist_kb", "construct",
              "n_samples", "status", "r", "pct", "mean_r", "n_genes",
              "cis_best_p", "cis_verdict", "source_file"]
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def test_pool_math_properties():
    rs, ns = [0.9, 0.6], [8, 20]
    pooled, lo, hi, i2 = card.fisher_z_pool(rs, ns)
    assert min(rs) <= pooled <= max(rs)
    assert lo < pooled < hi
    assert i2 is not None and i2 >= 0


def test_fisher_combine_dominance():
    p = card.fisher_combine_p([0.01, 0.02])
    assert p <= 0.01 + 1e-12
    assert card.fisher_combine_p([]) is None


def test_dirty_rows_excluded_from_pool_but_flag_sign_flip(tmp_path):
    rows = [
        {"gse": "GSE1", "allele": "A1", "driver": "Drv", "candidate": "G1",
         "n_samples": "10", "status": "HIGH_CONFIDENCE", "r": "0.9", "pct": "99",
         "mean_r": "0.01", "cis_best_p": "0.001", "source_file": "t.tsv"},
        {"gse": "GSE2", "allele": "A2", "driver": "Drv", "candidate": "G1",
         "n_samples": "8", "status": "BACKGROUND_TOO_HIGH", "r": "-0.85",
         "pct": "10", "mean_r": "0.7", "cis_best_p": "", "source_file": "t.tsv"},
    ]
    p = _store(tmp_path, rows)
    c = card.make_card("Drv", "G1", card.load_store(p))
    assert c["n_usable"] == 1                     # 脏行不入池
    assert c["pooled_r"] == 0.9
    assert c["sign_flip"] is True                 # 但符号翻转被观察并标记
    assert "needs-review" in c["labels"]
    assert any(d["note"] == "BACKGROUND_TOO_HIGH" for d in c["datasets"])


def test_replicated_label_requires_two_clean_same_direction(tmp_path):
    rows = [
        {"gse": "GSE1", "allele": "A1", "driver": "Drv", "candidate": "G1",
         "n_samples": "8", "status": "HIGH_CONFIDENCE", "r": "0.9", "pct": "99",
         "mean_r": "0.01", "cis_best_p": "0.005", "source_file": "t.tsv"},
        {"gse": "GSE2", "allele": "A2", "driver": "Drv", "candidate": "G1",
         "n_samples": "18", "status": "HIGH_CONFIDENCE", "r": "0.8", "pct": "98",
         "mean_r": "-0.1", "cis_best_p": "0.0009", "source_file": "t.tsv"},
    ]
    c = card.make_card("Drv", "G1", rows)
    assert c["n_usable"] == 2
    assert c["cis_combined_p"] < 1e-3
    assert "replicated" in c["labels"]
    assert not c["sign_flip"]


def test_l1_channel_graceful(tmp_path):
    p = _store(tmp_path, [])
    c = card.make_card("Prnp", "X", card.load_store(p),
                       ssot_row={"status": "confirmed", "evidence_level": "L1",
                                 "gene_in_bac_pct": "100", "gene_distance_kb": "0.2"})
    assert c["n_usable"] == 0
    assert any("L1-channel" in lab for lab in c["labels"])
    assert c["pooled_r"] is None


def test_n_lt4_not_pooled():
    rows = [{"gse": "GSE1", "allele": "A", "driver": "Drv", "candidate": "G",
             "n_samples": "3", "status": "HIGH_CONFIDENCE", "r": "0.99",
             "pct": "99", "mean_r": "0", "cis_best_p": "", "source_file": "t"}]
    c = card.make_card("Drv", "G", rows)
    assert c["n_usable"] == 0
    assert any(d["note"] == "n_lt4" for d in c["datasets"])


def test_score_components():
    from tgscan import score as S
    c = {"pooled_r": 0.9, "cis_combined_p": 1e-4, "n_usable": 2,
         "sign_flip": False, "I2_pct": 0.0, "datasets": [1], "structure": {"in_bac_pct": 100}}
    f = S.features_from_card(c)
    assert abs(f[0] - 0.9) < 1e-9 and f[3] == 1.0
    c2 = dict(c, sign_flip=True, I2_pct=95)
    f2 = S.features_from_card(c2)
    assert f2[4] == 1.0 and f2[5] == 0.95
    m = S.fit_logistic([f, f2], [1.0, 0.0])
    p_hi = S.apply_model(m, f); p_lo = S.apply_model(m, f2)
    assert p_hi > p_lo
