"""Validation of the van Elteren implementation in src/de_methods.py."""

import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from de_methods import van_elteren, wilcoxon_ranksum  # noqa: E402

SEED = 42


def test_single_stratum_matches_scipy():
    """(a) One stratum: p-value equals scipy mannwhitneyu asymptotic (tie corrected)."""
    rng = np.random.default_rng(SEED)
    n_genes, n = 500, 60
    is_tumor = np.zeros(n, bool)
    is_tumor[:23] = True
    # mix of continuous and heavily tied (count-like) genes
    x = np.vstack([rng.normal(size=(250, n)) + 0.5 * is_tumor,
                   rng.poisson(3, size=(250, n)).astype(float)])
    ve = van_elteren(x, is_tumor, np.zeros(n, int))
    ref = stats.mannwhitneyu(x[:, is_tumor], x[:, ~is_tumor], axis=1, alternative="two-sided",
                             method="asymptotic", use_continuity=False).pvalue
    np.testing.assert_allclose(ve["pvalue"], ref, rtol=0, atol=1e-8)
    # the project's plain Wilcoxon (P4) is the same test
    np.testing.assert_allclose(wilcoxon_ranksum(x, is_tumor)[0], ref, rtol=0, atol=1e-12)


def _confounded_null(rng, n_genes=4000):
    """Strata with very different levels, unbalanced tumor fraction, no condition effect."""
    sizes = [60, 60, 60, 60, 60]
    tumor_frac = [0.1, 0.3, 0.5, 0.7, 0.9]
    level = [0.0, 3.0, 6.0, 9.0, 12.0]
    strata, is_tumor, mu = [], [], []
    for k, (n, f, lv) in enumerate(zip(sizes, tumor_frac, level)):
        strata += [k] * n
        t = np.zeros(n, bool)
        t[: int(round(f * n))] = True
        is_tumor.append(t)
        mu += [lv] * n
    strata = np.array(strata)
    is_tumor = np.concatenate(is_tumor)
    x = rng.normal(size=(n_genes, len(strata))) + np.array(mu)
    return x, is_tumor, strata


def test_stratified_null_is_uniform_and_plain_wilcoxon_is_inflated():
    """(b) Under a strata-confounded null van Elteren keeps ~5% FPR; plain Wilcoxon does not."""
    rng = np.random.default_rng(SEED)
    x, is_tumor, strata = _confounded_null(rng)
    p_ve = van_elteren(x, is_tumor, strata)["pvalue"]
    p_w = wilcoxon_ranksum(x, is_tumor)[0]
    fpr_ve = float((p_ve < 0.05).mean())
    fpr_w = float((p_w < 0.05).mean())
    ks = stats.kstest(p_ve, "uniform")
    print(f"\nvan Elteren FPR={fpr_ve:.4f} KS p={ks.pvalue:.3f} | plain Wilcoxon FPR={fpr_w:.4f}")
    assert 0.04 <= fpr_ve <= 0.06
    assert ks.pvalue > 0.01
    assert fpr_w > 0.5


def test_one_condition_strata_are_dropped():
    rng = np.random.default_rng(SEED)
    x, is_tumor, strata = _confounded_null(rng, n_genes=200)
    extra = rng.normal(size=(200, 20)) + 50
    x2 = np.hstack([x, extra])
    t2 = np.concatenate([is_tumor, np.ones(20, bool)])
    s2 = np.concatenate([strata, np.full(20, 99)])
    a = van_elteren(x, is_tumor, strata)
    b = van_elteren(x2, t2, s2)
    np.testing.assert_allclose(a["pvalue"], b["pvalue"], atol=1e-12)
    assert b["n_used"] == a["n_used"] == 300
