"""The six vcc2026 metrics, implemented from Arc's published specification.

Ported from ``docs/vcc2026_metrics/vcc2026-metrics-brief.md`` in
`ArcInstitute/cell-eval2` (0.15.0, competition ``rule_version`` 3).  This is a
reimplementation rather than a wrapper, for two reasons: `cell-eval2`'s scoring
path wants the official reference bundles, which are not released, and the
point of having it here is to tune against the metrics on *public* data where
no bundle exists.

:mod:`virtualcell.metrics` implements the **2025** suite and is kept as-is for
comparability with this project's earlier results.  It is not the 2026
objective and should not be used for new tuning.  The differences that matter:

* discrimination is now **cosine**, so it is invariant to a global rescaling of
  the prediction.  The 2025 L1 version rose with prediction magnitude (Liu
  et al., arXiv 2511.16954), which is why ``Hyper.beta`` was tunable against
  it.  That lever is gone.
* **four of the six** metrics are now differential-expression metrics computed
  by a Wilcoxon rank-sum test **on the submitted cells**, so a submission's
  within-perturbation dispersion changes its score directly.
* ``overlap_at_k`` is gone.  Set agreement is a Jaccard index over
  BH-significant sets, which penalises an invented gene exactly as hard as a
  missed one -- and therefore caps at ``n_real/n_pred`` for a submission that
  over-calls.

Every public function here is checked against an invariant the specification
states explicitly (a constant submission scores exactly 0.5 on ``pds_cosine``,
a zero-fold-change submission scores exactly 1.0 on ``nmae``, and so on); see
``tests/test_vcc2026_metrics.py``.  Those invariants are the only available
ground truth without the reference bundles, and they catch the errors that
matter.

    python -m virtualcell.vcc2026_metrics --self-test
"""

from __future__ import annotations

import numpy as np

TARGET_SUM_BULK = 5.0e4      # pseudobulk profile normalisation, sections 1-2
TARGET_SUM_CELL = 1.0e6      # per-cell normalisation for the DE test
MIN_CPM = 5.0                # low-expression gate, control cells only
ALPHA = 0.05                 # BH threshold, per perturbation
EPS_LFC = 1.0e-9
REACH_PURITY_FLOOR = 0.9
MIN_GATE_SIZE = 10
NMAE_FLOOR = -6.0


# ---------------------------------------------------------------------------
# profiles
# ---------------------------------------------------------------------------

def bulk_profile(counts: np.ndarray, target_sum: float = TARGET_SUM_BULK
                 ) -> np.ndarray:
    """``log(1 + TS * P_g / sum_g P_g)`` over a group's summed raw counts.

    ``counts`` is (cells, genes) raw counts for one group; the sum is taken
    first and normalised once, which is not the same as normalising each cell
    and averaging.
    """
    total = np.asarray(counts, dtype=np.float64).sum(axis=0)
    denom = total.sum()
    if denom <= 0:
        return np.zeros_like(total)
    return np.log1p(target_sum * total / denom)


def jackknife_noise(counts: np.ndarray, target_sum: float = TARGET_SUM_BULK
                    ) -> float:
    """Delete-one jackknife estimate of a bulk profile's sampling variance.

    ``C_p = (n-1)/n * sum_g sum_i (v_ig - mean_g)^2`` with ``v_ig`` the profile
    recomputed with cell ``i`` removed.  Zero for fewer than two cells.
    """
    y = np.asarray(counts, dtype=np.float64)
    n = y.shape[0]
    if n < 2:
        return 0.0
    total = y.sum(axis=0)
    S = total.sum()
    ell = y.sum(axis=1)
    denom = S - ell                                    # (cells,)
    if np.any(denom <= 0):
        keep = denom > 0
        if keep.sum() < 2:
            return 0.0
        y, total, ell, denom = y[keep], total, ell[keep], denom[keep]
        n = y.shape[0]
    v = np.log1p(target_sum * (total[None, :] - y) / denom[:, None])
    dev = v - v.mean(axis=0, keepdims=True)
    return float((n - 1) / n * (dev ** 2).sum())


# ---------------------------------------------------------------------------
# 1. perturbation discrimination
# ---------------------------------------------------------------------------

def pds_cosine(pred_delta: np.ndarray, real_delta: np.ndarray) -> np.ndarray:
    """Per-perturbation cosine discrimination; 1 best, 0.5 uninformative.

    Both arrays are (n_pert, n_genes) effects against the *reference* control,
    already restricted to the common feature space with all panel target genes
    removed -- the specification removes all 300, not just each row's own, so
    that every distance is computed in one fixed space.  Row ``i`` of each is
    the same perturbation.

    Ties share the average of the positions they span, which is what makes a
    constant submission score exactly 0.5 rather than something seed-dependent.
    """
    P = np.ascontiguousarray(pred_delta, dtype=np.float64)
    R = np.ascontiguousarray(real_delta, dtype=np.float64)
    n = P.shape[0]
    if R.shape[0] != n:
        raise ValueError("pred and real must carry the same perturbations")

    pn = np.linalg.norm(P, axis=1)
    rn = np.linalg.norm(R, axis=1)
    # d = 1 - cos, and exactly 1 wherever either vector has zero norm
    with np.errstate(invalid="ignore", divide="ignore"):
        d = 1.0 - (P @ R.T) / np.outer(pn, rn)
    d[~np.isfinite(d)] = 1.0
    d[pn <= 0, :] = 1.0
    d[:, rn <= 0] = 1.0

    own = np.diagonal(d).copy()
    less = (d < own[:, None]).sum(axis=1)
    equal = (d == own[:, None]).sum(axis=1)
    k = less + 0.5 * (equal - 1)
    return 1.0 - k / (n - 1)


# ---------------------------------------------------------------------------
# 2. expression error
# ---------------------------------------------------------------------------

def expr_mse_unbiased_capped_norm(
        pred_b: np.ndarray, real_b: np.ndarray, ctrl_b: np.ndarray,
        pred_C: np.ndarray, real_C: np.ndarray, ctrl_C: float,
        keep: np.ndarray) -> float:
    """Panel-level ratio of noise-corrected prediction error to effect size.

    ``keep`` is a boolean (n_pert, n_genes) mask: True where the gene counts
    for that perturbation, i.e. everything except its own target gene.  Returns
    a single number for the panel, clamped to [0, 1]; there is no
    per-perturbation value.
    """
    Gp = keep.sum(axis=1).astype(np.float64)
    Gp = np.maximum(Gp, 1.0)

    sq_pred = (((pred_b - real_b) ** 2) * keep).sum(axis=1)
    sq_real = (((real_b - ctrl_b[None, :]) ** 2) * keep).sum(axis=1)
    credit = np.minimum(pred_C, real_C)

    # rho bounds the prediction's credited noise correction by the submission's
    # own across-perturbation spread, so a submission cannot claim a correction
    # larger than the variation it actually exhibits.
    w = keep.astype(np.float64) / Gp[:, None]
    wsum = w.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        wmean = np.where(wsum > 0, (w * pred_b).sum(axis=0) / wsum, 0.0)
    B = float((w * (pred_b - wmean[None, :]) ** 2).sum())
    pool = float((credit / Gp).sum())
    rho = min(1.0, B / pool) if pool > 0 else 1.0

    N = (sq_pred - rho * credit - real_C) / Gp
    D = (sq_real - real_C - ctrl_C) / Gp
    denom = float(D.sum())
    if denom == 0:
        return 1.0
    return float(np.clip(N.sum() / denom, 0.0, 1.0))


# ---------------------------------------------------------------------------
# 3. the differential-expression table
# ---------------------------------------------------------------------------

def cpm(counts: np.ndarray) -> np.ndarray:
    """Per-cell normalisation to a total of 1e6, as the DE test ranks."""
    y = np.asarray(counts, dtype=np.float64)
    tot = y.sum(axis=1, keepdims=True)
    return np.divide(y, tot, out=np.zeros_like(y), where=tot > 0) * TARGET_SUM_CELL


def expression_gate(ctrl_cpm: np.ndarray, min_cpm: float = MIN_CPM
                    ) -> np.ndarray:
    """Which genes are tested: mean CPM above the floor in *control* cells.

    The gate reads the reference's control group alone, so it admits the same
    gene set for every submission -- which is why it can be computed once and
    reused across arms of a comparison.
    """
    return ctrl_cpm.mean(axis=0) > min_cpm


def wilcoxon_p(group: np.ndarray, ctrl: np.ndarray, chunk: int = 1500
               ) -> np.ndarray:
    """Two-sided Mann-Whitney p per gene, normal approximation with tie correction.

    Vectorised over genes and chunked over them, because the intermediates are
    (cells x genes) and about eight of them are live at once: unchunked, a
    20,000-cell control pool against 8,000 genes wants ~10 GB and gets killed.
    ``chunk`` trades peak memory for a few percent of speed and changes no
    result.
    """
    if group.shape[1] > chunk:
        return np.concatenate([
            wilcoxon_p(group[:, i:i + chunk], ctrl[:, i:i + chunk], chunk)
            for i in range(0, group.shape[1], chunk)])

    n1, n2 = group.shape[0], ctrl.shape[0]
    stacked = np.concatenate([group, ctrl], axis=0)
    n = n1 + n2
    order = np.argsort(stacked, axis=0, kind="stable")
    srt = np.take_along_axis(stacked, order, axis=0)

    # Tied blocks, fully vectorised over genes.  For each element find the first
    # and last index of its tied run; the average rank of the run is then
    # (start + end)/2 + 1, and a run of length t contributes t^3 - t to the tie
    # correction -- which is the same as every one of its t elements
    # contributing t^2 - 1, so the correction needs no per-gene grouping either.
    # This matters for cost: the measurement below runs this hundreds of times
    # on (400 + 20,000) cells x ~8,000 genes, and a per-gene Python loop over
    # bincount dominated everything else.
    idx = np.arange(n)[:, None]
    newblock = np.empty(srt.shape, dtype=bool)
    newblock[0] = True
    np.not_equal(srt[1:], srt[:-1], out=newblock[1:])

    start = np.maximum.accumulate(np.where(newblock, idx, 0), axis=0)
    endmark = np.empty_like(newblock)
    endmark[:-1] = newblock[1:]
    endmark[-1] = True
    end = np.minimum.accumulate(np.where(endmark, idx, n - 1)[::-1], axis=0)[::-1]

    ranks = (start + end) * 0.5 + 1.0
    back = np.empty_like(ranks)
    np.put_along_axis(back, order, ranks, axis=0)

    R1 = back[:n1].sum(axis=0)
    U = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    t = (end - start + 1).astype(np.float64)
    tie = (t * t - 1.0).sum(axis=0)
    var = n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1)))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (U - mu) / np.sqrt(var)
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    from scipy.stats import norm
    return 2.0 * norm.sf(np.abs(z))


def log2fc(group_cpm: np.ndarray, ctrl_mean: np.ndarray) -> np.ndarray:
    """``log2((mean_pert + eps) / (mean_ctrl + eps))`` on CPM means."""
    return np.log2((group_cpm.mean(axis=0) + EPS_LFC) / (ctrl_mean + EPS_LFC))


def bh_reject(p: np.ndarray, alpha: float = ALPHA) -> np.ndarray:
    """Benjamini-Hochberg rejection mask over the genes handed in."""
    n = p.size
    if n == 0:
        return np.zeros(0, dtype=bool)
    order = np.argsort(p, kind="stable")
    passing = p[order] <= alpha * np.arange(1, n + 1) / n
    out = np.zeros(n, dtype=bool)
    if passing.any():
        out[order[: np.flatnonzero(passing)[-1] + 1]] = True
    return out


# ---------------------------------------------------------------------------
# 4-7. the DE metrics
# ---------------------------------------------------------------------------

def direction_fidelity(pred_sig: np.ndarray, pred_lfc: np.ndarray,
                       real_sig: np.ndarray, real_lfc: np.ndarray) -> float:
    """``k / max(n_pred, n_real)``: directional precision with a coverage cap.

    All four arrays are 1-D over the gated genes for a single perturbation,
    with that perturbation's own target gene already removed.  Returns NaN when
    the reference found nothing and the submission called nothing, which the
    specification drops from the mean.
    """
    adjudicable = np.isfinite(real_lfc) & (real_lfc != 0)
    call = pred_sig & adjudicable
    n_pred = int((pred_sig).sum())        # a called gene without a direction
    n_real = int(real_sig.sum())          # stays in n_pred and counts as a miss
    if n_pred == 0 and n_real == 0:
        return float("nan")
    if n_pred == 0:
        return 0.0
    k = int((np.sign(pred_lfc[call]) == np.sign(real_lfc[call])).sum())
    return k / max(n_pred, n_real)


def direction_reach(pred_sig: np.ndarray, pred_padj: np.ndarray,
                    pred_p: np.ndarray, pred_lfc: np.ndarray,
                    real_sig: np.ndarray, real_lfc: np.ndarray,
                    gene_names: np.ndarray,
                    purity_floor: float = REACH_PURITY_FLOOR) -> float:
    """How deep into the submission's own ranking its directions stay >=90% right.

    The pool is the reference's significant set.  Order is: significant calls
    first, then ascending predicted adjusted p, then ascending predicted p,
    then descending |predicted lfc|, then gene name.  ``k*`` is the *deepest*
    prefix clearing the floor, not the first -- purity is not monotone.
    """
    pool = np.flatnonzero(real_sig)
    if pool.size == 0:
        return float("nan")
    order = np.lexsort((gene_names[pool], -np.abs(pred_lfc[pool]),
                        pred_p[pool], pred_padj[pool], ~pred_sig[pool]))
    ranked = pool[order]

    adj = np.isfinite(real_lfc[ranked]) & (real_lfc[ranked] != 0)
    ranked = ranked[adj]
    if ranked.size == 0:
        return 0.0
    match = (np.sign(pred_lfc[ranked]) == np.sign(real_lfc[ranked]))
    purity = np.cumsum(match) / np.arange(1, match.size + 1)
    ok = np.flatnonzero(purity >= purity_floor)
    k_star = int(ok[-1] + 1) if ok.size else 0
    return k_star / float(pool.size)


def sig_jaccard(pred_sig: np.ndarray, real_sig: np.ndarray) -> float:
    """``|R n R^| / |R u R^|``, with an empty union defined as 1."""
    inter = int((pred_sig & real_sig).sum())
    union = int((pred_sig | real_sig).sum())
    return 1.0 if union == 0 else inter / union


def lfc_nmae(pred_lfc: np.ndarray, real_lfc: np.ndarray,
             real_sig: np.ndarray, min_gate: int = MIN_GATE_SIZE
             ) -> float:
    """Normalised absolute fold-change error on the reference's significant set.

    A submission predicting zero everywhere makes numerator equal denominator
    term by term and scores exactly 1 -- the specification's "no skill" point.
    Returns NaN when the gate holds fewer than ``min_gate`` genes.
    """
    gate = real_sig & np.isfinite(real_lfc)
    if int(gate.sum()) < min_gate:
        return float("nan")
    pred = np.where(np.isfinite(pred_lfc[gate]), pred_lfc[gate], 0.0)
    denom = float(np.abs(real_lfc[gate]).sum())
    if denom == 0 or not np.isfinite(denom):
        return float("nan")
    return float(np.abs(pred - real_lfc[gate]).sum() / denom)


# ---------------------------------------------------------------------------
# scaling
# ---------------------------------------------------------------------------

def scale(u: float, baseline: float, replicate: float) -> float:
    """``(u - b) / (r - b)``: 0 is the context mean, 1 is a replicate.

    Not clipped.  The specification reports values above 1 and below 0 rather
    than clamping them, and clamping while *tuning* is actively harmful: once a
    metric is below baseline, making it worse becomes free.
    """
    span = replicate - baseline
    return float("nan") if span == 0 else float((u - baseline) / span)


# ---------------------------------------------------------------------------
# self-test against the specification's stated invariants
# ---------------------------------------------------------------------------

def _self_test() -> None:
    rng = np.random.default_rng(0)
    n_pert, n_gene = 40, 300

    print("pds_cosine")
    real = rng.normal(size=(n_pert, n_gene))
    got = pds_cosine(real, real).mean()
    print(f"  exact reproduction                  {got:.4f}   expect 1.0000")
    assert abs(got - 1.0) < 1e-12

    shared = np.tile(rng.normal(size=(1, n_gene)), (n_pert, 1))
    got = pds_cosine(shared, real).mean()
    print(f"  one shared profile                  {got:.4f}   expect 0.5000")
    assert abs(got - 0.5) < 1e-12, got

    got = pds_cosine(np.zeros((n_pert, n_gene)), real).mean()
    print(f"  zero effect                         {got:.4f}   expect 0.5000")
    assert abs(got - 0.5) < 1e-12, got

    print("\nlfc_nmae")
    real_lfc = rng.normal(size=n_gene)
    real_sig = rng.random(n_gene) < 0.2
    got = lfc_nmae(np.zeros(n_gene), real_lfc, real_sig)
    print(f"  zero fold change everywhere         {got:.4f}   expect 1.0000")
    assert abs(got - 1.0) < 1e-12
    got = lfc_nmae(real_lfc, real_lfc, real_sig)
    print(f"  exact reproduction                  {got:.4f}   expect 0.0000")
    assert abs(got) < 1e-12
    got = lfc_nmae(real_lfc, real_lfc, np.zeros(n_gene, bool))
    print(f"  gate below {MIN_GATE_SIZE} genes                  "
          f"{got}      expect nan")
    assert np.isnan(got)

    print("\nsig_jaccard")
    a = rng.random(n_gene) < 0.1
    print(f"  identical sets                      {sig_jaccard(a, a):.4f}"
          f"   expect 1.0000")
    assert sig_jaccard(a, a) == 1.0
    print(f"  disjoint sets                       "
          f"{sig_jaccard(a, ~a & (rng.random(n_gene) < 0.5)):.4f}"
          f"   expect 0.0000")
    z = np.zeros(n_gene, bool)
    print(f"  empty union                         {sig_jaccard(z, z):.4f}"
          f"   expect 1.0000")
    assert sig_jaccard(z, z) == 1.0

    # the arithmetic cap that the wrap-up's over-calling finding runs into
    real_set = np.zeros(n_gene, bool); real_set[:30] = True
    pred_set = np.zeros(n_gene, bool); pred_set[:300] = True
    got = sig_jaccard(pred_set, real_set)
    print(f"  call 300 when truth is 30           {got:.4f}   expect "
          f"{30 / 300:.4f}  (= n_real/n_pred)")
    assert abs(got - 0.1) < 1e-12

    print("\ndirection_fidelity")
    real_lfc = rng.normal(size=n_gene)
    real_sig = np.zeros(n_gene, bool); real_sig[:50] = True
    perfect = direction_fidelity(real_sig, real_lfc, real_sig, real_lfc)
    print(f"  perfect calls and directions        {perfect:.4f}   expect 1.0000")
    assert abs(perfect - 1.0) < 1e-12
    got = direction_fidelity(np.zeros(n_gene, bool), real_lfc, real_sig, real_lfc)
    print(f"  calls nothing, reference found 50   {got:.4f}   expect 0.0000")
    assert got == 0.0
    # random directions over a call set at coverage 1 -> 0.5 in expectation
    big = np.zeros(n_gene, bool); big[:200] = True
    vals = [direction_fidelity(big, rng.normal(size=n_gene), real_sig, real_lfc)
            for _ in range(200)]
    print(f"  random directions, coverage >=1     {np.mean(vals):.4f}   "
          f"expect ~0.5000")
    assert abs(np.mean(vals) - 0.5) < 0.02

    print("\ndirection_reach")
    names = np.array([f"g{i:04d}" for i in range(n_gene)])
    padj = np.full(n_gene, 0.001)
    got = direction_reach(real_sig, padj, padj, real_lfc, real_sig, real_lfc,
                          names)
    print(f"  perfect directions                  {got:.4f}   expect 1.0000")
    assert abs(got - 1.0) < 1e-12
    # a single wrong call at the very first position, never recovering, is the
    # only way to drive k* to zero
    bad = real_lfc.copy(); bad[:] = -real_lfc
    got = direction_reach(real_sig, padj, padj, bad, real_sig, real_lfc, names)
    print(f"  every direction wrong               {got:.4f}   expect 0.0000")
    assert got == 0.0
    # The floor at P_0 = 0.9 has a sharp consequence the specification spells
    # out: a leading miss followed by nine matches gives P(10) = 9/10 and k*
    # reaches 10, but with only eight following matches P(9) = 8/9 < 0.9 and no
    # shorter prefix qualifies either, so k* collapses to 0.  Testing the pair
    # pins the floor; testing a 50-gene pool does not, because 49/50 = 0.98
    # clears the floor at full depth and k* is the *deepest* crossing, not the
    # first.
    rank_order = np.arange(n_gene, dtype=float)
    for pool_n, expect in ((10, 1.0), (9, 0.0)):
        sig = np.zeros(n_gene, bool); sig[:pool_n] = True
        pl = real_lfc.copy()
        pl[0] = -real_lfc[0]                     # one miss, at position 1
        got = direction_reach(sig, padj, rank_order, pl, sig, real_lfc, names)
        print(f"  leading miss, pool of {pool_n:<2}            {got:.4f}   "
              f"expect {expect:.4f}")
        assert abs(got - expect) < 1e-12, (pool_n, got)

    print("\nexpr_mse_unbiased_capped_norm")
    ctrl_b = rng.random(n_gene) * 2
    real_b = ctrl_b[None, :] + rng.normal(scale=0.3, size=(n_pert, n_gene))
    keep = np.ones((n_pert, n_gene), bool)
    for i in range(n_pert):
        keep[i, i % n_gene] = False
    z_C = np.zeros(n_pert)
    got = expr_mse_unbiased_capped_norm(real_b, real_b, ctrl_b, z_C, z_C, 0.0,
                                        keep)
    print(f"  exact reproduction                  {got:.4f}   expect 0.0000")
    assert abs(got) < 1e-12
    ctrl_rep = np.tile(ctrl_b[None, :], (n_pert, 1))
    got = expr_mse_unbiased_capped_norm(ctrl_rep, real_b, ctrl_b, z_C, z_C,
                                        0.0, keep)
    print(f"  emits the control unchanged         {got:.4f}   expect ~1.0000")
    assert abs(got - 1.0) < 0.05, got

    print("\nwilcoxon_p vs scipy")
    from scipy.stats import mannwhitneyu
    g = rng.poisson(5.0, size=(60, 12)).astype(float)
    c = rng.poisson(4.0, size=(90, 12)).astype(float)
    mine = wilcoxon_p(g, c)
    theirs = mannwhitneyu(g, c, axis=0, method="asymptotic",
                          use_continuity=False).pvalue
    print(f"  max |diff| over 12 genes            {np.abs(mine - theirs).max():.2e}"
          f"   expect ~0")
    assert np.allclose(mine, theirs, atol=1e-12)

    print("\nall invariants hold.")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.parse_args()
    _self_test()


if __name__ == "__main__":
    main()
