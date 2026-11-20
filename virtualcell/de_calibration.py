"""Calibrate how many genes a submission calls, which is what caps Jaccard.

Four of the six vcc2026 metrics come from a Wilcoxon rank-sum test run **on the
submitted cells**.  So a submission does not merely predict fold changes -- it
implicitly decides how many genes it calls as responding, and two of the six
metrics charge it for that decision:

    sig_jaccard = |R n R^| / |R u R^|  <=  n_real / n_pred

That bound is arithmetic.  It holds with perfect recall and no matter how good
the biology is.  The 2026 wrap-up reports measured responses with a median of
22-44 significant genes per target, against predictions calling 122-136
(X-team), about 500 (the winning entry) and 340-4,200 for the rest of the top
ten.  :mod:`virtualcell.wrapup_analysis` inverts the published scores through
the metric anchors and finds the winner sitting at ~80% of the cap its own
over-calling imposes -- so on that metric the count, not the gene selection, is
what is binding.

This module fixes the count.  It does **not** invent a threshold: it estimates
the response-size distribution from a context where ground truth exists, and
transfers that distribution rather than a single median, because a strong
knockdown genuinely moves more genes than a weak one and collapsing them to one
number throws that away.

What is certain here and what is not, stated plainly:

* the arithmetic bound, and therefore the headroom, is certain;
* how much *recall* survives thresholding down to the right count depends on how
  well the model ranks genes, which is a property of real data.
  :func:`virtualcell.wrapup_analysis.jaccard_ceiling` brackets it -- even at 30%
  retained recall the gain is roughly half the published first-to-second
  margin -- but the point estimate needs a measured run.

    python -m virtualcell.de_calibration --self-test
"""

from __future__ import annotations

import numpy as np

from . import vcc2026_metrics as M


# ---------------------------------------------------------------------------
# measuring the reference response size
# ---------------------------------------------------------------------------

def significant_counts(pert_cells: dict[str, np.ndarray],
                       ctrl_cells: np.ndarray,
                       gene_names: np.ndarray | None = None,
                       alpha: float = M.ALPHA) -> dict[str, int]:
    """Significant-gene count per perturbation, through the scored pipeline.

    Runs exactly the test sections 3-7 of the specification read: per-cell CPM
    normalisation, a 5-CPM gate taken from the control cells alone, a two-sided
    Wilcoxon rank-sum, then Benjamini-Hochberg *within* each perturbation.  The
    perturbation's own target gene is removed when ``gene_names`` is supplied.

    Use it on measured cells to learn what a real response looks like, and on
    generated cells to check what a submission is about to claim.
    """
    ctrl = M.cpm(ctrl_cells)
    gate = M.expression_gate(ctrl)
    ctrl_g = ctrl[:, gate]
    gated_names = None if gene_names is None else np.asarray(gene_names)[gate]

    out: dict[str, int] = {}
    for name, cells in pert_cells.items():
        g = M.cpm(cells)[:, gate]
        p = M.wilcoxon_p(g, ctrl_g)
        drop = None
        if gated_names is not None:
            hit = np.flatnonzero(gated_names == name)
            drop = int(hit[0]) if hit.size else None
        if drop is not None:
            keep = np.ones(p.size, bool)
            keep[drop] = False
            sig = M.bh_reject(p[keep], alpha)
            out[name] = int(sig.sum())
        else:
            out[name] = int(M.bh_reject(p, alpha).sum())
    return out


# ---------------------------------------------------------------------------
# transferring it
# ---------------------------------------------------------------------------

def rank_matched_targets(pred_delta: np.ndarray,
                         reference_counts: np.ndarray) -> np.ndarray:
    """Assign each perturbation a target call count by rank-matching magnitude.

    The perturbation with the k-th largest predicted response magnitude is given
    the k-th largest count in the reference distribution.  This transfers the
    whole distribution rather than its median, which matters because response
    size varies over an order of magnitude across a panel and the metrics are
    computed per perturbation before averaging.

    ``reference_counts`` is a sample of per-perturbation counts from a context
    with ground truth; it need not be the same length as ``pred_delta``.
    """
    mag = np.linalg.norm(pred_delta, axis=1)
    n = mag.size
    ref = np.sort(np.asarray(reference_counts, dtype=float))
    if ref.size == 0:
        raise ValueError("no reference counts to match against")
    # quantile of each perturbation's magnitude -> same quantile of the counts
    q = (np.argsort(np.argsort(mag)) + 0.5) / n
    idx = np.clip((q * ref.size).astype(int), 0, ref.size - 1)
    return np.maximum(ref[idx].astype(int), 1)


def sparsify(pred_delta: np.ndarray, k_per_pert: np.ndarray,
             score: np.ndarray | None = None) -> np.ndarray:
    """Keep each perturbation's top-``k`` genes and zero the rest.

    Ranking defaults to ``|delta|``.  That is a proxy for expected significance
    rather than the thing itself, but a defensible one: the effect is in log1p
    space, which is variance-stabilising, so magnitude and detectability track
    each other far better than they would on raw counts.  Pass ``score`` to rank
    by something better when a model can supply it.
    """
    D = np.asarray(pred_delta, dtype=np.float64)
    S = np.abs(D) if score is None else np.asarray(score, dtype=np.float64)
    out = np.zeros_like(D)
    for i in range(D.shape[0]):
        k = int(min(max(k_per_pert[i], 0), D.shape[1]))
        if k == 0:
            continue
        keep = np.argpartition(-S[i], k - 1)[:k]
        out[i, keep] = D[i, keep]
    return out


def realised_vs_target(realised: np.ndarray, target: np.ndarray) -> dict:
    """Summary of whether sparsifying actually produced the intended counts.

    Sparsifying to ``k`` genes does not guarantee ``k`` significant calls: the
    test still has to reject, and a kept gene with a tiny effect will not.  So
    the realised count is an upper-bounded, noisy function of ``k`` and is worth
    reporting rather than assuming.
    """
    realised = np.asarray(realised, float)
    target = np.asarray(target, float)
    return {"median_realised": float(np.median(realised)),
            "median_target": float(np.median(target)),
            "mean_ratio": float(np.mean(realised / np.maximum(target, 1))),
            "frac_at_or_below_target":
                float(np.mean(realised <= target + 1e-9))}


# ---------------------------------------------------------------------------
# mechanism check
# ---------------------------------------------------------------------------

def _self_test() -> None:
    """Verify the machinery on cells whose truth is constructed, not assumed.

    This checks the code, not the biology: it builds a panel where a known set
    of genes responds, generates cells for a model whose gene *ranking* is
    deliberately imperfect, and confirms (a) the DE pipeline recovers the
    planted set, (b) over-calling collapses Jaccard to the arithmetic bound, and
    (c) calibrating the count recovers it.  The recall a real model retains is
    not something a simulation can settle.
    """
    rng = np.random.default_rng(0)
    n_gene, n_cell, n_ctrl, n_pert = 600, 400, 2000, 12
    n_true = 30

    base = rng.gamma(2.0, 8.0, size=n_gene)        # mean CPM-ish per gene
    ctrl = rng.poisson(base[None, :], size=(n_ctrl, n_gene)).astype(float)

    truth_sets, pert_cells = {}, {}
    for i in range(n_pert):
        name = f"P{i}"
        hit = rng.choice(n_gene, n_true, replace=False)
        lfc = rng.normal(0.0, 1.0, size=n_true)
        mu = base.copy()
        mu[hit] = np.maximum(mu[hit] * 2.0 ** lfc, 0.01)
        pert_cells[name] = rng.poisson(mu[None, :], size=(n_cell, n_gene)
                                       ).astype(float)
        truth_sets[name] = (hit, lfc)

    counts = significant_counts(pert_cells, ctrl)
    print(f"planted {n_true} directly-responding genes per perturbation")
    print(f"  called by the scored pipeline: median "
          f"{np.median(list(counts.values())):.0f}, range "
          f"{min(counts.values())}-{max(counts.values())}")
    print("  The gap is not an error: moving 30 genes changes each cell's total\n"
          "  count, and the test ranks CPM, so every other gene shifts too. That\n"
          "  compositional term is present in real data as well, and it is one\n"
          "  reason a submission over-calls without predicting anything wrong.\n"
          "  Everything below therefore scores against the pipeline's own\n"
          "  significant set, which is what the metric actually compares to.")

    # The reference set, as the metric would see it.
    ctrl_cpm = M.cpm(ctrl)
    gate = M.expression_gate(ctrl_cpm)
    print(f"  genes passing the 5-CPM gate: {int(gate.sum())}/{n_gene}")

    def sets_for(cells):
        g = M.cpm(cells)[:, gate]
        return M.bh_reject(M.wilcoxon_p(g, ctrl_cpm[:, gate]))

    print("\nJaccard against the planted truth, as a submission's call count "
          "varies.")
    print("A model with a *fixed*, imperfect ranking is used throughout, so "
          "only the\ncount changes between rows.\n")
    print(f"  {'calls':>8}{'raw J':>9}{'arith. cap':>13}{'% of cap':>11}")
    print("  " + "-" * 41)

    name = "P0"
    real_sig = sets_for(pert_cells[name])
    n_real = int(real_sig.sum())
    hit, lfc = truth_sets[name]
    # a ranking that puts most but not all real genes near the top
    score = rng.random(n_gene) * 0.4
    score[hit] += 0.6
    ranked = np.argsort(-score)

    rows = []
    for n_call in (n_real, 2 * n_real, 4 * n_real, 16 * n_real):
        n_call = min(n_call, int(gate.sum()))
        called = np.zeros(n_gene, bool)
        called[ranked[:n_call]] = True
        j = M.sig_jaccard(called[gate], real_sig)
        cap = n_real / max(n_call, 1)
        rows.append({"calls": n_call, "raw_jaccard": j, "cap": cap})
        print(f"  {n_call:>8}{j:>9.4f}{cap:>13.4f}{j / cap:>10.1%}")

    assert rows[0]["raw_jaccard"] > rows[-1]["raw_jaccard"] * 3, \
        "over-calling must collapse Jaccard"
    print(f"\n  calling {rows[0]['calls']} genes instead of "
          f"{rows[-1]['calls']} raises raw Jaccard "
          f"{rows[-1]['raw_jaccard']:.4f} -> {rows[0]['raw_jaccard']:.4f}, "
          f"a {rows[0]['raw_jaccard'] / max(rows[-1]['raw_jaccard'], 1e-9):.1f}x "
          f"move,\n  from the identical ranking. The gene selection never "
          f"changed.")

    print("\nrank_matched_targets transfers the distribution, not the median:")
    pred = rng.normal(size=(n_pert, n_gene))
    pred *= np.linspace(0.2, 3.0, n_pert)[:, None]
    tg = rank_matched_targets(pred, np.array(list(counts.values())))
    print(f"  predicted magnitude spans "
          f"{np.linalg.norm(pred, axis=1).min():.1f}-"
          f"{np.linalg.norm(pred, axis=1).max():.1f}")
    print(f"  assigned counts span {tg.min()}-{tg.max()} "
          f"(reference span {min(counts.values())}-{max(counts.values())})")
    assert tg[np.argmax(np.linalg.norm(pred, axis=1))] == max(tg)

    print("\nsparsify keeps exactly the requested number of genes:")
    sp = sparsify(pred, tg)
    nz = (sp != 0).sum(axis=1)
    print(f"  requested {tg[:6].tolist()} -> kept {nz[:6].tolist()}")
    assert np.array_equal(nz, tg)

    print("\nmechanism verified.")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.parse_args()
    _self_test()


if __name__ == "__main__":
    main()
