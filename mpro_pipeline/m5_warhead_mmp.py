"""
Module 5 — covalent-warhead matched-pair test: what does opt_score actually see?

THE QUESTION
Boltz-2 treats binding as a noncovalent equilibrium, yet M3 measured
rho = +0.779 against enzymatic pIC50 on a chemotype-diverse set. Nirmatrelvir
(nitrile), GC373 (aldehyde) and PF-00835231 (hydroxymethylketone) all modify
Cys145 covalently, so where does that signal come from? An earlier hypothesis
in this work was that the model is warhead-blind and the correlation rides
entirely on P1/P2 occupancy -- which, if true, would mean covalent docking
could be skipped for ranking. This module tests it and REFUTES it.

DESIGN — matched molecular pairs mined from the measured data, not invented
m1's 3,629-record set is mined for compounds sharing a core after the warhead
atoms are deleted (11 warhead SMARTS classes). 29 such clusters exist; the 4
most informative were co-folded against the Mpro dimer (44 compounds, same
target/pocket/settings as M3):

  cluster 0  n=12  3 warheads   measured span 3.33 log
  cluster 1  n= 9  6 warheads   measured span 2.17 log  (nirmatrelvir scaffold)
  cluster 2  n=17  2 warheads   measured span 3.40 log  (15 aldehydes)
  cluster 3  n= 6  3 warheads   measured span 3.56 log  (PF-00835231 scaffold)

Cluster 1 varies the warhead on a fixed core; cluster 2 fixes the warhead
(aldehyde) and varies recognition. Together they separate the two axes.

WHAT HAPPENED, IN ORDER — the power trap
1. Pre-registered per-cluster criterion (rho >= 0.5, p < 0.05, CI above 0)
   returned BLIND for all four clusters: rho = +0.538 / +0.661 / +0.439 /
   -0.143, every bootstrap CI spanning 0. Read alone, that says the model has
   no within-series resolving power.
2. That reading would have been wrong. n = 6-17 per cluster cannot resolve a
   moderate effect. A blocked analysis -- z-scoring within each cluster, then
   pooling, which uses only within-cluster information -- gives
   rho = +0.536, p = 0.0002, and a within-cluster permutation test p = 0.0001
   against a null whose 95% upper bound is +0.262. Cluster-level bootstrap CI
   [+0.203, +0.615].
   ! The blocked test is POST-HOC: it was added after the per-cluster tests
     came back underpowered. It is a power fix rather than metric shopping
     (the permutation null guards the p-value), but it was not pre-registered
     and is labelled accordingly.
3. Blocked rho (+0.536) EXCEEDS the unstratified rho on the same 44 compounds
   (+0.360). 57% of opt_score variance is between-cluster, but that component
   is not aligned with potency -- it is a per-scaffold offset that adds noise
   when chemotypes are pooled. The model does better inside a series than
   across series.

THE ANSWER — pairwise concordance, which is how the score is actually used
(does it pick the better of two candidates?) Pairs with measured difference
below 0.5 log are excluded, that being the label noise floor.

  warhead swaps                  61/87  = 70.1%   p = 0.0001
  recognition swaps, same warhead 74/103 = 71.8%   p < 0.0001
  all within-cluster pairs       135/190 = 71.1%
  large effects (>= 1.5 log)      62/78  = 79.5%
    of those, warhead swaps       28/35  = 80.0%   p = 0.0003

So the model is NOT warhead-blind, and it is no better at recognition than at
warheads -- the two axes come out within 2 points of each other. The
hypothesis this module set out to test is refuted, and with it the hope of
justifying "skip covalent docking because the warhead does not register".

WHAT THIS DOES LICENSE
Noncovalent co-folding already calls ~70% of within-series pairs correctly,
~80% when the real gap is >= 1.5 log. That is a useful, cheap baseline and it
sets the bar a covalent/QM-MM method must clear to be worth its cost. It is
not reliable ranking: 3 pairs in 10 come out backwards, and some warhead
comparisons fail outright -- alpha-ketoamide vs chloroacetamide 17/17, but
aldehyde vs chloroacetamide 8/16, a coin flip (n per pair type is small; treat
the breakdown as a hypothesis).

DYNAMIC-RANGE COMPRESSION [limit]
  cluster 0  3.33 log -> 0.477 score units   ~7.0 log per unit
  cluster 1  2.17 log -> 0.314                ~6.9
  cluster 2  3.40 log -> 0.287               ~11.8
  cluster 3  3.56 log -> 0.090               ~39.4
Run-to-run precision is ~0.003-0.008 (the two nirmatrelvir stereoisomer
entries score 0.597 vs 0.599 here, 0.624 vs 0.616 in M3), so the spread is
real signal rather than noise -- but a 0.02 score gap can hide more than a log
of potency. Never read small opt_score differences as meaningful.
"""
import collections
import itertools
import json
import os

import numpy as np
from scipy import stats

METRIC = "optimization_score"
MIN_DY = 0.5        # label noise floor: ignore pairs closer than this
BIG_DY = 1.5        # "large effect" threshold
BOOT, SEED = 20000, 20260925


def load():
    here = os.path.join(os.path.dirname(__file__), "data")
    with open(os.path.join(here, "wh_meta.json")) as fh:
        meta = json.load(fh)
    with open(os.path.join(here, "wh_pred.json")) as fh:
        pred = json.load(fh)
    rows = [dict(id=k, cl=v["cluster"], wh=v["warhead"], y=v["pIC50"],
                 x=float(pred[k][METRIC]))
            for k, v in meta.items() if k in pred and METRIC in pred[k]]
    by = collections.defaultdict(list)
    for r in rows:
        by[r["cl"]].append(r)
    return rows, by


def blocked_rho(by):
    """Within-cluster z-score then pool: uses only within-series information."""
    xs, ys = [], []
    for v in by.values():
        x = np.array([r["x"] for r in v]); y = np.array([r["y"] for r in v])
        if x.std() == 0 or y.std() == 0:
            continue
        xs.append((x - x.mean()) / x.std())
        ys.append((y - y.mean()) / y.std())
    X, Y = np.concatenate(xs), np.concatenate(ys)
    rho, p = stats.spearmanr(X, Y)
    rng = np.random.default_rng(SEED)
    null = []
    for _ in range(BOOT):
        px = []
        for arr in xs:
            a = arr.copy(); rng.shuffle(a); px.append(a)
        null.append(stats.spearmanr(np.concatenate(px), Y)[0])
    null = np.asarray(null)
    return rho, p, float((null >= rho).mean()), float(np.percentile(null, 95))


def concordance(pairs):
    n = len(pairs)
    if n == 0:
        return None
    k = sum(1 for a, b in pairs if (a["x"] - b["x"]) * (a["y"] - b["y"]) > 0)
    p = stats.binomtest(k, n, 0.5, alternative="greater").pvalue
    return k, n, k / n, float(p)


def report():
    rows, by = load()
    print("=" * 72)
    print("MODULE 5 — covalent-warhead matched pairs: what does opt_score see?")
    print("=" * 72)
    print(f"{len(rows)} compounds in {len(by)} matched-core clusters\n")

    print("per-cluster (pre-registered criterion; all underpowered):")
    for cl in sorted(by):
        v = by[cl]
        x = np.array([r["x"] for r in v]); y = np.array([r["y"] for r in v])
        rho, p = stats.spearmanr(x, y)
        whs = len({r["wh"] for r in v})
        print(f"  cluster {cl}  n={len(v):<3} warheads={whs}  "
              f"rho={rho:+.3f} p={p:.4f}   measured span {y.max() - y.min():.2f} log"
              f" -> score span {x.max() - x.min():.3f}")

    rho, p, pperm, null95 = blocked_rho(by)
    print(f"\nblocked within-cluster analysis [post-hoc, power fix]:")
    print(f"  rho = {rho:+.3f}  p = {p:.4f}   permutation p = {pperm:.4f}"
          f"  (null 95% upper = {null95:+.3f})")
    gx = np.array([r["x"] for r in rows]); gy = np.array([r["y"] for r in rows])
    grho, _ = stats.spearmanr(gx, gy)
    print(f"  unstratified rho on the same 44 = {grho:+.3f}"
          f"   -> within-series is BETTER by {rho - grho:+.3f}")

    buckets = collections.defaultdict(list)
    for v in by.values():
        for a, b in itertools.combinations(v, 2):
            if abs(a["y"] - b["y"]) < MIN_DY:
                continue
            buckets["warhead swap" if a["wh"] != b["wh"]
                    else "recognition swap"].append((a, b))
            buckets["all"].append((a, b))
            if abs(a["y"] - b["y"]) >= BIG_DY:
                buckets["large effect"].append((a, b))
                if a["wh"] != b["wh"]:
                    buckets["large effect, warhead"].append((a, b))

    print(f"\npairwise concordance (pairs >= {MIN_DY} log apart):")
    for kind in ("warhead swap", "recognition swap", "all",
                 "large effect", "large effect, warhead"):
        r = concordance(buckets[kind])
        if r:
            k, n, f, pv = r
            print(f"  {kind:<24}{k:>4}/{n:<4} = {100 * f:>5.1f}%   p = {pv:.4f}")

    print("\nverdict: NOT warhead-blind. Warhead and recognition swaps are")
    print("  called at the same rate (70.1% vs 71.8%), so the hypothesis that")
    print("  covalent docking is skippable because warheads do not register is")
    print("  refuted. What is licensed: noncovalent co-folding is a ~70-80%")
    print("  pairwise baseline within a series, and that is the bar a")
    print("  covalent/QM-MM method must beat to justify its cost.")
    return dict(blocked_rho=float(rho), perm_p=pperm,
                warhead=concordance(buckets["warhead swap"]),
                recognition=concordance(buckets["recognition swap"]))


if __name__ == "__main__":
    report()
