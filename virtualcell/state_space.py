"""Is the reachable set of perturbation effects a low-dimensional cone?

The *cell functional state space* framework makes many proposals, most of which
need data this task does not have -- velocity fields need time series, fiber and
path-dependence need paired or repeated-stimulus designs, population flow needs
distributions at several timepoints.  The Virtual Cell Challenge gives one
steady-state timepoint per condition, so those parts are untestable here.

One proposal is testable, and the framework names this exact dataset for it
(§13.2): map every perturbation's displacement into local tangent coordinates
and ask what the set of displacements actually covers.

    Δ_u(z) = Log_z(z_u)

    "若 latent dimension 为 d，但数千个 perturbations 的 Δ_u 只稳定落在低秩
     子空间或少数角锥中，则说明细胞虽然测量维度很高，实际可控自由度可能很低。"

Replogle's genome-wide K562 arm is 9,866 displacements in 8,246 gene
coordinates, so the question is answerable directly rather than by analogy.

Two things are measured, and they are not the same claim:

**Low rank.**  How many directions carry the variance, against a null with the
same per-perturbation magnitudes but random directions.  If the answer is
"few", the measured dimension of a screen vastly overstates its controllable
degrees of freedom.

**Cone, not subspace.**  This is the sharper and more useful claim.  A subspace
is closed under negation: if a direction is reachable, so is its opposite.  A
*cone* is not.  Testing it means asking whether the coefficients along each
leading direction are one-sided.  It matters for modelling, because this
project already projects predictions onto a low-rank *subspace* (``rank_mix``)
and cross-validation chose rank 40-80 -- so if the reachable set is really a
cone, that projection is throwing away a constraint it could be enforcing.

    python -m virtualcell.state_space
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.utils.extmath import randomized_svd

from .gwps import load_gwps


def effective_rank(s: np.ndarray) -> dict:
    """Summaries of a spectrum that do not depend on where you cut it."""
    var = s ** 2
    p = var / var.sum()
    cum = np.cumsum(p)
    return {
        "participation_ratio": float(var.sum() ** 2 / (var ** 2).sum()),
        "entropy_rank": float(np.exp(-(p * np.log(p + 1e-300)).sum())),
        **{f"n_for_{int(q * 100)}pct": int(np.searchsorted(cum, q) + 1)
           for q in (0.5, 0.8, 0.9, 0.95)},
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rank", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/state_space.json"))
    args = ap.parse_args()

    src = load_gwps()
    D = src.delta                                   # (n_pert, n_gene)
    D = D - D.mean(axis=0)                          # tangent at the control state
    print(f"displacement set: {D.shape[0]:,} perturbations x {D.shape[1]:,} "
          f"gene coordinates")

    U, S, Vt = randomized_svd(D, n_components=args.rank, random_state=args.seed)
    real = effective_rank(S)

    # Null: same per-perturbation magnitudes, directions randomised.  This is the
    # right null -- it holds effect *size* fixed and destroys only the alignment
    # between perturbations, which is exactly the structure being claimed.
    rng = np.random.default_rng(args.seed)
    norms = np.linalg.norm(D, axis=1, keepdims=True)
    R = rng.normal(size=D.shape)
    R *= norms / (np.linalg.norm(R, axis=1, keepdims=True) + 1e-12)
    _, S_null, _ = randomized_svd(R, n_components=args.rank,
                                  random_state=args.seed)
    null = effective_rank(S_null)

    print(f"\n{'':<22}{'measured':>12}{'null':>12}{'ratio':>10}")
    print("-" * 56)
    for k in real:
        r = real[k] / max(null[k], 1e-9)
        print(f"{k:<22}{real[k]:>12.1f}{null[k]:>12.1f}{r:>10.3f}")

    # Cone test: along each leading direction, are the coefficients one-sided?
    # A subspace is closed under negation, a cone is not, so a leading direction
    # whose coefficients are strongly one-signed is evidence for a cone.
    coef = U[:, :20] * S[:20]
    frac = np.maximum((coef > 0).mean(0), (coef < 0).mean(0))
    # A symmetric direction sits at 0.5; the null gives the sampling floor.
    coef_null = rng.normal(size=coef.shape)
    frac_null = np.maximum((coef_null > 0).mean(0), (coef_null < 0).mean(0))
    print(f"\none-sidedness of the leading 20 directions "
          f"(0.5 = symmetric, 1.0 = a ray)")
    print(f"  measured: median {np.median(frac):.3f}, max {frac.max():.3f}, "
          f"{(frac > 0.65).sum()}/20 above 0.65")
    print(f"  null:     median {np.median(frac_null):.3f}, "
          f"max {frac_null.max():.3f}")
    print("  per-direction, strongest first:")
    print("   " + "  ".join(f"{f:.2f}" for f in frac[:10]))

    # What the model actually chose, for comparison.
    bench = Path("results/gwps_single_source.json")
    chosen = None
    if bench.exists():
        h = json.loads(bench.read_text())["hyper"]
        chosen = sorted({int(v["rank"]) for v in h.values()})

    out = {"measured": real, "null": null,
           "one_sidedness": frac.tolist(),
           "one_sidedness_null": frac_null.tolist(),
           "shape": list(D.shape), "cv_chosen_rank": chosen}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2))

    print(f"\n{real['n_for_90pct']} directions carry 90% of the variance across "
          f"{D.shape[0]:,} perturbations,")
    print(f"against {null['n_for_90pct']} for magnitude-matched random "
          f"directions.")
    if chosen:
        print(f"Cross-validation independently chose a program basis of rank "
              f"{chosen} for\nthis screen, without seeing any of the above.")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
