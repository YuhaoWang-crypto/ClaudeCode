"""A discrete connection between cell lines, and whether transporting helps.

The fiber framework asks for a local transport matrix ``R_AB`` between
neighbouring states, estimated from an ancestor-descendant coupling ``Γ_ij``.
In developmental data that coupling has to be *inferred* -- LineageOT, CellRank,
optimal transport -- and every inferred coupling adds noise to whatever is
computed on top of it.

This project has an unusual piece of luck: **the coupling is exact.**  Four
CRISPRi cell lines share 2,053 knockdowns, so knockdown *g* in line A and
knockdown *g* in line B are the same test force applied at two different base
states.  ``Γ`` is the identity on shared perturbations.  Nothing is inferred, so
whatever structure appears in ``R_AB`` is in the data rather than in a coupling
estimator.

The construction follows the proposal directly:

    coef_X(g) = U_X · Δ_X(g)                 local fiber coordinates
    M         = Σ_g coef_B(g) coef_A(g)ᵀ     cross-state covariance
    M = UΣVᵀ  →  R_AB = U Vᵀ                 orthogonal Procrustes

and the reporting obeys the gauge warning that comes with it.  A connection
matrix element is coordinate language: change the local fiber basis and it
changes.  So what is reported is the coordinate-independent part -- principal
angles between the local subspaces, the singular-value spectrum of the coupling,
path inconsistency ‖R_AC − R_BC R_AB‖, holonomy ‖H − I‖ around closed loops --
each against a null.

**And then the part that decides whether any of it matters.**  This model
currently transfers a source effect into a new line essentially as-is, modulated
by baseline expression.  If the connection is real, transporting first should
predict the target line's effect better than not transporting.  That is a
falsifiable claim about a model component, not an interpretation, and it is the
last thing this module computes.

    python -m virtualcell.connection
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import subspace_angles

from .data import load_all, shared_perturbations

K_DEFAULT = 30          # local fiber rank; the reachable set is low-rank at the top


def local_basis(delta: np.ndarray, k: int) -> np.ndarray:
    """Top-``k`` right singular vectors: a local frame for this line's fiber."""
    _, _, Vt = np.linalg.svd(delta - delta.mean(0), full_matrices=False)
    return Vt[:k]


def procrustes(coef_a: np.ndarray, coef_b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Orthogonal ``R`` minimising ‖R·coefA − coefB‖, plus the coupling spectrum.

    Rows of ``coef_*`` are perturbations, so the cross-covariance is
    ``coef_bᵀ coef_a`` and the Procrustes solution is ``U Vᵀ`` of its SVD.
    """
    M = coef_b.T @ coef_a
    U, S, Vt = np.linalg.svd(M)
    return U @ Vt, S


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=K_DEFAULT)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=Path("results/connection.json"))
    args = ap.parse_args()

    lines, _, symbols = load_all()
    common = np.array(shared_perturbations(lines))
    names = [cl.name for cl in lines]
    print(f"{len(names)} base states, {common.size:,} shared knockdowns "
          f"(the coupling is exact), local fiber rank {args.k}")

    delta, coef, basis = {}, {}, {}
    for cl in lines:
        idx = {n: i for i, n in enumerate(cl.names)}
        D = cl.delta[[idx[p] for p in common]]
        delta[cl.name] = D
        basis[cl.name] = local_basis(D, args.k)
        coef[cl.name] = D @ basis[cl.name].T

    # ---- are the local frames even comparable? -------------------------------
    print(f"\nprincipal angles between local fiber subspaces (degrees, "
          f"median over {args.k})")
    angles = {}
    for a, b in itertools.combinations(names, 2):
        ang = np.degrees(subspace_angles(basis[a].T, basis[b].T))
        angles[f"{a}|{b}"] = ang.tolist()
        print(f"  {a:>7} vs {b:<7} median {np.median(ang):5.1f}  "
              f"min {ang.min():5.1f}  max {ang.max():5.1f}")
    rng = np.random.default_rng(args.seed)
    rand = local_basis(rng.normal(size=(common.size, symbols.size)), args.k)
    null_ang = np.degrees(subspace_angles(basis[names[0]].T, rand.T))
    print(f"  {'null':>7} (random subspace)  median {np.median(null_ang):5.1f}")

    # ---- the connection ------------------------------------------------------
    R, spec = {}, {}
    for a, b in itertools.permutations(names, 2):
        R[(a, b)], s = procrustes(coef[a], coef[b])
        spec[f"{a}->{b}"] = s.tolist()

    # ---- path consistency: does going via a third state agree with going direct?
    print(f"\npath inconsistency  ‖R_AC − R_BC·R_AB‖_F   "
          f"(0 = flat, {np.sqrt(2 * args.k):.1f} = unrelated)")
    incon = {}
    for a, b, c in itertools.permutations(names, 3):
        d = float(np.linalg.norm(R[(a, c)] - R[(b, c)] @ R[(a, b)]))
        incon[f"{a}->{b}->{c}"] = d
    for k_, v in sorted(incon.items(), key=lambda kv: -kv[1])[:4]:
        print(f"  {k_:<26}{v:6.2f}")
    print(f"  {'mean over all 24 paths':<26}{np.mean(list(incon.values())):6.2f}")

    # ---- holonomy around closed loops ---------------------------------------
    print(f"\nholonomy  ‖H − I‖_F  around closed loops "
          f"(0 = fiber returns unchanged)")
    hol = {}
    for loop in itertools.permutations(names, 3):
        a, b, c = loop
        H = R[(c, a)] @ R[(b, c)] @ R[(a, b)]
        hol["->".join(loop) + f"->{a}"] = float(
            np.linalg.norm(H - np.eye(args.k)))
    for k_, v in sorted(hol.items(), key=lambda kv: -kv[1])[:4]:
        print(f"  {k_:<26}{v:6.2f}")
    print(f"  {'mean over all 24 loops':<26}{np.mean(list(hol.values())):6.2f}")
    # Null: the same statistic for random orthogonal matrices of this size.
    from scipy.stats import special_ortho_group
    null_h = [float(np.linalg.norm(
        special_ortho_group.rvs(args.k, random_state=args.seed + i)
        @ special_ortho_group.rvs(args.k, random_state=args.seed + 100 + i)
        @ special_ortho_group.rvs(args.k, random_state=args.seed + 200 + i)
        - np.eye(args.k))) for i in range(24)]
    print(f"  {'null (random rotations)':<26}{np.mean(null_h):6.2f}")

    # ---- the only question that changes the model ---------------------------
    # Predict line B's effect for knockdown g from line A's, with and without
    # transport, on knockdowns held out of the Procrustes fit.
    # Transporting also *projects onto the target line's own subspace*, which
    # denoises regardless of any rotation. Two controls separate the two: the
    # same projection with no rotation at all, and the same projection with a
    # random rotation. Only the gap above those is the connection's doing.
    from scipy.stats import special_ortho_group
    print(f"\nDoes transporting help? Median cosine to the target line's true "
          f"effect,\n5-fold over knockdowns, Procrustes fit on training folds "
          f"only.")
    folds = np.array_split(rng.permutation(common.size), 5)
    ARMS = ("as-is", "projected", "random rot.", "transported")
    rows: dict[str, dict[str, float]] = {}
    for a, b in itertools.permutations(names, 2):
        acc = {k_: [] for k_ in ARMS}
        for fi, f in enumerate(folds):
            tr = np.setdiff1d(np.arange(common.size), f)
            Rf, _ = procrustes(coef[a][tr], coef[b][tr])
            Q = special_ortho_group.rvs(args.k, random_state=args.seed + fi)
            true = delta[b][f]
            cand = {
                "as-is": delta[a][f],
                "projected": coef[a][f] @ basis[b],
                "random rot.": (Q @ coef[a][f].T).T @ basis[b],
                "transported": (Rf @ coef[a][f].T).T @ basis[b],
            }
            for k_, P in cand.items():
                num = (P * true).sum(1)
                den = np.linalg.norm(P, axis=1) * np.linalg.norm(true, axis=1)
                acc[k_].append(np.median(num / np.maximum(den, 1e-12)))
        rows[f"{a}->{b}"] = {k_: float(np.mean(v)) for k_, v in acc.items()}

    print("  " + f"{'source -> target':<20}" + "".join(f"{k_:>13}" for k_ in ARMS)
          + f"{'gain vs proj':>14}")
    for k_, r in rows.items():
        print(f"  {k_:<20}" + "".join(f"{r[a_]:>13.4f}" for a_ in ARMS)
              + f"{r['transported'] - r['projected']:>+14.4f}")
    mean = {a_: np.mean([r[a_] for r in rows.values()]) for a_ in ARMS}
    print(f"  {'mean':<20}" + "".join(f"{mean[a_]:>13.4f}" for a_ in ARMS)
          + f"{mean['transported'] - mean['projected']:>+14.4f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"principal_angles": angles, "coupling_spectrum": spec,
         "path_inconsistency": incon, "holonomy": hol,
         "holonomy_null_mean": float(np.mean(null_h)),
         "transport_gain": rows, "k": args.k,
         "arms": list(ARMS)}, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
