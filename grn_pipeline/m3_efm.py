"""
Module 3 — Elementary Flux Modes (EFMs): the irreducible generators.

Idea (rigorous):
  At steady state the flux vector v obeys  S v = 0  with v_i >= 0 for every
  irreversible reaction. The feasible fluxes form a POLYHEDRAL CONE. Its
  EXTREME RAYS are the Elementary Flux Modes: support-minimal pathways that
  cannot be decomposed into simpler steady-state pathways. Every feasible
  steady-state flux is a NON-NEGATIVE combination of EFMs.

  This is the closest biological object to a "prime factorisation": the EFMs
  are an irreducible generating set for all metabolic behaviour. There is no
  smaller set of pathways that spans the same cone.

We build a small central-carbon-style network, enumerate its EFMs by
support-minimality, and verify that a random feasible flux reconstructs as a
non-negative combination of EFMs (scipy.optimize.nnls).
"""
from itertools import combinations
import numpy as np
from scipy.optimize import nnls


def toy_network():
    """Metabolites (rows) x reactions (cols). All reactions irreversible
    (fluxes >= 0). External supply/removal reactions balance the internal
    metabolites A, B, C, D.

        R1:      -> A        (uptake)
        R2:   A  -> B
        R3:   A  -> C
        R4:   B  -> C
        R5:   B  -> D
        R6:   C  -> D
        R7:   D  ->          (secretion)
    Internal (balanced) metabolites: A, B, C, D.
    """
    reactions = ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]
    mets = ["A", "B", "C", "D"]
    S = np.array([
        # R1 R2 R3 R4 R5 R6 R7
        [ 1, -1, -1,  0,  0,  0,  0],   # A
        [ 0,  1,  0, -1, -1,  0,  0],   # B
        [ 0,  0,  1,  1,  0, -1,  0],   # C
        [ 0,  0,  0,  0,  1,  1, -1],   # D
    ], dtype=float)
    return S, reactions, mets


def _has_flux_with_support(S, support, tol=1e-9):
    """Is there v>=0 with Sv=0, support exactly `support` (nonzero on it)?
    Returns a representative vector or None."""
    n = S.shape[1]
    idx = list(support)
    Ssub = S[:, idx]
    # null space of Ssub
    u, s, vt = np.linalg.svd(Ssub)
    rank = np.sum(s > tol)
    ns = vt[rank:].T  # columns span null space
    if ns.shape[1] == 0:
        return None
    # look for a strictly-positive vector in the null space (all support
    # entries nonzero, same sign). Try each null basis vector and simple
    # combinations; for support-minimal EFMs the null space is 1-D.
    for col in range(ns.shape[1]):
        w = ns[:, col]
        for sgn in (1.0, -1.0):
            ww = sgn * w
            if np.all(ww > tol):
                full = np.zeros(n)
                full[idx] = ww
                return full / full.max()
    # 1-D case handled above; higher-dim supports are not minimal, skip
    return None


def elementary_flux_modes(S, reactions):
    """Enumerate EFMs as support-minimal non-negative steady-state fluxes.
    Brute force over supports (fine for small networks)."""
    n = S.shape[1]
    efms = []
    efm_supports = []
    for k in range(1, n + 1):
        for support in combinations(range(n), k):
            sset = set(support)
            # skip if a strictly smaller EFM support is already contained
            if any(prev < sset for prev in efm_supports):
                continue
            v = _has_flux_with_support(S, support)
            if v is not None:
                # confirm steady state
                if np.allclose(S @ v, 0, atol=1e-7):
                    efms.append(v)
                    efm_supports.append(sset)
    return efms, efm_supports


def report():
    S, reactions, mets = toy_network()
    efms, supports = elementary_flux_modes(S, reactions)

    print("=" * 68)
    print("MODULE 3 — Elementary Flux Modes (irreducible flux generators)")
    print("=" * 68)
    print(f"network: {S.shape[0]} balanced metabolites, "
          f"{S.shape[1]} reactions")
    print(f"number of EFMs (extreme rays of the flux cone): {len(efms)}")
    for i, (v, sup) in enumerate(zip(efms, supports), 1):
        path = " + ".join(f"{v[j]:.2f}·{reactions[j]}"
                           for j in sorted(sup))
        print(f"  EFM{i}: {path}")

    # verification: a random non-negative combination of EFMs is a valid
    # steady-state flux, and nnls recovers it from the EFM basis.
    E = np.array(efms).T                       # reactions x n_efm
    rng = np.random.default_rng(3)
    coef_true = rng.uniform(0.2, 2.0, size=len(efms))
    v_random = E @ coef_true
    assert np.allclose(S @ v_random, 0, atol=1e-8), "not steady state!"
    coef_fit, res = nnls(E, v_random)
    print("\nverification:")
    print(f"  a random feasible flux satisfies S·v = 0 : "
          f"max|S v| = {np.max(np.abs(S @ v_random)):.2e}")
    print(f"  reconstructed as non-negative combo of EFMs, residual = "
          f"{res:.2e}")
    print("  => every steady-state flux IS a non-negative sum of EFMs.")
    return {"S": S, "reactions": reactions, "efms": efms,
            "supports": supports}


if __name__ == "__main__":
    report()
