"""
Symmetry-indicator diagnosis of band inversion for the Wu-Hu crystal.

At Gamma the modes organise into C6 angular-momentum multiplets. The two
doublets straddling the gap are p-like (|l|=1) and d-like (|l|=2). Their
ordering is the RIGOROUS topological fingerprint, computed directly from the
PWE eigenvectors of the real photonic crystal (no effective model needed):

    trivial      :  p (|l|=1) below the gap, d (|l|=2) above
    topological  :  d (|l|=2) below the gap, p (|l|=1) above   <- INVERTED
"""
from __future__ import annotations
import numpy as np


def _c6_permutation(G):
    """Index map i -> j such that R60 . G_i == G_j (or -1 if rotated out of set)."""
    c, s = np.cos(np.pi / 3), np.sin(np.pi / 3)
    R = np.array([[c, -s], [s, c]])
    RG = G @ R.T
    perm = np.full(len(G), -1, int)
    for i, g in enumerate(RG):
        d = np.linalg.norm(G - g, axis=1)
        j = int(np.argmin(d))
        if d[j] < 1e-6:
            perm[i] = j
    return perm


def c6_angular_momentum(pc, doublet_bands):
    """Return the |l| values of a degenerate doublet at Gamma (e.g. [1,1] or [2,2])."""
    _, V = pc.solve_k([0, 0], max(doublet_bands) + 1, return_vecs=True)
    perm = _c6_permutation(pc.G)

    def C6(vec):
        out = np.zeros_like(vec)
        m = perm >= 0
        out[perm[m]] = vec[m]
        return out

    sub = V[:, doublet_bands]
    M = sub.conj().T @ np.column_stack([C6(sub[:, k]) for k in range(sub.shape[1])])
    ev = np.linalg.eigvals(M)
    ls = np.angle(ev) / (np.pi / 3)  # eigenvalue = exp(-i l pi/3)
    return np.sort(np.abs(np.round(ls)))


def diagnose(pc, lower=(1, 2), upper=(3, 4)):
    """Classify the crystal as 'topological' or 'trivial' from the doublet order."""
    l_low = c6_angular_momentum(pc, list(lower))
    l_high = c6_angular_momentum(pc, list(upper))
    inverted = np.mean(l_low) > np.mean(l_high)  # d (2) below p (1) => inverted
    return {
        "l_lower": l_low,
        "l_upper": l_high,
        "phase": "topological" if inverted else "trivial",
        "inverted": bool(inverted),
    }
