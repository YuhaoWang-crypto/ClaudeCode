"""Pure-numpy structural geometry: Kabsch superposition + RMSD."""

from __future__ import annotations

import numpy as np


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Minimal RMSD between two equal-length C-alpha traces after optimal
    rigid-body superposition (Kabsch algorithm).

    Parameters
    ----------
    P, Q : (N, 3) arrays of coordinates in the same residue order.

    Returns
    -------
    float : root-mean-square deviation in the same length units as the inputs.
    """
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    if P.shape != Q.shape:
        raise ValueError(f"shape mismatch {P.shape} vs {Q.shape}")
    if P.ndim != 2 or P.shape[1] != 3:
        raise ValueError("inputs must be (N, 3)")

    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)

    # Covariance and optimal rotation via SVD, with a reflection correction.
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T

    P_aligned = Pc @ R.T
    diff = P_aligned - Qc
    return float(np.sqrt((diff * diff).sum() / P.shape[0]))


def self_clash_count(coords: np.ndarray, contact: float = 3.0) -> int:
    """Count non-sequential C-alpha pairs closer than ``contact`` angstrom.

    A crude steric-clash proxy standing in for a real all-atom energy term.
    Neighbours |i-j| <= 1 are ignored (chain connectivity).
    """
    coords = np.asarray(coords, dtype=float)
    n = coords.shape[0]
    if n < 3:
        return 0
    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=-1))
    iu = np.triu_indices(n, k=2)          # skip self and immediate neighbour
    return int((dist[iu] < contact).sum())
