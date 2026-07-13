"""
Finite-Difference Frequency-Domain (FDFD) solver for the REAL rod crystal.

This replaces the effective BHZ ribbon of ``edge_states.py`` with an ab-initio
supercell calculation of the actual silicon-rod Wu-Hu crystal, so the edge
states are computed directly from Maxwell's equations (removes the ⚠️ label).

TM modes (E_z) on a rectangular Yee-less grid solve the generalized eigenproblem

        -∇² E_z = (ω/c)²  ε(r)  E_z

with a 5-point Laplacian, a Bloch phase exp(i k_x L_x) across the x-period, and
periodic wrapping in y. The triangular Wu-Hu lattice is represented on an
orthorhombic supercell of size a × √3·a (two clusters per cell); a strip is
built by stacking such cells in y with a topological | trivial domain wall, so
two helical edge channels appear inside the bulk gap.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

A = 1.0
SQ3 = np.sqrt(3.0)
EPS_ROD, EPS_BG, R_ROD = 11.7, 1.0, 0.12 * A


def _cluster_rods(cx, cy, R):
    return [(cx + R * np.cos(np.pi / 3 * k), cy + R * np.sin(np.pi / 3 * k))
            for k in range(6)]


def build_epsilon(Ncells_y, R_top, R_bot, Nx=54):
    """
    Orthorhombic strip: width Lx=a (periodic, Bloch in x), height Ly=Ncells_y*√3 a
    (periodic in y). Top half clusters use R_top, bottom half R_bot -> domain wall.
    Returns eps grid (Nx, Ny), spacing h, Lx, Ly.
    """
    h = A / Nx
    cell_h = SQ3 * A
    Ny_cell = int(round(cell_h / h))
    Ny = Ny_cell * Ncells_y
    Ly = Ny * h
    Lx = Nx * h
    xs = (np.arange(Nx) + 0.5) * h
    ys = (np.arange(Ny) + 0.5) * h
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    eps = np.full((Nx, Ny), EPS_BG)

    # lattice points of the triangular lattice inside each orthorhombic cell:
    # (0,0) and (a/2, √3 a/2) relative to the cell origin.
    for c in range(Ncells_y):
        y0 = c * cell_h
        R = R_top if c >= Ncells_y // 2 else R_bot
        centres = [(0.0, y0), (A / 2, y0 + cell_h / 2),
                   (A, y0), (0.0, y0 + cell_h), (A, y0 + cell_h)]  # +images at x=a
        rods = []
        for (cx, cy) in centres:
            rods += _cluster_rods(cx, cy, R)
        for (rx, ry) in rods:
            # nearest-image in x (period a) and y (period Ly)
            dx = X - rx
            dx -= Lx * np.round(dx / Lx)
            dy = Y - ry
            dy -= Ly * np.round(dy / Ly)
            eps[dx * dx + dy * dy < R_ROD * R_ROD] = EPS_ROD
    return eps, h, Lx, Ly


def _laplacian(Nx, Ny, h, kx, Lx):
    """-∇² with Bloch phase in x (period Lx) and periodic y. Complex Hermitian."""
    N = Nx * Ny
    idx = lambda i, j: (i % Nx) * Ny + (j % Ny)
    inv_h2 = 1.0 / (h * h)
    rows, cols, vals = [], [], []
    phase = np.exp(1j * kx * Lx)
    for i in range(Nx):
        for j in range(Ny):
            p = idx(i, j)
            rows.append(p); cols.append(p); vals.append(4.0 * inv_h2)
            # +x neighbour (wrap carries Bloch phase)
            ph = phase if i == Nx - 1 else 1.0
            rows.append(p); cols.append(idx(i + 1, j)); vals.append(-inv_h2 * ph)
            ph = np.conj(phase) if i == 0 else 1.0
            rows.append(p); cols.append(idx(i - 1, j)); vals.append(-inv_h2 * ph)
            rows.append(p); cols.append(idx(i, j + 1)); vals.append(-inv_h2)
            rows.append(p); cols.append(idx(i, j - 1)); vals.append(-inv_h2)
    L = sp.csr_matrix((vals, (rows, cols)), shape=(N, N), dtype=complex)
    return L


def strip_spectrum(kx, eps, h, Lx, nstates=40, sigma=None):
    """Lowest `nstates` frequencies (ωa/2πc) at wavevector kx for the strip."""
    Nx, Ny = eps.shape
    L = _laplacian(Nx, Ny, h, kx, Lx)
    Eps = sp.diags(eps.ravel())
    if sigma is None:
        sigma = (2 * np.pi * 0.46) ** 2  # near the gap
    w2, vecs = eigsh(L, k=nstates, M=Eps, sigma=sigma, which="LM")
    order = np.argsort(w2.real)
    w2 = np.clip(w2.real[order], 0, None)
    freqs = np.sqrt(w2) / (2 * np.pi)
    return freqs, vecs[:, order]


def edge_localization(vec, Nx, Ny):
    """Fraction of |E|² in the central quarter of y (domain-wall region)."""
    p = np.abs(vec.reshape(Nx, Ny)) ** 2
    py = p.sum(axis=0)
    py /= py.sum()
    q = Ny // 4
    return py[Ny // 2 - q // 2: Ny // 2 + q // 2].sum()
