"""
Real-space (spin) Bott index — the topological invariant that REMAINS valid when
disorder breaks periodicity (addresses the critique that k-space Berry curvature
no longer applies under disorder).

We put the two BHZ Kramers blocks on an L×L torus, add Anderson on-site disorder,
and compute the Loring-Hastings Bott index of each block from the occupied-state
projector and the exponentiated position operators:

    B = (1/2π) Im Tr log( U_y U_x U_y† U_x† ) ,   U_{x,y} = Ψ† e^{i2π r_{x,y}/L} Ψ

with Ψ the occupied eigenvectors. The spin Bott index B_s = (B_↑ − B_↓)/2 equals
the clean-limit spin-Chern number but is defined WITHOUT a Brillouin zone, so it
survives disorder. Sweeping the disorder strength W shows the quantized plateau
and the disorder-driven topological transition (and, near the phase boundary, a
topological-Anderson-insulator 0→1 jump, as in the cited ref [3]).
"""
from __future__ import annotations
import numpy as np
from numpy.linalg import eigh
from scipy.linalg import logm

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def _block_hamiltonian(L, M, chir, W, rng, A=1.0, B=1.0):
    """One 2-band block on an L×L torus. chir=+1/-1 sets the chirality (spin)."""
    N = L * L
    dim = 2 * N
    H = np.zeros((dim, dim), complex)
    # Kramers partner: H_dn(k) = H_up*(-k) flips ONLY the sigma_x (d_x) term,
    # which flips the Chern sign (flipping both d_x and d_y would not).
    onsite = (M - 4 * B) * SZ
    hop_x = B * SZ + chir * (A / (2j)) * SX
    hop_y = B * SZ + (A / (2j)) * SY

    def blk(i):
        return slice(2 * i, 2 * i + 2)

    for x in range(L):
        for y in range(L):
            i = x * L + y
            dis = rng.uniform(-W / 2, W / 2)
            H[blk(i), blk(i)] = onsite + dis * np.eye(2)
            xn = ((x + 1) % L) * L + y
            yn = x * L + (y + 1) % L
            H[blk(i), blk(xn)] += hop_x
            H[blk(xn), blk(i)] += hop_x.conj().T
            H[blk(i), blk(yn)] += hop_y
            H[blk(yn), blk(i)] += hop_y.conj().T
    return H


def _bott(H, L):
    """Bott index of the lower (occupied) half of spectrum."""
    N = L * L
    w, V = eigh(H)
    Psi = V[:, :N]  # half filling: lowest N states
    xs = np.repeat(np.arange(L), L)
    ys = np.tile(np.arange(L), L)
    px = np.repeat(np.exp(2j * np.pi * xs / L), 2)  # per orbital
    py = np.repeat(np.exp(2j * np.pi * ys / L), 2)
    Ux = Psi.conj().T @ (px[:, None] * Psi)
    Uy = Psi.conj().T @ (py[:, None] * Psi)
    Mmat = Uy @ Ux @ Uy.conj().T @ Ux.conj().T
    return np.imag(np.trace(logm(Mmat))) / (2 * np.pi)


def spin_bott(L, M, W, rng):
    Bup = _bott(_block_hamiltonian(L, M, +1, W, rng), L)
    Bdn = _bott(_block_hamiltonian(L, M, -1, W, rng), L)
    return Bup, Bdn, (Bup - Bdn) / 2


def sweep_disorder(L=14, M=1.0, Ws=None, nreal=6, seed=0):
    if Ws is None:
        Ws = np.linspace(0, 8, 17)
    rng = np.random.default_rng(seed)
    mean, std = [], []
    for W in Ws:
        vals = [spin_bott(L, M, W, rng)[2] for _ in range(nreal)]
        mean.append(np.mean(vals)); std.append(np.std(vals))
    return Ws, np.array(mean), np.array(std)


if __name__ == "__main__":
    r = np.random.default_rng(0)
    for M, tag in [(1.0, "topological"), (-1.0, "trivial")]:
        bu, bd, bs = spin_bott(12, M, 0.0, r)
        print(f"clean M={M:+.1f} ({tag}): B_up={bu:+.2f} B_dn={bd:+.2f} spin-Bott={bs:+.2f}")
