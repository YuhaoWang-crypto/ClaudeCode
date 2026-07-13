"""
Phononic (acoustic) twin of the Wu-Hu crystal — a sonic crystal on the SAME
hexagonal-cluster lattice, solved with a genuine acoustic PWE kernel.

Scalar pressure field p obeys   ∇·((1/ρ)∇p) + (ω²/B) p = 0 .
In a plane-wave basis this is the generalized eigenproblem

    Σ_G' (k+G)·(k+G') η(G−G') p_{G'} = ω² Σ_G' β(G−G') p_{G'} ,

with η = 1/ρ and β = 1/B. Rods have a density/modulus contrast with the
background (both η and β vary — the full acoustic problem, not a photonic
look-alike). The same geometry gives a double Dirac cone at R=a/3 and a
topological gap for R>a/3, providing the phononic system for the two-system
comparison of the complete package (cf. ref [4], GHz elastic waveguides).

Frequencies are reported as f = ωa/2πc0 with c0 the background sound speed.
"""
from __future__ import annotations
import numpy as np
from pwe import Lattice, eps_fourier
import geometry as geo

# nondimensional materials (background = air-like, c0=1): rods are denser & softer
RHO_BG, B_BG = 1.0, 1.0            # -> c0 = 1
RHO_ROD, B_ROD = 6.0, 2.0          # heavier, more compliant (slower, high contrast)


def _rods_with_value(cluster_R, value):
    rods = geo.wu_hu_rods(cluster_R)
    return [{"r": r["r"], "pos": r["pos"], "eps": value} for r in rods]


class AcousticCrystal:
    def __init__(self, cluster_R, nmax=9):
        self.lat = geo.wu_hu_lattice()
        self.G = self.lat.g_vectors(nmax)
        nG = len(self.G)
        Gdiff = (self.G[:, None, :] - self.G[None, :, :]).reshape(-1, 2)
        # eta = 1/rho , beta = 1/B  -> reuse rod-Fourier machinery
        eta = eps_fourier(Gdiff, _rods_with_value(cluster_R, 1 / RHO_ROD),
                          1 / RHO_BG, self.lat.area).reshape(nG, nG)
        beta = eps_fourier(Gdiff, _rods_with_value(cluster_R, 1 / B_ROD),
                           1 / B_BG, self.lat.area).reshape(nG, nG)
        self.eta = 0.5 * (eta + eta.conj().T)
        self.beta = 0.5 * (beta + beta.conj().T)

    def solve_k(self, k, nbands, return_vecs=False):
        kG = np.asarray(k) + self.G
        # A_{GG'} = (k+G).(k+G') eta_{G-G'}
        KK = kG @ kG.T
        A = KK * self.eta
        A = 0.5 * (A + A.conj().T)
        if return_vecs:
            from scipy.linalg import eigh
            w2, V = eigh(A, self.beta)
            order = np.argsort(w2.real)
            w2 = np.clip(w2.real[order], 0, None)
            return np.sqrt(w2[:nbands]) / (2 * np.pi), V[:, order[:nbands]]
        w2 = np.linalg.eigvals(np.linalg.solve(self.beta, A)).real
        w2 = np.sort(np.clip(w2, 0, None))
        return np.sqrt(w2[:nbands]) / (2 * np.pi)

    def bands(self, kpath, nbands=8):
        return np.array([self.solve_k(k, nbands) for k in kpath])


def gamma_gap(pc):
    f = pc.solve_k([0, 0], 6)
    return f[3] - f[2]
