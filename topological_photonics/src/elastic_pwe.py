"""
Full-vector 2D in-plane (P-SV) elastic plane-wave-expansion solver — the
phononic kernel for a ref-[4]-class case study (GHz elastic topological
waveguide of holes in silicon). Unlike the scalar-acoustic twin, this couples
the two displacement components u=(u_x,u_y) through the full elastic tensor.

Generalised eigenproblem (Kushwaha/Vasseur PWE for elastic phononic crystals):

   M_{αβ}(k) u = ω²  R  u ,
   M_{αβ}(G,G') = λ(G−G')(k+G)_α(k+G')_β
                + μ(G−G')[ (k+G)·(k+G') δ_{αβ} + (k+G)_β(k+G')_α ] ,
   R(G,G')      = ρ(G−G') δ_{αβ} .

Isotropic silicon matrix with cylindrical air holes (a honeycomb of two holes
per cell). Equal holes → Dirac cone at K; unequal holes (broken inversion) →
gap at K → valley-Hall phase (valley Chern ±1/2). Frequencies are reported as
f = ω a / (2π c_t), c_t = √(μ_Si/ρ_Si) the shear-wave speed.
"""
from __future__ import annotations
import numpy as np
from pwe import Lattice, eps_fourier

# --- silicon (isotropic approximation) and air-hole values (nondimensional) ---
RHO_SI, LAM_SI, MU_SI = 1.0, 1.28, 1.0          # units: rho_Si=1, mu_Si=1
CONTRAST = 1e-3                                   # air holes ~ vacuum
RHO_H, LAM_H, MU_H = (RHO_SI * CONTRAST, LAM_SI * CONTRAST, MU_SI * CONTRAST)


def honeycomb_lattice(a=1.0):
    a1 = np.array([a, 0.0])
    a2 = np.array([a / 2, a * np.sqrt(3) / 2])
    return Lattice(a1, a2)


def honeycomb_holes(rA, rB, a=1.0):
    lat = honeycomb_lattice(a)
    dA = (1 / 3) * lat.a1 + (1 / 3) * lat.a2
    dB = (2 / 3) * lat.a1 + (2 / 3) * lat.a2
    return [{"r": rA, "pos": dA}, {"r": rB, "pos": dB}], lat


class ElasticCrystal:
    def __init__(self, rA, rB, nmax=7, a=1.0):
        holes, self.lat = honeycomb_holes(rA, rB, a)
        self.G = self.lat.g_vectors(nmax)
        nG = len(self.G)
        Gdiff = (self.G[:, None, :] - self.G[None, :, :]).reshape(-1, 2)

        def field(val_hole, val_bg):
            rods = [{"r": h["r"], "pos": h["pos"], "eps": val_hole} for h in holes]
            return eps_fourier(Gdiff, rods, val_bg, self.lat.area).reshape(nG, nG)

        self.Lam = field(LAM_H, LAM_SI)
        self.Mu = field(MU_H, MU_SI)
        self.Rho = field(RHO_H, RHO_SI)
        self.nG = nG

    def _assemble(self, k):
        kG = k + self.G                      # (N,2)
        N = self.nG
        KK = kG @ kG.T                       # (k+G)·(k+G')
        kx, ky = kG[:, 0], kG[:, 1]
        oxx = np.outer(kx, kx); oxy = np.outer(kx, ky)
        oyx = np.outer(ky, kx); oyy = np.outer(ky, ky)
        Mxx = self.Lam * oxx + self.Mu * (KK + oxx)
        Myy = self.Lam * oyy + self.Mu * (KK + oyy)
        Mxy = self.Lam * oxy + self.Mu * oyx
        Myx = self.Lam * oyx + self.Mu * oxy
        M = np.zeros((2 * N, 2 * N), complex)
        M[0::2, 0::2] = Mxx; M[0::2, 1::2] = Mxy
        M[1::2, 0::2] = Myx; M[1::2, 1::2] = Myy
        R = np.zeros((2 * N, 2 * N), complex)
        R[0::2, 0::2] = self.Rho; R[1::2, 1::2] = self.Rho
        return 0.5 * (M + M.conj().T), 0.5 * (R + R.conj().T)

    def solve_k(self, k, nbands, return_vecs=False):
        from scipy.linalg import eigh
        M, R = self._assemble(np.asarray(k))
        w2, V = eigh(M, R)
        order = np.argsort(w2.real)
        w2 = np.clip(w2.real[order], 0, None)
        f = np.sqrt(w2) / (2 * np.pi)         # c_t=1 in these units
        if return_vecs:
            return f[:nbands], V[:, order[:nbands]], R
        return f[:nbands]

    def bands(self, kpath, nbands=8):
        return np.array([self.solve_k(k, nbands) for k in kpath])


def KMGamma_path(a=1.0, n=60):
    from pwe import high_symmetry_path
    lat = honeycomb_lattice(a)
    b1, b2 = lat.b1, lat.b2
    G0 = np.array([0.0, 0.0])
    K = (2 * b1 + b2) / 3.0
    M = 0.5 * b1
    pts = {"Γ": G0, "K": K, "M": M, "Γ ": G0}
    return high_symmetry_path(lat, pts, n)


if __name__ == "__main__":
    for tag, (rA, rB) in [("equal (Dirac)", (0.28, 0.28)),
                          ("broken (gap)", (0.34, 0.20))]:
        c = ElasticCrystal(rA, rB, nmax=7)
        path, ticks, keys = KMGamma_path(n=40)
        Kidx = ticks[1]
        fK = c.solve_k(path[Kidx], 8)
        print(f"{tag:16s} bands at K:", np.round(fK, 4))
