"""
Plane-Wave Expansion (PWE) engine for 2D photonic crystals — TM polarization.

Solves the Maxwell master equation for E_z (TM modes) of a 2D photonic crystal:

        |k+G|^2 c_G  =  (omega/c)^2  sum_G'  eps(G-G') c_G'

as a generalized eigenvalue problem in a plane-wave basis. Frequencies are
returned in the dimensionless unit  omega*a / (2*pi*c)  (a = lattice constant).

The crystal is defined by a Bravais lattice (a1, a2) and a list of circular
dielectric rods (radius, centre, permittivity) inside a background medium.
This is a general, reusable solver — the Wu-Hu topological crystal in
``geometry.py`` is just one particular rod list handed to it.

No proprietary code: pure NumPy/SciPy. Validated against the analytic empty-
lattice (eps=const) bands and against the known honeycomb Dirac cone.
"""
from __future__ import annotations
import numpy as np
from scipy.special import j1


class Lattice:
    """A 2D Bravais lattice with its reciprocal vectors and cell area."""

    def __init__(self, a1, a2):
        self.a1 = np.asarray(a1, float)
        self.a2 = np.asarray(a2, float)
        self.area = abs(self.a1[0] * self.a2[1] - self.a1[1] * self.a2[0])
        # 2D reciprocal vectors (b_i . a_j = 2*pi delta_ij)
        self.b1 = 2 * np.pi * np.array([self.a2[1], -self.a2[0]]) / self.area
        self.b2 = 2 * np.pi * np.array([-self.a1[1], self.a1[0]]) / self.area

    def frac_to_cart(self, u, v):
        return u * self.a1 + v * self.a2

    def g_vectors(self, nmax):
        """All reciprocal vectors m*b1 + n*b2 with |m|,|n| <= nmax."""
        idx = range(-nmax, nmax + 1)
        return np.array([m * self.b1 + n * self.b2 for m in idx for n in idx])


def _rod_form_factor(Gnorm, radius, area_cell):
    """Fourier transform of a filled cylinder (area fraction at G=0)."""
    area = np.pi * radius * radius
    out = np.empty_like(Gnorm)
    small = Gnorm < 1e-9
    out[small] = area / area_cell
    g = Gnorm[~small]
    out[~small] = (2 * area / area_cell) * j1(g * radius) / (g * radius)
    return out


def eps_fourier(Gdiff, rods, eps_bg, area_cell):
    """
    eps(G) for a set of rods on a background.

    rods: list of dicts {'r': radius, 'pos': (x,y), 'eps': eps_rod}
    Returns complex array len(Gdiff).
    """
    Gnorm = np.linalg.norm(Gdiff, axis=1)
    out = np.zeros(len(Gdiff), dtype=complex)
    for rod in rods:
        f = _rod_form_factor(Gnorm, rod["r"], area_cell)
        phase = np.exp(-1j * (Gdiff @ np.asarray(rod["pos"])))
        out += (rod["eps"] - eps_bg) * f * phase
    out[Gnorm < 1e-9] += eps_bg
    return out


class PhotonicCrystal2D:
    """2D photonic crystal; TM (E_z) band solver via plane-wave expansion."""

    def __init__(self, lattice: Lattice, rods, eps_bg=1.0, nmax=7):
        self.lat = lattice
        self.rods = rods
        self.eps_bg = eps_bg
        self.nmax = nmax
        self.G = lattice.g_vectors(nmax)
        # Precompute the eps(G-G') matrix once (k-independent).
        nG = len(self.G)
        Gdiff = (self.G[:, None, :] - self.G[None, :, :]).reshape(-1, 2)
        B = eps_fourier(Gdiff, rods, eps_bg, lattice.area).reshape(nG, nG)
        self.B = 0.5 * (B + B.conj().T)  # hermitize numerical noise

    def solve_k(self, k, nbands, return_vecs=False):
        """Return the lowest `nbands` frequencies (omega a/2pi c) at wavevector k."""
        kG = np.asarray(k) + self.G
        A = np.diag(np.sum(kG * kG, axis=1)).astype(complex)
        if return_vecs:
            # generalized hermitian eig:  A c = w2 B c
            from scipy.linalg import eigh
            w2, vecs = eigh(A, self.B)
            order = np.argsort(w2.real)
            w2 = np.clip(w2.real[order], 0, None)
            freqs = np.sqrt(w2) / (2 * np.pi)
            return freqs[:nbands], vecs[:, order[:nbands]]
        w2 = np.linalg.eigvals(np.linalg.solve(self.B, A)).real
        w2 = np.sort(np.clip(w2, 0, None))
        return np.sqrt(w2[:nbands]) / (2 * np.pi)

    def bands(self, kpath, nbands=8):
        return np.array([self.solve_k(k, nbands) for k in kpath])

    def field_Ez(self, k, band, ngrid=80, return_vecs=False):
        """Real-space |E_z| of a given band at k (for mode-symmetry inspection)."""
        _, vecs = self.solve_k(k, band + 1, return_vecs=True)
        c = vecs[:, band]
        xs = np.linspace(0, 1, ngrid)
        Ez = np.zeros((ngrid, ngrid), complex)
        for i, u in enumerate(xs):
            for j, v in enumerate(xs):
                r = self.lat.frac_to_cart(u, v)
                Ez[i, j] = np.sum(c * np.exp(1j * ((np.asarray(k) + self.G) @ r)))
        return Ez


def high_symmetry_path(lat: Lattice, points, n_per_seg=50):
    """Build a k-path through named high-symmetry points; return path + ticks."""
    path, ticks = [], [0]
    keys = list(points.keys())
    for p0, p1 in zip(keys[:-1], keys[1:]):
        A, Bp = points[p0], points[p1]
        for t in np.linspace(0, 1, n_per_seg, endpoint=False):
            path.append(A + t * (Bp - A))
        ticks.append(len(path) + 1)
    path.append(points[keys[-1]])
    return np.array(path), ticks, keys
