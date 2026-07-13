"""
Effective BHZ (Bernevig-Hughes-Zhang) model for the Wu-Hu photonic pseudospin.

The C6 band inversion diagnosed in ``symmetry.py`` maps, near Gamma, onto a
two-copy (Kramers pair) BHZ Hamiltonian in the (p+-, d+-) basis. This is the
standard effective description of the photonic quantum-spin-Hall crystal
(Wu & Hu 2015). We use it to compute the spin-Chern number (Berry curvature +
Fukui-Hatsugai-Suzuki) and to obtain the helical edge spectrum of a ribbon.

The mass M is calibrated FROM the PWE gap: M = +gap/2 in the topological
(expanded) phase, M = -gap/2 in the trivial (shrunken) phase, so the effective
model is tied to the real crystal rather than carrying arbitrary numbers.

    Lattice-regularised single spin block:
        d_x = A sin kx
        d_y = A sin ky
        d_z = M - 2B (2 - cos kx - cos ky)
        H_up(k) = d . sigma ;   H_dn(k) = [H_up(-k)]*
"""
from __future__ import annotations
import numpy as np
from numpy.linalg import eigh

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def h_up(kx, ky, M, A=1.0, B=1.0):
    dx = A * np.sin(kx)
    dy = A * np.sin(ky)
    dz = M - 2 * B * (2 - np.cos(kx) - np.cos(ky))
    return dx * SX + dy * SY + dz * SZ


def h_dn(kx, ky, M, A=1.0, B=1.0):
    return h_up(-kx, -ky, M, A, B).conj()


# ----------------------------------------------------------------------
# Berry curvature + Chern via Fukui-Hatsugai-Suzuki (per spin block)
# ----------------------------------------------------------------------
def chern_fhs(hfun, M, N=48, **kw):
    ks = np.linspace(-np.pi, np.pi, N, endpoint=False)
    vecs = np.empty((N, N, 2), complex)
    for i, kx in enumerate(ks):
        for j, ky in enumerate(ks):
            w, v = eigh(hfun(kx, ky, M, **kw))
            vecs[i, j] = v[:, 0]  # lower band
    F = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ip, jp = (i + 1) % N, (j + 1) % N
            U1 = np.vdot(vecs[i, j], vecs[ip, j])
            U2 = np.vdot(vecs[ip, j], vecs[ip, jp])
            U3 = np.vdot(vecs[ip, jp], vecs[i, jp])
            U4 = np.vdot(vecs[i, jp], vecs[i, j])
            prod = (U1 / abs(U1)) * (U2 / abs(U2)) * (U3 / abs(U3)) * (U4 / abs(U4))
            F[i, j] = np.angle(prod)
    return np.sum(F) / (2 * np.pi), F, ks


def spin_chern(M, N=48, **kw):
    C_up, F_up, ks = chern_fhs(h_up, M, N, **kw)
    C_dn, F_dn, _ = chern_fhs(h_dn, M, N, **kw)
    return {
        "C_up": C_up,
        "C_dn": C_dn,
        "C_total": C_up + C_dn,          # = 0 (time-reversal symmetric)
        "C_spin": (C_up - C_dn) / 2,     # spin-Chern: +-1 topological, 0 trivial
        "F_up": F_up,
        "ks": ks,
    }


# ----------------------------------------------------------------------
# Ribbon (strip) spectrum: helical edge states
# ----------------------------------------------------------------------
def ribbon_spectrum(kx, M, Ny=40, A=1.0, B=1.0, disorder=0.0, rng=None, both_spins=True):
    """
    4-band (or 2-band) ribbon: periodic in x (wavevector kx), Ny sites in y.
    Returns sorted eigenvalues. Optional on-site disorder amplitude.
    """
    blocks = [(+1,)] if not both_spins else [(+1,), (-1,)]
    all_E = []
    for spin in blocks:
        s = spin[0]
        # on-site (same y) 2x2 block in x-momentum space
        # H = dx(kx) sx + dy sy-part(from y hopping) + dz(kx, y-hopping) sz
        dim = 2 * Ny
        H = np.zeros((dim, dim), complex)
        for n in range(Ny):
            kxx = kx if s == +1 else -kx
            dx = A * np.sin(kxx)
            # dz(kx,ky) = M - 2B(2 - cos kx) + 2B cos ky ; cos ky -> y-hopping.
            dz0 = M - 2 * B * (2 - np.cos(kxx))  # ky-independent on-site part
            on = dx * SX + dz0 * SZ
            if s == -1:
                on = on.conj()
            H[2 * n:2 * n + 2, 2 * n:2 * n + 2] = on
            if n + 1 < Ny:
                # y-hopping: from -B cos ky and A sin ky terms
                # sin ky -> (T - T^dag)/2i on sy ; cos ky -> (T + T^dag)/2 on sz
                hop = (B) * SZ + (A / (2j)) * SY
                if s == -1:
                    hop = ((B) * SZ + (A / (2j)) * SY).conj()
                H[2 * n:2 * n + 2, 2 * (n + 1):2 * (n + 1) + 2] = hop
                H[2 * (n + 1):2 * (n + 1) + 2, 2 * n:2 * n + 2] = hop.conj().T
        if disorder > 0 and rng is not None:
            diag = rng.uniform(-disorder, disorder, dim)
            H += np.diag(diag.astype(complex))
        all_E.append(np.linalg.eigvalsh(H))
    return np.sort(np.concatenate(all_E))
