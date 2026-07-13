"""
Wu-Hu topological photonic crystal geometry (Wu & Hu, PRL 114, 223901, 2015).

A triangular Bravais lattice with a hexagonal cluster of 6 dielectric rods per
unit cell. The single tuning parameter is R = cluster radius (distance of each
rod from the cell centre), in units of the lattice constant a:

    R = a/3   ->  the 6 rods complete a perfect honeycomb lattice; the two
                  graphene Dirac cones at K/K' fold onto Gamma and merge into a
                  DOUBLE DIRAC CONE (four-fold degeneracy at Gamma).
    R < a/3   ->  "shrunken" cluster -> trivial gap opens at Gamma.
    R > a/3   ->  "expanded" cluster -> band-inverted (topological) gap: the
                  p-like (l=+-1) and d-like (l=+-2) doublets swap order.

This band inversion is the photonic analogue of the quantum spin Hall effect;
the pseudospin comes from the C6 rotational symmetry (p+-, d+- orbitals).
"""
from __future__ import annotations
import numpy as np
from pwe import Lattice, PhotonicCrystal2D, high_symmetry_path

A = 1.0  # lattice constant (all lengths in units of a)
EPS_ROD = 11.7  # silicon rods (matches ref [4] unsuspended-Si platform)
EPS_BG = 1.0  # air
R_ROD = 0.12 * A  # individual rod radius


def wu_hu_lattice() -> Lattice:
    a1 = np.array([A, 0.0])
    a2 = np.array([A / 2, A * np.sqrt(3) / 2])
    return Lattice(a1, a2)


def wu_hu_rods(cluster_R, centre=(0.0, 0.0), r_rod=R_ROD, eps_rod=EPS_ROD):
    """6 rods on a hexagon of radius `cluster_R` about `centre`."""
    cx, cy = centre
    rods = []
    for k in range(6):
        th = np.pi / 3 * k
        rods.append(
            {
                "r": r_rod,
                "pos": (cx + cluster_R * np.cos(th), cy + cluster_R * np.sin(th)),
                "eps": eps_rod,
            }
        )
    return rods


def wu_hu_crystal(cluster_R, nmax=9) -> PhotonicCrystal2D:
    lat = wu_hu_lattice()
    return PhotonicCrystal2D(lat, wu_hu_rods(cluster_R), eps_bg=EPS_BG, nmax=nmax)


def gamma_M_K_path(n_per_seg=60):
    """Gamma - M - K - Gamma path for the triangular lattice."""
    lat = wu_hu_lattice()
    b1, b2 = lat.b1, lat.b2
    G = np.array([0.0, 0.0])
    M = 0.5 * b1
    K = (2 * b1 + b2) / 3.0
    pts = {"Γ": G, "M": M, "K": K, "Γ ": G}
    return high_symmetry_path(lat, pts, n_per_seg)


# Critical / representative cluster radii
R_CRITICAL = A / 3.0
R_SHRUNK = 0.30 * A   # trivial
R_EXPAND = 0.36 * A   # topological
