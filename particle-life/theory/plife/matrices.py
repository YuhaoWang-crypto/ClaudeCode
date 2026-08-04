"""
Interaction-matrix generators.
==============================

Ports of the *named* matrix generators shipped by sandbox-science.com/particle-life
(ids 0-30, category "Experimental"), recovered verbatim from the site bundle, plus
a few analytically clean matrices used by the theory experiments.

Every generator returns an ``S x S`` array ``A`` with ``A[a, b] in [-1, 1]``:
the force felt by a particle of type ``a`` due to a particle of type ``b``.
``A[a,b] > 0`` = ``a`` is attracted to ``b``.

The *structure* of A is the genome of the universe.  What matters physically is
its decomposition

        A = A_S + A_N,   A_S = (A + A^T)/2,   A_N = (A - A^T)/2

``A_S`` (reciprocal) can be derived from a potential and can only build *static*
structures; ``A_N`` (non-reciprocal) breaks Newton's third law and is the sole
source of self-propulsion, rotation and travelling patterns.  See
``nonreciprocity()`` below and ``THEORY.md`` sections 3 and 6.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "GENERATORS",
    "generate",
    "random_radii",
    "symmetrise",
    "antisymmetrise",
    "nonreciprocity",
    "chase_pair",
]


# --------------------------------------------------------------------------- #
#  the site's generators (1:1 ports)
# --------------------------------------------------------------------------- #
def g_random(S, rng=None):
    """id 0 - Random.  A[a,b] ~ U(-1, 1), fully non-reciprocal on average."""
    rng = np.random.default_rng() if rng is None else rng
    return rng.uniform(-1.0, 1.0, size=(S, S))


def g_symmetric(S, rng=None):
    """id 1 - Symmetric.  Random, then mirrored: A = A^T.  Newton's third law holds."""
    rng = np.random.default_rng() if rng is None else rng
    A = rng.uniform(-1.0, 1.0, size=(S, S))
    iu = np.triu_indices(S, 1)
    A[(iu[1], iu[0])] = A[iu]
    return A


def g_snake(S, rng=None):
    """id 2 - Snake.  Self-attraction 1, weak forward link 0.2, nothing backwards.

    Maximally non-reciprocal cyclic motif: type i is pulled towards i+1 but i+1
    ignores i.  Produces crawling filaments.
    """
    A = np.zeros((S, S))
    for i in range(S):
        A[i, i] = 1.0
        A[i, (i + 1) % S] = 0.2
    return A


def g_chains1(S, rng=None):
    """id 3 - Chains 1.  Neighbours on the type-cycle attract (1), all else repel (-1)."""
    A = -np.ones((S, S))
    for i in range(S):
        for j in (i, (i + 1) % S, (i - 1) % S):
            A[i, j] = 1.0
    return A


def g_chains2(S, rng=None):
    """id 4 - Chains 2.  Self 1, cycle-neighbours 0.2, else -1."""
    A = -np.ones((S, S))
    for i in range(S):
        A[i, i] = 1.0
        A[i, (i + 1) % S] = 0.2
        A[i, (i - 1) % S] = 0.2
    return A


def g_chains3(S, rng=None):
    """id 5 - Chains 3.  Self 1, cycle-neighbours 0.2, else 0 (no long-range repulsion)."""
    A = np.zeros((S, S))
    for i in range(S):
        A[i, i] = 1.0
        A[i, (i + 1) % S] = 0.2
        A[i, (i - 1) % S] = 0.2
    return A


def g_rps(S, rng=None):
    """id 6 - Rock-Paper-Scissors.  Cyclic pursuit: +0.9 forward, -0.7 backward.

    The archetypal non-reciprocal matrix: A_{i,i+1} = +0.9 while A_{i+1,i} = -0.7,
    so every ordered pair is a chase.  Generates rotating clusters and travelling
    bands (net internal torque, see twobody.cluster_torque).
    """
    A = np.zeros((S, S))
    for i in range(S):
        A[i, i] = -0.1
        A[i, (i + 1) % S] = 0.9
        A[i, (i - 1) % S] = -0.7
    return A


def g_bipartite(S, rng=None):
    """id 7 - Bipartite Alliances.  Same parity attract (0.8), different repel (-0.8)."""
    A = np.empty((S, S))
    for i in range(S):
        for j in range(S):
            A[i, j] = 0.2 if i == j else (0.8 if (i % 2) == (j % 2) else -0.8)
    return A


def g_hub(S, rng=None):
    """id 8 - Hub and Spokes.  Type 0 is a hub: attracts everything, is attracted back."""
    A = np.zeros((S, S))
    for i in range(S):
        for j in range(S):
            if i == j:
                A[i, j] = 0.0
            elif i == 0:
                A[i, j] = 1.0
            elif j == 0:
                A[i, j] = 0.6
    return A


def g_shells(S, rng=None):
    """id 9 - Concentric Shells.  Self 0.9, next 0.3, everything else -0.6.

    Builds core-shell "cells": a condensed nucleus wrapped in a membrane of the
    next species.
    """
    A = np.full((S, S), -0.6)
    for i in range(S):
        A[i, i] = 0.9
        A[i, (i + 1) % S] = 0.3
    return A


def g_checker(S, rng=None):
    """id 11 - Checker Offsets.  Nearest type-neighbours repel, next-nearest attract."""
    A = np.zeros((S, S))
    for i in range(S):
        for j in range(S):
            if i == j:
                A[i, j] = 0.2
            elif j in ((i + 1) % S, (i - 1) % S):
                A[i, j] = -0.8
            elif j in ((i + 2) % S, (i - 2) % S):
                A[i, j] = 0.8
    return A


def g_dimers(S, rng=None):
    """id 12 - Dimers & Chains.  Partner = i XOR 1: strictly paired, everything else repels."""
    A = np.full((S, S), -0.9)
    for i in range(S):
        A[i, i] = 0.0
        A[i, i ^ 1 if (i ^ 1) < S else i] = 1.0
    return A


def g_triads(S, rng=None):
    """id 13 - Triad Flocks.  Types clustered into groups of three."""
    A = np.empty((S, S))
    grp = lambda i: i // 3
    for i in range(S):
        for j in range(S):
            A[i, j] = 0.1 if i == j else (0.9 if grp(i) == grp(j) else -0.7)
    return A


def g_spiral_conveyor(S, rng=None):
    """id 14 - Spiral Conveyor.  Asymmetric forward attraction 0.7 / 0.3, backward -0.6."""
    A = np.full((S, S), -0.6)
    for i in range(S):
        A[i, i] = -0.1
        A[i, (i + 1) % S] = 0.7
        A[i, (i + 2) % S] = 0.3
    return A


def g_patchwork(S, rng=None):
    """id 15 - Patchwork.  Deterministic hash -> sparse +-0.9 blocks."""
    A = np.zeros((S, S))
    for i in range(S):
        for j in range(S):
            if i == j:
                continue
            h = abs(np.sin(float((i * 73856093) ^ (j * 19349663))) * 43758.5453 % 1.0)
            A[i, j] = 0.9 if h < 0.35 else (-0.9 if h < 0.7 else 0.0)
    return A


def g_wavefield(S, rng=None):
    """id 16 - Wavefield.  A[i,j] = 0.9 sin(2 pi (j-i)/S).

    Purely antisymmetric off-diagonal (sin is odd) - the cleanest possible
    non-reciprocal matrix.  A_S ~ 0, so *all* of the dynamics is non-variational.
    """
    A = np.empty((S, S))
    for i in range(S):
        for j in range(S):
            A[i, j] = -0.05 if i == j else 0.9 * np.sin(2 * np.pi * (j - i) / S)
    return A


def g_chiral_bandpass(S, rng=None):
    """id 17 - Chiral Bandpass.  sin + 2nd harmonic + a sign-of-offset chirality term."""
    A = np.empty((S, S))
    for i in range(S):
        for j in range(S):
            if i == j:
                A[i, j] = -0.05
                continue
            d = (j - i) % S
            c = 2 * np.pi * d / S
            A[i, j] = 0.7 * np.sin(c) + 0.35 * np.sin(2 * c) + 0.15 * (1 if d > 0 else -1)
    return A


def g_rotating_conveyor(S, rng=None):
    """id 18 - Rotating Conveyor.  Wavefield with a per-row phase ramp (breaks cyclic symmetry)."""
    off = 2 * np.pi / (S * 1.5)
    A = np.empty((S, S))
    for i in range(S):
        for j in range(S):
            if i == j:
                A[i, j] = -0.05
                continue
            A[i, j] = 0.8 * np.sin(2 * np.pi * (j - i) / S + off * i)
    return A


def g_parity_vortex(S, rng=None):
    """id 20 - Parity Vortex.  even->odd 0.8, odd->even 0.6, like-parity repel."""
    A = np.empty((S, S))
    ev = lambda i: i % 2 == 0
    for i in range(S):
        for j in range(S):
            if i == j:
                A[i, j] = -0.05
            elif ev(i) and not ev(j):
                A[i, j] = 0.8
            elif not ev(i) and ev(j):
                A[i, j] = 0.6
            else:
                A[i, j] = -0.4 if ev(i) else -0.6
    return A


# --------------------------------------------------------------------------- #
#  analytically clean matrices used by the theory
# --------------------------------------------------------------------------- #
def g_gas(S, rng=None):
    """All-repulsive.  Uniform phase, used as the negative control."""
    return np.full((S, S), -0.6)


def g_condensate(S, rng=None):
    """Purely reciprocal self-attraction.  Static blobs; Lyapunov functional exists."""
    A = np.full((S, S), -0.3)
    np.fill_diagonal(A, 0.9)
    return A


def chase_pair(a12=0.9, a21=-0.9, a11=-0.1, a22=-0.1):
    """The minimal non-reciprocal universe: 2 species, 1 chases, 1 flees.

    ``a12 > 0 > a21``: species 0 is attracted to species 1, species 1 is repelled
    by species 0.  This 2x2 matrix is the building block for every "chaser" and
    the exactly solvable case in ``twobody.py`` / ``stability.exceptional_point``.
    """
    return np.array([[a11, a12], [a21, a22]], dtype=float)


GENERATORS = {
    "random": g_random,
    "symmetric": g_symmetric,
    "snake": g_snake,
    "chains1": g_chains1,
    "chains2": g_chains2,
    "chains3": g_chains3,
    "rps": g_rps,
    "bipartite": g_bipartite,
    "hub": g_hub,
    "shells": g_shells,
    "checker": g_checker,
    "dimers": g_dimers,
    "triads": g_triads,
    "spiral_conveyor": g_spiral_conveyor,
    "patchwork": g_patchwork,
    "wavefield": g_wavefield,
    "chiral_bandpass": g_chiral_bandpass,
    "rotating_conveyor": g_rotating_conveyor,
    "parity_vortex": g_parity_vortex,
    "gas": g_gas,
    "condensate": g_condensate,
}

# human-readable names, matching the site's preset menu where applicable
DISPLAY_NAMES = {
    "random": "Random",
    "symmetric": "Symmetric",
    "snake": "Snake",
    "chains1": "Chains 1",
    "chains2": "Chains 2",
    "chains3": "Chains 3",
    "rps": "Rock-Paper-Scissors",
    "bipartite": "Bipartite Alliances",
    "hub": "Hub and Spokes",
    "shells": "Concentric Shells",
    "checker": "Checker Offsets",
    "dimers": "Dimers & Chains",
    "triads": "Triad Flocks",
    "spiral_conveyor": "Spiral Conveyor",
    "patchwork": "Patchwork",
    "wavefield": "Wavefield",
    "chiral_bandpass": "Chiral Bandpass",
    "rotating_conveyor": "Rotating Conveyor",
    "parity_vortex": "Parity Vortex",
    "gas": "All-repulsive gas (control)",
    "condensate": "Reciprocal condensate (control)",
}


def generate(name, S, rng=None):
    """Build the named interaction matrix for ``S`` species."""
    if name not in GENERATORS:
        raise KeyError(f"unknown generator {name!r}; have {sorted(GENERATORS)}")
    A = np.asarray(GENERATORS[name](S, rng), dtype=float)
    return np.clip(np.round(A, 2), -1.0, 1.0)


# --------------------------------------------------------------------------- #
#  radii
# --------------------------------------------------------------------------- #
def random_radii(S, min_range=(12.0, 24.0), max_range=(32.0, 64.0), rng=None,
                 symmetric=False):
    """Per-ordered-pair min/max radii, sampled the way the site's randomiser does.

    Defaults are the platform's GPU-engine defaults ``minRadiusRange=[12,24]``,
    ``maxRadiusRange=[32,64]``.  ``beta`` is clamped above ``alpha`` so the tent
    never inverts.
    """
    rng = np.random.default_rng() if rng is None else rng
    alpha = rng.uniform(*min_range, size=(S, S))
    beta = rng.uniform(*max_range, size=(S, S))
    if symmetric:
        alpha = 0.5 * (alpha + alpha.T)
        beta = 0.5 * (beta + beta.T)
    beta = np.maximum(beta, alpha * 1.05 + 1e-6)
    return alpha, beta


def uniform_radii(S, alpha=15.0, beta=45.0):
    """Constant radii for every pair - the clean setting for analytic work."""
    return np.full((S, S), float(alpha)), np.full((S, S), float(beta))


# --------------------------------------------------------------------------- #
#  reciprocity diagnostics
# --------------------------------------------------------------------------- #
def symmetrise(A):
    """Reciprocal (variational) part ``A_S = (A + A^T)/2``."""
    A = np.asarray(A, dtype=float)
    return 0.5 * (A + A.T)


def antisymmetrise(A):
    """Non-reciprocal (non-variational) part ``A_N = (A - A^T)/2``."""
    A = np.asarray(A, dtype=float)
    return 0.5 * (A - A.T)


def nonreciprocity(A):
    r"""Scalar non-reciprocity index

    .. math::  \nu = \frac{\|A - A^{T}\|_F}{\|A + A^{T}\|_F + \|A-A^{T}\|_F} \in [0,1].

    ``nu = 0``: Newton's third law holds exactly -> gradient dynamics -> only static
    structures.  ``nu -> 1``: purely antisymmetric -> no energy function at all.
    """
    A = np.asarray(A, dtype=float)
    s = np.linalg.norm(A + A.T)
    n = np.linalg.norm(A - A.T)
    return float(n / (s + n)) if (s + n) > 0 else 0.0
