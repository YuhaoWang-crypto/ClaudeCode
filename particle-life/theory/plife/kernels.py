"""
Pair-interaction kernels of Particle Life.
=========================================

This module contains the *exact* force law used by sandbox-science.com/particle-life
(reverse-engineered from its WebGPU compute shader, see docs/REVERSE_ENGINEERING.md),
together with the analytic pair potential, its 2-D Hankel transform, and the
scalar diagnostics (bond length, well depth, integrated interaction) that the
continuum theory in ``stability.py`` needs.

Sign convention (identical to the shader)
-----------------------------------------
``f_ab(r) > 0``  means particle of type ``a`` is **attracted** towards a particle
of type ``b`` that sits a distance ``r`` away.  The force *vector* on particle i is

        F_i = kappa * sum_j  f_{s_i s_j}(r_ij) * rhat_ij ,
        rhat_ij = (x_j - x_i)/|x_j - x_i| .

Force law
---------
With ``alpha = minRadius[a][b]``, ``beta = maxRadius[a][b]``, ``A = force[a][b]``,
``R = repel`` (a single global scalar) and ``m = (alpha+beta)/2``, ``w = (beta-alpha)/2``:

        f(r) =  R * (r/alpha - 1)                      0 <= r < alpha     (core, always repulsive)
                A * (1 - |r - m| / w)                  alpha <= r < beta  (tent, sign = sign(A))
                0                                      r >= beta

Note two structural facts that drive *everything* downstream:

1.  ``f(alpha) = 0`` with ``f'(alpha^-) > 0`` and ``f'(alpha^+) = A/w``.  For ``A > 0``
    the distance ``r = alpha`` is therefore a **stable mechanical bond length**.
    The min-radius matrix is literally the bond-length matrix of the "chemistry".
2.  The force depends on the *ordered* pair (a, b).  Newton's third law is only
    satisfied if ``A`` (and alpha, beta) are symmetric.  The antisymmetric part is
    what makes Particle Life alive.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import j0


@lru_cache(maxsize=16)
def _leggauss(n):
    """Gauss-Legendre nodes/weights.  Cached: building them is O(n^2)."""
    return np.polynomial.legendre.leggauss(n)


__all__ = [
    "pair_force",
    "pair_potential",
    "hankel_transform",
    "hankel_matrix",
    "well_depth",
    "bond_length",
    "integrated_interaction",
]


# --------------------------------------------------------------------------- #
#  force
# --------------------------------------------------------------------------- #
def pair_force(r, A, alpha, beta, repel=1.0):
    """Scalar radial force f_ab(r).  Positive = attractive.

    All arguments broadcast; ``r`` may be an array of any shape.
    """
    r = np.asarray(r, dtype=float)
    A = np.asarray(A, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    beta = np.asarray(beta, dtype=float)

    m = 0.5 * (alpha + beta)
    w = 0.5 * (beta - alpha)

    core = r < alpha
    tent = (r >= alpha) & (r < beta)

    out = np.zeros(np.broadcast(r, A, alpha, beta).shape, dtype=float)
    # np.where keeps everything vectorised and avoids division warnings
    with np.errstate(divide="ignore", invalid="ignore"):
        f_core = repel * (r / alpha - 1.0)
        f_tent = A * (1.0 - np.abs(r - m) / w)
    out = np.where(core, f_core, np.where(tent, f_tent, 0.0))
    return out


# --------------------------------------------------------------------------- #
#  potential  (only meaningful for the reciprocal part, but always well defined)
# --------------------------------------------------------------------------- #
def pair_potential(r, A, alpha, beta, repel=1.0):
    r"""Pair potential ``U_ab(r)`` with ``U'(r) = f(r)`` and ``U(\infty) = 0``.

    Closed form (piecewise quadratic), derived by integrating ``pair_force``:

    .. math::
        U(r) = \begin{cases}
          \dfrac{R}{2\alpha}(r-\alpha)^2 - \dfrac{A(\beta-\alpha)}{2}, & 0\le r\le\alpha\\[2mm]
          -\dfrac{A}{2w}\bigl[2w^2-(r-\alpha)^2\bigr],                 & \alpha\le r\le m\\[2mm]
          -\dfrac{A}{2w}(\beta-r)^2,                                   & m\le r\le\beta\\[2mm]
          0,                                                           & r\ge\beta
        \end{cases}

    The well minimum sits exactly at ``r = alpha`` with depth ``-A(beta-alpha)/2``.
    """
    r = np.asarray(r, dtype=float)
    A = np.asarray(A, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    beta = np.asarray(beta, dtype=float)

    m = 0.5 * (alpha + beta)
    w = 0.5 * (beta - alpha)

    with np.errstate(divide="ignore", invalid="ignore"):
        u_core = repel / (2 * alpha) * (r - alpha) ** 2 - A * w
        u_in = -A / (2 * w) * (2 * w**2 - (r - alpha) ** 2)
        u_out = -A / (2 * w) * (beta - r) ** 2

    out = np.where(
        r < alpha,
        u_core,
        np.where(r < m, u_in, np.where(r < beta, u_out, 0.0)),
    )
    return out


def well_depth(A, alpha, beta):
    """Depth of the attractive well, ``U(alpha) = -A (beta-alpha)/2``.

    This is the natural *energy* scale of the model (for the reciprocal part):
    the binding energy of one bond.
    """
    return -np.asarray(A) * (np.asarray(beta) - np.asarray(alpha)) / 2.0


def bond_length(A, alpha, beta):
    """Stable pair separation.  ``alpha`` when ``A>0``, ``nan`` (unbound) otherwise."""
    A = np.asarray(A, dtype=float)
    return np.where(A > 0, np.asarray(alpha, dtype=float), np.nan)


# --------------------------------------------------------------------------- #
#  Hankel (2-D Fourier) transform of the kernel  -> continuum theory
# --------------------------------------------------------------------------- #
def hankel_transform(k, A, alpha, beta, repel=1.0, n_quad=4000):
    r"""2-D Fourier transform of the isotropic kernel ``U_ab``:

    .. math::  \hat U(k) = \int \mathrm{d}^2 r\, U(|r|) e^{-i k\cdot r}
                        = 2\pi\int_0^{\beta} r\,J_0(kr)\,U(r)\,\mathrm{d}r .

    Evaluated by Gauss-Legendre quadrature on ``[0, beta]`` (``U`` vanishes beyond).
    ``k`` may be an array; ``A, alpha, beta`` must be scalars.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    nodes, weights = _leggauss(n_quad)
    r = 0.5 * beta * (nodes + 1.0)          # map [-1,1] -> [0,beta]
    jac = 0.5 * beta
    u = pair_potential(r, A, alpha, beta, repel)
    integrand = (r * u)[None, :] * j0(np.outer(k, r))
    return 2 * np.pi * jac * (integrand * weights[None, :]).sum(axis=1)


def hankel_matrix(k, Amat, alpha_mat, beta_mat, repel=1.0, n_quad=3000):
    r"""``\hat U_{ab}(k)`` for every ordered species pair.

    Returns an array of shape ``(len(k), S, S)``.

    Uses one shared quadrature grid on ``[0, max(beta)]`` so the (expensive)
    Bessel matrix ``J_0(k_i r_j)`` is built once instead of ``S^2`` times.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    S = Amat.shape[0]
    rmax = float(np.max(beta_mat))

    nodes, weights = _leggauss(n_quad)
    r = 0.5 * rmax * (nodes + 1.0)
    jac = 0.5 * rmax
    B = j0(np.outer(k, r)) * (r * weights)[None, :] * (2 * np.pi * jac)   # (K, nq)

    Uall = pair_potential(
        r[None, None, :],
        Amat[:, :, None], alpha_mat[:, :, None], beta_mat[:, :, None], repel,
    )                                                                     # (S,S,nq)
    out = np.einsum("kq,abq->kab", B, Uall)
    return out


def hankel_basis(k, alpha_mat, beta_mat, n_quad=3000):
    r"""Decompose the kernel into its two *linear* channels.

    The pair potential is exactly linear in the two strength parameters::

        U(r; A, R, alpha, beta) = R * C(r; alpha) + A * T(r; alpha, beta)

    so the transform splits as :math:`\hat U_{ab}(k) = R\,\hat C_{ab}(k)
    + A_{ab}\,\hat T_{ab}(k)`.  Precomputing :math:`(\hat C,\hat T)` for a fixed
    pair of radius matrices lets an entire phase diagram over ``A`` be evaluated
    without touching a Bessel function again.

    Returns ``(Chat, That)``, each of shape ``(len(k), S, S)``.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    alpha_mat = np.asarray(alpha_mat, float)
    beta_mat = np.asarray(beta_mat, float)
    S = alpha_mat.shape[0]
    rmax = float(np.max(beta_mat))

    nodes, weights = _leggauss(n_quad)
    r = 0.5 * rmax * (nodes + 1.0)
    jac = 0.5 * rmax
    B = j0(np.outer(k, r)) * (r * weights)[None, :] * (2 * np.pi * jac)

    zeros = np.zeros((S, S))
    ones = np.ones((S, S))
    Ucore = pair_potential(r[None, None, :], zeros[:, :, None],
                           alpha_mat[:, :, None], beta_mat[:, :, None], 1.0)
    Uatt = pair_potential(r[None, None, :], ones[:, :, None],
                          alpha_mat[:, :, None], beta_mat[:, :, None], 0.0)
    return np.einsum("kq,abq->kab", B, Ucore), np.einsum("kq,abq->kab", B, Uatt)


def integrated_interaction(A, alpha, beta, repel=1.0):
    r"""``\hat U(0) = 2\pi \int_0^\infty r\,U(r)\,dr`` in closed form.

    Splitting into core and tent contributions:

    .. math::
        \hat U(0) = \underbrace{\frac{\pi R\,\alpha^3}{12}}_{\text{core (repulsive, }>0)}
                  \;-\;\underbrace{2\pi A\Bigl[\tfrac{\alpha(\beta-\alpha)}{2}\cdot\alpha
                  + \dots\Bigr]}_{\text{tent (}<0\text{ if }A>0)} .

    The sign of :math:`\hat U(0)` decides between *macroscopic* phase separation
    (:math:`\hat U(0)<0`, unbounded coarsening) and *microphase* separation
    (:math:`\hat U(0)>0` while :math:`\min_k\hat U(k)<0`: finite-size blobs).
    Computed here by quadrature to keep the piecewise algebra honest.
    """
    return hankel_transform(np.array([0.0]), A, alpha, beta, repel)[0]
