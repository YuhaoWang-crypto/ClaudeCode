r"""
Continuum theory of Particle Life: linear stability of the uniform state.
=========================================================================

Coarse-graining
---------------
Let :math:`\rho_a(x,t)` be the density of species *a*.  Mass conservation plus the
momentum balance of the microscopic model (``model.py``) give, to linear order in
the fluctuations,

.. math::
    \partial_t \rho_a + \nabla\!\cdot(\rho_a u_a) = 0, \qquad
    \partial_t u_a = -60\kappa\,\nabla\Phi_a - \gamma u_a, \qquad
    \Phi_a = \sum_b U_{ab}\star\rho_b .

The convolution kernel is the pair potential of ``kernels.pair_potential`` -- note
that *no* free energy exists unless :math:`U_{ab}=U_{ba}`, yet the *linear operator*
is perfectly well defined either way.  That asymmetry is the whole story.

Dispersion relation
-------------------
Writing :math:`\delta\rho_a \propto e^{i k\cdot x + \lambda t}` and letting
:math:`\sigma_n(k)` be the eigenvalues of the (generally non-symmetric) matrix

.. math::  \mathcal{M}(k) = \operatorname{diag}(\bar\rho)\,\hat U(k),

each branch obeys the quadratic

.. math::  \lambda^2 + \gamma\lambda + 60\kappa k^2 \sigma_n(k) = 0
           \;\Longrightarrow\;
           \lambda_\pm = \tfrac12\Bigl(-\gamma \pm \sqrt{\gamma^2 - 240\kappa k^2\sigma_n}\Bigr).

Three exact consequences
------------------------
**(i) Reciprocal matrices can only make static patterns.**
If :math:`\hat U = \hat U^{T}` and the mean densities are equal, :math:`\mathcal M`
is symmetric, so every :math:`\sigma_n` is real.  Then
:math:`\operatorname{Re}\lambda_+>0` requires :math:`\sigma_n<0`, which makes the
square root real: the growing mode has :math:`\operatorname{Im}\lambda=0`.
Growing modes never oscillate -> **stationary structures only**.

**(ii) Non-reciprocity is the sole route to travelling structures.**
If :math:`\sigma = \sigma_r + i\sigma_i` with :math:`\sigma_i \neq 0`, then
:math:`\operatorname{Re}\lambda_+>0` iff

.. math::  \boxed{\;\sigma_i^{\,2} \;>\; \frac{\gamma^{2}}{60\,\kappa\,k^{2}}\;\sigma_r \;}

(automatically satisfied when :math:`\sigma_r<0`).  When :math:`\sigma_r>0` --
i.e. when the *reciprocal* part alone would be stable -- the instability exists
**only because of non-reciprocity**, and the unstable mode necessarily has
:math:`\operatorname{Im}\lambda\neq0`: it is a travelling wave with phase speed
:math:`c=|\operatorname{Im}\lambda|/k`.  That is the linear-theory fingerprint of
a "chaser", a crawling snake, a rotating cell.

**(iii) Structure size is set by** :math:`k^\*=\arg\max_k\operatorname{Re}\lambda(k)`,
and the maximum is at a **finite** :math:`k` for a Cahn-Hilliard reason: the
dynamics conserves particle number, which puts an explicit :math:`k^2` in front of
:math:`\sigma`, killing long wavelengths, while the repulsive core makes
:math:`\sigma(k)>0` at large :math:`k`, killing short ones.  The competition
selects :math:`\ell^\*=2\pi/k^\*`.

Note what does *not* happen: for a **single** species the kernel
:math:`\hat U(k)` never has an interior minimum -- over 282 random single-species
parameter sets with :math:`\hat U(0)<0`, the minimum sat at :math:`k=0` every
single time.  A lone species can therefore only coarsen without bound.  Finite,
repeating structures always come from **several species competing**: short-range
self-attraction plus longer-range cross-repulsion is a SALR (short-attraction /
long-repulsion) kernel realised in *species space*, and eliminating the fast
species gives the mediated effective interaction

.. math::  \hat U^{\rm eff}_{aa}(k)=\hat U_{aa}(k)
           -\frac{\hat U_{ab}(k)\hat U_{ba}(k)}{\hat U_{bb}(k)},

whose interior minimum is what fixes the cell size.  ``ell*`` scales linearly with
``minRadius`` -- the quantitative version of the platform's own rule of thumb.
Being a mean-field theory it gets the *scaling* right and the *prefactor* only to
within a constant (measured: :math:`\ell_{\rm meas}\approx0.63\,\ell^\*`); see
``experiments/e04_linear_stability.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .kernels import hankel_basis, hankel_matrix, hankel_transform

__all__ = [
    "DispersionResult",
    "LinearTheory",
    "sigma_matrix",
    "dispersion",
    "scan",
    "most_unstable",
    "classify_linear",
    "exceptional_point_2species",
]


# --------------------------------------------------------------------------- #
def sigma_matrix(k, params, rho_bar):
    r"""Eigenvalues :math:`\sigma_n(k)` of :math:`\operatorname{diag}(\bar\rho)\hat U(k)`.

    Returns ``(K, S)`` complex array (eigenvalues are complex whenever the kernel
    matrix is non-symmetric).
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    Uhat = hankel_matrix(k, params.A, params.alpha, params.beta, params.repel)
    rho = np.asarray(rho_bar, dtype=float)
    M = rho[None, :, None] * Uhat                    # diag(rho) @ Uhat
    return np.linalg.eigvals(M)                      # (K, S)


def dispersion(k, params, rho_bar):
    r"""Growth rates :math:`\lambda_\pm(k)` for every branch.

    Returns ``(K, S, 2)`` complex array; ``[..., 0]`` is the :math:`\lambda_+` root.
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    sig = sigma_matrix(k, params, rho_bar)           # (K,S)
    g = params.gamma
    disc = g**2 - 240.0 * params.force_factor * (k**2)[:, None] * sig
    root = np.sqrt(disc.astype(complex))
    lam_p = 0.5 * (-g + root)
    lam_m = 0.5 * (-g - root)
    return np.stack([lam_p, lam_m], axis=-1)


def _summarise(k, lam, Uhat0, sigma, tol=1e-9):
    """Reduce a (K, S, 2) spectrum to the dominant static and oscillatory modes."""
    flat = lam.reshape(len(k), -1)
    re, im = flat.real, np.abs(flat.imag)

    i_all = np.unravel_index(int(np.argmax(re)), re.shape)
    lam_star = flat[i_all]

    osc = im > tol
    if osc.any():
        masked = np.where(osc, re, -np.inf)
        i_osc = np.unravel_index(int(np.argmax(masked)), re.shape)
        osc_growth = float(re[i_osc])
        osc_k = float(k[i_osc[0]])
        osc_omega = float(im[i_osc])
    else:
        osc_growth, osc_k, osc_omega = -np.inf, np.nan, 0.0

    growth = float(lam_star.real)
    omega = float(abs(lam_star.imag))
    return DispersionResult(
        k=k, lam=lam, sigma=sigma,
        k_star=float(k[i_all[0]]),
        growth=growth, omega=omega,
        length=float(2 * np.pi / k[i_all[0]]),
        speed=float(omega / k[i_all[0]]),
        travelling=bool(growth > 0 and omega > tol),
        osc_growth=osc_growth, osc_k=osc_k, osc_omega=osc_omega,
        osc_speed=float(osc_omega / osc_k) if osc_k == osc_k and osc_k > 0 else 0.0,
        Uhat0=Uhat0,
    )


@dataclass
class DispersionResult:
    """Outcome of the linear analysis.

    Two branches are reported because they answer different questions.

    ``growth`` / ``k_star`` -- the *fastest* growing mode of any kind.  This sets
    the length scale :math:`\\ell^\\*=2\\pi/k^\\*` of the first structures to appear.

    ``osc_growth`` / ``osc_k`` / ``osc_speed`` -- the fastest growing mode that is
    genuinely **oscillatory** (:math:`\\operatorname{Im}\\lambda\\neq0`).  In almost
    every Particle-Life universe the dominant instability is a *real* condensation
    mode -- clumps form first -- and the oscillatory branch grows more slowly.  It
    is the oscillatory branch that says whether the resulting structures will
    *move*, so it must be tracked separately rather than being hidden behind the
    global maximum.
    """

    k: np.ndarray
    lam: np.ndarray            # (K,S,2) complex
    sigma: np.ndarray          # (K,S)   complex
    k_star: float
    growth: float              # max Re lambda over k and branches
    omega: float               # |Im lambda| of that mode
    length: float              # 2 pi / k_star
    speed: float               # |Im lambda| / k_star  (phase velocity)
    travelling: bool
    osc_growth: float          # max Re lambda among modes with Im != 0
    osc_k: float
    osc_omega: float
    osc_speed: float
    Uhat0: np.ndarray          # (S,S) kernel at k -> 0

    @property
    def has_oscillatory(self) -> bool:
        return bool(self.osc_growth > 0)

    def regime(self):
        if self.growth <= 0:
            return "uniform"
        if self.travelling:
            return "travelling"
        return "condense+move" if self.has_oscillatory else "static"

    def summary(self):
        return (f"{self.regime():>14s}: k*={self.k_star:.4f}  l*={self.length:7.1f}  "
                f"Reλ={self.growth:+.4f}  |  osc: Reλ={self.osc_growth:+.4f} "
                f"k={self.osc_k:.4f} c={self.osc_speed:.3f}")


def scan(params, rho_bar, kmin=1e-3, kmax=None, nk=400):
    """Full linear analysis over a band of wavenumbers."""
    if kmax is None:
        kmax = 6.0 * np.pi / float(np.min(params.alpha))
    k = np.linspace(kmin, kmax, nk)
    lam = dispersion(k, params, rho_bar)
    sig = sigma_matrix(k, params, rho_bar)
    Uhat0 = hankel_matrix(np.array([1e-6]), params.A, params.alpha,
                          params.beta, params.repel)[0]
    return _summarise(k, lam, Uhat0, sig)


def most_unstable(params, rho_bar, **kw):
    r = scan(params, rho_bar, **kw)
    return r.k_star, r.growth, r.travelling


def classify_linear(params, rho_bar, **kw):
    """``'uniform'`` | ``'static'`` | ``'condense+move'`` | ``'travelling'``."""
    return scan(params, rho_bar, **kw).regime()


# --------------------------------------------------------------------------- #
#  the exactly solvable two-species case
# --------------------------------------------------------------------------- #
def exceptional_point_2species(k, params, rho_bar):
    r"""Discriminant of the 2x2 problem; zero <=> **exceptional point**.

    For two species with :math:`\hat u_{ab}(k)`,

    .. math::
        \sigma_\pm = \tfrac12(\rho_1\hat u_{11}+\rho_2\hat u_{22})
                     \pm\sqrt{\Delta^2+\rho_1\rho_2\hat u_{12}\hat u_{21}},\quad
        \Delta=\tfrac12(\rho_1\hat u_{11}-\rho_2\hat u_{22}).

    :math:`\sigma_\pm` become a complex-conjugate pair exactly when

    .. math::  D \;\equiv\; \Delta^2 + \rho_1\rho_2\,\hat u_{12}\hat u_{21} \;<\;0,

    which requires :math:`\hat u_{12}\hat u_{21}<0` -- i.e. **one species chases
    while the other flees**.  The locus :math:`D=0` is the exceptional point where
    the two eigenvectors coalesce; crossing it turns a pair of static demixing
    modes into a single travelling mode.  This is the Particle-Life instance of
    the non-reciprocal phase transition of Fruchart et al., *Nature* **592**, 363
    (2021).
    """
    k = np.atleast_1d(np.asarray(k, dtype=float))
    U = hankel_matrix(k, params.A, params.alpha, params.beta, params.repel)
    r1, r2 = rho_bar
    delta = 0.5 * (r1 * U[:, 0, 0] - r2 * U[:, 1, 1])
    cross = r1 * r2 * U[:, 0, 1] * U[:, 1, 0]
    return delta**2 + cross


# --------------------------------------------------------------------------- #
#  cached engine: the kernel is linear in A, so a whole phase diagram is cheap
# --------------------------------------------------------------------------- #
class LinearTheory:
    r"""Reusable linear-stability engine for a **fixed pair of radius matrices**.

    Because :math:`\hat U = R\hat C + A\odot\hat T` (see ``kernels.hankel_basis``),
    the Bessel quadrature only has to be done once.  Scanning thousands of force
    matrices then costs nothing but ``S x S`` eigenvalue problems.
    """

    def __init__(self, alpha, beta, gamma, force_factor=1.0, repel=1.0,
                 kmin=1e-3, kmax=None, nk=400, n_quad=3000):
        self.alpha = np.asarray(alpha, float)
        self.beta = np.asarray(beta, float)
        self.gamma = float(gamma)
        self.kappa = float(force_factor)
        self.repel = float(repel)
        if kmax is None:
            kmax = 6.0 * np.pi / float(self.alpha.min())
        self.k = np.linspace(kmin, kmax, nk)
        self.Chat, self.That = hankel_basis(self.k, self.alpha, self.beta, n_quad)

    @classmethod
    def from_params(cls, params, **kw):
        return cls(params.alpha, params.beta, params.gamma,
                   params.force_factor, params.repel, **kw)

    def Uhat(self, A):
        return self.repel * self.Chat + np.asarray(A, float)[None, :, :] * self.That

    def sigma(self, A, rho_bar):
        M = np.asarray(rho_bar, float)[None, :, None] * self.Uhat(A)
        return np.linalg.eigvals(M)

    def analyse(self, A, rho_bar):
        sig = self.sigma(A, rho_bar)
        disc = self.gamma**2 - 240.0 * self.kappa * (self.k**2)[:, None] * sig
        root = np.sqrt(disc.astype(complex))
        lam = np.stack([0.5 * (-self.gamma + root), 0.5 * (-self.gamma - root)], -1)
        return _summarise(self.k, lam, self.Uhat(A)[0], sig)


def two_species_phase_boundary(a12_grid, a21_grid, base_params, rho_bar, **kw):
    """Maps of ``growth``, ``osc_growth`` and ``osc_speed`` over the ``(A12, A21)`` plane."""
    lt = LinearTheory.from_params(base_params, **kw)
    shape = (len(a21_grid), len(a12_grid))
    growth = np.zeros(shape)
    osc = np.full(shape, -np.inf)
    speed = np.zeros(shape)
    A = base_params.A.copy()
    for jj, a21 in enumerate(a21_grid):
        for ii, a12 in enumerate(a12_grid):
            A[0, 1], A[1, 0] = a12, a21
            r = lt.analyse(A, rho_bar)
            growth[jj, ii] = r.growth
            osc[jj, ii] = r.osc_growth
            speed[jj, ii] = r.osc_speed
    return growth, osc, speed
