r"""
Phase behaviour of a protein mixture from measured Kd's, sizes and valences.
===========================================================================

Free-energy density (units: kT per nm^3), the standard stickers-and-spacers
thermodynamics:

.. math::
    f(\{\rho\}) = \underbrace{\sum_a \rho_a[\ln\rho_a - 1]}_{\text{ideal}}
    \;+\; \underbrace{f^{\rm ex}_{\rm HS}(\{\rho\})}_{\text{excluded volume, BMCSL}}
    \;+\; \underbrace{\sum_a \rho_a n_a\Bigl[\ln X_a - \tfrac{X_a}{2} + \tfrac12\Bigr]}_{\text{association, Wertheim TPT1}}

with the unbonded-site fractions solving the mass-action closure

.. math::  X_a = \Bigl[1 + \sum_b \rho_b\,n_b\,X_b\,\Delta_{ab}\Bigr]^{-1},
           \qquad \Delta_{ab} = \frac{g^{HS}_{ab}(\sigma_{ab})}{K_{d,ab}^{(n)}} .

Why Wertheim and not a mean-field well
--------------------------------------
Because **valence is not optional**.  An isotropic attractive well lets a
particle make ~12 bonds; a protein with one binding site makes one.  Wertheim
counts bonds correctly, so it reproduces the qualitative facts that any usable
model must get right:

* ``n = 1`` -> only dimers.  **No phase separation at any concentration**, no
  matter how tight the Kd.  (Verified numerically in ``test_bio.py``.)
* ``n = 2`` -> linear chains; still no bulk condensation on its own.
* ``n >= 3`` -> branching, gelation and liquid-liquid phase separation.

and because :math:`\Delta_{ab}` *is* :math:`1/K_d` -- the measured quantity
enters with no fitted well width.

What this does and does not give you
------------------------------------
``✅`` binodal (the c_sat you measure), tie lines, partition coefficients,
spinodal, which species end up in which phase.

``⚠️`` interfacial tensions and therefore the core-shell/separate topology:
these need a gradient coefficient that is *estimated* from the interaction
range, so read the *ordering* of the tensions, not their absolute values.

``✗`` anything kinetic (nucleation times, coarsening, gel arrest), anything
driven (ATP-dependent transport, active droplets), and sequence-specific
effects beyond what the Kd/valence summary carries.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import ConvexHull

from .units import M_TO_NM3, conc_to_density, density_to_conc

__all__ = ["Mixture", "free_energy", "f_excess", "chemical_potentials", "osmotic_pressure",
           "spinodal_eig", "spinodal_modes", "mode_coherence", "hessian", "coexistence",
           "csat_scan", "partition_coefficients"]


# --------------------------------------------------------------------------- #
@dataclass
class Mixture:
    """A protein mixture specified the way a bench experiment specifies it."""

    names: list                       # species labels
    diameters: np.ndarray             # nm, = 2 * radius
    valence: np.ndarray               # number of equivalent binding sites
    kd: np.ndarray                    # (S,S) dissociation constants in M; inf = no binding
    kT: float = 1.0                   # energies are already in kT

    def __post_init__(self):
        self.diameters = np.asarray(self.diameters, float)
        self.valence = np.asarray(self.valence, float)
        self.kd = np.asarray(self.kd, float)
        S = len(self.names)
        assert self.diameters.shape == (S,) and self.valence.shape == (S,)
        assert self.kd.shape == (S, S)
        if not np.allclose(self.kd, self.kd.T, equal_nan=True):
            raise ValueError(
                "Kd must be symmetric: proteins obey Newton's third law. "
                "A non-symmetric interaction matrix is only meaningful for an "
                "ATP-driven (active) system, which this module does not model.")

    @property
    def S(self):
        return len(self.names)

    @property
    def sigma(self):
        """Contact distances sigma_ab = (d_a + d_b)/2."""
        d = self.diameters
        return 0.5 * (d[:, None] + d[None, :])

    def copy(self, **kw):
        d = dict(names=list(self.names), diameters=self.diameters.copy(),
                 valence=self.valence.copy(), kd=self.kd.copy(), kT=self.kT)
        d.update(kw)
        return Mixture(**d)

    def assoc_f(self, rho):
        r"""Wertheim association free-energy density (kT/nm^3).

        Factored out as a method so that everything built on top of ``f`` --
        :func:`coexistence`, :func:`hessian`, :func:`csat_scan`,
        ``phases.phase_diagram`` -- works unchanged for the multi-site-class
        generalisation in :mod:`plife.bio.sites`, which supplies its own.
        """
        X = unbonded_fractions(rho, self)
        return (rho * self.valence[None, :]
                * (np.log(np.maximum(X, 1e-300)) - 0.5 * X + 0.5)).sum(axis=1)


# --------------------------------------------------------------------------- #
#  hard-sphere reference: BMCSL mixture
# --------------------------------------------------------------------------- #
def _xi(rho, d):
    """Scaled moments xi_m = (pi/6) sum_a rho_a d_a^m, m = 0..3."""
    rho = np.atleast_2d(np.asarray(rho, float))
    return [(np.pi / 6.0) * (rho * d[None, :] ** m).sum(axis=1) for m in range(4)]


def f_hs_excess(rho, d):
    """BMCSL excess free-energy density (kT / nm^3) of a hard-sphere mixture."""
    x0, x1, x2, x3 = _xi(rho, d)
    x3 = np.clip(x3, 0.0, 0.99)
    one = 1.0 - x3
    with np.errstate(divide="ignore", invalid="ignore"):
        term = ((x2**3 / np.maximum(x3, 1e-30) ** 2 - x0) * np.log(one)
                + 3.0 * x1 * x2 / one
                + x2**3 / (np.maximum(x3, 1e-30) * one**2))
    return np.where(x3 > 1e-14, term, 0.0)


def g_contact(rho, d):
    """Boublik contact values ``g_ab(sigma_ab)`` -> shape (K, S, S)."""
    _, _, x2, x3 = _xi(rho, d)
    x3 = np.clip(x3, 0.0, 0.99)
    one = (1.0 - x3)[:, None, None]
    x2 = x2[:, None, None]
    dd = (d[:, None] * d[None, :]) / (d[:, None] + d[None, :])
    dd = dd[None, :, :]
    return 1.0 / one + 3.0 * x2 * dd / one**2 + 2.0 * x2**2 * dd**2 / one**3


# --------------------------------------------------------------------------- #
#  Wertheim TPT1 association
# --------------------------------------------------------------------------- #
def _delta(rho, mix):
    r"""Association volumes :math:`\Delta_{ab}` (nm^3) with the crowding factor.

    .. math::  \Delta_{ab} = s_{ab}\,\frac{g^{HS}_{ab}(\sigma_{ab})}{K_{d,ab}^{(n)}},
               \qquad s_{ab} = \begin{cases}2 & a=b\\ 1 & a\neq b\end{cases}

    The symmetry factor :math:`s_{aa}=2` is **not** a fudge: laboratory Kd for a
    homodimer is defined by :math:`K_d = [A]^2/[A_2]`, and because the two
    partners are the same molecule that reaction has twice the site-pair
    multiplicity of the heterodimer convention :math:`K_d=[A][B]/[AB]`.  Without
    it, Wertheim returns exactly half the dimers required by the law of mass
    action -- which is how the factor was found (see ``test_bio.py``, where both
    limits are now reproduced to 1e-4).
    """
    g = g_contact(rho, mix.diameters)                        # (K,S,S)
    kd_n = mix.kd * M_TO_NM3                                 # nm^-3
    with np.errstate(divide="ignore", invalid="ignore"):
        base = np.where(np.isfinite(kd_n) & (kd_n > 0), 1.0 / kd_n, 0.0)
    sym = np.ones_like(base) + np.eye(base.shape[0])         # 2 on the diagonal
    return g * (sym * base)[None, :, :]


def unbonded_fractions(rho, mix, tol=1e-12, itmax=500):
    r"""Solve the Wertheim closure for ``X_a`` by damped fixed-point iteration.

    ``X_a`` is the fraction of species-*a* sites that are **not** bonded; it is
    the single most informative intermediate quantity in the whole model, since
    ``1 - X_a`` is the bond saturation you would compare against a pull-down or
    an FCS-derived bound fraction.
    """
    rho = np.atleast_2d(np.asarray(rho, float))
    K, S = rho.shape
    D = _delta(rho, mix)                                     # (K,S,S)
    n = mix.valence
    X = np.ones((K, S))
    for _ in range(itmax):
        # sum_b rho_b n_b X_b Delta_ab
        acc = np.einsum("kab,kb->ka", D, rho * n[None, :] * X)
        Xnew = 1.0 / (1.0 + acc)
        if np.max(np.abs(Xnew - X)) < tol:
            X = Xnew
            break
        X = 0.5 * (X + Xnew)                                 # damping
    return X


def free_energy(rho, mix):
    """Total free-energy density ``f`` in kT/nm^3.  ``rho`` is (K, S) or (S,)."""
    rho = np.atleast_2d(np.asarray(rho, float))
    with np.errstate(divide="ignore", invalid="ignore"):
        f_id = np.where(rho > 0, rho * (np.log(np.maximum(rho, 1e-300)) - 1.0), 0.0
                        ).sum(axis=1)
    return f_id + f_hs_excess(rho, mix.diameters) + mix.assoc_f(rho)


def f_excess(rho, mix):
    """Non-ideal part of ``f``: hard-sphere + association (kT/nm^3).

    Kept separate because the ideal term ``sum rho(ln rho - 1)`` has an exactly
    known derivative that spans many decades; differentiating it numerically
    together with the rest destroys the low-density behaviour.
    """
    rho = np.atleast_2d(np.asarray(rho, float))
    return f_hs_excess(rho, mix.diameters) + mix.assoc_f(rho)


def _rel_step(rho, a, h):
    """Per-row finite-difference step, relative to that row's own density."""
    S = rho.shape[1]
    st = np.zeros_like(rho)
    col = rho[:, a]
    st[:, a] = h * np.maximum(col, 1e-12 * np.maximum(rho.max(axis=1), 1e-30))
    st[:, a] = np.maximum(st[:, a], 1e-300)
    return st


def chemical_potentials(rho, mix, h=1e-5):
    r"""``mu_a = ln rho_a + d f_ex / d rho_a``  (kT).

    The ideal part is analytic; only the excess part is differentiated, with a
    step proportional to each point's own density so the result stays accurate
    over the eight decades of concentration a phase diagram spans.
    """
    rho = np.atleast_2d(np.asarray(rho, float))
    S = rho.shape[1]
    mu = np.empty_like(rho)
    for a in range(S):
        st = _rel_step(rho, a, h)
        d = 2 * st[:, a]
        mu[:, a] = (f_excess(rho + st, mix) - f_excess(np.maximum(rho - st, 0.0), mix)) / d
    with np.errstate(divide="ignore"):
        mu += np.log(np.maximum(rho, 1e-300))
    return mu


def hessian(rho, mix, h=1e-4):
    r"""Stability matrix :math:`\partial^2 f/\partial\rho_a\partial\rho_b` -> (K,S,S).

    Written as :math:`\delta_{ab}/\rho_a + \partial^2 f_{\rm ex}/\partial\rho_a\partial\rho_b`
    so the ideal (stabilising, :math:`1/\rho`) contribution is exact.  That term
    diverges as :math:`\rho\to0`, which is precisely why a dilute mixture is
    always stable and why there is a **threshold concentration** at all -- the
    feature that athermal Particle Life does not have.
    """
    rho = np.atleast_2d(np.asarray(rho, float))
    K, S = rho.shape
    H = np.zeros((K, S, S))
    for a in range(S):
        st = _rel_step(rho, a, h)
        d = 2 * st[:, a]
        gp = np.empty((K, S))
        gm = np.empty((K, S))
        for b in range(S):
            stb = _rel_step(rho, b, h)
            db = 2 * stb[:, b]
            gp[:, b] = (f_excess(rho + st + stb, mix)
                        - f_excess(np.maximum(rho + st - stb, 0.0), mix)) / db
            gm[:, b] = (f_excess(np.maximum(rho - st, 0.0) + stb, mix)
                        - f_excess(np.maximum(np.maximum(rho - st, 0.0) - stb, 0.0), mix)) / db
        H[:, a, :] = (gp - gm) / d[:, None]
    H = 0.5 * (H + np.transpose(H, (0, 2, 1)))
    for a in range(S):
        H[:, a, a] += 1.0 / np.maximum(rho[:, a], 1e-300)
    return H


def spinodal_eig(rho, mix, h=1e-4):
    r"""Smallest eigenvalue of :func:`hessian`.

    Positive => locally stable.  Negative => inside the spinodal.
    """
    return np.linalg.eigvalsh(hessian(rho, mix, h))[:, 0]


def spinodal_modes(rho, mix, h=1e-4, scale=True):
    r"""Eigenvalues **and eigenvectors** of the stability matrix.

    The eigenvector of every *negative* eigenvalue is the composition
    fluctuation that grows, and its sign pattern says what kind of structure
    that growth produces -- the protein analogue of the mode-coherence test used
    on the Particle Life dispersion relation:

    * all components the **same sign**  -> those species move into the dense
      phase *together*: co-condensation, one shared body;
    * **mixed signs**                   -> the species with one sign concentrate
      while the others are expelled: demixing, and the body excludes them.

    The number of negative eigenvalues bounds how many independent demixing
    directions the mixture has, hence roughly ``n_bodies - 1``.

    With ``scale=True`` the Hessian is symmetrically rescaled by
    :math:`\sqrt{\rho_a\rho_b}` first.  Without that, the exact :math:`1/\rho_a`
    diagonal makes the eigenvectors of a mixture spanning several decades of
    concentration point almost entirely along the most dilute species, which
    says nothing about the physics.  Rescaling asks the scale-free question
    "which *relative* composition fluctuation grows?", and leaves the *sign* of
    every eigenvalue unchanged (it is a congruence transform).

    Returns ``(vals, vecs)`` with shapes ``(K,S)`` and ``(K,S,S)``; ``vecs[k,:,i]``
    is the eigenvector for ``vals[k,i]``, ascending.
    """
    rho = np.atleast_2d(np.asarray(rho, float))
    H = hessian(rho, mix, h)
    if scale:
        s = np.sqrt(np.maximum(rho, 1e-300))
        H = H * s[:, :, None] * s[:, None, :]
    return np.linalg.eigh(H)


def mode_coherence(vec, tol=1e-12):
    r"""``|sum_a v_a| / sum_a |v_a|`` for a composition eigenvector.

    1 => every species moves the same way (condensation);
    0 => the species split evenly into enriched and expelled groups (demixing).
    The same order parameter that separated condensation from demixing in the
    Particle Life dispersion analysis, reused verbatim.
    """
    v = np.asarray(vec, float)
    den = np.abs(v).sum(axis=-1)
    return np.where(den > tol, np.abs(v.sum(axis=-1)) / np.maximum(den, tol), 1.0)


# --------------------------------------------------------------------------- #
#  coexistence by common-tangent (lower convex hull) construction
# --------------------------------------------------------------------------- #
def coexistence(mix, grid, tol_rel=1e-9):
    r"""Global phase behaviour by the **lower convex hull** of ``f``.

    Sampling ``f`` on a composition grid and taking the lower convex hull is the
    common-tangent construction done globally: any composition lying *above* the
    hull lowers its free energy by splitting into the phases at the vertices of
    the facet beneath it.  This is far more robust than Newton-solving the
    equal-chemical-potential equations, and it finds three-phase coexistence
    automatically.

    Returns a dict with ``on_hull`` (bool per grid point: single phase),
    ``gap`` (``f - f_hull``, the driving force for demixing) and ``facets``
    (the coexisting compositions).
    """
    rho = np.asarray(grid, float)
    f = free_energy(rho, mix)

    # Drop composition axes that do not vary over the sampled set (e.g. a species
    # held at zero, or a 1-D ray through a 2-D space).  Qhull cannot build a hull
    # in a degenerate embedding, and the common-tangent construction restricted to
    # the sampled affine subspace is exactly what those axes mean.
    span = rho.max(axis=0) - rho.min(axis=0)
    keep = span > 1e-14 * max(rho.max(), 1e-300)
    if not keep.any():
        keep = np.zeros(rho.shape[1], bool)
        keep[0] = True
    red = rho[:, keep]
    S = int(keep.sum())

    pts = np.column_stack([red, f])
    try:
        hull = ConvexHull(pts, qhull_options="Qt")
    except Exception:
        # Nearly-degenerate point sets (two composition columns almost equal, or
        # a very flat free-energy surface) make Qhull's exact predicates fail.
        # Joggling perturbs the input below the resolution of the phase boundary
        # and always produces a simplicial hull.
        hull = ConvexHull(pts, qhull_options="QJ")

    # lower hull = facets whose outward normal has a negative f-component
    lower = [i for i, eq in enumerate(hull.equations) if eq[S] < -1e-12]
    facets = []
    for i in lower:
        v = hull.simplices[i]
        facets.append(dict(vertices=v.tolist(), rho=rho[v].tolist(),
                           f=f[v].tolist(), equation=hull.equations[i].tolist()))

    # For every point, the supporting plane is the HIGHEST lower-hull plane
    # beneath it; that plane is the common tangent, and its facet names the
    # coexisting phases.  Vectorised: (n_lower, K) then argmax.
    if lower:
        eqs = hull.equations[lower]                       # (L, S+2)
        vals = -(red @ eqs[:, :S].T + eqs[:, -1][None, :]) / eqs[:, S][None, :]
        best = np.argmax(vals, axis=1)
        f_hull = vals[np.arange(len(rho)), best]
        facet_of = np.asarray(lower)[best]
    else:
        f_hull = f.copy()
        facet_of = np.full(len(rho), -1)

    gap = f - f_hull
    scale = np.abs(f).max() + 1e-300
    on_hull = gap <= tol_rel * scale
    return dict(rho=rho, f=f, f_hull=f_hull, gap=gap, on_hull=on_hull,
                facets=facets, lower=lower, facet_of=facet_of, keep=keep,
                simplices=hull.simplices, n_lower=len(lower))


def _log_grid(lo_uM, hi_uM, n):
    return conc_to_density(np.geomspace(lo_uM, hi_uM, n))


def csat_scan(mix, species=0, partner_uM=0.0, lo_uM=1e-4, hi_uM=1e5, n=900):
    r"""Saturation concentration of ``species`` at a fixed partner concentration.

    Scans along a ray in composition space and returns the lowest concentration
    at which the mixture stops being a single phase (the **binodal**, i.e. the
    experimentally measured c_sat), plus the spinodal for reference.

    Returns ``dict(csat_uM, cspin_uM, conc_uM, gap, stable)``; ``csat_uM`` is
    ``inf`` when the mixture never phase separates over the scanned range --
    which is the correct answer for, e.g., a monovalent binder.
    """
    S = mix.S
    c = np.geomspace(lo_uM, hi_uM, n)
    rho = np.zeros((n, S))
    rho[:, species] = conc_to_density(c)
    for b in range(S):
        if b != species:
            rho[:, b] = conc_to_density(np.atleast_1d(partner_uM)[0]
                                        if np.isscalar(partner_uM) else partner_uM[b])
    res = coexistence(mix, rho)
    unstable = ~res["on_hull"]
    csat = float(c[unstable][0]) if unstable.any() else np.inf

    eig = spinodal_eig(rho, mix)
    neg = eig < 0
    cspin = float(c[neg][0]) if neg.any() else np.inf
    return dict(csat_uM=csat, cspin_uM=cspin, conc_uM=c, gap=res["gap"],
                stable=res["on_hull"], spinodal_eig=eig)


def partition_coefficients(mix, rho_dilute, rho_dense):
    """``K_a = rho_a^dense / rho_a^dilute`` for each species."""
    d = np.atleast_2d(np.asarray(rho_dense, float))
    l = np.atleast_2d(np.asarray(rho_dilute, float))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(l > 0, d / l, np.inf)
