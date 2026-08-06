r"""
Proteins with more than one kind of binding site.
=================================================

:class:`~plife.bio.thermo.Mixture` gives every protein **one class of equivalent
sites**.  That is fine for a homotypic sticker but wrong for essentially every
real signalling protein, and the error is not small.  Grb2 is ``SH3-SH2-SH3``:
one SH2 that binds phospho-LAT and two SH3s that bind SOS1.  Collapse that to
"valence 3, binds LAT at 1 uM" and the model lets a single Grb2 bridge *two*
LAT molecules through two SH2 domains it does not have -- which manufactures a
percolating network out of a mixture that experimentally stays a soup.

So here the site index is promoted to its own dimension.  Species *a* carries
site classes :math:`i` with multiplicities :math:`n_{a i}`, and the Wertheim
closure runs over the flattened list :math:`\alpha=(a,i)`:

.. math::
    X_\alpha = \Bigl[1+\sum_\beta \rho_{s(\beta)}\,n_\beta X_\beta
               \Delta_{\alpha\beta}\Bigr]^{-1},
    \qquad
    \Delta_{\alpha\beta} = s_{\alpha\beta}\,
        \frac{g^{HS}_{s(\alpha)s(\beta)}(\sigma)}{K^{(n)}_{d,\alpha\beta}},

.. math::
    f_{\rm assoc} = \sum_\alpha \rho_{s(\alpha)} n_\alpha
                    \Bigl[\ln X_\alpha-\tfrac{X_\alpha}2+\tfrac12\Bigr],

with :math:`s(\alpha)` the species owning site class :math:`\alpha` and the
symmetry factor :math:`s_{\alpha\beta}=2` only when :math:`\alpha=\beta`, i.e.
for a genuinely self-complementary site class.  Two *different* classes on the
same molecule binding each other (head-to-tail polymerisation, DIX/DAX) get 1.

Percolation with heterogeneous sites
------------------------------------
The branching argument has to be redone, and this is where the physics actually
changes.  Arriving at a molecule of species *a* through a site of class *i*
consumes that site, leaving :math:`n_{ak}-\delta_{ik}` sites of each class *k*:

.. math::  M_{(ai),(bj)} \;=\; \sum_k \bigl(n_{ak}-\delta_{ik}\bigr)\,
                              \pi_{(ak),(bj)},
           \qquad \pi_{\alpha\beta}=X_\alpha\Delta_{\alpha\beta}
                                    \rho_{s(\beta)}n_\beta X_\beta ,

and the network percolates when :math:`\varrho(M)\ge1`.  Enter Grb2 by its SH2
and only the two SH3s remain, so it can pass the walk on to SOS1 but never back
to a second LAT -- exactly the constraint the single-class model threw away.
Everything reduces to :math:`(n_A-1)(n_B-1)p_Ap_B\ge1` when each species has one
class.

Everything downstream is unchanged: :class:`SiteMixture` exposes the same
``diameters``/``assoc_f`` interface as :class:`~plife.bio.thermo.Mixture`, so
``coexistence``, ``hessian``, ``csat_scan`` and ``phase_diagram`` all work on it
as they stand.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .thermo import g_contact
from .units import M_TO_NM3, conc_to_density

__all__ = ["Site", "SiteMixture", "unbonded_fractions", "branching_matrix", "percolates"]


@dataclass
class Site:
    """One class of equivalent binding sites on one protein."""

    species: str
    name: str
    n: float

    @property
    def label(self):
        return f"{self.species}.{self.name}"


class SiteMixture:
    """A mixture whose proteins may each carry several distinct site classes."""

    def __init__(self, names, diameters, sites, kd, kT=1.0):
        self.names = list(names)
        self.diameters = np.asarray(diameters, float)
        self.sites = list(sites)
        self.kd = np.asarray(kd, float)
        self.kT = kT
        S, A = len(self.names), len(self.sites)
        assert self.diameters.shape == (S,)
        assert self.kd.shape == (A, A), f"kd must be {A}x{A}, got {self.kd.shape}"
        if not np.allclose(self.kd, self.kd.T, equal_nan=True):
            raise ValueError("Kd must be symmetric: proteins obey Newton's third law. "
                             "A non-symmetric interaction matrix is only meaningful "
                             "for an ATP-driven (active) system.")
        idx = {n: i for i, n in enumerate(self.names)}
        self.owner = np.array([idx[s.species] for s in self.sites])
        self.n = np.array([s.n for s in self.sites], float)

    @classmethod
    def from_valence(cls, names, diameters, valence, kd, kT=1.0):
        """One site class per species -- the :class:`~plife.bio.thermo.Mixture` case.

        Reproduces ``Mixture`` bit-for-bit (checked in ``test_pathways.py``), so
        it is the migration path for code that has only a total valence.
        """
        sites = [Site(n, "site", float(v)) for n, v in zip(names, valence)]
        return cls(names, diameters, sites, kd, kT)

    # -- interface shared with thermo.Mixture -------------------------------- #
    @property
    def S(self):
        return len(self.names)

    @property
    def A(self):
        return len(self.sites)

    @property
    def valence(self):
        """Total sites per species, for reporting only."""
        v = np.zeros(self.S)
        np.add.at(v, self.owner, self.n)
        return v

    @property
    def sigma(self):
        d = self.diameters
        return 0.5 * (d[:, None] + d[None, :])

    @property
    def labels(self):
        return [s.label for s in self.sites]

    def delta(self, rho):
        r"""Association volumes :math:`\Delta_{\alpha\beta}` -> (K, A, A)."""
        g = g_contact(rho, self.diameters)                    # (K,S,S)
        g_sites = g[:, self.owner[:, None], self.owner[None, :]]
        kd_n = self.kd * M_TO_NM3
        with np.errstate(divide="ignore", invalid="ignore"):
            base = np.where(np.isfinite(kd_n) & (kd_n > 0), 1.0 / kd_n, 0.0)
        sym = np.ones_like(base) + np.eye(self.A)             # 2 only on the diagonal
        return g_sites * (sym * base)[None, :, :]

    def assoc_f(self, rho):
        """Wertheim association free-energy density (kT/nm^3), summed over sites."""
        rho = np.atleast_2d(np.asarray(rho, float))
        X = unbonded_fractions(rho, self)
        rho_site = rho[:, self.owner]                         # (K,A)
        return (rho_site * self.n[None, :]
                * (np.log(np.maximum(X, 1e-300)) - 0.5 * X + 0.5)).sum(axis=1)


def unbonded_fractions(rho, sm: SiteMixture, tol=1e-12, itmax=500):
    """Solve the site-resolved closure for ``X_alpha`` -> (K, A)."""
    rho = np.atleast_2d(np.asarray(rho, float))
    K = rho.shape[0]
    D = sm.delta(rho)                                         # (K,A,A)
    rho_site = rho[:, sm.owner]                               # (K,A)
    X = np.ones((K, sm.A))
    for _ in range(itmax):
        acc = np.einsum("kab,kb->ka", D, rho_site * sm.n[None, :] * X)
        Xnew = 1.0 / (1.0 + acc)
        if np.max(np.abs(Xnew - X)) < tol:
            X = Xnew
            break
        X = 0.5 * (X + Xnew)
    return X


def branching_matrix(rho, sm: SiteMixture):
    r"""Mean-offspring matrix :math:`M_{(ai),(bj)}` of the bond-network walk."""
    rho = np.atleast_2d(np.asarray(rho, float))
    X = unbonded_fractions(rho, sm)[0]
    D = sm.delta(rho)[0]
    rho_site = rho[0][sm.owner]
    pi = X[:, None] * D * (rho_site * sm.n * X)[None, :]      # (A,A), rows sum to 1-X
    # same_mol[i,k] = 1 when site classes i and k sit on the same species
    same = (sm.owner[:, None] == sm.owner[None, :]).astype(float)
    remaining = same * sm.n[None, :] - np.eye(sm.A)           # n_ak - delta_ik
    return np.maximum(remaining, 0.0) @ pi


def percolates(conc_uM, sm: SiteMixture):
    """``(bool, spectral_radius, bonded_site_fractions)`` at one composition."""
    rho = conc_to_density(np.atleast_2d(np.asarray(conc_uM, float)))
    M = branching_matrix(rho, sm)
    lam = float(np.max(np.abs(np.linalg.eigvals(M))))
    return lam >= 1.0, lam, 1 - unbonded_fractions(rho, sm)[0]
