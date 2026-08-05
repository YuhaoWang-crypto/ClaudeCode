r"""
Two-protein phase diagrams, partitioning, and the morphology call.
==================================================================

Everything here is built on the lower-convex-hull (common-tangent) construction
in ``thermo.coexistence``.  For a two-component mixture that gives, at every
composition, either

* a point **on** the hull -> one homogeneous phase, or
* a point **above** it    -> the hull facet underneath names the coexisting
  phases and the lever rule gives how much of each.

From the facet we read off everything an experiment measures: the saturation
concentration, the tie line, the partition coefficient of each protein, and how
many distinct dense phases there are.

Morphology
----------
The number and composition of coexisting phases already fixes most of the
picture::

    1 phase                      -> homogeneous (sub-saturation)
    2 phases, both enriched      -> co-condensation (one droplet type, both proteins)
    2 phases, one depleted       -> condensate that EXCLUDES the partner
    3 phases                     -> two immiscible condensates

For three phases the *arrangement* (separate droplets vs one engulfing the
other) is set by interfacial tensions through the spreading coefficient
:math:`S = \gamma_{13}-\gamma_{12}-\gamma_{23}`; that part is ⚠️
semi-quantitative here because the square-gradient coefficient is estimated
from the interaction range rather than measured.  Read the *ordering*, not the
absolute tensions.
"""

from __future__ import annotations

import numpy as np

from .thermo import (Mixture, coexistence, free_energy, spinodal_eig,
                     unbonded_fractions)
from .units import conc_to_density, density_to_conc, volume_fraction

__all__ = ["grid2d", "phase_diagram", "partner_window", "describe_point"]


# --------------------------------------------------------------------------- #
def grid2d(cA_uM, cB_uM):
    """Cartesian product of two concentration axes -> (K,2) density array."""
    A, B = np.meshgrid(conc_to_density(np.asarray(cA_uM, float)),
                       conc_to_density(np.asarray(cB_uM, float)), indexing="ij")
    return np.column_stack([A.ravel(), B.ravel()]), A.shape


def _facet_split(res, idx, rho, sliver_tol=0.03):
    r"""Coexisting compositions at grid point ``idx``, or ``None`` if single phase.

    Qhull triangulates the lower hull, so a **two-phase tie-line region** — which
    is a ruled (developable) surface — comes back as long thin triangles whose
    three vertices are nearly collinear in composition space.  Counting facet
    vertices therefore reports "3 phases" everywhere unless those slivers are
    detected.  Here a facet is called degenerate when its area is a negligible
    fraction of its longest edge squared, in which case only the two extreme
    vertices are kept and the region is correctly a two-phase tie line.
    """
    if res["on_hull"][idx]:
        return None
    fi = res["facet_of"][idx]
    if fi < 0:
        return None
    v = rho[res["simplices"][fi]]                    # candidate coexisting phases
    x = rho[idx]
    if len(v) < 3:                                   # reduced (1-D) hull: a tie line
        seg = v[1] - v[0]
        t = float(np.dot(x - v[0], seg) / max(np.dot(seg, seg), 1e-300))
        t = min(max(t, 0.0), 1.0)
        return v, np.array([1 - t, t])

    # collinearity test in composition space
    e1, e2 = v[1] - v[0], v[2] - v[0]
    area = abs(e1[0] * e2[1] - e1[1] * e2[0])
    longest = max(np.linalg.norm(v[i] - v[j]) for i in range(3) for j in range(i + 1, 3))
    degenerate = area <= sliver_tol * longest**2

    if degenerate:                                   # 2-phase tie line
        d = [np.linalg.norm(v[i] - v[j]) for i, j in ((0, 1), (0, 2), (1, 2))]
        i, j = [(0, 1), (0, 2), (1, 2)][int(np.argmax(d))]
        v2 = v[[i, j]]
        seg = v2[1] - v2[0]
        t = float(np.dot(x - v2[0], seg) / max(np.dot(seg, seg), 1e-300))
        t = min(max(t, 0.0), 1.0)
        return v2, np.array([1 - t, t])

    Aeq = np.vstack([v.T, np.ones(len(v))])
    lam, *_ = np.linalg.lstsq(Aeq, np.append(x, 1.0), rcond=None)
    return v, lam


def phase_diagram(mix: Mixture, cA_uM, cB_uM, tol_rel=1e-9):
    r"""Full two-protein phase diagram over a concentration grid.

    Returns a dict with, on the ``(nA, nB)`` grid:

    ``single``      True where the mixture stays homogeneous
    ``n_phases``    1, 2 or 3 coexisting phases
    ``gap``         :math:`f - f_{\rm hull}`, the demixing driving force
    ``spinodal``    True inside the spinodal (barrierless demixing)
    ``K_A``, ``K_B``  partition coefficients (dense/dilute) where two phases coexist
    ``bound_A``, ``bound_B``  bond saturation ``1 - X_a`` (comparable to a pull-down)
    ``label``       integer morphology code, see ``LABELS``
    """
    rho, shape = grid2d(cA_uM, cB_uM)
    res = coexistence(mix, rho, tol_rel=tol_rel)
    spin = spinodal_eig(rho, mix) < 0
    X = unbonded_fractions(rho, mix)

    n_ph = np.ones(len(rho), dtype=int)
    KA = np.full(len(rho), np.nan)
    KB = np.full(len(rho), np.nan)
    label = np.zeros(len(rho), dtype=int)

    for i in range(len(rho)):
        sp = _facet_split(res, i, rho)
        if sp is None:
            label[i] = 0                                  # homogeneous
            continue
        v, lam = sp
        keep = lam > 1e-3
        if keep.sum() >= 2:
            v, lam = v[keep], lam[keep]
        n_ph[i] = len(v)
        order = np.argsort(v.sum(axis=1))                 # dilute -> dense
        v, lam = v[order], lam[order]
        dil, dense = v[0], v[-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            KA[i] = dense[0] / dil[0] if dil[0] > 0 else np.inf
            KB[i] = dense[1] / dil[1] if dil[1] > 0 else np.inf
        if n_ph[i] >= 3:
            label[i] = 3                                  # two immiscible condensates
        elif KB[i] >= 2.0 and KA[i] >= 2.0:
            label[i] = 1                                  # co-condensation
        else:
            label[i] = 2                                  # condensate excluding the partner

    rs = lambda a: a.reshape(shape)
    return dict(
        cA=np.asarray(cA_uM, float), cB=np.asarray(cB_uM, float),
        rho=rho, single=rs(res["on_hull"]), gap=rs(res["gap"]),
        spinodal=rs(spin), n_phases=rs(n_ph), K_A=rs(KA), K_B=rs(KB),
        label=rs(label), bound_A=rs(1 - X[:, 0]), bound_B=rs(1 - X[:, 1]),
        phi=rs(volume_fraction(rho.sum(axis=1), mix.diameters.mean() / 2)),
    )


LABELS = {0: "homogeneous", 1: "co-condensation", 2: "condensate excludes partner",
          3: "two immiscible condensates"}


# --------------------------------------------------------------------------- #
def partner_window(pd, cA_query_uM, want=(1, 2, 3)):
    r"""**Given [A], what range of [B] gives the observed structure?**

    Reads the phase diagram along the ``cA = cA_query`` column and returns the
    contiguous ``[B]`` intervals whose label is in ``want``.  This is the direct
    answer to "I measured protein A; what does that tell me about B?" -- and the
    honest form of that answer is an *interval*, because the phase boundary is a
    curve in the plane, not a point.

    Returns a list of ``(lo_uM, hi_uM)`` intervals (``hi = inf`` if the region
    runs off the top of the scanned grid).
    """
    ia = int(np.argmin(np.abs(pd["cA"] - cA_query_uM)))
    lab = pd["label"][ia]
    cB = pd["cB"]
    inside = np.isin(lab, list(want))
    out, start = [], None
    for j, v in enumerate(inside):
        if v and start is None:
            start = j
        elif not v and start is not None:
            out.append((float(cB[start]), float(cB[j - 1])))
            start = None
    if start is not None:
        out.append((float(cB[start]), np.inf if start < len(cB) - 1 else float(cB[-1])))
    return out


def describe_point(mix, cA_uM, cB_uM, pd=None):
    """Human-readable read-out at one composition."""
    rho = conc_to_density(np.array([[cA_uM, cB_uM]]))
    X = unbonded_fractions(rho, mix)[0]
    out = dict(cA_uM=cA_uM, cB_uM=cB_uM,
               bound_A=float(1 - X[0]), bound_B=float(1 - X[1]),
               phi=float(volume_fraction(rho.sum(), mix.diameters.mean() / 2)),
               spinodal=bool(spinodal_eig(rho, mix)[0] < 0))
    if pd is not None:
        ia = int(np.argmin(np.abs(pd["cA"] - cA_uM)))
        ib = int(np.argmin(np.abs(pd["cB"] - cB_uM)))
        out["label"] = LABELS[int(pd["label"][ia, ib])]
        out["K_A"] = float(pd["K_A"][ia, ib])
        out["K_B"] = float(pd["K_B"][ia, ib])
    return out
