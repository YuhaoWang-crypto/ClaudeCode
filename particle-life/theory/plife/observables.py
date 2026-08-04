r"""
Order parameters that turn "it looks like snakes" into numbers.
===============================================================

Everything here is computed from a snapshot ``(pos, vel, types)`` on a torus.
The seven quantities below are enough to separate every morphology the platform
produces, and each one is tied to a specific term of the theory:

======================  =========================================================
observable              what it measures / which theory object it tests
======================  =========================================================
``speed``               mean :math:`|v|`; the *motility* order parameter.
                        Zero (up to the drag scale) iff the dynamics is
                        variational -> tests the reciprocity theorem.
``cluster_speed``       speed of cluster centres of mass; isolates *coherent*
                        self-propulsion from internal jiggle.  Predicted by
                        ``twobody.cluster_internal_force``.
``spin``                mean cluster angular velocity; predicted by
                        ``twobody.cluster_internal_torque``.
``k_peak``              peak of the static structure factor :math:`S(k)`;
                        compared against :math:`k^\*` from ``stability.scan``.
``psi6``                hexatic bond-orientational order; 1 = crystal, 0 = fluid.
``dim2``                correlation (box-counting) dimension of the point set on
                        scales :math:`[\alpha, 3\alpha]`.  ~1 = filaments/snakes,
                        ~2 = blobs/sheets.
``n_clusters``,         cluster count and size distribution -> is the structure
``mean_cluster``        microphase separated (many finite blobs) or coarsening?
======================  =========================================================
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

__all__ = [
    "clusters", "cluster_kinematics", "structure_factor", "pair_correlation",
    "psi6", "correlation_dimension", "observables", "classify",
]


# --------------------------------------------------------------------------- #
def _tree(pos, box):
    p = np.mod(pos, box)
    return cKDTree(p, boxsize=box)


def clusters(pos, box, cutoff):
    """Connected components of the ``cutoff``-proximity graph.  Returns labels."""
    tree = _tree(pos, box)
    pairs = tree.query_pairs(cutoff, output_type="ndarray")
    n = len(pos)
    if pairs.size == 0:
        return np.arange(n), n
    g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    ncomp, labels = connected_components(g, directed=False)
    return labels, ncomp


def cluster_kinematics(pos, vel, box, cutoff, min_size=4):
    r"""Per-cluster centre-of-mass speed and angular velocity.

    Positions are unwrapped relative to each cluster's first member before
    averaging so that clusters straddling the periodic seam are handled.
    """
    labels, ncomp = clusters(pos, box, cutoff)
    box = np.asarray(box, float)
    speeds, spins, sizes = [], [], []
    order = np.argsort(labels, kind="stable")
    lab_sorted = labels[order]
    bounds = np.searchsorted(lab_sorted, np.arange(ncomp + 1))
    for c in range(ncomp):
        idx = order[bounds[c]:bounds[c + 1]]
        m = len(idx)
        if m < min_size:
            continue
        p = pos[idx]
        ref = p[0]
        d = p - ref
        d -= box * np.round(d / box)                  # unwrap
        p = ref + d
        X = p.mean(axis=0)
        V = vel[idx].mean(axis=0)
        rel = p - X
        vrel = vel[idx] - V
        Lz = (rel[:, 0] * vrel[:, 1] - rel[:, 1] * vrel[:, 0]).sum()
        I = (rel ** 2).sum()
        speeds.append(np.hypot(*V))
        spins.append(Lz / I if I > 1e-12 else 0.0)
        sizes.append(m)
    if not sizes:
        return dict(n_clusters=0, mean_cluster=0.0, cluster_speed=0.0,
                    spin=0.0, abs_spin=0.0, max_cluster=0)
    sizes = np.array(sizes, float)
    return dict(
        n_clusters=len(sizes),
        mean_cluster=float(sizes.mean()),
        max_cluster=int(sizes.max()),
        cluster_speed=float(np.average(speeds, weights=sizes)),
        spin=float(np.average(spins, weights=sizes)),
        abs_spin=float(np.average(np.abs(spins), weights=sizes)),
    )


# --------------------------------------------------------------------------- #
def structure_factor(pos, box, ngrid=256, kmax=None, types=None, which=None):
    r"""Radially averaged static structure factor :math:`S(k)` via a density FFT.

    ``which=(a,b)`` gives the partial structure factor
    :math:`S_{ab}(k)=\langle \delta\hat\rho_a \delta\hat\rho_b^*\rangle`, whose
    *sign* distinguishes condensation (in-phase, :math:`S_{ab}>0`) from demixing
    (anti-phase, :math:`S_{ab}<0`).
    """
    box = np.asarray(box, float)
    p = np.mod(pos, box)

    def grid_of(sel):
        H, _, _ = np.histogram2d(p[sel, 0], p[sel, 1], bins=ngrid,
                                 range=[[0, box[0]], [0, box[1]]])
        return H - H.mean()

    if which is None:
        f1 = np.fft.rfft2(grid_of(slice(None)))
        power = np.abs(f1) ** 2
        norm = len(pos)
    else:
        a, b = which
        f1 = np.fft.rfft2(grid_of(types == a))
        f2 = np.fft.rfft2(grid_of(types == b))
        power = (f1 * np.conj(f2)).real
        norm = np.sqrt(max(1, (types == a).sum()) * max(1, (types == b).sum()))

    kx = 2 * np.pi * np.fft.fftfreq(ngrid, d=box[0] / ngrid)
    ky = 2 * np.pi * np.fft.rfftfreq(ngrid, d=box[1] / ngrid)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    K = np.hypot(KX, KY)

    if kmax is None:
        kmax = 0.5 * K.max()
    nb = 120
    edges = np.linspace(0, kmax, nb + 1)
    which_bin = np.digitize(K.ravel(), edges) - 1
    valid = (which_bin >= 0) & (which_bin < nb)
    s = np.bincount(which_bin[valid], power.ravel()[valid], minlength=nb)
    c = np.bincount(which_bin[valid], minlength=nb)
    with np.errstate(invalid="ignore", divide="ignore"):
        S = np.where(c > 0, s / np.maximum(c, 1) / norm, np.nan)
    kcen = 0.5 * (edges[1:] + edges[:-1])
    return kcen, S


def peak_wavenumber(pos, box, kmin=0.005, **kw):
    k, S = structure_factor(pos, box, **kw)
    m = (k > kmin) & np.isfinite(S)
    if not m.any():
        return np.nan
    return float(k[m][np.nanargmax(S[m])])


# --------------------------------------------------------------------------- #
def pair_correlation(pos, box, rmax, nbins=120, types=None, which=None):
    r"""Radial distribution function :math:`g(r)` (or partial :math:`g_{ab}`)."""
    box = np.asarray(box, float)
    area = box[0] * box[1]
    if which is None:
        pa = pb = np.mod(pos, box)
        na = nb = len(pos)
        same = True
    else:
        a, b = which
        pa = np.mod(pos[types == a], box)
        pb = np.mod(pos[types == b], box)
        na, nb = len(pa), len(pb)
        same = a == b
    if na == 0 or nb == 0:
        return np.linspace(0, rmax, nbins), np.zeros(nbins)
    ta, tb = cKDTree(pa, boxsize=box), cKDTree(pb, boxsize=box)
    edges = np.linspace(0, rmax, nbins + 1)
    counts = ta.count_neighbors(tb, edges, cumulative=False)[1:]
    if same:
        counts = counts.astype(float)
        counts[0] -= na                      # remove self pairs from the first bin
    shell = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    ideal = na * nb / area * shell
    r = 0.5 * (edges[1:] + edges[:-1])
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(ideal > 0, counts / ideal, np.nan)
    return r, g


# --------------------------------------------------------------------------- #
def psi6(pos, box, cutoff, kneigh=6):
    r"""Hexatic bond-orientational order :math:`|\langle\psi_6\rangle|`,
    :math:`\psi_6^{(i)} = \frac1{n_i}\sum_{j} e^{6 i \theta_{ij}}`."""
    box = np.asarray(box, float)
    p = np.mod(pos, box)
    tree = cKDTree(p, boxsize=box)
    dd, ii = tree.query(p, k=kneigh + 1)
    psis = []
    for a in range(len(p)):
        nb = [(d, j) for d, j in zip(dd[a][1:], ii[a][1:]) if d < cutoff]
        if len(nb) < 3:
            continue
        z = 0j
        for _, j in nb:
            d = p[j] - p[a]
            d -= box * np.round(d / box)
            z += np.exp(6j * np.arctan2(d[1], d[0]))
        psis.append(z / len(nb))
    return float(abs(np.mean(psis))) if psis else 0.0


def correlation_dimension(pos, box, r_lo, r_hi, nr=14):
    r"""Slope of :math:`\log C(r)` vs :math:`\log r`, i.e. the correlation dimension.

    ~1 for filaments (snakes, chains), ~2 for space-filling blobs/sheets,
    between 1 and 2 for branched networks and membranes.
    """
    box = np.asarray(box, float)
    tree = _tree(pos, box)
    rs = np.logspace(np.log10(r_lo), np.log10(r_hi), nr)
    counts = np.array([tree.count_neighbors(tree, r) for r in rs], float)
    counts = counts - len(pos)                       # drop self pairs
    m = counts > 0
    if m.sum() < 4:
        return np.nan
    A = np.polyfit(np.log(rs[m]), np.log(counts[m]), 1)
    return float(A[0])


# --------------------------------------------------------------------------- #
def observables(sim, cutoff_scale=1.35):
    """The full diagnostic vector of a live simulation."""
    p = sim.p
    box = np.asarray(p.L, float)
    alpha_typ = float(np.median(p.alpha))
    cutoff = cutoff_scale * alpha_typ

    out = dict(
        t=sim.t,
        speed=sim.mean_speed(),
        potential=sim.potential_energy() / sim.N,
        alpha_typ=alpha_typ,
    )
    out.update(cluster_kinematics(sim.pos, sim.vel, box, cutoff))
    out["k_peak"] = peak_wavenumber(sim.pos, box, kmax=4 * np.pi / alpha_typ)
    out["length"] = 2 * np.pi / out["k_peak"] if out["k_peak"] > 0 else np.nan
    out["psi6"] = psi6(sim.pos, box, cutoff)
    out["dim2"] = correlation_dimension(sim.pos, box, alpha_typ, 3.5 * alpha_typ)
    # normalised motility: speed measured in "bond lengths per second"
    out["motility"] = out["cluster_speed"] / alpha_typ
    return out


def classify(obs, motility_thresh=0.02, dim_filament=1.55, psi6_cryst=0.42,
             gas_cluster=6.0):
    """Coarse morphology label from the observable vector.

    The thresholds are empirical dividing lines, but each axis is physical:
    ``motility`` separates variational from non-variational dynamics, ``dim2``
    separates filaments from blobs, ``psi6`` separates crystal from fluid, and
    ``mean_cluster`` separates a gas from condensed matter.
    """
    if obs["n_clusters"] == 0 or obs["mean_cluster"] < gas_cluster:
        return "gas"
    moving = obs["motility"] > motility_thresh
    spinning = abs(obs["spin"]) > 0.15
    if obs["dim2"] < dim_filament:
        return "crawling filaments" if moving else "static chains / network"
    if moving:
        return "rotating cells" if spinning else "migrating cells / chasers"
    return "crystalline clusters" if obs["psi6"] > psi6_cryst else "static blobs"
