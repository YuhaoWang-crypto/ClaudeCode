"""Markov-state model over the accumulated Gen-COMPAS sampling.

The whole point of Gen-COMPAS is that every trajectory it collects is
*unbiased*: short runs started from many different places.  That is exactly
the input a Markov-state model needs, so the equilibrium distribution -- and
therefore the free-energy landscape -- can be recovered from it without any
reweighting of a biasing potential.

Clustering is done in the same CV-free feature space the networks use; the
projection onto (phi, psi) happens only at the plotting stage.
"""

import numpy as np
from scipy.cluster.vq import kmeans2


def cluster(features, n_clusters=150, seed=0, iters=60):
    rng = np.random.default_rng(seed)
    init = features[rng.choice(len(features), n_clusters, replace=False)]
    centers, labels = kmeans2(features, init, minit="matrix", iter=iters)
    # drop empty clusters and relabel
    used, labels = np.unique(labels, return_inverse=True)
    return centers[used], labels


def count_matrix(labels, traj_lengths, lag):
    n = labels.max() + 1
    c = np.zeros((n, n))
    start = 0
    for L in traj_lengths:
        seg = labels[start:start + L]
        start += L
        if L > lag:
            np.add.at(c, (seg[:-lag], seg[lag:]), 1.0)
    return c


def largest_connected_set(c):
    """Indices of the largest strongly connected component of the count graph."""
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import csr_matrix

    adj = csr_matrix((c > 0).astype(int))
    n_comp, lab = connected_components(adj, directed=True, connection="strong")
    if n_comp == 1:
        return np.arange(c.shape[0])
    sizes = np.bincount(lab)
    return np.where(lab == sizes.argmax())[0]


def reversible_mle(c, max_iter=3000, tol=1e-10):
    """Maximum-likelihood reversible transition matrix (Prinz et al. 2011).

    Returns (T, pi).  The reversible estimator is used rather than simple
    count symmetrisation because the trajectories start far from equilibrium.
    """
    c = np.asarray(c, dtype=float)
    n = c.shape[0]
    ci = c.sum(1)
    x = c + c.T                      # current estimate of unnormalised x_ij
    x /= x.sum()
    for _ in range(max_iter):
        xi = x.sum(1)
        x_new = np.zeros_like(x)
        denom = ci[:, None] / np.where(xi[:, None] == 0, 1, xi[:, None]) + \
            ci[None, :] / np.where(xi[None, :] == 0, 1, xi[None, :])
        num = c + c.T
        with np.errstate(divide="ignore", invalid="ignore"):
            x_new = np.where(denom > 0, num / denom, 0.0)
        s = x_new.sum()
        if s <= 0:
            break
        x_new /= s
        if np.abs(x_new - x).max() < tol:
            x = x_new
            break
        x = x_new
    pi = x.sum(1)
    pi = pi / pi.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(x.sum(1, keepdims=True) > 0,
                     x / x.sum(1, keepdims=True), 0.0)
    return t, pi


def frame_weights(labels, pi):
    """Equilibrium weight of each frame: pi(cluster) / n_frames_in_cluster."""
    counts = np.bincount(labels, minlength=len(pi)).astype(float)
    w = np.zeros(len(labels))
    good = counts > 0
    per = np.zeros(len(pi))
    per[good] = pi[good] / counts[good]
    w = per[labels]
    return w / w.sum()


def mfpt(t, pi, source, target, lag_ps):
    """Mean first passage time from a set of states to a target set.

    Solves (I - T_rest) m = tau on the states outside the target, then
    averages over the source set with equilibrium weights.
    """
    n = t.shape[0]
    rest = np.setdiff1d(np.arange(n), target)
    if len(rest) == 0 or len(target) == 0:
        return float("nan")
    a = np.eye(len(rest)) - t[np.ix_(rest, rest)]
    try:
        m = np.linalg.solve(a, np.full(len(rest), lag_ps))
    except np.linalg.LinAlgError:
        return float("nan")
    full = np.zeros(n)
    full[rest] = m
    src = np.intersect1d(source, rest)
    if len(src) == 0:
        return float("nan")
    w = pi[src]
    if w.sum() <= 0:
        return float(full[src].mean())
    return float((full[src] * w).sum() / w.sum())


def implied_timescales(c, lags, labels, traj_lengths, n_ts=3):
    out = []
    for lag in lags:
        cm = count_matrix(labels, traj_lengths, lag)
        keep = largest_connected_set(cm)
        cm = cm[np.ix_(keep, keep)]
        t, pi = reversible_mle(cm)
        ev = np.linalg.eigvals(t)
        ev = np.sort(np.abs(ev))[::-1][1:1 + n_ts]
        ev = np.clip(ev, 1e-12, 1 - 1e-12)
        out.append(-lag / np.log(ev))
    return np.array(out)
