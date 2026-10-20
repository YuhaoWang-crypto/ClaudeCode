"""Select a gene-pair panel per cell type (numpy; used at build time only).

Greedy k-TSP-style selection with disjoint genes (each gene used once, so a
k-pair panel needs exactly 2k genes). A pair's value for separating classes u
and v is the total-variation distance between its ternary symbol distributions
under u and v (see pair_code.py), averaged over response strengths 0.5 and 1.0
so pairs that only work at full strength are not favoured. For two stimuli the
second one may take any strength, so a pair must differ in shape, not only size.
The panel maximises sum(log(0.1 + separation)) over class pairs: a soft minimum
that keeps every class pair covered without stalling on an inseparable one.
"""
import itertools
import numpy as np

from .pair_code import STRENGTH_GRID

SELECT_STRENGTHS = (0.5, 1.0)


def phi_np(x):
    """Normal CDF via erf (Abramowitz-Stegun 7.1.26, |error| < 1.5e-7)."""
    z = np.asarray(x, float) / np.sqrt(2)
    s = np.sign(z)
    z = np.abs(z)
    t = 1 / (1 + 0.3275911 * z)
    y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * np.exp(-z * z)
    return 0.5 * (1 + s * y)


def probs_np(d, sigma, tau):
    s2 = sigma * np.sqrt(2)
    p_pos = 1 - phi_np((tau - d) / s2)
    p_neg = phi_np((-tau - d) / s2)
    p_zero = np.clip(1 - p_pos - p_neg, 0, None)
    return p_pos, p_zero, p_neg


def _tv(du, dv, sigma, tau):
    pu, pv = probs_np(du, sigma, tau), probs_np(dv, sigma, tau)
    return 0.5 * sum(np.abs(x - y) for x, y in zip(pu, pv))


def select_panel(genes, a, g, classes, k, sigma, tau, strengths=SELECT_STRENGTHS):
    """classes: identifiable stimuli, e.g. ["IFNA", "IFNG"] or ["IFNA"]."""
    genes = list(genes)
    a = np.asarray(a, float)
    g = np.asarray(g, float)
    ii, jj = np.triu_indices(len(genes), 1)
    da, dg = a[ii] - a[jj], g[ii] - g[jj]
    d = {"none": np.zeros_like(da), "IFNA": da, "IFNG": dg}
    labels = ["none"] + list(classes)
    pairs = list(itertools.combinations(labels, 2))
    if not pairs:
        return []
    def sep(u, v):
        if u == "none":
            return np.mean([_tv(0 * d[v], s * d[v], sigma, tau) for s in strengths], axis=0)
        # stimulus vs stimulus: v may appear at any strength, so take the worst case over the grid.
        # This favours pairs that differ in shape, not just in size.
        return np.mean([np.min([_tv(s * d[u], t * d[v], sigma, tau) for t in STRENGTH_GRID], axis=0)
                        for s in strengths], axis=0)
    tv = np.stack([np.minimum(sep(u, v), sep(v, u)) if u != "none" else sep(u, v) for u, v in pairs])
    acc = np.zeros(len(pairs))
    used = np.zeros(len(genes), bool)
    panel = []
    for rank in range(k):
        free = ~(used[ii] | used[jj])
        if not free.any():
            break
        tot = acc[:, None] + tv
        # soft minimum: balanced across class pairs, but one inseparable class pair
        # (value stuck near 0) cannot freeze the choice the way a hard minimum does
        obj = np.log(0.1 + tot).sum(0)
        obj[~free] = -np.inf
        p = int(np.argmax(obj))
        acc += tv[:, p]
        used[ii[p]] = used[jj[p]] = True
        i, j = ii[p], jj[p]
        if da[p] < 0:  # orient so the IFNa difference is non-negative (cosmetic)
            i, j = j, i
        panel.append({"rank": rank + 1, "gene_i": genes[i], "gene_j": genes[j],
                      "d_ifna": float(a[i] - a[j]), "d_ifng": float(g[i] - g[j])})
    return panel


def loglik_np(panel, B, gene_index, sigma, tau, d_key, strengths):
    """Vectorised twin of pair_code.class_loglik for many samples.

    B: samples x genes (NaN = not measured). Returns best log-likelihood per sample.
    """
    i = np.array([gene_index[p["gene_i"]] for p in panel])
    j = np.array([gene_index[p["gene_j"]] for p in panel])
    diff = B[:, i] - B[:, j]
    sym = np.where(diff > tau, 0, np.where(diff < -tau, 2, 1))  # 0:+, 1:0, 2:-
    measured = ~np.isnan(diff)
    d = np.array([p[d_key] for p in panel])
    best = np.full(B.shape[0], -np.inf)
    for s in strengths:
        pp, pz, pn = probs_np(s * d, sigma, tau)
        logp = np.log(np.stack([pp, pz, pn]) + 1e-9)  # 3 x k
        ll = np.where(measured, logp[sym, np.arange(len(panel))], 0.0).sum(1)
        best = np.maximum(best, ll)
    return best
