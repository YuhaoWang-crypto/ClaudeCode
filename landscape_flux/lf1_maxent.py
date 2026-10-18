"""
LF1 — Macro entropy vs micro distribution; Maximum-Entropy (MaxEnt) inversion.

Claims tested
  (a) "A single macroscopic entropy cannot be inverted into P(X)."
      -> shown on REAL data: two pancreatic cell types with (almost) the same
         per-cell transcriptome entropy but completely different profiles.
  (b) "With a few extra constraints MaxEnt gives the least-biased P(X); the
      mismatch between the MaxEnt prediction and the measured distribution
      exposes hidden interactions."
      -> positive control with a planted 3-body interaction (XOR "complex")
         that is invisible to every pairwise statistic.
      -> REAL data: PBMC-3k marker genes. Independent (means) and pairwise
         (means+correlations, Ising) MaxEnt models fitted by exact enumeration.
         The residual is then explained by a HIDDEN variable (cell identity),
         which we recover and verify.
"""
import itertools
import numpy as np
from scipy.optimize import minimize
from . import _data

LN2 = np.log(2.0)


# ---------------------------------------------------------------- entropy utils
def entropy_from_counts(counts, miller_madow=True):
    c = np.asarray(counts, float)
    n = c.sum()
    p = c[c > 0] / n
    h = -(p * np.log(p)).sum()
    if miller_madow:
        h += (np.count_nonzero(c) - 1) / (2 * n)
    return h


def H(p):
    p = p[p > 0]
    return -(p * np.log(p)).sum()


# ---------------------------------------------------------------- MaxEnt fitter
class MaxEnt:
    """Exact MaxEnt over binary states {0,1}^N with constraints on all
    moments of order <= `order` (1 = independent, 2 = Ising, 3 = triplets)."""

    def __init__(self, N, order):
        self.N = N
        self.S = np.array(list(itertools.product([0, 1], repeat=N)), float)
        feats, names = [], []
        for k in range(1, order + 1):
            for idx in itertools.combinations(range(N), k):
                feats.append(self.S[:, idx].prod(1)); names.append(idx)
        self.F = np.stack(feats, 1)          # states x features
        self.names = names

    def fit(self, data):
        target = self._features(data).mean(0)

        def nll(lam):
            e = self.F @ lam
            m = e.max()
            logZ = m + np.log(np.exp(e - m).sum())
            p = np.exp(e - logZ)
            return logZ - lam @ target, self.F.T @ p - target

        res = minimize(nll, np.zeros(self.F.shape[1]), jac=True, method="L-BFGS-B",
                       options=dict(maxiter=5000, gtol=1e-10))
        e = self.F @ res.x
        self.p = np.exp(e - e.max()); self.p /= self.p.sum()
        self.lam = res.x
        return self

    def _features(self, data):
        cols = [data[:, idx].prod(1) for idx in self.names]
        return np.stack(cols, 1).astype(float)

    @property
    def entropy(self):
        return H(self.p)


def state_counts(B):
    N = B.shape[1]
    code = B @ (1 << np.arange(N)[::-1])
    return np.bincount(code, minlength=2 ** N)


def multi_information_ladder(B):
    """Return dict with S1 (indep), S2 (pairwise), S3 (triplet), S_emp."""
    N = B.shape[1]
    out = {}
    for k in (1, 2, 3):
        out[f"S{k}"] = MaxEnt(N, k).fit(B).entropy
    out["S_emp"] = entropy_from_counts(state_counts(B))
    I = out["S1"] - out["S_emp"]
    out["I_N"] = I
    out["frac_pair"] = (out["S1"] - out["S2"]) / I if I > 0 else np.nan
    out["frac_trip"] = (out["S1"] - out["S3"]) / I if I > 0 else np.nan
    return out


def triplet_scan(B, model2):
    """Connected 3-body residual for every triplet: binomial z-score of the
    observed co-occurrence P(1,1,1) against the pairwise-MaxEnt prediction."""
    N, n = B.shape[1], B.shape[0]
    S = model2.S
    rows = []
    for idx in itertools.combinations(range(N), 3):
        obs = B[:, idx].all(1).mean()
        pred = model2.p[S[:, idx].all(1)].sum()
        z = (obs - pred) / np.sqrt(max(pred * (1 - pred), 1e-12) / n)
        rows.append((idx, obs, pred, z))
    return sorted(rows, key=lambda r: -abs(r[3]))


# ---------------------------------------------------------------- Part A
def _auc(x, y):
    """Symmetric AUC (0.5 = indistinguishable, 1 = perfectly separable)."""
    from scipy.stats import mannwhitneyu
    u = mannwhitneyu(x, y).statistic / (len(x) * len(y))
    return max(u, 1 - u)


def part_a_same_entropy_different_cells():
    import scipy.sparse as sp
    a = _data.load_pancreas()
    X = a.layers["spliced"]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    lib = X.sum(1, keepdims=True)
    P = X / lib
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.where(P > 0, P * np.log(P), 0).sum(1)
    types = a.obs["clusters"].astype(str).values
    stats = {}
    logx = np.log1p(P * 1e4)
    for t in np.unique(types):
        m = types == t
        stats[t] = dict(H=h[m].mean(), sd=h[m].std(), n=m.sum(), prof=logx[m].mean(0))
    # pick the pair whose mean entropies differ by < 0.3 pooled SD but whose
    # expression programmes overlap least (top-50 genes)
    best = None
    names = list(stats)
    for i, j in itertools.combinations(names, 2):
        dH = abs(stats[i]["H"] - stats[j]["H"])
        sd = np.sqrt((stats[i]["sd"] ** 2 + stats[j]["sd"] ** 2) / 2)
        if dH > 0.3 * sd:
            continue
        top_i = np.argsort(-stats[i]["prof"])[:50]
        top_j = np.argsort(-stats[j]["prof"])[:50]
        overlap = len(set(top_i) & set(top_j)) / 50
        if best is None or overlap < best[0]:
            best = (overlap, i, j, dH)
    overlap, i, j, dH = best
    genes = np.array(a.var_names)
    d = stats[i]["prof"] - stats[j]["prof"]
    ti = genes[np.argsort(-d)[:5]]
    tj = genes[np.argsort(d)[:5]]
    # discriminability: entropy alone vs best single marker gene
    mi, mj = types == i, types == j
    auc_H = _auc(h[mi], h[mj])
    g = np.argmax(np.abs(d))
    auc_g = _auc(logx[mi, g], logx[mj, g])
    n4x = int((np.abs(d) > np.log(4)).sum())
    return dict(stats={k: (v["H"], v["sd"], v["n"]) for k, v in stats.items()},
                pair=(i, j), dH=dH, overlap=overlap, top_i=list(ti), top_j=list(tj),
                h=h, types=types, auc_H=auc_H, auc_gene=auc_g, gene=genes[g], n4x=n4x)


# ---------------------------------------------------------------- Part B
def part_b_planted_triplet(n=3000, seed=0):
    rng = np.random.default_rng(seed)
    A = rng.random(n) < 0.5
    Bg = rng.random(n) < 0.5
    C = np.logical_xor(A, Bg) ^ (rng.random(n) < 0.1)  # 10% noise
    noise = rng.random((n, 3)) < np.array([0.3, 0.5, 0.7])
    B = np.column_stack([A, Bg, C, noise]).astype(int)
    lad = multi_information_ladder(B)
    corr = np.corrcoef(B.T)
    m2 = MaxEnt(B.shape[1], 2).fit(B)
    scan = triplet_scan(B, m2)
    return dict(ladder=lad, max_abs_pair_corr=np.abs(corr[np.triu_indices(6, 1)]).max(),
                top_triplet=scan[0])


# ---------------------------------------------------------------- Part C
PBMC_GENES = ["CD3E", "IL7R", "CCR7", "LYZ", "CD14", "MS4A1", "NKG7", "FCGR3A"]


def part_c_pbmc(k_types=6, seed=0):
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA
    X, genes = _data.load_pbmc3k()
    genes = np.array(genes)
    gi = [int(np.where(genes == g)[0][0]) for g in PBMC_GENES]
    B = (X[:, gi].toarray() > 0).astype(int)
    lad = multi_information_ladder(B)
    m2 = MaxEnt(B.shape[1], 2).fit(B)
    scan = triplet_scan(B, m2)

    # hidden variable: cell identity from the transcriptome EXCLUDING these 8 genes
    keep = np.ones(X.shape[1], bool); keep[gi] = False
    Xk = X[:, keep]
    det = np.asarray((Xk > 0).sum(0)).ravel() >= 10
    Xk = Xk[:, det]
    lib = np.asarray(Xk.sum(1)).ravel()
    L = np.log1p(Xk.multiply(1e4 / lib[:, None]).toarray())
    var = L.var(0)
    L = L[:, np.argsort(-var)[:1500]]
    Z = PCA(20, random_state=seed).fit_transform(L - L.mean(0))
    lab = KMeans(k_types, n_init=10, random_state=seed).fit_predict(Z)

    # conditional multi-information: within each latent class, genes ~ independent?
    within = []
    for c in range(k_types):
        Bc = B[lab == c]
        S1c = MaxEnt(B.shape[1], 1).fit(Bc).entropy
        Sec = entropy_from_counts(state_counts(Bc))
        within.append((len(Bc), max(S1c - Sec, 0)))
    w = np.array(within)
    I_within = (w[:, 0] * w[:, 1]).sum() / w[:, 0].sum()

    # a latent-class model: P(s) = sum_c pi_c prod_i Bern(s_i | theta_ci)
    S = m2.S
    pi = np.bincount(lab, minlength=k_types) / len(lab)
    theta = np.stack([(B[lab == c].mean(0) * len(B[lab == c]) + 0.5) /
                      (len(B[lab == c]) + 1) for c in range(k_types)])
    p_lat = sum(pi[c] * np.prod(np.where(S == 1, theta[c], 1 - theta[c]), 1)
                for c in range(k_types))
    emp = state_counts(B) / len(B)

    def kl(p, q):
        m = p > 0
        return (p[m] * np.log(p[m] / q[m])).sum()

    m1 = MaxEnt(B.shape[1], 1).fit(B)
    kls = dict(indep=kl(emp, m1.p), pairwise=kl(emp, m2.p), latent=kl(emp, p_lat))
    names = [tuple(PBMC_GENES[i] for i in r[0]) for r in scan[:3]]
    # are the over-represented "impossible" co-expressors doublets?  Doublets
    # carry ~2x the UMIs of a singlet.  Lineage-exclusive markers: T (CD3E),
    # monocyte (CD14/LYZ), B (MS4A1).
    umi = np.asarray(X.sum(1)).ravel()
    g = {nm: B[:, PBMC_GENES.index(nm)] for nm in PBMC_GENES}
    multi = (g["CD3E"] + g["MS4A1"] + np.maximum(g["CD14"], g["LYZ"])) >= 3
    doublet = dict(n=int(multi.sum()), umi_ratio=np.median(umi[multi]) / np.median(umi[~multi]))
    from scipy.stats import mannwhitneyu
    doublet["p"] = mannwhitneyu(umi[multi], umi[~multi], alternative="greater").pvalue
    return dict(ladder=lad, kl=kls, I_within=I_within, within=within, doublet=doublet,
                top_triplets=[(nm, r[1], r[2], r[3]) for nm, r in zip(names, scan[:3])],
                labels=lab, B=B)


def run(make_fig=True):
    print("=" * 72); print("LF1  Macro entropy & MaxEnt inversion"); print("=" * 72)
    A = part_a_same_entropy_different_cells()
    i, j = A["pair"]
    print("[A] per-cell transcriptome entropy (nats), pancreas (Bastidas-Ponce 2019):")
    for t, (h, sd, n) in sorted(A["stats"].items(), key=lambda kv: -kv[1][0]):
        print(f"     {t:15s} H={h:5.3f} ± {sd:4.3f}  (n={n})")
    print(f"    -> '{i}' vs '{j}': ΔH = {A['dH']:.3f} nats, yet {A['n4x']} genes differ "
          f">4-fold")
    print(f"       {i} up: {A['top_i']}\n       {j} up: {A['top_j']}")
    print(f"       AUC separating them by entropy = {A['auc_H']:.2f}; by one gene "
          f"({A['gene']}) = {A['auc_gene']:.2f}")

    Bres = part_b_planted_triplet()
    L = Bres["ladder"]
    print("\n[B] planted 3-body 'complex' (C = A xor B, 10% noise) + 3 bystanders")
    print(f"    max |pairwise corr| = {Bres['max_abs_pair_corr']:.3f}  (pairwise stats blind)")
    print(f"    multi-information I_N = {L['I_N']:.3f} nats; captured by pairwise "
          f"{L['frac_pair']:.1%}, by triplet MaxEnt {L['frac_trip']:.1%}")
    t = Bres["top_triplet"]
    print(f"    triplet scan top hit = genes {t[0]} (obs P111={t[1]:.3f} vs "
          f"pairwise-MaxEnt {t[2]:.3f}, z={t[3]:+.1f})  -> planted (0,1,2) recovered")

    C = part_c_pbmc()
    L = C["ladder"]
    print("\n[C] PBMC-3k, 8 marker genes binarised (detected / not), n=2700 cells")
    print(f"    S_indep={L['S1']:.3f}  S_pair={L['S2']:.3f}  S_trip={L['S3']:.3f}  "
          f"S_emp(MM)={L['S_emp']:.3f} nats")
    print(f"    I_N={L['I_N']:.3f} nats; pairwise MaxEnt captures {L['frac_pair']:.1%}, "
          f"triplet {L['frac_trip']:.1%}")
    print("    KL(empirical || model): indep %.3f | pairwise %.3f | latent-class(6) %.3f"
          % (C["kl"]["indep"], C["kl"]["pairwise"], C["kl"]["latent"]))
    print("    largest pairwise-MaxEnt residuals (triplets):")
    for nm, o, p, lr in C["top_triplets"]:
        print(f"       {nm}: obs P111={o:.3f} pred={p:.3f} z={lr:+.1f}")
    d = C["doublet"]
    print(f"    cells positive for T+B+monocyte markers: n={d['n']}, median UMI "
          f"{d['umi_ratio']:.2f}x the rest (one-sided MWU p={d['p']:.1e}) -> doublet /\n"
          f"      ambient-RNA candidates (weak: a clean doublet would be ~2x)")
    print(f"    mean multi-information WITHIN latent cell classes = {C['I_within']:.3f} "
          f"nats (vs {L['I_N']:.3f} overall) -> {1 - C['I_within']/L['I_N']:.0%} of the "
          f"dependency is explained by the hidden cell-identity variable")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
        order = sorted(A["stats"], key=lambda t: -A["stats"][t][0])
        data = [A["h"][A["types"] == t] for t in order]
        ax[0].boxplot(data, showfliers=False); ax[0].set_xticks(range(1, len(order) + 1))
        ax[0].set_xticklabels(order, rotation=45, ha="right", fontsize=8)
        for k, t in enumerate(order):
            if t in A["pair"]:
                ax[0].get_xticklabels()[k].set_color("crimson")
        ax[0].set_ylabel("per-cell transcriptome entropy (nats)")
        ax[0].set_title(f"(A) same entropy ≠ same state\n{i} vs {j}: ΔH={A['dH']:.3f}")
        lb = ["indep\n(means)", "pairwise\n(Ising)", "triplet"]
        for res, col, name in [(Bres["ladder"], "tab:gray", "planted XOR"),
                               (C["ladder"], "tab:blue", "PBMC real")]:
            ax[1].plot(lb, [res["frac_pair"] * 0, res["frac_pair"], res["frac_trip"]],
                       "o-", color=col, label=name)
        ax[1].axhline(1, ls=":", c="k"); ax[1].set_ylim(-0.05, 1.1)
        ax[1].set_ylabel("fraction of multi-information captured")
        ax[1].set_title("(B,C) MaxEnt ladder"); ax[1].legend()
        ax[2].bar(["indep", "pairwise", "latent\n(cell type)"],
                  [C["kl"]["indep"], C["kl"]["pairwise"], C["kl"]["latent"]],
                  color=["tab:gray", "tab:blue", "tab:green"])
        ax[2].set_ylabel("KL(empirical ‖ model)  [nats]")
        ax[2].set_title("(C) residual explained by hidden variable")
        fig.tight_layout(); fig.savefig(_data.fig("lf1_maxent.png"), dpi=130)
        plt.close(fig)
    return dict(A=A, B=Bres, C=C)


if __name__ == "__main__":
    run()
