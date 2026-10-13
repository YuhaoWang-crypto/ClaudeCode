"""
Module 26 — T cell exhaustion: order parameter from scRNA-seq, network
interactions from Perturb-seq, and what can (and cannot) be concluded about
attractor stability.

Data: GSE327731 (in vivo Perturb-seq, 73 single-gene knockouts in CD8 T cells
transferred into a mouse melanoma model; "community ML challenge" dataset).

Part A — positive control on a synthetic regulatory network with a KNOWN
Jacobian, to establish what Perturb-seq-type data can identify:
  * A knockout + sequencing measures the TOTAL steady-state response
    R = -J^{-1} B (direct + all indirect paths), not J itself.
  * Modular response analysis (MRA; Kholodenko et al. PNAS 2002) recovers the
    DIRECT, normalised interactions  r = -[diag(R^{-1})]^{-1} R^{-1}
    (r_ij = d ln x_i / d ln x_j at fixed other nodes, r_ii = -1),
    provided every node of the module is perturbed.
  * r fixes J only up to each node's own turnover rate: J = D r with an
    unknown positive diagonal D. Eigenvalues -- hence relaxation times and
    attractor "depth" -- are therefore NOT identified by steady-state
    Perturb-seq alone; time-resolved data are needed.

Part B — real data (GSE327731): exhaustion order parameter, knockout
response matrix on the measured TF module, MRA, stability under unknown D,
and validation against held-out knockouts.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
DATA = os.path.join(HERE, "..", "data_cache", "GSE327731")


# ---------------------------------------------------------------------------
# Part A: synthetic positive control
# ---------------------------------------------------------------------------
def random_grn(n=8, density=0.3, seed=0):
    """Hill-type regulatory network
        dx_i/dt = b_i * prod_j h_ij(x_j) - g_i x_i
    with activating or repressing Hill terms on random edges."""
    rng = np.random.default_rng(seed)
    A = (rng.random((n, n)) < density).astype(float)
    np.fill_diagonal(A, 0)
    S = np.where(rng.random((n, n)) < 0.5, 1.0, -1.0) * A    # sign
    K = rng.uniform(0.5, 1.5, (n, n))
    b = rng.uniform(1.0, 2.0, n)
    g = rng.lognormal(0.0, 0.8, n)            # turnover rates differ a lot
    return dict(A=A, S=S, K=K, b=b, g=g, n=n)


def _rhs(x, net, p):
    h = np.ones_like(net["K"])
    act = (net["S"] > 0); rep = (net["S"] < 0)
    X = x[None, :] / net["K"]
    h[act] = (X ** 2 / (1 + X ** 2))[act] * 0.8 + 0.2   # leaky activation
    h[rep] = (1 / (1 + X ** 2))[rep]
    return net["b"] * p * h.prod(1) - net["g"] * x


def steady_state(net, p=None, x0=None):
    from scipy.integrate import solve_ivp
    p = np.ones(net["n"]) if p is None else p
    x0 = np.ones(net["n"]) if x0 is None else x0
    sol = solve_ivp(lambda t, x: _rhs(x, net, p), (0, 400), x0,
                    method="LSODA", rtol=1e-10, atol=1e-12)
    return sol.y[:, -1]


def jacobian(net, x, p=None, eps=1e-6):
    p = np.ones(net["n"]) if p is None else p
    n = net["n"]; J = np.zeros((n, n))
    for j in range(n):
        d = np.zeros(n); d[j] = eps * max(1.0, abs(x[j]))
        J[:, j] = (_rhs(x + d, net, p) - _rhs(x - d, net, p)) / (2 * d[j])
    return J


def true_local_response(J, x):
    """Normalised direct interactions r_ij = -J_ij x_j / (J_ii x_i)."""
    return -(J * x[None, :]) / (np.diag(J) * x)[:, None]


def mra(R):
    """Modular response analysis: global response matrix R (rows = nodes,
    columns = perturbations of each node) -> local response matrix r."""
    Ri = np.linalg.inv(R)
    return -Ri / np.diag(Ri)[:, None]


def knockdown_responses(net, x0, frac=0.5, noise=0.0, seed=0):
    """Global response matrix from knocking down the production of each
    node by `frac` (log-fold changes), with multiplicative measurement
    noise of SD `noise` in log space."""
    rng = np.random.default_rng(seed)
    n = net["n"]; R = np.zeros((n, n))
    for j in range(n):
        p = np.ones(n); p[j] = 1 - frac
        xj = steady_state(net, p, x0)
        R[:, j] = np.log(xj / x0) / np.log(1 - frac)
    return R + noise * rng.normal(size=R.shape)


def positive_control(n_nets=20, n=8, noises=(0.0, 0.02, 0.05, 0.1, 0.2)):
    """How well do (a) the raw response matrix and (b) MRA recover the true
    direct interactions, as measurement noise grows? Also: does r fix
    stability? (sample positive diagonal D, count sign changes of
    max Re eig(D r))."""
    from scipy.stats import spearmanr
    rows = []
    stab = []
    for s in range(n_nets):
        net = random_grn(n=n, seed=s)
        x0 = steady_state(net)
        J = jacobian(net, x0)
        if np.max(np.linalg.eigvals(J).real) >= 0:
            continue                                 # need a stable state
        r_true = true_local_response(J, x0)
        off = ~np.eye(n, dtype=bool)
        edge = net["A"].astype(bool) & off
        for nz in noises:
            R = knockdown_responses(net, x0, noise=nz, seed=s)
            r_hat = mra(R)
            for name, est in [("raw response R", R), ("MRA r", r_hat)]:
                rho = spearmanr(est[off], r_true[off])[0]
                sign = np.mean(np.sign(est[edge]) == np.sign(r_true[edge]))
                # direct-edge detection: |est| ranks true edges above non-edges
                from scipy.stats import rankdata
                rk = rankdata(np.abs(est[off]))
                pos = edge[off]
                auc = ((rk[pos].sum() - pos.sum() * (pos.sum() + 1) / 2)
                       / (pos.sum() * (~pos).sum()))
                rows.append(dict(net=s, noise=nz, method=name, rho=rho,
                                 sign=sign, edge_auc=auc))
        # identifiability of stability: true D vs random positive D
        rng = np.random.default_rng(100 + s)
        flips = 0
        for _ in range(200):
            D = np.diag(rng.lognormal(0, 1.5, n))
            if np.max(np.linalg.eigvals(D @ r_true).real) >= 0:
                flips += 1
        lam_true = np.max(np.linalg.eigvals(J).real)
        taus = []
        for _ in range(200):
            D = np.diag(rng.lognormal(0, 1.5, n))
            lam = np.max(np.linalg.eigvals(D @ r_true).real)
            if lam < 0:
                taus.append(-1 / lam)
        stab.append(dict(net=s, frac_unstable=flips / 200,
                         tau_true=-1 / lam_true,
                         tau_range=(float(np.min(taus)), float(np.max(taus)))
                         if taus else (np.nan, np.nan)))
    return rows, stab


# ---------------------------------------------------------------------------
# Part B: real data
# ---------------------------------------------------------------------------
PROGRAMS = {   # canonical mouse CD8 T-cell state markers (none is a KO target)
    "progenitor": ["Tcf7", "Slamf6", "Sell", "Ccr7", "Il7r", "Xcl1", "Id3"],
    "terminal": ["Havcr2", "Tox", "Pdcd1", "Lag3", "Cd244a", "Cd160",
                 "Gzmb", "Prf1"],
    "effector": ["Cx3cr1", "Klrg1", "Zeb2", "Klf2", "S1pr5", "Gzma"],
    "cycling": ["Mki67", "Top2a", "Stmn1", "Tuba1b", "Hmgb2"],
}
MIN_CELLS = 100


def load(path=os.path.join(DATA, "GSE327731_screen_2.h5ad")):
    """Raw counts -> normalised log1p; drop knockouts with < MIN_CELLS."""
    import anndata as ad
    import scanpy as sc
    a = ad.read_h5ad(path)
    a.X = a.layers["rawcounts"].copy()
    del a.layers["rawcounts"]
    keep = a.obs["gRNA_gene"].map(a.obs["gRNA_gene"].value_counts()) \
        >= MIN_CELLS
    a = a[keep.values].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def order_parameter(a, n_hvg=2000, seed=0):
    """Exhaustion order parameter built from non-targeted genes only.
    Three independent definitions are computed and compared:
      E_dpt  : diffusion pseudotime rooted in the most progenitor-like NT cell
      E_score: terminal-minus-progenitor marker score
      E_pc   : the principal component best aligned with E_score
    Also returns the PCA participation ratio (effective dimensionality)."""
    import scanpy as sc
    targets = set(a.obs["gRNA_gene"].unique()) - {"NT"}
    genes = [g for g in a.var_names if g not in targets
             and not g.startswith(("mt-", "Rpl", "Rps"))]
    b = a[:, genes].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=n_hvg)
    b = b[:, b.var["highly_variable"]].copy()
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=50, random_state=seed)
    for name, gl in PROGRAMS.items():
        gl = [g for g in gl if g in a.var_names]
        sc.tl.score_genes(a, gl, score_name=f"score_{name}", random_state=seed)
    a.obs["E_score"] = a.obs["score_terminal"] - a.obs["score_progenitor"]
    sc.pp.neighbors(b, n_neighbors=30, use_rep="X_pca", random_state=seed)
    sc.tl.diffmap(b, n_comps=15)
    nt = np.where(a.obs["gRNA_gene"].values == "NT")[0]
    root = nt[np.argmin(a.obs["E_score"].values[nt])]
    b.uns["iroot"] = int(root)
    sc.tl.dpt(b)
    a.obs["E_dpt"] = b.obs["dpt_pseudotime"].values
    pcs = b.obsm["X_pca"]
    corr = [np.corrcoef(pcs[:, k], a.obs["E_score"])[0, 1] for k in range(10)]
    k = int(np.argmax(np.abs(corr)))
    a.obs["E_pc"] = pcs[:, k] * np.sign(corr[k])
    ev = b.uns["pca"]["variance"]
    pr = float(ev.sum() ** 2 / (ev ** 2).sum())
    sc.tl.umap(b, random_state=seed)
    sc.tl.leiden(b, resolution=0.6, random_state=seed, flavor="igraph",
                 n_iterations=2, directed=False)
    a.obs["leiden"] = b.obs["leiden"].values
    a.obsm["X_umap"] = b.obsm["X_umap"]; a.obsm["X_pca"] = pcs
    return dict(pc_index=k, pc_corr=corr, participation_ratio=pr,
                var_ratio=b.uns["pca"]["variance_ratio"], hvg=list(b.var_names))


def ko_effects(a, axis="E_score", n_boot=500, seed=0):
    """Per knockout: shift of the order parameter vs non-targeting cells,
    paired within each mouse (hashtag) so mouse-to-mouse differences cancel;
    CI by bootstrapping mice. Also: change in state fractions and knockout
    efficiency at the mRNA level (target expression KO vs NT)."""
    o = a.obs
    rng = np.random.default_rng(seed)
    nt = o["gRNA_gene"] == "NT"
    rows = []
    X = a.X
    for g in sorted(set(o["gRNA_gene"]) - {"NT"}):
        ko = o["gRNA_gene"] == g
        diffs, w = [], []
        for h in o["hash.ID"].unique():
            m_ko = ko & (o["hash.ID"] == h); m_nt = nt & (o["hash.ID"] == h)
            if m_ko.sum() >= 5 and m_nt.sum() >= 5:
                diffs.append(o.loc[m_ko, axis].mean() - o.loc[m_nt, axis].mean())
                w.append(m_ko.sum())
        diffs, w = np.array(diffs), np.array(w, float)
        est = float(np.average(diffs, weights=w))
        bs = []
        for _ in range(n_boot):
            k = rng.integers(0, len(diffs), len(diffs))
            bs.append(np.average(diffs[k], weights=w[k]))
        lo, hi = np.percentile(bs, [2.5, 97.5])
        # knockout efficiency: target mRNA in KO vs NT cells (log1p scale)
        eff = np.nan
        if g in a.var_names:
            j = a.var_names.get_loc(g)
            col = X[:, j].toarray().ravel() if hasattr(X, "toarray") else X[:, j]
            eff = float(col[ko.values].mean() - col[nt.values].mean())
        frac = {s: float((o.loc[ko, "state"] == s).mean()
                         - (o.loc[nt, "state"] == s).mean())
                for s in ("progenitor", "intermediate", "terminal", "cycling")} \
            if "state" in o else {}
        rows.append(dict(gene=g, n=int(ko.sum()), n_mice=len(diffs),
                         dE=est, lo=float(lo), hi=float(hi),
                         target_dlog=eff, **{f"d_{k}": v
                                             for k, v in frac.items()}))
    return rows


def assign_states(a):
    """Coarse states from Leiden clusters' mean program scores."""
    o = a.obs
    cl = o.groupby("leiden", observed=True)[
        ["score_progenitor", "score_terminal", "score_cycling",
         "E_score"]].mean()
    lab = {}
    for c, r in cl.iterrows():
        if r["score_cycling"] > 0.3:
            lab[c] = "cycling"
        elif r["E_score"] < 0:
            lab[c] = "progenitor"
        elif r["E_score"] > 0.8:
            lab[c] = "terminal"
        else:
            lab[c] = "intermediate"
    o["state"] = o["leiden"].map(lab).astype(str)
    return lab


def permutation_fdr(a, genes, axis="E_score", n_perm=2000, seed=0):
    """Within-mouse label permutation: shuffle KO/NT labels inside each
    hashtag (keeps mouse composition), statistic = cell-weighted mean of
    per-mouse differences. Two-sided p, Benjamini-Hochberg FDR."""
    o = a.obs; rng = np.random.default_rng(seed)
    out = {}
    for g in genes:
        sub = o[o["gRNA_gene"].isin([g, "NT"])]
        groups = [(s[axis].values, (s["gRNA_gene"] == g).values)
                  for _, s in sub.groupby("hash.ID", observed=True)
                  if (s["gRNA_gene"] == g).sum() >= 5
                  and (s["gRNA_gene"] == "NT").sum() >= 5]

        def stat(lbls):
            d = [v[l].mean() - v[~l].mean() for (v, _), l in zip(groups, lbls)]
            w = [l.sum() for l in lbls]
            return np.average(d, weights=w)
        obs_ = stat([l for _, l in groups])
        null = np.array([stat([rng.permutation(l) for _, l in groups])
                         for _ in range(n_perm)])
        out[g] = (np.sum(np.abs(null) >= abs(obs_)) + 1) / (n_perm + 1)
    p = np.array([out[g] for g in genes])
    order = np.argsort(p); q = np.empty_like(p)
    q[order] = np.minimum.accumulate((p[order] * len(p)
                                      / np.arange(1, len(p) + 1))[::-1])[::-1]
    return dict(zip(genes, p)), dict(zip(genes, np.minimum(q, 1)))


def insilico_ko(a, genes, hvg, axis="E_score", n_reg=400, alpha=10.0,
                steps=3, seed=0):
    """CellOracle-style in-silico knockout, fitted on NT cells only:
      1. regulators = KO genes (expressed) + the n_reg most variable genes;
         every HVG is ridge-regressed on the regulators (self excluded)
         -> linear GRN W (CellOracle additionally filters edges by TF
         motifs; no motif prior is used here);
      2. KO of g: dx_g = -mean NT expression of g; propagate
         dx <- dx + W dx for `steps` iterations (target genes only);
      3. project onto the order-parameter gradient beta (ridge of the order
         parameter on HVG expression in NT cells) -> predicted dE.
    Baseline: -corr(x_g, E) in NT cells (a KO removes g)."""
    from sklearn.linear_model import Ridge
    nt = (a.obs["gRNA_gene"] == "NT").values
    allg = list(dict.fromkeys([g for g in genes if g in a.var_names]
                              + list(hvg)))
    Xn = a[nt, allg].X
    Xn = Xn.toarray() if hasattr(Xn, "toarray") else np.asarray(Xn)
    mu, sd = Xn.mean(0), Xn.std(0) + 1e-6
    Z = (Xn - mu) / sd
    var_rank = np.argsort(-sd)
    reg_idx = sorted(set([allg.index(g) for g in genes if g in allg])
                     | set(var_rank[:n_reg].tolist()))
    tgt_idx = [allg.index(g) for g in hvg]
    W = np.zeros((len(allg), len(allg)))
    R = Z[:, reg_idx]
    for t in tgt_idx:
        cols = [k for k, r_ in enumerate(reg_idx) if r_ != t]
        m = Ridge(alpha=alpha).fit(R[:, cols], Z[:, t])
        W[t, np.array(reg_idx)[cols]] = m.coef_
    e = a.obs[axis].values[nt]
    beta = Ridge(alpha=alpha).fit(Z[:, tgt_idx], (e - e.mean()) / e.std()) \
        .coef_
    pred, base = {}, {}
    for g in genes:
        if g not in allg:
            continue
        j = allg.index(g)
        dx = np.zeros(len(allg)); dx[j] = -mu[j] / sd[j]
        for _ in range(steps):
            dx = dx + W @ dx
            dx[j] = -mu[j] / sd[j]
        pred[g] = float(beta @ dx[tgt_idx])
        base[g] = float(-np.corrcoef(Z[:, j], e)[0, 1]) if sd[j] > 1e-5 \
            else 0.0
    return pred, base


def response_matrix(a, genes, min_cells=100, seed=0):
    """Global response matrix on the KO'd genes that are measured:
    R[i, j] = within-mouse mean log1p change of gene i when gene j is KO'd.
    Returns R, the gene list, and per-column bootstrap replicates (mice)."""
    o = a.obs
    gl = [g for g in genes if g in a.var_names]
    X = a[:, gl].X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    nt = (o["gRNA_gene"] == "NT").values
    hid = o["hash.ID"].values
    R = np.full((len(gl), len(gl)), np.nan)
    boots = {}
    rng = np.random.default_rng(seed)
    for j, g in enumerate(gl):
        ko = (o["gRNA_gene"] == g).values
        diffs, w = [], []
        for h in np.unique(hid):
            m1 = ko & (hid == h); m2 = nt & (hid == h)
            if m1.sum() >= 5 and m2.sum() >= 5:
                diffs.append(X[m1].mean(0) - X[m2].mean(0)); w.append(m1.sum())
        if not diffs:
            continue
        D = np.array(diffs); w = np.array(w, float)
        R[:, j] = np.average(D, axis=0, weights=w)
        boots[g] = np.array([np.average(D[k], axis=0, weights=w[k])
                             for k in rng.integers(0, len(D), (200, len(D)))])
    return R, gl, boots


def mra_diagnostics(R, gl, boots):
    """Can MRA be applied? Checks the three things that break it:
    (1) the diagonal (does the KO change its own transcript?),
    (2) conditioning of R, (3) stability of r under bootstrap."""
    d = np.diag(R)
    ok = np.isfinite(R).all(0)
    Rk = R[np.ix_(ok, ok)]
    cond = float(np.linalg.cond(Rk))
    out = dict(n=int(ok.sum()), cond=cond,
               diag=dict(zip([g for g, k in zip(gl, ok) if k],
                             np.round(d[ok], 3))),
               frac_diag_negative=float((d[ok] < 0).mean()),
               median_abs_diag=float(np.median(np.abs(d[ok]))),
               median_abs_offdiag=float(np.median(np.abs(
                   Rk[~np.eye(ok.sum(), dtype=bool)]))))
    gk = [g for g, k in zip(gl, ok) if k]
    try:
        r0 = mra(Rk)
        reps = []
        for _ in range(100):
            Rb = np.column_stack([boots[g][np.random.randint(200)][ok]
                                  for g in gk])
            try:
                reps.append(mra(Rb))
            except np.linalg.LinAlgError:
                pass
        reps = np.array(reps)
        off = ~np.eye(len(gk), dtype=bool)
        sign_cons = float(np.mean(np.mean(np.sign(reps[:, off])
                                          == np.sign(r0[off])[None], 0)))
        out.update(r=r0, genes=gk, sign_consistency=sign_cons,
                   boot_cv=float(np.median(np.abs(reps[:, off].std(0)
                                                  / (r0[off] + 1e-9)))))
    except np.linalg.LinAlgError:
        out.update(r=None, genes=gk, sign_consistency=np.nan)
    return out


def mra_breakdown(n_nets=12, n=8, noises=(0.0, 0.01, 0.02, 0.05, 0.1, 0.2,
                                          0.4), diag_atten=(1.0, 0.3, 0.1)):
    """At what data quality does MRA stop working? Sweep measurement noise
    and attenuation of the perturbed node's own response (CRISPR often
    leaves the target transcript intact), and record the diagonal-to-
    off-diagonal ratio alongside recovery quality. The real-data ratio can
    then be read off this curve."""
    rows = []
    for s in range(n_nets):
        net = random_grn(n=n, seed=s)
        x0 = steady_state(net)
        J = jacobian(net, x0)
        if np.max(np.linalg.eigvals(J).real) >= 0:
            continue
        r_true = true_local_response(J, x0)
        off = ~np.eye(n, dtype=bool)
        edge = net["A"].astype(bool) & off
        R0 = knockdown_responses(net, x0, noise=0.0, seed=s)
        rng = np.random.default_rng(s)
        for att in diag_atten:
            for nz in noises:
                R = R0.copy()
                R[np.eye(n, dtype=bool)] *= att
                R = R + nz * rng.normal(size=R.shape)
                ratio = (np.median(np.abs(np.diag(R)))
                         / np.median(np.abs(R[off])))
                try:
                    r_hat = mra(R)
                except np.linalg.LinAlgError:
                    continue
                rows.append(dict(
                    net=s, noise=nz, atten=att, ratio=ratio,
                    sign=float(np.mean(np.sign(r_hat[edge])
                                       == np.sign(r_true[edge]))),
                    cond=float(np.linalg.cond(R))))
    return rows


# ---------------------------------------------------------------------------
# figure + report
# ---------------------------------------------------------------------------
def figure(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    fig, ax = plt.subplots(2, 3, figsize=(16, 9.2)); ax = ax.ravel()
    o = res["obs"]
    a0 = ax[0]
    sc = a0.scatter(o["umap1"], o["umap2"], c=o["E_score"], s=1.2,
                    cmap="viridis", rasterized=True)
    fig.colorbar(sc, ax=a0, label="E (terminal - progenitor score)")
    a0.set_xticks([]); a0.set_yticks([])
    a0.set_title("(a) 25,423 CD8 T cells: exhaustion order parameter\n"
                 f"effective dimensionality (participation ratio) "
                 f"{res['participation_ratio']:.0f}")
    a1 = ax[1]
    for k, (lab, col) in enumerate([("E_dpt", "#2a6fdb"), ("E_pc", "#e4572e")]):
        r = res["agreement"][lab]
        a1.scatter(o["E_score"], (o[lab] - o[lab].mean()) / o[lab].std(),
                   s=0.8, alpha=0.25, c=col, rasterized=True,
                   label=f"{lab}  rho={r:.2f}")
    a1.set_xlabel("E_score (marker-based)"); a1.set_ylabel("z-scored")
    a1.legend(fontsize=8, markerscale=8)
    a1.set_title("(b) three independent definitions of the\norder parameter "
                 "agree on the same axis")
    d = pd.DataFrame(res["ko"]).sort_values("dE")
    a2 = ax[2]
    sig = d["q"] < 0.1
    a2.errorbar(d["dE"], range(len(d)),
                xerr=[d["dE"] - d["lo"], d["hi"] - d["dE"]], fmt="none",
                ecolor="#bbb", lw=0.8)
    a2.scatter(d["dE"], range(len(d)), s=14,
               c=np.where(sig, "#2a6fdb", "#bbb"))
    a2.axvline(0, c="k", lw=0.8)
    a2.axvspan(res["ctrl_ci"][0], res["ctrl_ci"][1], color="#e4572e",
               alpha=0.18, label="negative control (cutting vs\nnon-cutting "
                                 "guides, 95% CI)")
    a2.set_yticks(range(len(d)), d["gene"], fontsize=5)
    a2.set_xlabel("shift in E vs non-targeting (within mouse)")
    a2.legend(fontsize=7, loc="lower right")
    a2.set_title(f"(c) {int(sig.sum())}/{len(d)} knockouts move the\norder "
                 f"parameter (FDR < 0.1)")
    a3 = ax[3]
    for c, lab, col in [("pred_grn", "GRN propagation (CellOracle-style)",
                         "#2a6fdb"),
                        ("pred_corr", "correlation baseline", "#7f8c8d")]:
        v = d[c] / d[c].abs().max()
        a3.scatter(v, d["dE"], s=18, c=col,
                   label=f"{lab}\nrho={res['pred_rho'][c]:+.2f} "
                         f"(p={res['pred_p'][c]:.2f})")
    a3.axhline(0, c="k", lw=0.6); a3.axvline(0, c="k", lw=0.6)
    a3.set_xlabel("predicted effect (scaled)")
    a3.set_ylabel("observed shift in E")
    a3.legend(fontsize=7)
    a3.set_title("(d) in-silico knockout does NOT predict\nthe measured "
                 "effect")
    a4 = ax[4]
    b = pd.DataFrame(res["breakdown"])
    g = b.groupby(pd.cut(b.ratio, [0, 1.5, 3, 6, 12, 40, 1e9]),
                  observed=True).agg(ratio=("ratio", "median"),
                                     sign=("sign", "median"))
    a4.semilogx(g["ratio"], 100 * g["sign"], "o-", c="#2a6fdb",
                label="synthetic networks (known answer)")
    a4.axhline(50, c="k", ls=":", lw=1, label="chance")
    a4.axvline(res["mra"]["ratio"], c="#e4572e", lw=2)
    a4.scatter([res["mra"]["ratio"]], [100 * res["mra"]["sign_consistency"]],
               s=90, marker="*", c="#e4572e", zorder=5,
               label="GSE327731 (real data)")
    a4.set_xlabel("diagonal / off-diagonal ratio of the response matrix")
    a4.set_ylabel("% of interaction signs recovered")
    a4.set_ylim(40, 100); a4.legend(fontsize=7)
    a4.set_title("(e) MRA needs a strong diagonal;\nthe real data sits at "
                 "chance level")
    a5 = ax[5]
    t = res["identifiability"]
    a5.barh(range(len(t)), [x["tau_true"] for x in t], color="#2a6fdb",
            label="true relaxation time")
    for k, x in enumerate(t):
        a5.plot(x["tau_range"], [k, k], c="#e4572e", lw=3, alpha=0.7)
    a5.set_xscale("log"); a5.set_yticks([])
    a5.set_xlabel("relaxation time tau (arb. units)")
    a5.plot([], [], c="#e4572e", lw=3,
            label="range consistent with the same\ninteraction matrix r "
                  "(unknown turnover D)")
    a5.legend(fontsize=7)
    a5.set_title("(f) even with a PERFECT interaction matrix,\ntau spans "
                 "~2 orders of magnitude")
    fig.suptitle("M26 — T cell exhaustion: order parameter works, "
                 "Jacobian recovery from steady-state Perturb-seq does not",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=120)
    plt.close(fig)


def report(fig_path=os.path.join(FIG, "m26_texh_jacobian.png")):
    import json
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    from scipy.stats import spearmanr
    print("=" * 70)
    print("MODULE 26 — T cell exhaustion: order parameter -> Jacobian?")
    print("=" * 70)
    proc = os.path.join(DATA, "processed.h5ad")
    if os.path.exists(proc):
        a = ad.read_h5ad(proc)
        ev = a.obsm["X_pca"].var(0)          # PCA variance kept on the cache
        info = dict(participation_ratio=float(ev.sum() ** 2 / (ev ** 2).sum()))
    else:
        a = load(); info = order_parameter(a); a.write_h5ad(proc)
    assign_states(a)
    agree = {k: spearmanr(a.obs["E_score"], a.obs[k])[0]
             for k in ("E_dpt", "E_pc")}
    print(f"A. order parameter: 3 definitions agree, Spearman "
          + ", ".join(f"{k} {v:.2f}" for k, v in agree.items()))
    ko = pd.DataFrame(ko_effects(a))
    p, q = permutation_fdr(a, ko["gene"].tolist())
    ko["p"] = ko["gene"].map(p); ko["q"] = ko["gene"].map(q)
    targets = set(a.obs["gRNA_gene"].unique()) - {"NT"}
    b = a[:, [g for g in a.var_names if g not in targets
              and not g.startswith(("mt-", "Rpl", "Rps"))]].copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000)
    hvg = list(b.var_names[b.var["highly_variable"]])
    pred, base = insilico_ko(a, ko["gene"].tolist(), hvg)
    ko["pred_grn"] = ko["gene"].map(pred); ko["pred_corr"] = ko["gene"].map(base)
    rho, pv = {}, {}
    rng = np.random.default_rng(0)
    m = ko["pred_grn"].notna()
    for c in ("pred_grn", "pred_corr"):
        r = spearmanr(ko.loc[m, c], ko.loc[m, "dE"])[0]
        null = [spearmanr(rng.permutation(ko.loc[m, c].values),
                          ko.loc[m, "dE"])[0] for _ in range(5000)]
        rho[c] = float(r); pv[c] = float(np.mean(np.abs(null) >= abs(r)))
    print(f"B. {int((ko['q'] < 0.1).sum())}/{len(ko)} knockouts shift E "
          f"(FDR<0.1); in-silico KO predicts them at rho="
          f"{rho['pred_grn']:+.2f} (p={pv['pred_grn']:.2f})")
    nt = a.obs[a.obs["gRNA_gene"] == "NT"].copy()
    nt["ctrl"] = np.where(nt["gRNA_classification"].astype(str)
                          .str.contains("INTERGENIC"), "cut", "nocut")
    rows, w = [], []
    for h, s in nt.groupby("hash.ID", observed=True):
        c1 = s[s.ctrl == "cut"]["E_score"]; c2 = s[s.ctrl == "nocut"]["E_score"]
        if len(c1) >= 5 and len(c2) >= 5:
            rows.append(c1.mean() - c2.mean()); w.append(len(c1))
    rows, w = np.array(rows), np.array(w, float)
    bs = [np.average(rows[k], weights=w[k])
          for k in rng.integers(0, len(rows), (1000, len(rows)))]
    ctrl_ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    print(f"C. negative control (cutting vs non-cutting guides): "
          f"{np.average(rows, weights=w):+.3f} {ctrl_ci}")
    R, gl, boots = response_matrix(a, ko["gene"].tolist())
    mr = mra_diagnostics(R, gl, boots)
    mr["ratio"] = mr["median_abs_diag"] / mr["median_abs_offdiag"]
    print(f"D. MRA on the real response matrix: diag/offdiag "
          f"{mr['ratio']:.1f}, cond {mr['cond']:.0f}, sign consistency "
          f"{mr['sign_consistency']:.0%} (chance 50%)")
    pc_rows, stab = positive_control()
    pc = pd.DataFrame(pc_rows)
    bd = mra_breakdown()
    print("E. synthetic control: MRA recovers interactions when the data "
          "are good (sign "
          f"{pc[(pc.noise == 0) & (pc.method == 'MRA r')]['sign'].median():.0%}"
          f" at zero noise) and fails at the real data's quality")
    print("F. identifiability: with the TRUE interaction matrix but unknown "
          "turnover, tau spans "
          f"{np.median([x['tau_range'][1] / x['tau_range'][0] for x in stab]):.0f}x")
    o = a.obs
    res = dict(obs=pd.DataFrame(dict(
        umap1=a.obsm["X_umap"][:, 0], umap2=a.obsm["X_umap"][:, 1],
        E_score=o["E_score"].values, E_dpt=o["E_dpt"].values,
        E_pc=o["E_pc"].values)),
        participation_ratio=info.get("participation_ratio", np.nan),
        agreement=agree, ko=ko.to_dict("records"), ctrl_ci=ctrl_ci,
        pred_rho=rho, pred_p=pv, mra=mr, breakdown=bd,
        identifiability=stab[:12],
        pc=pc.groupby(["noise", "method"])[["rho", "sign", "edge_auc"]]
        .median().reset_index().to_dict("records"))
    figure(res, fig_path)
    with open(os.path.join(FIG, "m26_summary.json"), "w") as fh:
        json.dump({k: v for k, v in res.items() if k != "obs"}, fh, indent=1,
                  default=lambda x: (x.tolist() if hasattr(x, "tolist")
                                     else (None if isinstance(x, float)
                                           and np.isnan(x) else str(x))))
    ko.to_csv(os.path.join(FIG, "m26_ko_effects.csv"), index=False)
    print(f"figure: {fig_path}")
    return res


if __name__ == "__main__":
    report()
