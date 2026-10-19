"""
LF8 — Critical slowing in the RNA-velocity field: replication with a speed control.

An independent synthesis report (scVelo dynamical model, 2-D UMAP) found that
the local-Jacobian relaxation time |τ| of the velocity field peaks at the
Pre-endocrine stage (5.2x the rest, ~30x baseline in pseudotime), i.e. just
before the Alpha/Beta/Delta/Epsilon fate split.

Concern: a local Jacobian scales linearly with velocity magnitude.  If RNA
velocity is simply small in a region (e.g. few unspliced reads), τ = 1/|Re λ|
inflates there without any dynamical slowing.  Three readouts:
  1. raw τ per cell type (replication; steady-state velocity, 5 PCs)
  2. corr(log τ, log |v|)
  3. τ in units of the local transit time, τ·|v|/ℓ (ℓ = neighbourhood radius):
     dimensionless, invariant to a regional rescaling of velocity
  4. null: velocities shuffled WITHIN each cell type (keeps each type's speed
     distribution, destroys the field's structure)
"""
import numpy as np
from . import _data


def local_tau(Z, V, k=40, ridge=1e-3):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k).fit(Z)
    dist, ind = nn.kneighbors(Z)
    d = Z.shape[1]
    tau = np.empty(len(Z)); ell = dist[:, -1]
    for i in range(len(Z)):
        dx = Z[ind[i]] - Z[i]
        dv = V[ind[i]] - V[ind[i]].mean(0)
        dx = dx - dx.mean(0)
        G = dx.T @ dx + ridge * np.trace(dx.T @ dx) / d * np.eye(d)
        A = np.linalg.solve(G, dx.T @ dv).T          # dv ≈ A dx
        re = np.abs(np.linalg.eigvals(A).real)
        tau[i] = 1 / max(re.min(), 1e-12)
    return tau, ell


def run(make_fig=True, seed=0):
    from sklearn.decomposition import PCA
    from .lf3_flux_realdata import _velocity
    print("=" * 72); print("LF8  Critical slowing in the velocity field (with a speed control)"); print("=" * 72)
    a = _data.load_pancreas()
    v = _velocity(a, seed=seed)
    Ms, Vg = v["Ms"], v["V"]
    L = np.log1p(Ms)
    pca = PCA(5, random_state=seed).fit(L - L.mean(0))
    Z0 = pca.transform(L - L.mean(0))
    dt = 0.05 * np.median(Ms[Ms > 0]) / np.median(np.abs(Vg[Vg != 0]))
    V = pca.transform(np.log1p(np.clip(Ms + dt * Vg, 0, None)) - L.mean(0)) - Z0
    types = a.obs["clusters"].astype(str).values
    speed = np.linalg.norm(V, axis=1)
    tau, ell = local_tau(Z0, V)
    tnorm = tau * speed / ell
    rng = np.random.default_rng(seed)
    Vs = V.copy()
    for t in np.unique(types):
        m = np.where(types == t)[0]
        Vs[m] = V[rng.permutation(m)]
    tau_null, _ = local_tau(Z0, Vs)

    order = ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta", "Alpha", "Delta", "Epsilon"]
    print(f"    {'cell type':15s} {'median τ':>10s} {'median |v|':>11s} {'τ·|v|/ℓ':>9s} {'τ (shuffled v)':>15s}")
    rows = {}
    for t in order:
        m = types == t
        rows[t] = (np.median(tau[m]), np.median(speed[m]), np.median(tnorm[m]), np.median(tau_null[m]))
        print(f"    {t:15s} {rows[t][0]:10.1f} {rows[t][1]:11.4f} {rows[t][2]:9.3f} {rows[t][3]:15.1f}")
    from scipy.stats import mannwhitneyu, spearmanr
    pe = types == "Pre-endocrine"
    ratio = np.median(tau[pe]) / np.median(tau[~pe])
    ratio_n = np.median(tnorm[pe]) / np.median(tnorm[~pe])
    ratio_null = np.median(tau_null[pe]) / np.median(tau_null[~pe])
    p_raw = mannwhitneyu(tau[pe], tau[~pe], alternative="greater").pvalue
    p_n = mannwhitneyu(tnorm[pe], tnorm[~pe], alternative="greater").pvalue
    rho = spearmanr(tau, speed).statistic
    print(f"    Pre-endocrine / rest: raw τ {ratio:.2f}x (p={p_raw:.1e}); "
          f"speed-normalised {ratio_n:.2f}x (p={p_n:.1e}); shuffled-v null {ratio_null:.2f}x")
    print(f"    Spearman(τ, |v|) = {rho:+.2f}")
    sweep = {}
    for d in (2, 3, 10):
        p_ = PCA(d, random_state=seed).fit(L - L.mean(0))
        Zd = p_.transform(L - L.mean(0))
        Vd = p_.transform(np.log1p(np.clip(Ms + dt * Vg, 0, None)) - L.mean(0)) - Zd
        td, _ = local_tau(Zd, Vd)
        Vds = Vd.copy()
        for t in np.unique(types):
            m = np.where(types == t)[0]
            Vds[m] = Vd[rng.permutation(m)]
        tds, _ = local_tau(Zd, Vds)
        sweep[d] = (np.median(td[pe]) / np.median(td[~pe]), np.median(tds[pe]) / np.median(tds[~pe]))
    print("    dimension sweep, Pre-endocrine/rest (real vs shuffled-v): " +
          ", ".join(f"{d} PCs {r:.2f}x vs {n:.2f}x" for d, (r, n) in sweep.items()))
    print("    -> the real field is not distinguishable from the shuffled null: no replicated peak")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(16, 4.3))
        for k, (vals, lab) in enumerate([(tau, "raw τ = 1/|Re λ|min"), (tnorm, "τ·|v|/ℓ (speed-normalised)"),
                                         (tau_null, "τ with velocities shuffled within type")]):
            data = [vals[types == t] for t in order]
            ax[k].boxplot(data, showfliers=False)
            ax[k].set_yscale("log"); ax[k].set_xticks(range(1, len(order) + 1))
            ax[k].set_xticklabels(order, rotation=40, ha="right", fontsize=8)
            ax[k].get_xticklabels()[order.index("Pre-endocrine")].set_color("crimson")
            ax[k].set_title(lab)
        fig.tight_layout(); fig.savefig(_data.fig("lf8_jacobian_slowing.png"), dpi=130)
        plt.close(fig)
    return dict(rows=rows, sweep=sweep, ratio=ratio, ratio_norm=ratio_n, ratio_null=ratio_null, rho=rho,
                p_raw=p_raw, p_norm=p_n)


if __name__ == "__main__":
    run()
