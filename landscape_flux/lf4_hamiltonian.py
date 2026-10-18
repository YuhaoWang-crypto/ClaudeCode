"""
LF4 — Is there a biological "Hamiltonian"?  Testing U + K = const.

Proposal: deviation (potential U) + flux rate (kinetic K) = conserved total
energy H, e.g. "differentiation potential + metabolic activity = const".

Tests
 (a) THEORY, on the LF2 activator–inhibitor circuit.
     * Sample states from P_ss; U = −ln P_ss, K = ½|F|².  Conservation would
       require corr(U, K) ≈ −1.
     * What IS conserved: with F = −D∇U + V and ∇·J = 0, one gets
       V·∇U = ∇·V, and in the weak-noise limit V ⟂ ∇Φ (Φ = quasi-potential).
       So the ROTATIONAL part of the dynamics moves along iso-potential
       contours (Hamiltonian-like), while the gradient part dissipates.
       We measure the |cos| between V and ∇U as D → 0.
 (b) REAL DATA, endocrine pancreas (RNA velocity, see LF3):
       U = −ln (kernel density in PCA space), K = ½|velocity|².
 (c) REAL DATA, the proposal's own example:
       potency = number of genes detected per cell (the CytoTRACE principle,
       Gulati et al. 2020 Science), metabolic activity = mean z-score of
       OXPHOS + glycolysis genes; pancreas and Paul et al. 2015 haematopoiesis.
     A linear conservation law aU + bK = const exists iff |corr| → 1 (the
     best combination has variance 1 − |corr| in z-units).
"""
import re
import numpy as np
import scipy.sparse as sp
from . import _data
from .lf2_landscape_flux import Grid2D, F_activator_inhibitor

GLY = ["Hk1", "Hk2", "Gpi1", "Pfkl", "Pfkp", "Pfkm", "Aldoa", "Tpi1", "Gapdh", "Pgk1",
       "Pgam1", "Eno1", "Pkm", "Ldha", "Slc2a1"]
OXPHOS = re.compile(r"^(Nduf|Sdh[abcd]$|Uqcr|Cox[4-8]|Atp5)")


def conservation_stats(u, k):
    zu = (u - u.mean()) / u.std(); zk = (k - k.mean()) / k.std()
    r = np.corrcoef(zu, zk)[0, 1]
    return dict(r=r, var_sum=np.var(zu + zk) / 2,  # 1 + r ; = 0 if U+K const
                best_resid=1 - abs(r))


def theory(seed=0):
    rng = np.random.default_rng(seed)
    out = {}
    g = Grid2D(F_activator_inhibitor, (0, 0), (3.2, 3.2), 110, 0.05).stationary()
    p = g.pi / g.pi.sum()
    idx = rng.choice(len(p), 5000, p=p)
    Fx, Fy = F_activator_inhibitor(g.X.ravel()[idx], g.Y.ravel()[idx])
    U = g.U.ravel()[idx]; K = 0.5 * (Fx ** 2 + Fy ** 2)
    out["UK"] = conservation_stats(U, K)
    cos = []
    for D in (0.2, 0.1, 0.05, 0.025, 0.0125):
        n = 110 if D >= 0.025 else 160
        gD = Grid2D(F_activator_inhibitor, (0, 0), (3.2, 3.2), n, D).stationary()
        c, ratio = gD.flux_vs_gradient()
        cos.append((D, c, ratio))
    out["cos"] = cos
    return out


def _counts(a, layer=None):
    X = a.layers[layer] if layer else a.X
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def potency_metabolism(a, layer=None):
    X = _counts(a, layer)
    genes = np.array(a.var_names)
    lib = np.asarray(X.sum(1)).ravel()
    ngenes = np.asarray((X > 0).sum(1)).ravel().astype(float)
    L = (sp.diags(1e4 / lib) @ X).toarray()
    L = np.log1p(L)
    Z = (L - L.mean(0)) / (L.std(0) + 1e-9)
    met_idx = np.isin(genes, GLY) | np.array([bool(OXPHOS.match(x)) for x in genes])
    met = Z[:, met_idx].mean(1)
    st = conservation_stats(ngenes, met)
    # partial correlation controlling for library size (both scale with UMIs)
    ll = np.log(lib)

    def resid(v):
        b = np.polyfit(ll, v, 1); return v - np.polyval(b, ll)
    st["r_partial"] = np.corrcoef(resid(ngenes), resid(met))[0, 1]
    st["n_met_genes"] = int(met_idx.sum())
    return st, ngenes, met


def pancreas_UK(seed=0):
    from sklearn.decomposition import PCA
    from sklearn.neighbors import KernelDensity
    from .lf3_flux_realdata import _velocity
    a = _data.load_pancreas()
    v = _velocity(a, seed=seed)
    Ms, V = v["Ms"], v["V"]
    L = np.log1p(Ms)
    pca = PCA(10, random_state=seed).fit(L - L.mean(0))
    Z0 = pca.transform(L - L.mean(0))
    dt = 0.05 * np.median(Ms[Ms > 0]) / np.median(np.abs(V[V != 0]))
    Z1 = pca.transform(np.log1p(np.clip(Ms + dt * V, 0, None)) - L.mean(0))
    kde = KernelDensity(bandwidth="scott").fit(Z0)
    U = -kde.score_samples(Z0)
    K = 0.5 * ((Z1 - Z0) ** 2).sum(1)
    types = a.obs["clusters"].astype(str).values
    per = {}
    for t in np.unique(types):
        m = types == t
        if m.sum() > 100:
            per[t] = conservation_stats(U[m], K[m])["r"]
    return conservation_stats(U, K), per, U, K, types


def run(make_fig=True):
    print("=" * 72); print("LF4  Biological 'Hamiltonian': is U + K conserved?"); print("=" * 72)
    T = theory()
    print("[a] theory (activator–inhibitor, states sampled from P_ss):")
    print(f"    corr(U, K) = {T['UK']['r']:+.2f}   Var(zU+zK)/2 = {T['UK']['var_sum']:.2f} "
          f"(0 would mean U+K conserved)")
    print("    flux velocity V vs ∇U, weighted |cos| (0 = V moves along iso-U contours):")
    for D, c, ratio in T["cos"]:
        print(f"       D={D:<7} |cos(V,∇U)| = {c:.3f}   |V|/|D∇U| = {ratio:.2f}")

    UK, per, U, K, types = pancreas_UK()
    print("\n[b] pancreas, U = −ln density (10 PCs), K = ½|RNA velocity|²")
    print(f"    corr(U, K) = {UK['r']:+.2f}   best linear invariant leaves "
          f"{UK['best_resid']:.0%} of variance")
    print("    within cell types: " + ", ".join(f"{t} {r:+.2f}" for t, r in per.items()))

    rows = {}
    pa = _data.load_pancreas()
    rows["pancreas"] = potency_metabolism(pa, "spliced")
    rows["Paul15 haematopoiesis"] = potency_metabolism(_data.load_paul15())
    print("\n[c] 'differentiation potential + metabolic activity = const'?")
    for k, (st, *_rest) in rows.items():
        print(f"    {k:22s} corr(potency, metabolism) = {st['r']:+.2f}  "
              f"(library-size-controlled {st['r_partial']:+.2f}; {st['n_met_genes']} "
              f"metabolic genes); sum would be const only if corr = −1")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 4, figsize=(19, 4.3))
        Ds = [c[0] for c in T["cos"]]; cs = [c[1] for c in T["cos"]]
        ax[0].semilogx(Ds, cs, "o-"); ax[0].set_xlabel("noise D"); ax[0].invert_xaxis()
        ax[0].minorticks_off(); ax[0].set_xticks(Ds); ax[0].set_xticklabels([f"{d:g}" for d in Ds])
        ax[0].set_ylabel("|cos(V, ∇U)|")
        ax[0].set_title("(a) the rotational flux runs along\niso-potential contours as D→0")
        cmap = plt.get_cmap("tab10")
        for k_, t_ in enumerate(np.unique(types)):
            m_ = types == t_
            ax[1].scatter(U[m_], K[m_], s=3, color=cmap(k_), label=t_)
        ax[1].legend(fontsize=6, markerscale=3, ncol=2)
        ax[1].set_xlabel("U = −ln density"); ax[1].set_ylabel("K = ½|v|²")
        ax[1].set_yscale("log")
        ax[1].set_title(f"(b) pancreas: corr(U,K) = {UK['r']:+.2f}\n(U+K const needs −1)")
        for k, (name, (st, ng, met)) in enumerate(rows.items()):
            ax[2 + k].scatter(ng, met, s=3, alpha=0.5)
            ax[2 + k].set_xlabel("genes detected (potency proxy)")
            ax[2 + k].set_ylabel("OXPHOS+glycolysis score")
            ax[2 + k].set_title(f"(c) {name}\ncorr = {st['r']:+.2f}")
        fig.tight_layout(); fig.savefig(_data.fig("lf4_hamiltonian.png"), dpi=130)
        plt.close(fig)
    return dict(theory=T, pancreas_UK=(UK, per), potency=rows)


if __name__ == "__main__":
    run()
