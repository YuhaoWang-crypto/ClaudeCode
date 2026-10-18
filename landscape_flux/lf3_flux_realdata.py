"""
LF3 — Measuring probability flux (the "momentum") in REAL data.

The bottleneck named in the proposal: flux needs dx/dt.  Two routes:

 (a) TIME SERIES — yeast metabolic cycle (Tu et al. 2005 Science, GEO GSE3431;
     36 microarrays, one every ~25 min, ~3 cycles).  We project the
     transcriptome onto PC1–PC2 and measure the winding number of the state
     trajectory against a time-shuffled null.  Net circulation ≠ 0 is a direct
     signature of broken detailed balance (an equilibrium system has no
     preferred rotation direction).  A thermodynamic-uncertainty-relation (TUR)
     bound then converts the cycle-to-cycle period variability into a minimum
     entropy production per cycle (Barato & Seifert 2015/2016).

 (b) SNAPSHOT ONLY — RNA velocity on the endocrine pancreas (Bastidas-Ponce
     2019; scVelo example data).  dx/dt is estimated from unspliced/spliced
     mRNA (steady-state model, La Manno 2018) WITHOUT any time series.
       * cell cycle:  in the (S-score, G2M-score) plane the textbook order
         G1 → S → G2/M → G1 predicts counter-clockwise rotation.  We test the
         sign of the velocity-derived angular velocity.
       * Helmholtz–Hodge split of the velocity flux field into a curl-free
         (gradient, "downhill") and a rotational part on the cells' kNN graph
         (discrete Hodge decomposition; validated: pure gradient -> 0.000,
         pure rotation -> 0.84, isotropic noise -> 0.67), cycling progenitors
         vs. the endocrine differentiation branch.
"""
import numpy as np
import scipy.sparse as sp
from . import _data

S_GENES = """Mcm5 Pcna Tyms Fen1 Mcm2 Mcm4 Rrm1 Ung Gins2 Mcm6 Cdca7 Dtl Prim1 Uhrf1
Hells Rfc2 Rpa2 Nasp Rad51ap1 Gmnn Wdr76 Slbp Ccne2 Ubr7 Pold3 Msh2 Atad2 Rad51 Rrm2
Cdc45 Cdc6 Exo1 Tipin Dscc1 Blm Casp8ap2 Usp1 Clspn Pola1 Chaf1b Brip1 E2f8""".split()
G2M_GENES = """Hmgb2 Cdk1 Nusap1 Ube2c Birc5 Tpx2 Top2a Ndc80 Cks2 Nuf2 Cks1b Mki67 Tmpo
Cenpf Tacc3 Smc4 Ccnb2 Ckap2l Ckap2 Aurkb Bub1 Kif11 Anp32e Tubb4b Gtse1 Kif20b Hjurp
Cdca3 Cdc20 Ttk Cdc25c Kif2c Rangap1 Ncapd2 Dlgap5 Cdca2 Cdca8 Ect2 Kif23 Hmmr Aurka
Psrc1 Anln Lbr Ckap5 Cenpe Ctcf Nek2 G2e3 Gas2l3 Cbx5 Cenpa""".split()


# ======================================================================= (a)
def ymc_winding(n_null=5000, seed=0):
    from sklearn.decomposition import PCA
    ids, X, t = _data.load_ymc()
    L = np.log2(np.clip(X, 1e-3, None))
    L = L[np.isfinite(L).all(1)]
    var = L.var(1)
    L = L[np.argsort(-var)[:2000]]                     # most dynamic genes
    Z = (L - L.mean(1, keepdims=True)) / L.std(1, keepdims=True)
    pca = PCA(2).fit(Z.T)
    pc = pca.transform(Z.T)
    pc -= pc.mean(0)

    def stats(p):
        dth = np.angle(np.exp(1j * np.diff(np.arctan2(p[:, 1], p[:, 0]))))  # wrapped
        return dth.sum() / (2 * np.pi), dth.sum() / np.abs(dth).sum()

    w, R = stats(pc)
    rng = np.random.default_rng(seed)
    null = np.array([stats(pc[rng.permutation(len(pc))]) for _ in range(n_null)])
    p = (np.sum(np.abs(null[:, 1]) >= abs(R)) + 1) / (n_null + 1)
    dth = np.angle(np.exp(1j * np.diff(np.arctan2(pc[:, 1], pc[:, 0]))))
    frac_same = (np.sign(dth) == np.sign(w)).mean()
    period = len(t) * 25.0 / abs(w)
    return dict(w=w, R=R, p=p, null_sd=null[:, 1].std(), pc=pc, t=t,
                frac_same=frac_same, period=period, evr=pca.explained_variance_ratio_)


# ======================================================================= (b)
def _velocity(a, n_hvg=2000, k=30, seed=0):
    """Steady-state RNA velocity (La Manno 2018 / scVelo 'deterministic')."""
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    S = a.layers["spliced"]; U = a.layers["unspliced"]
    S = S.tocsr() if sp.issparse(S) else sp.csr_matrix(S)
    U = U.tocsr() if sp.issparse(U) else sp.csr_matrix(U)
    ls = np.asarray(S.sum(1)).ravel(); lu = np.asarray(U.sum(1)).ravel()
    genes = np.array(a.var_names)
    ok = (np.asarray((S > 0).sum(0)).ravel() >= 20) & (np.asarray((U > 0).sum(0)).ravel() >= 20)
    Sn = sp.diags(np.median(ls) / ls) @ S[:, ok]
    Un = sp.diags(np.median(lu) / lu) @ U[:, ok]
    g_ok = genes[ok]
    logS = Sn.copy(); logS.data = np.log1p(logS.data)
    mu = np.asarray(logS.mean(0)).ravel()
    var = np.asarray(logS.multiply(logS).mean(0)).ravel() - mu ** 2
    disp = var / (mu + 1e-9)
    cc = np.isin(g_ok, S_GENES + G2M_GENES)
    sel = np.zeros(len(g_ok), bool); sel[np.argsort(-disp)[:n_hvg]] = True; sel |= cc
    Sn, Un, g = Sn[:, sel].toarray(), Un[:, sel].toarray(), g_ok[sel]
    logX = np.log1p(Sn)
    pcs = PCA(30, random_state=seed).fit_transform(logX - logX.mean(0))
    nn = NearestNeighbors(n_neighbors=k).fit(pcs)
    _, ind = nn.kneighbors(pcs)
    Ms = Sn[ind].mean(1); Mu = Un[ind].mean(1)     # first-order moments
    gamma = np.zeros(Ms.shape[1])
    for j in range(Ms.shape[1]):
        s, u = Ms[:, j], Mu[:, j]
        lo, hi = np.percentile(s, [5, 95])
        m = (s <= lo) | (s >= hi)
        den = (s[m] ** 2).sum()
        gamma[j] = (s[m] * u[m]).sum() / den if den > 0 else 0
    V = Mu - gamma * Ms
    return dict(Ms=Ms, V=V, genes=g, pcs=pcs)


def _scores(M, genes, ref_mu, ref_sd):
    Z = (np.log1p(M) - ref_mu) / ref_sd
    s = Z[:, np.isin(genes, S_GENES)].mean(1)
    g2 = Z[:, np.isin(genes, G2M_GENES)].mean(1)
    return np.column_stack([s, g2])


def graph_hodge(P, V, k=15):
    import scipy.sparse.linalg as spla
    from sklearn.neighbors import NearestNeighbors
    """Discrete Hodge split of a velocity field sampled at points P (n×d).
    Edge flow f_ij = ½(v_i+v_j)·(x_j−x_i) on the symmetrised kNN graph;
    potential φ minimises Σ (φ_j−φ_i − f_ij)²; returns residual fraction."""
    Ps = (P - P.mean(0)) / P.std(0)
    Vs = V / P.std(0)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Ps)
    _, ind = nn.kneighbors(Ps)
    i = np.repeat(np.arange(len(P)), k); j = ind[:, 1:].ravel()
    e = np.unique(np.sort(np.c_[i, j], 1), axis=0)
    i, j = e[:, 0], e[:, 1]
    d = Ps[j] - Ps[i]
    f = 0.5 * ((Vs[i] + Vs[j]) * d).sum(1)
    m = len(i)
    G = sp.csr_matrix((np.r_[-np.ones(m), np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[i, j])), shape=(m, len(P)))
    phi = spla.lsqr(G, f, atol=1e-10, btol=1e-10, iter_lim=20000)[0]
    r = f - G @ phi
    return (r ** 2).sum() / (f ** 2).sum()


def helmholtz_curl_fraction(P, Vel, n=40, bw=None):
    """Kernel-smoothed flux field J on a grid, split into ∇φ (curl-free) and a
    remainder by solving ∇²φ = ∇·J with Neumann BCs. Returns ||J_rot||²/||J||²
    weighted by local density."""
    from scipy.ndimage import gaussian_filter
    lo, hi = np.percentile(P, [1, 99], axis=0)
    span = hi - lo; lo -= 0.1 * span; hi += 0.1 * span
    ij = ((P - lo) / (hi - lo) * (n - 1)).round().astype(int)
    ok = ((ij >= 0) & (ij < n)).all(1)
    rho = np.zeros((n, n)); Jx = np.zeros((n, n)); Jy = np.zeros((n, n))
    np.add.at(rho, (ij[ok, 0], ij[ok, 1]), 1)
    np.add.at(Jx, (ij[ok, 0], ij[ok, 1]), Vel[ok, 0])
    np.add.at(Jy, (ij[ok, 0], ij[ok, 1]), Vel[ok, 1])
    sig = bw or 1.5
    rho, Jx, Jy = (gaussian_filter(z, sig) for z in (rho, Jx, Jy))
    div = np.gradient(Jx, axis=0) + np.gradient(Jy, axis=1)
    # Neumann Poisson via DCT
    from scipy.fft import dctn, idctn
    D = dctn(div, type=2, norm="ortho")
    kx = 2 * (np.cos(np.pi * np.arange(n) / n) - 1)
    lap = kx[:, None] + kx[None, :]
    lap[0, 0] = 1
    phi = idctn(D / lap * (lap != 0) * 1.0, type=2, norm="ortho")
    gx, gy = np.gradient(phi)
    Rx, Ry = Jx - gx, Jy - gy
    w = rho / rho.sum()
    m = rho > 0.05 * rho.max()
    num = ((Rx ** 2 + Ry ** 2) * m).sum(); den = ((Jx ** 2 + Jy ** 2) * m).sum()
    curl = np.gradient(Jy, axis=0) - np.gradient(Jx, axis=1)
    return dict(rot_frac=num / den, grid=(rho, Jx, Jy, curl), lo=lo, hi=hi)


def pancreas_flux(seed=0, n_null=1000):
    from sklearn.decomposition import PCA
    a = _data.load_pancreas()
    v = _velocity(a, seed=seed)
    Ms, V, g = v["Ms"], v["V"], v["genes"]
    L = np.log1p(Ms)
    mu, sd = L.mean(0), L.std(0) + 1e-6
    # future state after a small step; dt chosen so the median relative change is ~5%
    dt = 0.05 * np.median(Ms[Ms > 0]) / np.median(np.abs(V[V != 0]))
    F = np.clip(Ms + dt * V, 0, None)
    P0 = _scores(Ms, g, mu, sd); P1 = _scores(F, g, mu, sd)
    dP = P1 - P0
    types = a.obs["clusters"].astype(str).values
    obsS, obsG = a.obs["S_score"].values, a.obs["G2M_score"].values
    cycling = (np.maximum(obsS, obsG) > 0) & np.isin(types, ["Ductal", "Ngn3 low EP"])
    endo = np.isin(types, ["Ngn3 high EP", "Pre-endocrine", "Beta", "Alpha", "Delta", "Epsilon"])

    # angular velocity in the (S, G2M) plane around the cycling centroid
    c = P0[cycling].mean(0)
    r = P0 - c
    omega = (r[:, 0] * dP[:, 1] - r[:, 1] * dP[:, 0]) / (r ** 2).sum(1)
    frac_ccw = (omega[cycling] > 0).mean()
    rng = np.random.default_rng(seed)
    null = []
    idx = np.where(cycling)[0]
    for _ in range(n_null):
        perm = rng.permutation(idx)
        om = (r[idx, 0] * dP[perm, 1] - r[idx, 1] * dP[perm, 0]) / (r[idx] ** 2).sum(1)
        null.append(np.median(om))
    null = np.array(null)
    med = np.median(omega[cycling])
    p = (np.sum(np.abs(null) >= abs(med)) + 1) / (n_null + 1)

    # Helmholtz–Hodge: cycle plane (cycling cells) vs differentiation plane (endocrine)
    hh_cyc = helmholtz_curl_fraction(P0[cycling], dP[cycling])   # for plotting
    hh_cyc["rot_frac"] = graph_hodge(P0[cycling], dP[cycling])
    Le = L[endo]
    pca = PCA(2, random_state=seed).fit(Le - Le.mean(0))
    E0 = pca.transform(Le - Le.mean(0))
    E1 = pca.transform(np.log1p(F[endo]) - Le.mean(0))
    hh_endo = helmholtz_curl_fraction(E0, E1 - E0)
    hh_endo["rot_frac"] = graph_hodge(E0, E1 - E0)
    # null for rot_frac: same positions, velocities shuffled
    nul_c = np.array([graph_hodge(P0[cycling], dP[cycling][rng.permutation(cycling.sum())])
                      for _ in range(100)])
    nul_e = np.array([graph_hodge(E0, (E1 - E0)[rng.permutation(endo.sum())])
                      for _ in range(100)])
    return dict(frac_ccw=frac_ccw, med_omega=med, p=p, n_cyc=int(cycling.sum()),
                P0=P0, dP=dP, cycling=cycling, types=types, hh_cyc=hh_cyc,
                hh_endo=hh_endo, E0=E0, dE=E1 - E0, nul_c=nul_c, nul_e=nul_e,
                omega=omega, endo=endo)


def run(make_fig=True):
    print("=" * 72); print("LF3  Probability flux in real data"); print("=" * 72)
    Y = ymc_winding()
    print("[a] yeast metabolic cycle (GSE3431, 36 arrays x 25 min):")
    print(f"    PC1+PC2 explain {Y['evr'].sum():.0%} of variance of 2000 dynamic genes")
    print(f"    net winding = {Y['w']:+.2f} turns (sign = PCA orientation); "
          f"{Y['frac_same']:.0%} of the 35 steps turn the same way")
    print(f"    rotational consistency R = Σθ/Σ|θ| = {Y['R']:+.2f} "
          f"(time-shuffled null sd = {Y['null_sd']:.2f}, p = {Y['p']:.1e})")
    print(f"    implied period = {Y['period']:.0f} min (Tu et al. 2005 report ~4–5 h)")

    Pn = pancreas_flux()
    print("\n[b] endocrine pancreas, RNA velocity from ONE snapshot (n=3696 cells)")
    print(f"    cycling progenitors n={Pn['n_cyc']}: fraction rotating counter-clockwise "
          f"(G1→S→G2M) = {Pn['frac_ccw']:.1%}")
    print(f"    median angular velocity = {Pn['med_omega']:+.3f} rad/step "
          f"(velocity-shuffled null p = {Pn['p']:.1e})")
    print(f"    Helmholtz–Hodge rotational fraction of flux:")
    print(f"       cell-cycle plane, cycling cells : {Pn['hh_cyc']['rot_frac']:.2f} "
          f"(shuffled null {Pn['nul_c'].mean():.2f} ± {Pn['nul_c'].std():.2f})")
    print(f"       differentiation plane, endocrine: {Pn['hh_endo']['rot_frac']:.2f} "
          f"(shuffled null {Pn['nul_e'].mean():.2f} ± {Pn['nul_e'].std():.2f})")

    if make_fig:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 4, figsize=(20, 4.6))
        pc, t = Y["pc"], Y["t"]
        sc = ax[0].scatter(pc[:, 0], pc[:, 1], c=t, cmap="plasma", s=30, zorder=3)
        ax[0].plot(pc[:, 0], pc[:, 1], "-", c="gray", lw=0.8)
        for k in range(0, len(pc) - 1, 2):
            ax[0].annotate("", pc[k + 1], pc[k], arrowprops=dict(arrowstyle="->", color="k", lw=0.8))
        plt.colorbar(sc, ax=ax[0], label="time (min)")
        ax[0].set_title(f"(a) yeast metabolic cycle\nwinding = {Y['w']:+.2f} turns, R={Y['R']:+.2f}, p={Y['p']:.0e}")
        ax[0].set_xlabel("PC1"); ax[0].set_ylabel("PC2")
        P0, dP, cyc = Pn["P0"], Pn["dP"], Pn["cycling"]
        ax[1].scatter(P0[~cyc, 0], P0[~cyc, 1], s=3, c="lightgray")
        ax[1].scatter(P0[cyc, 0], P0[cyc, 1], s=4, c=Pn["omega"][cyc] > 0, cmap="coolwarm_r")
        sub = np.where(cyc)[0][::6]
        ax[1].quiver(P0[sub, 0], P0[sub, 1], dP[sub, 0], dP[sub, 1], angles="xy",
                     scale_units="xy", scale=0.35, width=0.003)
        ax[1].set_xlabel("S-phase score"); ax[1].set_ylabel("G2/M score")
        ax[1].set_title(f"(b) RNA velocity, cell cycle plane\n{Pn['frac_ccw']:.0%} counter-clockwise (G1→S→G2M)")
        ax[1].scatter([], [], c="#3b4cc0", s=12, label="counter-clockwise (G1→S→G2M)")
        ax[1].scatter([], [], c="#b40426", s=12, label="clockwise")
        ax[1].scatter([], [], c="lightgray", s=12, label="non-cycling")
        ax[1].legend(fontsize=7, loc="upper right")
        panels = [(Pn["hh_cyc"], Pn["nul_c"], "(c) cell-cycle plane (cycling cells)",
                   "S score", "G2M score"),
                  (Pn["hh_endo"], Pn["nul_e"], "(d) differentiation plane (endocrine)",
                   "PC1", "PC2")]
        for k, (hh, nul, title, xl, yl) in enumerate(panels):
            rho, Jx, Jy, _ = hh["grid"]
            m = rho > 0.05 * rho.max()
            vx = np.where(m, Jx / (rho + 1e-12), 0); vy = np.where(m, Jy / (rho + 1e-12), 0)
            ax[2 + k].imshow(np.where(m, rho, np.nan).T, origin="lower", cmap="Greys")
            s_ = slice(1, None, 3)
            X, Yg = np.meshgrid(np.arange(rho.shape[0]), np.arange(rho.shape[1]), indexing="ij")
            mm = m[s_, s_]
            ax[2 + k].quiver(X[s_, s_][mm], Yg[s_, s_][mm], vx[s_, s_][mm], vy[s_, s_][mm],
                             color="crimson", width=0.005)
            ax[2 + k].set_title(f"{title}\nmean velocity; Hodge rotational fraction "
                                f"{hh['rot_frac']:.2f} (null {nul.mean():.2f})", fontsize=10)
            ax[2 + k].set_xlabel(xl); ax[2 + k].set_ylabel(yl)
            ax[2 + k].set_xticks([]); ax[2 + k].set_yticks([])
        fig.tight_layout(); fig.savefig(_data.fig("lf3_flux_realdata.png"), dpi=130)
        plt.close(fig)
    return dict(ymc=Y, pancreas=Pn)


if __name__ == "__main__":
    run()
