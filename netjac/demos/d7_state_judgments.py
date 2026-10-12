"""
D7 — State judgments beyond eigenvalues, on the user's macro model
   dx/dt = 0.5 y - gamma x ,  dy/dt = -0.8 x - delta y
(x = unemployment gap, y = inflation gap), then the same matrix estimated
from REAL US data (FRED: UNRATE - NROU, CPI YoY - 2) era by era.

  A  phase map in (gamma, delta): node / focus / centre / saddle / source
  B  energy function: gamma, delta are the ONLY dissipation channels; the
     off-diagonal loop exchanges energy (sets frequency, not decay)
  C  reactivity / transient amplification of a stable system
  D  random shocks: quasi-cycles from a STABLE focus, lead-lag, and why the
     contemporaneous Phillips correlation can vanish
  E  delay: the policy lag that turns a soft landing into a boom-bust cycle
  F  controllability / observability, pole placement (feedback design)
  G  coupling sensitivity: which entry sets damping vs frequency
  H  REAL DATA (FRED, quarterly): era Jacobians 1960-83, 1984-2007,
     2008-19, 2020-; eigenvalues, reactivity, sign pattern, quasi-cycle
     period, lead-lag; residual-bootstrap 90 % intervals.
"""
import io
import os
import subprocess

import numpy as np

from netjac import FIG, DATA_CACHE
from netjac import statespace as ss

Jf = lambda g, d: np.array([[-g, 0.5], [-0.8, -d]])
SOFT = (0.2, 0.1)


def part_model():
    g, d = SOFT
    J = Jf(g, d)
    out = {"J": J, "eig": np.linalg.eigvals(J)}
    gs = np.linspace(-1.5, 1.5, 241)
    out["map_axes"] = gs
    out["map"] = ss.phase_map(Jf, gs, gs)
    out["energy"] = ss.quadratic_energy(J)
    out["energy_centre"] = ss.quadratic_energy(Jf(0.0, 0.0))
    # reactivity: symmetric part [[-g, -0.15], [-0.15, -d]] -> reactive iff
    # g*d < 0.15^2 (both damping terms positive)
    out["react"] = ss.reactivity(J)
    out["amp"] = ss.max_amplification(J)
    out["react_threshold"] = ((0.5 - 0.8) / 2) ** 2
    out["noise"] = ss.noise_response(J)
    out["period"] = 2 * np.pi / abs(out["eig"][0].imag)
    cd = ss.critical_delay(J, 0, 1)
    out["tau_star"], out["w_star"] = cd
    t1, X1 = ss.simulate_dde(J, 0, 1, 0.9 * cd[0], [0, 2], 300)
    t2, X2 = ss.simulate_dde(J, 0, 1, 1.1 * cd[0], [0, 2], 300)
    out["dde"] = (t1, X1, t2, X2)
    # delay threshold over the (gamma, delta) plane
    out["tau_grid"] = [(gg, ss.critical_delay(Jf(gg, gg / 2), 0, 1))
                       for gg in (0.05, 0.1, 0.2, 0.4, 0.8)]
    out["gram_x"] = ss.gramians(J, B=[1, 0], Cm=[1, 0])
    out["gram_y"] = ss.gramians(J, B=[0, 1], Cm=[0, 1])
    out["place"] = ss.pole_placement(J, [1, 0], [-0.5, -0.6])
    out["lqr"] = ss.lqr(J, [1, 0], np.eye(2), np.array([[1.0]]))
    out["sens"], _ = ss.coupling_sensitivity(J)
    return out


# ----------------------------------------------------------- real data
def _fred(series):
    import pandas as pd
    os.makedirs(os.path.join(DATA_CACHE, "fred"), exist_ok=True)
    p = os.path.join(DATA_CACHE, "fred", f"{series}.csv")
    if not os.path.exists(p):
        r = subprocess.run(["curl", "-sf", "--max-time", "60",
                            "https://fred.stlouisfed.org/graph/fredgraph.csv"
                            f"?id={series}"], capture_output=True, check=True)
        open(p, "wb").write(r.stdout)
    d = pd.read_csv(p)
    d.columns = ["date", series]
    d["date"] = pd.to_datetime(d["date"])
    d[series] = pd.to_numeric(d[series], errors="coerce")
    return d.set_index("date")[series]


def macro_data():
    """Quarterly x = UNRATE - NROU (CBO natural rate), y = CPI YoY - 2."""
    import pandas as pd
    u = _fred("UNRATE").resample("QS").mean()
    n = _fred("NROU").resample("QS").mean()
    cpi = _fred("CPIAUCSL").resample("QS").mean()
    infl = 100 * (cpi / cpi.shift(4) - 1)
    df = pd.concat({"x": u - n, "y": infl - 2.0}, axis=1).dropna()
    return df[df.index <= pd.Timestamp.today()]


ERAS = [("1960–1983 Great Inflation", "1960", "1983"),
        ("1984–2007 Great Moderation", "1984", "2007"),
        ("2008–2019 post-crisis", "2008", "2019"),
        ("2020– pandemic era", "2020", "2030")]


def era_analysis(Z, dt=0.25, n_boot=1000, seed=0):
    rng = np.random.default_rng(seed)
    f = ss.var_jacobian(Z, dt)
    J = f["J"]
    lam = np.linalg.eigvals(J)
    nr = ss.noise_response(J, D=f["Sigma"] / dt)
    res = {"n": len(Z), "J": J, "eig": lam, "kind": ss.classify_2x2(J),
           "react": ss.reactivity(J),
           "period": (2 * np.pi / abs(lam[0].imag)) if abs(lam[0].imag)
           > 1e-9 else None,
           "lag": nr["lag"], "lag_corr": nr["lag_corr"],
           "contemp": np.corrcoef(Z[:, 0], Z[:, 1])[0, 1]}
    # residual bootstrap (fixed design) for J entries and lead Re(lambda)
    A, E = f["A"], f["resid"]
    c = np.linalg.lstsq(np.hstack([np.ones((len(Z) - 1, 1)), Z[:-1]]),
                        Z[1:], rcond=None)[0][0]
    boots = []
    for _ in range(n_boot):
        Eb = E[rng.integers(0, len(E), len(E))]
        Zb = np.zeros_like(Z)
        Zb[0] = Z[0]
        for t in range(1, len(Z)):
            Zb[t] = c + A @ Zb[t - 1] + Eb[t - 1]
        try:
            Jb = ss.var_jacobian(Zb, dt)["J"]
            lb = np.linalg.eigvals(Jb)
            boots.append([Jb[0, 1], Jb[1, 0], Jb[0, 0], Jb[1, 1],
                          lb.real.max()])
        except Exception:                               # noqa: BLE001
            pass
    B = np.array(boots)
    q = lambda k: tuple(np.nanpercentile(B[:, k], [5, 95]))
    res["ci"] = {"J12": q(0), "J21": q(1), "J11": q(2), "J22": q(3),
                 "re_lead": q(4)}
    return res


def part_data():
    df = macro_data()
    out = {"df": df, "eras": {}}
    for lab, a, b in ERAS:
        sub = df.loc[a:b]
        if len(sub) < 20:
            continue
        out["eras"][lab] = era_analysis(sub[["x", "y"]].values)
    return out


# --------------------------------------------------------------- report
def report(make_fig=True, with_data=True):
    M = part_model()
    D = None
    if with_data:
        try:
            D = part_data()
        except Exception as e:                          # noqa: BLE001
            print(f"D7-H real data skipped ({type(e).__name__}: {e})")
    print("=" * 68)
    print("D7  STATE JUDGMENTS BEYOND EIGENVALUES — user's macro model")
    print("=" * 68)
    e = M["energy"]["diag"]
    print(f"[B] energy V = x² + {e['p'][1]:.3f} y²  ->  dV/dt = "
          f"-{e['dissipation'][0]:.3f} x² - {e['dissipation'][1]:.3f} y² "
          f"(cross terms cancel: {e['cross_terms_cancel']})  ✅")
    ec = M["energy_centre"]["diag"]
    print(f"    centre (γ=δ=0): dissipation {ec['dissipation']} -> V "
          f"conserved -> closed ellipses, size fixed by the shock")
    print(f"[C] reactivity {M['react']:+.4f} (reactive iff γδ < "
          f"{M['react_threshold']:.4f}; here γδ = {SOFT[0]*SOFT[1]:.3f}); "
          f"max ||e^Jt|| = {M['amp']['amp']:.4f} at t = {M['amp']['t']:.2f}"
          f"; worst shock direction {np.round(M['amp']['worst_shock'], 2)}")
    nz = M["noise"]
    print(f"[D] random shocks: spectral peak ω = {nz['peak_freq'][0]:.3f} "
          f"(eigen-frequency {abs(M['eig'][0].imag):.3f}; period "
          f"{M['period']:.2f}) = quasi-cycle from a STABLE focus")
    print(f"    contemporaneous corr(x, y) = {nz['contemp_corr']:+.3f}, but "
          f"corr peaks at {nz['lag_corr']:+.2f} with lag {nz['lag']:+.2f} "
          f"(inflation leads unemployment)")
    print(f"[E] critical policy delay on the y→x channel: τ* = "
          f"{M['tau_star']:.3f} (ω* = {M['w_star']:.3f}); DDE check: late "
          f"|y| at 0.9τ* = {np.abs(M['dde'][1][-3000:, 1]).max():.3f}, at "
          f"1.1τ* = {np.abs(M['dde'][3][-3000:, 1]).max():.3f}  ✅")
    print("    τ* vs damping (δ = γ/2): " + ", ".join(
        f"γ={g}: {c[0]:.2f}" for g, c in M["tau_grid"] if c))
    print(f"[F] controllable from either channel / observable from either "
          f"readout: {M['gram_x']['controllable']}, {M['gram_x']['observable']}"
          f", {M['gram_y']['controllable']}, {M['gram_y']['observable']}")
    K, cl = M["place"]
    fb = lambda K: f"u = {-K[0,0]:+.3f}·x {-K[0,1]:+.3f}·y"
    print(f"    feedback into the x-equation {fb(K)} -> eigenvalues "
          f"{np.round(cl, 3)} (no oscillation)")
    K2, cl2 = M["lqr"]
    print(f"    LQR (Q=I, R=1): {fb(K2)} -> {np.round(cl2, 3)}")
    S = M["sens"]
    print(f"[G] dλ/dJ: diagonal → Re {S[0,0].real:+.2f}, {S[1,1].real:+.2f} "
          f"(damping); off-diagonal → Im {S[0,1].imag:+.3f}, "
          f"{S[1,0].imag:+.3f} (frequency only)")
    if D:
        print("[H] REAL US DATA (quarterly; x = UNRATE−NROU, y = CPI YoY − 2)")
        for lab, r in D["eras"].items():
            ci = r["ci"]
            per = f"{r['period']:.1f} y" if r["period"] else "—"
            print(f"    {lab} (n={r['n']}): J = {np.round(r['J'], 3).tolist()}"
                  f"; λ = {np.round(r['eig'], 3)}; {r['kind']}")
            print(f"        J12 (y→x, 'policy') {r['J'][0,1]:+.3f} "
                  f"[{ci['J12'][0]:+.2f},{ci['J12'][1]:+.2f}]; J21 (x→y, "
                  f"'Phillips') {r['J'][1,0]:+.3f} [{ci['J21'][0]:+.2f},"
                  f"{ci['J21'][1]:+.2f}]; lead Re λ 90% "
                  f"[{ci['re_lead'][0]:+.3f},{ci['re_lead'][1]:+.3f}]")
            print(f"        period {per}; reactivity {r['react']:+.3f}; "
                  f"contemp corr {r['contemp']:+.2f}; model lag "
                  f"{r['lag']:+.2f} y (corr {r['lag_corr']:+.2f})")
    if make_fig:
        _figure(M, D)
    return {"model": M, "data": D}


def _figure(M, D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(16, 10))
    a = ax[0, 0]
    cats = ["stable node", "stable focus", "centre", "unstable focus",
            "unstable node", "saddle"]
    cols = ["#7fb7a4", "#2f8f78", "#d9a441", "#e07a7a", "#b33a48", "#7a6fb0"]
    idx = np.vectorize(lambda s: cats.index(s))(M["map"])
    g = M["map_axes"]
    a.imshow(idx, origin="lower", extent=[g[0], g[-1], g[0], g[-1]],
             cmap=ListedColormap(cols), vmin=0, vmax=5, aspect="auto")
    a.plot(*SOFT, "k*", ms=12)
    a.plot(0, 0, "ko", mfc="w")
    if D:     # real-era damping terms (their own off-diagonals differ)
        for k, (lab, r) in enumerate(D["eras"].items()):
            gx, dx = -r["J"][0, 0], -r["J"][1, 1]
            if abs(gx) < 1.5 and abs(dx) < 1.5:
                a.plot(gx, dx, "D", color="w", mec="k", ms=6)
                a.annotate(lab.split(" ")[0], (gx, dx), fontsize=7,
                           xytext=(4, 4), textcoords="offset points")
    for k, c in enumerate(cats):
        a.plot([], [], "s", color=cols[k], label=c)
    a.legend(fontsize=7, loc="lower left")
    a.set_xlabel("γ (self-repair of unemployment)")
    a.set_ylabel("δ (inflation anchoring)")
    a.set_title("A  phase map: ★ soft landing, ○ centre, ◇ US eras")
    a = ax[0, 1]
    am = M["amp"]
    a.plot(am["ts"], am["norms"])
    a.axhline(1, color="k", lw=.5)
    a.set_xlabel("time"); a.set_ylabel("max ||e^{Jt}|| over shocks")
    a.set_title(f"C  reactivity {M['react']:+.4f}: tiny transient growth")
    a = ax[0, 2]
    nz = M["noise"]
    a.plot(nz["taus"], nz["xcorr"])
    a.axvline(0, color="k", lw=.5)
    a.set_xlabel("lag τ (y leads x for τ<0)")
    a.set_ylabel("corr( x(t), y(t+τ) )")
    a.set_title(f"D  contemporaneous {nz['contemp_corr']:+.2f}, at lag "
                f"{nz['lag']:+.1f}: {nz['lag_corr']:+.2f}")
    a = ax[1, 0]
    t1, X1, t2, X2 = M["dde"]
    a.plot(t1, X1[:, 1], label="lag 0.9 τ*")
    a.plot(t2, X2[:, 1], label="lag 1.1 τ*", alpha=.8)
    a.set_ylim(-4, 4)
    a.set_xlabel("time"); a.set_ylabel("inflation gap y")
    a.set_title(f"E  policy lag τ* = {M['tau_star']:.2f} turns soft landing "
                "into cycles")
    a.legend(fontsize=8)
    if D:
        a = ax[1, 1]
        df = D["df"]
        for lab, s0, s1 in ERAS:
            sub = df.loc[s0:s1]
            a.plot(sub.x, sub.y, lw=.9, label=lab.split(" ")[0])
        a.set_xlabel("unemployment gap x (pp)")
        a.set_ylabel("inflation gap y (pp)")
        a.set_title("H  real US state space, by era")
        a.legend(fontsize=7)
        a = ax[1, 2]
        for k, (lab, r) in enumerate(D["eras"].items()):
            lam = r["eig"]
            a.scatter(lam.real, lam.imag, s=40, label=lab.split(" ")[0],
                      color=f"C{k}")
            lo, hi = r["ci"]["re_lead"]
            a.plot([lo, hi], [k * 0.05 - 0.1] * 2, color=f"C{k}", lw=3,
                   alpha=.5)
        a.axvline(0, color="k", lw=.6)
        a.set_xlabel("Re λ (1/year)"); a.set_ylabel("Im λ (rad/year)")
        a.set_title("H  era Jacobian eigenvalues (bars: 90% CI of lead Re λ)")
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d7_state_judgments.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    report()
