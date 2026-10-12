"""
D8 — The state-space modules applied to concrete physiological cases.

  1 ENERGY / SENSITIVITY — "which buffer to strengthen":
    Topp 2000 glucose-insulin-beta-cell model (BIOMD0000000341). For each
    intervention lever, compare (a) the linear view: total sensitivity of the
    slowest eigenvalue (resilience) and (b) the nonlinear view: the largest
    fraction of beta-cell mass that can be lost and still recover (distance to
    the diabetic attractor).
  2 DELAY / REACTIVITY — "when does it oscillate / overshoot":
    Sturis 1991 ultradian insulin model (BIOMD0000000382): critical hepatic
    delay for the ultradian oscillation; reactivity of the stable states.
  3 OBSERVABILITY / CONTROLLABILITY — "which readout suffices, which drug
    steers what": Topp readouts (glucose / insulin / beta-cell mass);
    Lemaire 2004 bone remodelling (BIOMD0000000278) drug inputs PTH, RANKL,
    OPG: steady-state gain matrix -J^-1 B and its rank.
  4 POLE PLACEMENT / LQR — "what feedback gain":
    artificial pancreas on the Topp model in type-1 state (beta-cell mass 0)
    with a subcutaneous insulin-absorption compartment (1 h, ⚠️ assumed): glucose-only feedback vs LQR
    (glucose + insulin-on-board), critical gain (delay-loop Hopf), meal test
    with the physical constraint u >= 0.
  5 SEGMENTED VAR — "let the state space grow out of data":
    PhysioNet Fantasia (20 young / 20 elderly, 2 h ECG + respiration):
    cardiorespiratory 2x2 Jacobian per 25-min segment, respiratory-sinus-
    arrhythmia gain (positive control: falls with age), within-subject
    reliability across segments.
"""
import os

import numpy as np

from netjac import FIG
from netjac import loops as L
from netjac import statespace as ss

LEVERS = {   # Topp 2000 parameter -> what it stands for clinically
    "si": "insulin sensitivity (exercise, TZD)",
    "R0": "hepatic glucose output (metformin target)",
    "Eg0": "glucose effectiveness",
    "k": "insulin clearance",
    "sigma": "max insulin secretion per beta-cell",
    "alpha": "secretion half-saturation (glucose sensing)",
    "d0": "beta-cell death rate",
    "r1": "glucose-driven beta-cell replication",
    "r2": "glucotoxic beta-cell loss",
}


# ================================================================ case 1
def case1(step=1.1):
    from scipy.optimize import root
    S, _ = L.sbml_system("BIOMD0000000341")
    xh = np.array([100.0, 10.0, 300.0])

    def eq(p):
        return root(lambda z: S.rhs(z, p), xh, method="hybr", tol=1e-12).x

    def slow(p):
        x = eq(p)
        return np.linalg.eigvals(S.jacobian(x, p)).real.max(), x

    def recovers(f, p):
        x = eq(p).copy()
        x[2] *= 1 - f
        _, X = S.simulate(x, 4000, 400, p)
        return X[-1, 0] < 150

    def margin(p):
        lo, hi = 0.0, 0.999
        for _ in range(22):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if recovers(mid, p) else (lo, mid)
        return 0.5 * (lo + hi)
    l0, x0 = slow({})
    f0 = margin({})
    J = S.jacobian(x0)
    rows = []
    for p, lab in LEVERS.items():
        v = S.defaults[p]
        lu, xu = slow({p: v * step})
        ld, _ = slow({p: v / step})
        fu, fd = margin({p: v * step}), margin({p: v / step})
        rows.append({"lever": p, "meaning": lab,
                     "dlam": (lu - ld) / (2 * np.log(step)),
                     "dmargin": (fu - fd) / (2 * np.log(step)),
                     "m_up": fu, "m_dn": fd, "B_up": xu[2]})
    return {"x": x0, "J": J, "eig": np.linalg.eigvals(J), "lam": l0,
            "margin": f0, "rows": rows,
            "tip_G": 0.00084 / (2 * 2.4e-6)}


# ================================================================ case 2
def case2():
    from scipy.optimize import root
    S, _ = L.sbml_system("BIOMD0000000382")
    att = S.attractor(3000)
    xs = root(lambda z: S.rhs(z), att["x_end"], method="hybr").x
    rows, ev = L.bifurcation_scan(S, "t3", np.linspace(5, 60, 221), xs)
    crit = [r for r in rows if r[2] is not None]
    pairs = []
    for v, x, e in crit:
        cp = e[e.imag > 1e-9]
        if len(cp):
            c = cp[np.argmax(cp.real)]
            pairs.append((v, c.real, 2 * np.pi / c.imag))
    # reactivity of a stable state (t3 = 15 min) in relative coordinates
    s = root(lambda z: S.rhs(z, {"t3": 15.0}), xs, method="hybr")
    J = S.jacobian(s.x, {"t3": 15.0})
    Jr = np.diag(1 / s.x) @ J @ np.diag(s.x)
    amp = ss.max_amplification(Jr, t_max=600)
    # Topp healthy state for comparison
    T, _ = L.sbml_system("BIOMD0000000341")
    xt = np.array([100.0, 10.0, 300.0])
    Jt = np.diag(1 / xt) @ T.jacobian(xt) @ np.diag(xt)
    ampt = ss.max_amplification(Jt, t_max=5)
    return {"events": ev, "pairs": pairs, "t3_default": S.defaults["t3"],
            "names": S.names, "react": ss.reactivity(Jr), "amp": amp,
            "topp_react": ss.reactivity(Jt), "topp_amp": ampt,
            "sim": S.attractor(3000)}


# ================================================================ case 3
def case3():
    from scipy.optimize import root
    T, _ = L.sbml_system("BIOMD0000000341")
    xt = np.array([100.0, 10.0, 300.0])
    Jr = np.diag(1 / xt) @ T.jacobian(xt) @ np.diag(xt)     # relative units
    obs = {}
    for lab, C in [("glucose", [1, 0, 0]), ("insulin", [0, 1, 0]),
                   ("β-cell mass", [0, 0, 1]),
                   ("glucose + insulin", [[1, 0, 0], [0, 1, 0]])]:
        g = ss.gramians(Jr, Cm=C)
        w = g["Wo_eigs"]
        obs[lab] = {"rank_ok": g["observable"], "w_min": w.min(),
                    "cond": w.max() / max(w.min(), 1e-300)}
    S, _ = L.sbml_system("BIOMD0000000278")
    att = S.attractor(500)
    x = root(lambda z: S.rhs(z), att["x_end"], method="hybr").x
    J = S.jacobian(x)
    drugs = {"I_P": "PTH (teriparatide)", "I_L": "RANKL",
             "I_O": "OPG (≈ denosumab: blocks RANKL)"}
    Bm, gains = [], {}
    for k in drugs:
        h = 1e-3 * max(1.0, abs(S.defaults[k]))
        b = (S.rhs(x, {k: S.defaults[k] + h}) -
             S.rhs(x, {k: S.defaults[k] - h})) / (2 * h)
        Bm.append(b)
        gains[k] = -np.linalg.solve(J, b) / x            # relative change
    Bm = np.array(Bm).T
    G = np.array([gains[k] for k in drugs]).T
    sv = np.linalg.svd(G / np.linalg.norm(G, axis=0), compute_uv=False)
    ctrl = {k: ss.gramians(np.diag(1 / x) @ J @ np.diag(x),
                           B=Bm[:, i] / x)["controllable"]
            for i, k in enumerate(drugs)}
    return {"obs": obs, "bone_x": x, "bone_names": S.names, "drugs": drugs,
            "B_nonzero_rows": [S.names[i] for i in
                               np.where(np.abs(Bm).sum(1) > 0)[0]],
            "gains": gains, "gain_rank_sv": sv, "ctrl": ctrl}


# ================================================================ case 4
def case4(tau_a=1 / 24, meal=450.0, t_meal=0.03):
    from scipy.integrate import solve_ivp
    S, _ = L.sbml_system("BIOMD0000000341")
    p = S.defaults
    G0 = 100.0
    I0 = (p["R0"] / G0 - p["Eg0"]) / p["si"]
    u0 = p["k"] * I0
    xs = np.array([G0, I0, u0 * tau_a])
    Jt = S.jacobian([G0, I0, 0.0])[:2, :2]
    J = np.zeros((3, 3))
    J[:2, :2] = Jt
    J[1, 2] = 1 / tau_a
    J[2, 2] = -1 / tau_a
    Bv = np.array([0.0, 0.0, 1.0])
    Ra = lambda t: meal * t / t_meal ** 2 * np.exp(-t / t_meal)

    def f(t, x, K):
        d = S.rhs([x[0], x[1], 0.0])[:2]
        u = max(0.0, u0 - K @ (x - xs))
        return [d[0] + Ra(t), d[1] + x[2] / tau_a, u - x[2] / tau_a]

    def sim(K, days=2.0):
        s = solve_ivp(f, (0, days), xs, args=(K,), max_step=0.001, rtol=1e-7)
        G = s.y[0]
        u = np.array([max(0.0, u0 - K @ (x - xs)) for x in s.y.T])
        late = s.t > days - 0.5                       # last 12 h
        return {"t": s.t * 24, "G": G, "peak": G.max(), "min": G.min(),
                "h_above180": np.trapezoid(G > 180, s.t) * 24,
                "h_below70": np.trapezoid(G < 70, s.t) * 24,
                "insulin": np.trapezoid(u, s.t) / u0,    # in basal-days
                "late_swing": G[late].max() - G[late].min(),
                "u_zero_frac": float(np.mean(u[s.t < 1] == 0))}
    res = {"open": sim(np.zeros(3)), "J": J, "basal": u0}
    # glucose-only proportional feedback: closed-loop eigenvalues vs gain
    ks = np.geomspace(10, 1e5, 400)
    lead = [np.linalg.eigvals(J - np.outer(Bv, [-k, 0, 0])).real.max()
            for k in ks]
    k_crit = ks[np.argmax(np.array(lead) > 0)] if max(lead) > 0 else None
    res["gain_scan"] = (ks, np.array(lead))
    res["k_crit"] = k_crit
    res["prop"] = {k: sim(np.array([-k, 0.0, 0.0]))
                   for k in (200.0, 3200.0, 30000.0)}
    res["lqr"] = {}
    for q in (50.0, 20.0, 5.0):
        K, cl = ss.lqr(J, Bv, np.diag([1 / q ** 2, 0, 0]),
                       np.array([[1 / u0 ** 2]]))
        res["lqr"][q] = {"K": K[0], "poles": cl, "sim": sim(K[0])}
    return res


# ================================================================ case 5
FS = 4.0


def _prep_fantasia(npz, seg_min=25, n_seg=4):
    from scipy.interpolate import interp1d
    from scipy.signal import butter, filtfilt
    d = np.load(npz, allow_pickle=True)
    b = d["beats"]
    rr = np.diff(b) * 1000
    tb = b[1:]
    med = np.convolve(rr, np.ones(11) / 11, mode="same")
    ok = (rr > 300) & (rr < 2000) & (np.abs(rr - med) < 0.2 * med)
    resp = d["resp"]
    tr = np.arange(len(resp)) / FS
    t0, t1 = max(tb[ok][0], 0), min(tb[ok][-1], tr[-1])
    grid = np.arange(t0, t1, 1 / FS)
    R = interp1d(tb[ok], rr[ok], kind="cubic", bounds_error=False,
                 fill_value="extrapolate")(grid)
    P = np.interp(grid, tr, resp)
    bb, aa = butter(2, [0.03 / (FS / 2), 0.8 / (FS / 2)], "band")
    Z = np.column_stack([filtfilt(bb, aa, R), filtfilt(bb, aa, P)])
    n = int(seg_min * 60 * FS)
    segs = [Z[i * n:(i + 1) * n] for i in range(n_seg) if (i + 1) * n <= len(Z)]
    return segs, float(d["age"]), str(d["group"]), float(np.median(rr[ok]))


def _segment_metrics(Z):
    from scipy.signal import csd, welch, coherence
    from netjac.demos.d6_baroreflex_jacobian import fit_var
    f, Pr = welch(Z[:, 1], FS, nperseg=512)
    band = (f > 0.1) & (f < 0.5)
    fr = f[band][np.argmax(Pr[band])]                 # breathing frequency
    _, Pxy = csd(Z[:, 1], Z[:, 0], FS, nperseg=512)
    _, Coh = coherence(Z[:, 1], Z[:, 0], FS, nperseg=512)
    k = np.argmin(np.abs(f - fr))
    rsa = np.abs(Pxy[k]) / Pr[k] * np.std(Z[:, 1])     # ms per SD of resp.
    Zs = Z / Z.std(0)
    J = ss.var_jacobian(Zs, 1 / FS)["J"]
    fits = [(p,) + fit_var(Zs, p) for p in range(1, 17)]
    p, A, S, _ = min(fits, key=lambda r: r[3])
    from netjac.demos.d6_baroreflex_jacobian import var_spectrum
    # directional gains at the breathing frequency from the VAR spectrum
    H = var_spectrum(A, S, [fr], 1 / FS)[0]
    return {"f_resp": fr, "rsa": rsa, "coh": Coh[k], "J": J, "p": p,
            # respiration is uncalibrated (polarity arbitrary): magnitudes only
            "J_RR_resp": abs(J[0, 1]), "J_resp_RR": abs(J[1, 0]),
            "J_RR_RR": J[0, 0], "J_resp_resp": J[1, 1]}


def case5():
    import glob
    from scipy.stats import mannwhitneyu, spearmanr
    from netjac.data import fantasia
    fantasia.fetch_all()
    subs = []
    for p in sorted(glob.glob(os.path.join(fantasia.DIR, "*.npz"))):
        segs, age, group, rr = _prep_fantasia(p)
        ms = [_segment_metrics(z) for z in segs]
        if not ms:
            continue
        subs.append({"rec": os.path.basename(p)[:-4], "age": age,
                     "group": group, "rr": rr, "segs": ms})
    keys = ["rsa", "f_resp", "coh", "J_RR_resp", "J_resp_RR", "J_RR_RR",
            "J_resp_resp"]
    tab = {k: np.array([np.median([m[k] for m in s["segs"]]) for s in subs])
           for k in keys}
    grp = np.array([s["group"] for s in subs])
    age = np.array([s["age"] for s in subs])
    comp = {}
    for k in keys:
        y, o = tab[k][grp == "young"], tab[k][grp == "elderly"]
        comp[k] = (np.median(y), np.median(o),
                   mannwhitneyu(y, o).pvalue, spearmanr(tab[k], age).statistic)
    # within-subject reliability: correlation between segment 1 and 2..4
    rel = {}
    for k in ("rsa", "J_RR_resp", "f_resp"):
        a = np.array([s["segs"][0][k] for s in subs if len(s["segs"]) > 1])
        b = np.array([np.mean([m[k] for m in s["segs"][1:]])
                      for s in subs if len(s["segs"]) > 1])
        rel[k] = spearmanr(a, b).statistic
    Jy = np.median(np.stack([np.median([m["J"] for m in s["segs"]], axis=0)
                             for s in subs if s["group"] == "young"]), axis=0)
    Jo = np.median(np.stack([np.median([m["J"] for m in s["segs"]], axis=0)
                             for s in subs if s["group"] == "elderly"]), axis=0)
    return {"n": len(subs), "comp": comp, "rel": rel, "Jy": Jy, "Jo": Jo,
            "tab": tab, "grp": grp, "age": age, "subs": subs}


# ================================================================ report
def report(make_fig=True, with_data=True):
    C1, C2, C3, C4 = case1(), case2(), case3(), case4()
    C5 = None
    if with_data:
        try:
            C5 = case5()
        except Exception as e:                       # noqa: BLE001
            print(f"case 5 skipped ({type(e).__name__}: {e})")
    print("=" * 68)
    print("D8  MODULE CASE STUDIES")
    print("=" * 68)
    print(f"[1] Topp healthy state {C1['x']}: eigenvalues "
          f"{np.round(C1['eig'], 3)} (/day); slowest = β-cell mass, "
          f"τ = {1/abs(C1['lam']):.0f} d; β-cell loss that still recovers: "
          f"{C1['margin']*100:.1f} %")
    print("    lever                          dReλ_slow/dln p   dmargin/dln p"
          "   margin at ×1.1 / ÷1.1")
    for r in sorted(C1["rows"], key=lambda r: -abs(r["dmargin"])):
        print(f"    {r['lever']:6s} {r['meaning'][:26]:26s} {r['dlam']:+9.5f}"
              f"       {r['dmargin']:+7.3f}      {r['m_up']*100:5.1f} % / "
              f"{r['m_dn']*100:5.1f} %")
    print(f"    tipping glucose r1/(2 r2) = {C1['tip_G']:.0f} mg/dL depends "
          f"only on β-cell turnover terms")
    print(f"[2] Sturis: {C2['events']} (hepatic delay t3; model default "
          f"{C2['t3_default']:.0f} min)")
    pr = [p for p in C2["pairs"] if p[0] in (10.0, 20.0, 25.0, 30.0, 36.0)]
    print("    critical pair Re / period by t3: " + "; ".join(
        f"t3={v:.0f}: {re:+.4f}, {per:.0f} min" for v, re, per in pr))
    a = C2["amp"]
    print(f"    stable state (t3=15): reactivity {C2['react']:+.3f}/min, max "
          f"amplification {a['amp']:.2f}× at {a['t']:.1f} min; worst shock "
          f"{dict(zip(C2['names'], np.round(a['worst_shock'], 2)))}")
    ta = C2["topp_amp"]
    print(f"    Topp healthy: reactivity {C2['topp_react']:+.1f}/day, max "
          f"amplification {ta['amp']:.2f}× at {ta['t']*24:.1f} h; worst "
          f"shock {np.round(ta['worst_shock'], 3)} (G, I, B)")
    print("[3] observability of the Topp state from one readout "
          "(relative units):")
    for k, v in C3["obs"].items():
        print(f"    {k:18s} rank-observable {v['rank_ok']}; weakest direction "
              f"Gramian {v['w_min']:.2e}; condition {v['cond']:.1e}")
    print(f"    Lemaire bone drugs enter only: {C3['B_nonzero_rows']}")
    for k, lab in C3["drugs"].items():
        g = C3["gains"][k]
        print(f"    {lab:32s} steady-state relative gain (R, B, C) = "
              f"{np.round(g / np.abs(g).max(), 3)}; controllable "
              f"{C3['ctrl'][k]}")
    sv = C3["gain_rank_sv"]
    print(f"    normalised gain matrix singular values {np.round(sv, 4)} -> "
          f"effective rank {int((sv > 1e-6 * sv[0]).sum())}")
    print(f"[4] artificial pancreas (Topp T1D + 1 h absorption, meal "
          f"450 mg/dL): open loop peak {C4['open']['peak']:.0f}, "
          f"{C4['open']['h_above180']:.2f} h > 180, insulin "
          f"{C4['open']['insulin']:.2f} basal-days")
    print(f"    glucose-only gain: Hopf (oscillation) above k ≈ "
          f"{C4['k_crit']:.0f}")
    for k, s in C4["prop"].items():
        print(f"      k={k:7.0f}: peak {s['peak']:.0f}, min {s['min']:.0f}, "
              f"{s['h_above180']:.2f} h > 180, {s['h_below70']:.2f} h < 70; "
              f"insulin {s['insulin']:.2f} basal-days; pump off "
              f"{s['u_zero_frac']*100:.0f} % of day 1; glucose swing in "
              f"last 12 h {s['late_swing']:.1f} mg/dL")
    for q, d in C4["lqr"].items():
        s = d["sim"]
        print(f"    LQR (q={q:.0f}): K = {np.round(d['K'], 1)} (G, I, "
              f"insulin-on-board); poles {np.round(d['poles'], 1)}; peak "
              f"{s['peak']:.0f}, min {s['min']:.0f}, {s['h_above180']:.2f} h"
              f" > 180; insulin {s['insulin']:.2f}; late swing "
              f"{s['late_swing']:.1f}")
    if C5:
        print(f"[5] Fantasia cardiorespiratory (n={C5['n']}; per subject "
              "median of 4 × 25-min segments)")
        lab = {"rsa": "RSA gain (ms/SD resp) [positive control]",
               "f_resp": "breathing frequency (Hz)",
               "coh": "resp–RR coherence at breathing f",
               "J_RR_resp": "|J| resp → RR (coupling)",
               "J_resp_RR": "|J| RR → resp (feedback)",
               "J_RR_RR": "J: RR self-damping",
               "J_resp_resp": "J: resp self-damping"}
        for k, (y, o, p, rho) in C5["comp"].items():
            print(f"    {lab[k]:42s} young {y:+.3f} elderly {o:+.3f} "
                  f"(p={p:.1g}; ρ_age {rho:+.2f})")
        print("    within-subject reliability (segment 1 vs 2-4, Spearman): " +
              ", ".join(f"{k} {v:.2f}" for k, v in C5["rel"].items()))
        print(f"    median J young {np.round(C5['Jy'], 3).tolist()}; elderly "
              f"{np.round(C5['Jo'], 3).tolist()}")
    if make_fig:
        _figure(C1, C2, C3, C4, C5)
    return {"c1": C1, "c2": C2, "c3": C3, "c4": C4, "c5": C5}


def _figure(C1, C2, C3, C4, C5):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    a = ax[0, 0]
    rows = sorted(C1["rows"], key=lambda r: r["dmargin"])
    y = np.arange(len(rows))
    a.barh(y, [r["dmargin"] for r in rows], color="C3", alpha=.8,
           label="β-cell-loss margin")
    a.set_yticks(y)
    a.set_yticklabels([f"{r['lever']}: {r['meaning'][:22]}" for r in rows],
                      fontsize=7)
    a.axvline(0, color="k", lw=.5)
    a.set_xlabel("d margin / d ln(lever)")
    a.set_title("1  Topp: which buffer protects against diabetes")
    a = ax[0, 1]
    v = [p[0] for p in C2["pairs"]]
    a.plot(v, [p[1] for p in C2["pairs"]])
    a.axhline(0, color="k", lw=.5)
    a.axvline(C2["t3_default"], color="gray", ls=":", label="model default")
    for t3, kind in C2["events"]:
        a.axvline(t3, color="r", ls="--", label=f"Hopf t3={t3:.1f}")
    a.set_xlabel("hepatic delay t3 (min)")
    a.set_ylabel("Re λ of oscillatory pair (1/min)")
    a.set_title("2  Sturis: ultradian rhythm needs delay > t3*")
    a.legend(fontsize=7)
    a = ax[0, 2]
    labs = list(C3["obs"])
    a.barh(range(len(labs)), [np.log10(C3["obs"][k]["cond"]) for k in labs],
           color="C0")
    a.set_yticks(range(len(labs))); a.set_yticklabels(labs, fontsize=8)
    a.set_xlabel("log10 observability condition number (lower = better)")
    a.set_title("3  which single readout reveals the Topp state")
    a = ax[1, 0]
    ks, lead = C4["gain_scan"]
    a.semilogx(ks, lead)
    a.axhline(0, color="k", lw=.5)
    if C4["k_crit"]:
        a.axvline(C4["k_crit"], color="r", ls="--")
    a.set_xlabel("glucose-only feedback gain k")
    a.set_ylabel("Re λ_lead closed loop (1/day)")
    a.set_title("4  insulin feedback through absorption lag → Hopf")
    a = ax[1, 1]
    o = C4["open"]
    a.plot(o["t"], o["G"], label="basal only", color="gray")
    for k, s in C4["prop"].items():
        a.plot(s["t"], s["G"], lw=.9, label=f"glucose-only k={k:.0f}")
    q, d = list(C4["lqr"].items())[-1]
    a.plot(d["sim"]["t"], d["sim"]["G"], lw=2, color="C3",
           label="LQR (G + insulin on board)")
    a.axhline(180, color="k", ls=":", lw=.6); a.axhline(70, color="k", ls=":",
                                                         lw=.6)
    a.set_xlim(0, 12); a.set_xlabel("hours after meal")
    a.set_ylabel("glucose (mg/dL)")
    a.set_title("4  meal test with u ≥ 0")
    a.legend(fontsize=7)
    a = ax[1, 2]
    if C5:
        g = C5["grp"]
        for k, (lab, col) in enumerate([("young", "C0"), ("elderly", "C1")]):
            vals = C5["tab"]["rsa"][g == lab]
            a.scatter(np.full(len(vals), k) + np.random.default_rng(k)
                      .uniform(-.1, .1, len(vals)), vals, color=col)
            a.plot([k - .2, k + .2], [np.median(vals)] * 2, color="k")
        a.set_xticks([0, 1]); a.set_xticklabels(["young (21–34)",
                                                 "elderly (68–85)"])
        a.set_yscale("log")
        a.set_ylabel("RSA gain (ms per SD respiration)")
        a.set_title("5  Fantasia: data-grown cardiorespiratory coupling")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d8_module_cases.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    report()
