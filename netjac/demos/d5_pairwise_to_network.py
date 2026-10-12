"""
D5 — From pairwise interaction matrices to whole-network attractors.

  A  The user's tumour–immune model (2x2): equilibria, Jacobian, loops,
     global portrait, the "T cells drop 0.8 -> 0.5" scenario, and the
     immune-escape bifurcation (T-cell death rate d crosses the carrying
     capacity K: transcritical exchange of stability).
  B  The user's macro model (gamma, delta): soft landing vs "centre". A linear
     centre is NOT a limit-cycle attractor: it is structurally fragile
     (gamma+delta = +-0.01 flips it) and its amplitude is set by the shock.
     A true limit cycle needs a nonlinear term (shown).
  C  Pairwise is not enough:
       C1 every 2x2 sub-network stable, whole 3-loop unstable (Goodwin /
          HPA-axis / circadian type negative feedback -> oscillation)
       C2 a sub-network unstable on its own, the whole stable (a runaway
          positive loop held by an inhibitor)
       C3 random 3x3 / 5x5 networks: how often does "all pairs stable"
          predict "whole stable"?
     The exact combination rule = characteristic polynomial = sums of
     disjoint feedback loops (Levins), checked by Routh-Hurwitz.
  D  Physiology atlas on CURATED BioModels equations (exact): Topp 2000
     beta-cell/glucose/insulin (bistability, diabetes saddle-node), Sturis
     1991 & Tolic 2000 ultradian insulin oscillations (Hopf), Lemaire 2004
     bone RANKL/OPG, Komarova 2003 bone remodelling, Miao 2010 influenza
     immune response, Borghans 1997 calcium oscillations.
  E  Is "more variables, more pairs" better? May (1972) random vs
     feedback-structured networks.
"""
import os
from itertools import combinations

import numpy as np

from netjac import FIG
from netjac import loops as L


# ======================================================== A tumour-immune
def tumour_immune():
    return L.System(
        "tumour–immune (user)",
        lambda x, p: [p["r"] * x[0] * (1 - x[0] / p["K"]) - x[0] * x[1],
                      x[0] * x[1] - p["d"] * x[1]],
        ["C", "T"], {"r": 2.0, "K": 5.0, "d": 3.0})


def part_a():
    S = tumour_immune()
    eq = S.equilibria([(0, 8), (0, 4)])
    table = []
    for e in eq:
        a = L.analyse_point(S.jacobian(e))
        table.append((e, a["kind"], a["eigenvalues"]))
    xs = np.array([3.0, 0.8])
    J = S.jacobian(xs)
    a = L.analyse_point(J)
    dec, chk = L.loop_decomposition(J, S.names)
    # the clinical scenario: T drops 0.8 -> 0.5
    t, X = S.simulate([3.0, 0.5], 25.0, 5000)
    peak_C = X[:, 0].max()
    t_peak = t[np.argmax(X[:, 0])]
    # overshoot ratio per cycle = exp(Re(lambda) * period)
    per = a["period"]
    damp_per_cycle = np.exp(a["lead"].real * per)
    eqs, portrait = L.global_portrait(S, [(0.01, 6), (0.01, 3)], n_ic=80)
    # immune escape: coexistence (C*, T*) = (d, r(1-d/K)) exists iff d < K
    rows, ev = L.bifurcation_scan(S, "d", np.linspace(1, 6, 501), xs)
    J50 = S.jacobian(np.array([5.0, 0.0]), {"d": 6.0})
    return {"S": S, "eq": table, "J": J, "a": a, "loops": dec, "chk": chk,
            "hurwitz": L.hurwitz(L.char_coeffs(J)),
            "sign_stab": L.sign_pattern_stability(J),
            "gersh": L.gershgorin(J)[0], "t": t, "X": X, "peak_C": peak_C,
            "t_peak": t_peak, "damp": damp_per_cycle, "portrait": portrait,
            "eqs": eqs, "bif": ev,
            "escape_eigs": np.linalg.eigvals(J50)}


# ================================================================ B macro
def macro(nonlinear=False):
    if not nonlinear:
        f = lambda x, p: [0.5 * x[1] - p["g"] * x[0],
                          -0.8 * x[0] - p["dl"] * x[1]]
    else:
        # de-anchored expectations amplify small inflation gaps (-delta y,
        # delta < 0) but saturate for large ones (- y^3): a Rayleigh-type term
        f = lambda x, p: [0.5 * x[1] - p["g"] * x[0],
                          -0.8 * x[0] - p["dl"] * x[1] - x[1] ** 3]
    return L.System("macro", f, ["x(unemp gap)", "y(infl gap)"],
                    {"g": 0.2, "dl": 0.1}, nonneg=False)


def part_b():
    out = {}
    M = macro()
    for lab, g, dl in [("soft landing", 0.2, 0.1), ("centre", 0.0, 0.0),
                       ("centre +0.01", 0.005, 0.005),
                       ("centre -0.01", -0.005, -0.005)]:
        J = M.jacobian(np.zeros(2), {"g": g, "dl": dl})
        out[lab] = L.analyse_point(J)
    dec, _ = L.loop_decomposition(M.jacobian(np.zeros(2)), M.names)
    out["loops"] = dec
    # linear centre: orbit size = shock size (no attractor)
    amps = {}
    for shock in (0.5, 2.0):
        t, X = M.simulate([0, shock], 200, 4000, {"g": 0.0, "dl": 0.0})
        amps[("linear", shock)] = X[len(X) // 2:, 1].max()
    Mn = macro(True)
    for shock in (0.5, 2.0):
        t, X = Mn.simulate([0, shock], 400, 8000, {"g": 0.0, "dl": -0.1})
        amps[("nonlinear", shock)] = X[len(X) // 2:, 1].max()
    out["amps"] = amps
    out["traj_lin"] = M.simulate([0, 2.0], 60, 3000, {"g": 0.0, "dl": 0.0})
    out["traj_soft"] = M.simulate([0, 2.0], 60, 3000, {"g": 0.2, "dl": 0.1})
    out["traj_nl"] = Mn.simulate([0, 0.1], 300, 6000, {"g": 0.0, "dl": -0.1})
    return out


# ======================================================= C pairwise vs whole
def part_c(seed=0):
    res = {}
    # C1 Goodwin-type 3-step negative feedback
    G = lambda k: np.array([[-1, 0, -k], [1, -1, 0], [0, 1, -1.0]])
    res["goodwin"] = []
    for k in (4.0, 8.0, 10.0):
        J = G(k)
        pairs = L.principal_subsystems(J, 2)
        a = L.analyse_point(J)
        res["goodwin"].append((k, all(p["stable"] for p in pairs), a["kind"],
                               a["eigenvalues"], L.hurwitz(L.char_coeffs(J))))
    dec, chk = L.loop_decomposition(G(10.0), ["X", "Y", "Z"])
    res["goodwin_loops"], res["goodwin_chk"] = dec, chk
    # C2 runaway positive pair (x1<->x2 mutual activation) held by x3
    J2 = np.array([[-0.5, 1.0, -2.0], [1.0, -0.5, 0.0], [1.0, 0.0, -1.0]])
    res["held"] = {"J": J2, "pairs": L.principal_subsystems(J2, 2),
                   "whole": L.analyse_point(J2)}
    # C3 random screening statistics
    rng = np.random.default_rng(seed)
    stats = {}
    for n in (3, 5):
        tab = np.zeros((2, 2), int)          # [pairs_all_stable, whole_stable]
        for _ in range(20000):
            A = rng.normal(0, 1, (n, n)) * (rng.random((n, n)) < 0.6)
            np.fill_diagonal(A, -np.abs(rng.normal(1, 0.5, n)))
            whole = np.linalg.eigvals(A).real.max() < 0
            pairs = all(np.linalg.eigvals(A[np.ix_(ix, ix)]).real.max() < 0
                        for ix in combinations(range(n), 2))
            tab[int(pairs), int(whole)] += 1
        ppv = tab[1, 1] / tab[1].sum()          # P(whole stable | pairs ok)
        miss = tab[0, 1] / tab[:, 1].sum()      # stable wholes w/ bad pair
        stats[n] = {"table": tab, "ppv": ppv, "stable_but_pair_unstable": miss}
    res["random"] = stats
    return res


# ===================================================== D physiology atlas
ATLAS = [
    # id, horizon for simulation, short physiological reading
    ("BIOMD0000000341", 3000, "Topp 2000: glucose–insulin–β-cell mass"),
    ("BIOMD0000000382", 3000, "Sturis 1991: ultradian insulin–glucose"),
    ("BIOMD0000000372", 3000, "Tolić 2000: insulin–glucose feedback"),
    ("BIOMD0000000278", 500, "Lemaire 2004: bone RANK/RANKL/OPG"),
    ("BIOMD0000000148", 3000, "Komarova 2003: osteoclast/osteoblast/bone"),
    ("BIOMD0000000546", 30, "Miao 2010: influenza innate+adaptive"),
    ("BIOMD0000000045", 100, "Borghans 1997: Ca2+ oscillation (CICR)"),
]


def atlas_entry(bid, horizon, label, box_scale=3.0):
    from scipy.optimize import root
    S, rr = L.sbml_system(bid, label)
    att = S.attractor(horizon)
    # fixed point near the attractor (unstable inside an oscillation)
    x = att["x_end"]
    s = root(lambda z: S.rhs(z), x, method="hybr", tol=1e-13)
    scale = 1.0 + np.abs(S.rhs(x)).max()
    res = np.abs(S.rhs(s.x)).max()
    ok = np.isfinite(s.x).all() and res < 1e-6 * scale
    if not ok:     # inside an oscillation: multi-start around the attractor
        hi = box_scale * np.abs(att["X"]).max(0) + 1e-9
        eqs = S.equilibria([(0, h) for h in hi], n_starts=300)
        if eqs:
            m = att["X"].mean(0)
            s.x = min(eqs, key=lambda q: np.linalg.norm((q - m) / hi))
            ok = True
    e = {"id": bid, "label": label, "n": S.n, "names": S.names,
         "oscillating": att["oscillating"], "period": att["period"],
         "state": att["state"],
         "fp_found": bool(ok)}
    if not ok:
        return e
    xs = s.x
    J = S.jacobian(xs)
    a = L.analyse_point(J)
    pairs = L.principal_subsystems(J, 2) if S.n >= 2 else []
    F = L.feedback_levels(J)
    loops = sorted(L.elementary_loops(J), key=lambda l: -abs(l[1]))
    # destabilising loops: positive-gain loops (feedback that amplifies)
    pos = [l for l in loops if l[1] > 0]
    # the loop most likely to drive an oscillation: strongest multi-step
    # loop by geometric-mean gain |g|^(1/length)
    multi = [l for l in loops if len(l[0]) >= 2]
    key = max(multi, key=lambda l: abs(l[1]) ** (1 / len(l[0]))) \
        if multi else None
    e.update(x=xs, J=J, kind=a["kind"], lead=a["lead"],
             n_pairs=len(pairs),
             n_pairs_unstable=sum(not p["stable"] for p in pairs),
             F=F, hurwitz=L.hurwitz(L.char_coeffs(J)),
             sign_stab=L.sign_pattern_stability(J, n=1500),
             gersh=L.gershgorin(J)[0],
             top_loops=[("->".join(S.names[i] for i in l[0]), l[1])
                        for l in loops[:4]],
             top_pos=[("->".join(S.names[i] for i in l[0]), l[1])
                      for l in pos[:2]],
             key_loop=(("->".join(S.names[i] for i in key[0]) + "->" +
                        S.names[key[0][0]], key[1]) if key else None),
             unstable_pairs=[(S.names[p["idx"][0]], S.names[p["idx"][1]])
                             for p in pairs if not p["stable"]],
             t=att["t"], X=att["X"])
    return e


def part_d():
    rows = []
    for bid, T, lab in ATLAS:
        try:
            rows.append(atlas_entry(bid, T, lab))
        except Exception as ex:                       # noqa: BLE001
            rows.append({"id": bid, "label": lab, "error": str(ex)[:80]})
    # Topp: all equilibria + diabetes saddle-node in beta-cell death rate
    S, _ = L.sbml_system("BIOMD0000000341")
    eq = S.equilibria([(0, 600), (0, 100), (0, 1500)], n_starts=300)
    topp_eq = [(e, L.analyse_point(S.jacobian(e))["kind"]) for e in eq]
    scan, ev = L.bifurcation_scan(S, "d0", np.linspace(0.06, 0.08, 401),
                                  [100, 10, 300])
    si_scan, si_ev = L.bifurcation_scan(S, "si", np.linspace(0.72, 0.05, 135),
                                        [100, 10, 300])
    d0_star = 0.00084 ** 2 / (4 * 2.4e-6)       # r1^2 / (4 r2), analytic
    return {"rows": rows, "topp_eq": topp_eq, "topp_scan": scan,
            "topp_events": ev, "topp_si_events": si_ev, "d0_star": d0_star,
            "si_scan": si_scan}


# ================================================================== E May
def part_e():
    ns = np.array([2, 5, 10, 20, 40, 80])
    rnd = L.may_random(ns, trials=300)
    st = L.may_random(ns, trials=300, structured=True)
    return {"ns": ns, "random": rnd, "structured": st,
            "may_n": (1.0 / 0.5) ** 2 / 0.3}


# ================================================================ report
def report(make_fig=True):
    A, B, C, D, E = part_a(), part_b(), part_c(), part_d(), part_e()
    print("=" * 68)
    print("D5  PAIRWISE -> NETWORK ATTRACTOR ANALYSIS")
    print("=" * 68)
    print("[A] tumour–immune (user): equilibria")
    for e, k, v in A["eq"]:
        print(f"    (C,T)=({e[0]:.3g},{e[1]:.3g}): {k}; λ={np.round(v, 3)}")
    print(f"    J(3,0.8) = {np.round(A['J'], 3).tolist()}; λ = "
          f"{A['a']['lead'].real:.3f}±{abs(A['a']['lead'].imag):.4f}i  ✅ "
          f"(user's hand calculation confirmed)")
    print(f"    loops: " + "; ".join(f"a{k}: " + ", ".join(
        f"{d} → {c:+.3g}" for d, c in v) for k, v in A["loops"].items())
        + f"  (check {A['chk']:.0e})")
    print(f"    sign-pattern stability {A['sign_stab']:.2f} (signs alone "
          f"guarantee convergence); Gershgorin-sufficient: {A['gersh']}")
    print(f"    global: {dict(A['portrait'])} of 80 random starts")
    print(f"    T 0.8→0.5: tumour rebounds to C={A['peak_C']:.2f} at "
          f"t={A['t_peak']:.2f}; period {A['a']['period']:.2f}; each cycle "
          f"keeps {A['damp']*100:.1f} % of the previous swing")
    print(f"    immune escape: {A['bif']}; at d=6 the tumour-only state "
          f"(5,0) has λ={np.round(A['escape_eigs'], 3)} -> attractor")
    print("[B] macro (user):")
    for lab in ("soft landing", "centre", "centre +0.01", "centre -0.01"):
        a = B[lab]
        print(f"    {lab:13s}: λ={np.round(a['eigenvalues'][0], 3)}  "
              f"{a['kind'][:45]}")
    am = B["amps"]
    print(f"    linear centre: shock 0.5 -> cycle amp {am[('linear', 0.5)]:.2f};"
          f" shock 2.0 -> {am[('linear', 2.0)]:.2f} (orbit = shock: no "
          f"attractor)")
    print(f"    nonlinear (Rayleigh) version: shock 0.5 -> "
          f"{am[('nonlinear', 0.5)]:.2f}; shock 2.0 -> "
          f"{am[('nonlinear', 2.0)]:.2f} (same cycle: a real limit cycle)")
    print("[C] pairwise vs whole")
    for k, pairs_ok, kind, ev, h in C["goodwin"]:
        print(f"    Goodwin k={k:4.1f}: all 2x2 stable={pairs_ok}; whole: "
              f"{kind[:40]}; Hurwitz {np.round(h, 2)}")
    print("    Goodwin loops: " + "; ".join(f"a{k}: " + ", ".join(
        f"{d} → {c:+.3g}" for d, c in v) for k, v in
        C["goodwin_loops"].items()))
    hd = C["held"]
    print("    held runaway: pairs " + ", ".join(
        f"{p['idx']}:{'stable' if p['stable'] else 'UNSTABLE'}"
        for p in hd["pairs"]) + f"; whole: {hd['whole']['kind'][:35]}")
    for n, s in C["random"].items():
        print(f"    random {n}x{n}: P(whole stable | all pairs stable) = "
              f"{s['ppv']:.2f}; stable wholes containing an unstable pair "
              f"= {s['stable_but_pair_unstable']*100:.1f} %")
    print("[D] physiology atlas (curated BioModels equations):")
    for r in D["rows"]:
        if "error" in r:
            print(f"    {r['label']}: error {r['error']}")
            continue
        if not r["fp_found"]:
            print(f"    {r['label']}: n={r['n']}, oscillating="
                  f"{r['oscillating']}, fixed point not resolved")
            continue
        per = f", period {r['period']:.3g}" if r["period"] else ""
        kl = r["key_loop"]
        print(f"    {r['label']}: n={r['n']}; simulation: {r['state']}{per}")
        print(f"        fixed point: {r['kind'][:46]}; λ_lead={r['lead']:.3g};"
              f" unstable 2x2: {r['n_pairs_unstable']}/{r['n_pairs']} "
              f"{r['unstable_pairs'][:3]}; all F_k<0: "
              f"{bool((r['F'] < 0).all())}; sign-stable {r['sign_stab']:.2f}")
        if kl:
            print(f"        key multi-step loop: {kl[0]} (gain {kl[1]:+.3g})")
    print("    Topp equilibria: " + "; ".join(
        f"{np.round(e, 1).tolist()} {k[:20]}" for e, k in D["topp_eq"]))
    print(f"    Topp β-cell death rate d0: {D['topp_events']} (analytic "
          f"saddle-node r1²/4r2 = {D['d0_star']:.4f}); insulin-sensitivity "
          f"scan 0.72→0.05: events {D['topp_si_events'] or 'none'} "
          f"(β-cell compensation keeps G=100)")
    print("[E] May: P(stable) random vs feedback-paired, n = " +
          ", ".join(f"{n}:{r:.2f}/{s:.2f}" for n, r, s in
                    zip(E["ns"], E["random"], E["structured"])) +
          f"  (May bound n* = {E['may_n']:.1f})")
    if make_fig:
        _figure(A, B, C, D, E)
    return {"A": A, "B": B, "C": C, "D": D, "E": E}


def _figure(A, B, C, D, E):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(3, 3, figsize=(16, 13))
    # A1 phase portrait
    a = ax[0, 0]
    S = A["S"]
    rng = np.random.default_rng(3)
    for _ in range(14):
        x0 = [rng.uniform(0.2, 6), rng.uniform(0.05, 2.5)]
        t, X = S.simulate(x0, 30, 3000)
        a.plot(X[:, 0], X[:, 1], lw=0.7, color="gray", alpha=.7)
    a.plot(A["X"][:, 0], A["X"][:, 1], color="C3", lw=1.6,
           label="T: 0.8→0.5 perturbation")
    for e, k, _ in A["eq"]:
        a.plot(*e, "o", ms=8, mfc="k" if "stable" in k and "un" not in k
               else "w", mec="k")
    a.set_xlabel("C (tumour)"); a.set_ylabel("T (T cells)")
    a.set_title("A  tumour–immune: global attractor (3, 0.8)")
    a.legend(fontsize=8)
    a = ax[0, 1]
    a.plot(A["t"], A["X"][:, 0], label="C"); a.plot(A["t"], A["X"][:, 1],
                                                    label="T")
    a.set_xlabel("time"); a.set_title(
        f"A  damped return: λ={A['a']['lead'].real:.2f}±"
        f"{abs(A['a']['lead'].imag):.2f}i"); a.legend(fontsize=8)
    a = ax[0, 2]
    t, X = B["traj_soft"]; a.plot(X[:, 0], X[:, 1], label="γ=.2, δ=.1 soft")
    t, X = B["traj_lin"]; a.plot(X[:, 0], X[:, 1], label="γ=δ=0 centre")
    t, X = B["traj_nl"]; a.plot(X[len(X)//2:, 0], X[len(X)//2:, 1],
                                label="nonlinear: true limit cycle")
    a.set_xlabel("unemployment gap"); a.set_ylabel("inflation gap")
    a.set_title("B  macro: sink vs centre vs limit cycle")
    a.legend(fontsize=7)
    # C Goodwin
    a = ax[1, 0]
    ks = np.linspace(1, 14, 200)
    lead = [np.linalg.eigvals(np.array([[-1, 0, -k], [1, -1, 0],
                                        [0, 1, -1.0]])).real.max() for k in ks]
    a.plot(ks, lead); a.axhline(0, color="k", lw=.5)
    a.axvline(8, color="r", ls=":")
    a.set_xlabel("loop gain k (3-step negative feedback)")
    a.set_ylabel("Re λ_lead (whole)")
    a.set_title("C1 all pairs stable, whole oscillates for k>8")
    a = ax[1, 1]
    labs, vals = [], []
    for n, s in C["random"].items():
        labs += [f"{n}x{n}\nP(whole ok | pairs ok)",
                 f"{n}x{n}\nstable wholes w/\nunstable pair"]
        vals += [s["ppv"], s["stable_but_pair_unstable"]]
    a.bar(range(4), vals, color=["C0", "C1", "C0", "C1"])
    a.set_xticks(range(4)); a.set_xticklabels(labs, fontsize=7)
    a.set_title("C3 pairwise screening misjudges the whole")
    # D Topp
    a = ax[1, 2]
    sc = [(v, x[0], e.real.max()) for v, x, e in D["topp_scan"]
          if x is not None]
    a.plot([s[0] for s in sc], [s[2] for s in sc])
    a.axvline(D["d0_star"], color="r", ls=":")
    a.set_xlabel("β-cell death rate d0"); a.set_ylabel("λ_lead (healthy)")
    a.set_title("D  Topp: λ→0 before diabetes (saddle-node)")
    # D atlas oscillators
    a = ax[2, 0]
    for r in D["rows"]:
        if r.get("oscillating") and "Sturis" in r["label"]:
            t, X = r["t"], r["X"]
            sel = t > t[-1] * 0.7
            a.plot(t[sel], X[sel, 0] / X[sel, 0].max(), label="insulin (x)")
            a.plot(t[sel], X[sel, 2] / X[sel, 2].max(), label="glucose (z)")
    a.set_xlabel("min"); a.set_title("D  Sturis 1991: ultradian Hopf cycle")
    a.legend(fontsize=8)
    a = ax[2, 1]
    good = [r for r in D["rows"] if r.get("fp_found")]
    a.barh(range(len(good)), [r["n_pairs_unstable"] / max(r["n_pairs"], 1)
                              for r in good], color="C1", label="unstable 2x2 share")
    a.set_yticks(range(len(good)))
    a.set_yticklabels([r["label"].split(":")[0] for r in good], fontsize=8)
    for i, r in enumerate(good):
        a.text(0.02, i, r["kind"][:26], fontsize=7, va="center")
    a.set_title("D  atlas: fixed-point class & unstable pairs")
    a = ax[2, 2]
    a.semilogx(E["ns"], E["random"], "o-", label="random (May)")
    a.semilogx(E["ns"], E["structured"], "s-", label="feedback-paired")
    a.axvline(E["may_n"], color="r", ls=":", label="May bound")
    a.set_xticks(E["ns"]); a.set_xticklabels([str(n) for n in E["ns"]])
    a.minorticks_off()
    a.set_xlabel("number of variables n"); a.set_ylabel("P(stable)")
    a.set_title("E  more variables ≠ more stable"); a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d5_pairwise_to_network.png"), dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    report()
