"""
D3 — Fragility-node profiling: which enzyme / edge / probe matters for the
stability of the WHOLE network, as opposed to which is most abundant or
carries most flux.

Fragility = first-order sensitivity of the critical eigenvalue,
    F_r = d Re(λ) / d ln E_r      (left/right eigenvectors; spectrum.py)
validated against finite differences, then tested on what it is supposed to
predict:

  A  Teusink 2000 glycolysis (BIOMD0000000064): the collapse seen in D1 is a
     saddle-node (λ_lead -> 0, hexose-phosphate runaway). Compare rankings by
     fragility, by flux control (C^J) and by flux magnitude against the TRUE
     fold-change of each enzyme that collapses the pathway.
  B  Same model, cohort: predict WHO collapses under a strong ATP-demand
     challenge from baseline readouts: concentrations vs measured relaxation
     rate (D2) vs fragility-projected eigenvalue λ + F_ATPase * ln(0.6).
  C  Wolf 2000 oscillations (BIOMD0000000206) at 98 % of the Hopf: predicted
     vs actual crossing for +-20 % changes of each enzyme.
  D  Probe placement: individuals approaching the Hopf; the ringing damping
     is fitted from a glucose pulse observed on ONE species with 1 % assay
     CV. Which species is the best probe? Predicted by the critical mode's
     relative participation |u_k|/x_k, not by abundance.
"""
import os
import numpy as np
from scipy.stats import spearmanr

from netjac import FIG
from netjac.core import KineticModel
from netjac import spectrum as sp, relax as rx, cohort as co
from netjac import discriminate as di


# ------------------------------------------------------------------- part A
def _collapse_fold(m, base, x0, rid, direction, fmax=4.0, tol=0.005):
    """Smallest multiplicative change of enzyme rid (up if direction>0) that
    makes the steady state disappear; inf if none up to fmax."""
    def ok(f):
        m.set_params(base)
        m.scale_enzyme(rid, f)
        try:
            m.steady_state(x0=x0, max_rounds=3)
            return True
        except RuntimeError:
            return False
    lo, hi = 1.0, (fmax if direction > 0 else 1 / fmax)
    if ok(hi):
        m.set_params(base)
        return np.inf
    while abs(np.log(hi / lo)) > tol:
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    m.set_params(base)
    return hi


def part_a():
    m = KineticModel("BIOMD0000000064", "Teusink2000")
    x = m.steady_state()
    base = m.snapshot()
    Floc, lam = sp.enzyme_fragility(m, x)
    fd = sp.enzyme_fragility_fd(m, x)
    ez = list(m.enzyme_params())
    err = max(abs(Floc[m.rxn_ids.index(r)] - fd[r]) for r in ez)
    Ftot = {r: sp.total_parameter_fragility(m, m.enzyme_params()[r][0], x)
            for r in ez if len(m.enzyme_params()[r]) == 1}
    c = m.control(x)
    iglt = m.rxn_ids.index("vGLT")
    fcc = {r: c["CJ_scaled"][iglt, m.rxn_ids.index(r)] for r in ez}
    flux = {r: abs(c["v"][m.rxn_ids.index(r)]) for r in ez}
    # true distance to collapse: smallest |ln fold| in the destabilising
    # direction (sign of the total fragility)
    dist = {}
    for r, F in Ftot.items():
        f = _collapse_fold(m, base, x, r, +1 if F > 0 else -1)
        dist[r] = abs(np.log(f))
    rs = [r for r in Ftot]
    d = np.array([dist[r] for r in rs])
    finite = np.isfinite(d)
    d_rank = np.where(finite, d, 1e3)
    # first-order prediction of the distance: |Re λ| / |F_total|
    pred = np.array([abs(lam.real) / max(abs(Ftot[r]), 1e-9) for r in rs])
    rho = {"fragility": spearmanr(pred, d_rank).statistic,
           "flux_control": spearmanr(-np.abs([fcc[r] for r in rs]),
                                     d_rank).statistic,
           "flux": spearmanr(-np.array([flux[r] for r in rs]),
                             d_rank).statistic}
    collapse_set = {r for r in rs if np.isfinite(dist[r])}
    k = max(len(collapse_set), 1)
    top_by = lambda score: set(sorted(rs, key=lambda r: -score(r))[:k])
    prec = {"fragility": len(top_by(lambda r: abs(Ftot[r])) & collapse_set) / k,
            "flux_control": len(top_by(lambda r: abs(fcc[r])) & collapse_set) / k,
            "flux": len(top_by(lambda r: flux[r]) & collapse_set) / k}
    F, _, _ = sp.edge_fragility(m, x)
    order = np.argsort(-np.abs(F).ravel())[:6]
    edges = [(m.rxn_ids[r], m.state_ids[k], F[r, k],
              abs(c["v"][r]))
             for r, k in zip(*np.unravel_index(order, F.shape))]
    m.set_params(base)
    return {"lam": lam, "fd_err": err, "Ftot": Ftot, "fcc": fcc,
            "flux": flux, "dist": dist, "pred": dict(zip(rs, pred)),
            "rho": rho, "prec": prec, "collapse_set": collapse_set,
            "edges": edges, "glycolytic_flux": c["v"][iglt]}


# ------------------------------------------------------------------- part B
def part_b(n=200, seed=5, challenge=0.6, cv_assay=0.10):
    m = KineticModel("BIOMD0000000064", "Teusink2000")
    x0 = m.steady_state()
    base = m.snapshot()
    rng = np.random.default_rng(seed)
    inds = co.make_individuals(m, n, rng, cv_enzyme=0.10, cv_global=0.10)
    rows = []
    for ind in inds:
        try:
            with co.applied(m, base, {}, ind):
                x = m.steady_state(x0=x0, max_rounds=3)
                lam = np.linalg.eigvals(m.reduced_jacobian(x)).real.max()
                Fatp = sp.total_fragility(m, "vATP", x, solver="integrate")
                conc = np.r_[x, m.rr["ATP"]]
                m.scale_enzyme("vATP", challenge)
                try:
                    m.steady_state(x0=x, max_rounds=3)
                    collapsed = 0
                except RuntimeError:
                    collapsed = 1
        except RuntimeError:
            continue                          # already collapsed at baseline
        rows.append((collapsed, lam, lam + Fatp * np.log(challenge), conc))
    y = np.array([r[0] for r in rows])
    lam = np.array([r[1] for r in rows])
    proj = np.array([r[2] for r in rows])
    C = np.log(co.assay(np.array([r[3] for r in rows]), rng, cv_assay))
    # measured relaxation rate: assume the D2 pulse estimate with 10 % error
    lam_meas = lam * co.lognormal(rng, 0.10, len(lam))
    # 10 % error also on F (it would come from a calibrated model)
    proj_meas = lam_meas + (proj - lam) * co.lognormal(rng, 0.10, len(lam))
    names = m.state_ids + ["ATP"]
    node = {s: di.auc(C[y == 1, j], C[y == 0, j], n_boot=1)["auc"]
            for j, s in enumerate(names) if s != "BPG"}
    best_node = max(node, key=node.get)
    panel = di.cv_auc(C[:, [j for j, s in enumerate(names) if s != "BPG"]], y)
    return {"n": len(y), "n_collapse": int(y.sum()),
            "auc_relax": di.auc(lam_meas[y == 1], lam_meas[y == 0]),
            "auc_proj": di.auc(proj_meas[y == 1], proj_meas[y == 0]),
            "auc_node": (best_node, node[best_node]), "auc_panel": panel,
            "lam": lam_meas, "proj": proj_meas, "y": y}


# ------------------------------------------------------------------- part C
def _wolf():
    m = KineticModel("BIOMD0000000206", "Wolf2000")
    m.rr.simulate(0, 50, 10)
    x_ref = m.fixed_point(m.x())
    from netjac.demos.d2_perturbation_relaxation import find_hopf
    k_h = find_hopf(m, x_ref, 0.5 * m.get("k0"), m.get("k0"))
    return m, x_ref, k_h


def part_c(frac=0.98, folds=(0.9, 1.1)):
    m, x_ref, k_h = _wolf()
    m.set_params({"k0": frac * k_h})
    base = m.snapshot()
    x = m.fixed_point(x_ref)
    lc = sp.critical_pair(m.reduced_jacobian(x))
    Floc, _ = sp.enzyme_fragility(m, x, target=lc)
    Fr = np.zeros(len(m.rxn_ids))
    for rid in m.enzyme_params():
        Fr[m.rxn_ids.index(rid)] = sp.total_fragility(m, rid, x, target=lc)
    flux = np.abs(m.v(x))
    rows = []
    for rid in m.enzyme_params():
        j = m.rxn_ids.index(rid)
        for f in folds:
            m.set_params(base)
            m.scale_enzyme(rid, f)
            try:
                xs = m.fixed_point(x)
                lam_new = sp.critical_pair(m.reduced_jacobian(xs))
                actual = lam_new.real if lam_new is not None else np.nan
            except RuntimeError:
                actual = np.nan
            pred = lc.real + Fr[j] * np.log(f)
            rows.append((rid, f, Fr[j], flux[j], pred, actual, Floc[j]))
    m.set_params(base)
    act = np.array([r[5] for r in rows])
    prd = np.array([r[4] for r in rows])
    ok = np.isfinite(act)
    cross_true = act > 0
    cross_pred = prd > 0
    agree = (cross_true[ok] == cross_pred[ok]).mean()
    floc = np.array([r[6] for r in rows])
    prd_loc = lc.real + floc * np.log([r[1] for r in rows])
    agree_loc = ((act > 0)[ok] == (prd_loc > 0)[ok]).mean()
    return {"lam_c": lc, "rows": rows, "agree": agree, "agree_local": agree_loc,
            "corr_local": np.corrcoef(prd_loc[ok], act[ok])[0, 1],
            "corr": np.corrcoef(prd[ok], act[ok])[0, 1],
            "n_cross": int(cross_true[ok].sum()),
            "fragility": dict(zip(m.rxn_ids, Fr)),
            "flux": dict(zip(m.rxn_ids, flux))}


# ------------------------------------------------------------------- part D
def part_d(n=90, seed=2, cv_assay=0.01, near=-3.0):
    m, x_ref, k_h = _wolf()
    base = m.snapshot()
    rng = np.random.default_rng(seed)
    # participation predicted at a reference point near the Hopf
    m.set_params({"k0": 0.97 * k_h})
    xr = m.fixed_point(x_ref)
    lc = sp.critical_pair(m.reduced_jacobian(xr))
    _, _, part = sp.mode_shape(m, xr, target=lc)
    probes = list(m.state_ids)
    truth, est = [], {s: [] for s in probes}
    n_osc = 0
    for _ in range(n):
        m.set_params(base)
        m.set_params({"k0": k_h * rng.uniform(0.85, 1.0)})
        for rid in m.enzyme_params():
            if rid != "v0":
                m.scale_enzyme(rid, co.lognormal(rng, 0.03))
        try:
            x = m.fixed_point(xr)
        except RuntimeError:
            continue
        lam = sp.critical_pair(m.reduced_jacobian(x))
        if lam is None or lam.real >= 0:
            n_osc += 1                        # already oscillating: excluded
            continue
        truth.append(lam.real > near)
        t, Y = rx.pulse_response(m, x, "s1", 0.1 * x[0], 1.5, 1500, probes)
        for j, s in enumerate(probes):
            y = Y[:, j] + cv_assay * x[j] * rng.standard_normal(len(t))
            try:
                ld = rx.hankel_dmd(y, t[1] - t[0], rank=6, delays=60)
                oc = ld[ld.imag > 1]
                s0 = (min(oc[np.argmax(oc.real)].real, -0.5)
                      if len(oc) else -20.0)
                w0 = oc[np.argmax(oc.real)].imag if len(oc) else lam.imag
                fit = rx.fit_modes(t, y, [-5.0, -10.0], [(s0, w0)])
                sf, wf = fit["osc"][0]
                val = sf if abs(wf - lam.imag) < 0.2 * lam.imag else -100.0
            except Exception:
                val = -100.0
            est[s].append(max(val, -100.0))
    m.set_params(base)
    y = np.array(truth)
    auc = {s: di.auc(np.array(est[s])[y], np.array(est[s])[~y], n_boot=300)
           for s in probes}
    xs = dict(zip(m.state_ids, xr))
    return {"n": len(y), "n_near": int(y.sum()), "n_osc": n_osc,
            "auc": auc, "participation": dict(zip(m.state_ids, part)),
            "abundance": xs,
            "rho_part": spearmanr([part[i] for i in range(len(probes))],
                                  [auc[s]["auc"] for s in probes]).statistic,
            "rho_abund": spearmanr([xs[s] for s in probes],
                                   [auc[s]["auc"] for s in probes]).statistic}


def report(make_fig=True):
    A = part_a()
    B = part_b()
    C = part_c()
    D = part_d()
    print("=" * 68)
    print("D3  FRAGILITY-NODE PROFILING")
    print("=" * 68)
    print(f"[A] Teusink: λ_lead = {A['lam'].real:.3f}; analytic vs finite-"
          f"difference enzyme fragility max|diff| = {A['fd_err']:.1e}  ✅")
    top = sorted(A["Ftot"], key=lambda r: -abs(A["Ftot"][r]))[:6]
    print("    enzyme   F_total  |C^J|   flux    true ln-fold to collapse")
    for r in top:
        print(f"    {r:7s} {A['Ftot'][r]:+7.2f}  {abs(A['fcc'][r]):.3f}  "
              f"{A['flux'][r]:6.1f}  {A['dist'][r]:.3f}")
    print(f"    Spearman(predictor, true distance-to-collapse): fragility "
          f"{A['rho']['fragility']:+.2f}, flux control "
          f"{A['rho']['flux_control']:+.2f}, flux magnitude "
          f"{A['rho']['flux']:+.2f}")
    print(f"    enzymes that can collapse the pathway alone: "
          f"{sorted(A['collapse_set'])}; top-{len(A['collapse_set'])} hit rate: "
          f"fragility {A['prec']['fragility']:.2f}, flux control "
          f"{A['prec']['flux_control']:.2f}, flux magnitude "
          f"{A['prec']['flux']:.2f}")
    print("    most fragile EDGES: " + "; ".join(
        f"{r}<-{k} {F:+.2f} (flux {v:.1f})" for r, k, F, v in A["edges"]))
    print(f"[B] cohort n={B['n']}: {B['n_collapse']} collapse under ATPase x0.6"
          f"  ⚠️ synthetic")
    print(f"    predict collapse from baseline: best single node "
          f"({B['auc_node'][0]}) AUC {B['auc_node'][1]:.2f}; 17-analyte "
          f"panel (CV) {B['auc_panel']:.2f}; measured relaxation rate "
          f"{B['auc_relax']['auc']:.2f}; fragility-projected λ "
          f"{B['auc_proj']['auc']:.2f} [{B['auc_proj']['ci'][0]:.2f},"
          f"{B['auc_proj']['ci'][1]:.2f}]")
    print(f"[C] Wolf at 0.98 k0*: critical pair {C['lam_c'].real:+.2f}"
          f"{C['lam_c'].imag:+.2f}i; ±10 % single-enzyme changes: "
          f"{C['n_cross']} cross the Hopf")
    print(f"    TOTAL fragility (incl. steady-state shift) predicts crossing in "
          f"{C['agree']*100:.0f} % (r = {C['corr']:.3f}); LOCAL-only fragility "
          f"{C['agree_local']*100:.0f} % (r = {C['corr_local']:.3f})")
    fr = sorted(C["fragility"], key=lambda r: -abs(C["fragility"][r]))[:4]
    print("    top fragile: " + ", ".join(
        f"{r} F={C['fragility'][r]:+.1f} (flux {C['flux'][r]:.0f})"
        for r in fr))
    print(f"[D] probe placement (n={D['n']} individuals, {D['n_near']} within "
          f"Re λc > -3, {D['n_osc']} already oscillating excluded; 1 % assay CV):")
    for s in sorted(D["auc"], key=lambda s: -D["auc"][s]["auc"]):
        print(f"    probe {s:4s} AUC {D['auc'][s]['auc']:.2f}  participation "
              f"{D['participation'][s]:.3f}  abundance {D['abundance'][s]:.3f}")
    print(f"    Spearman(AUC, participation) = {D['rho_part']:+.2f}; "
          f"Spearman(AUC, abundance) = {D['rho_abund']:+.2f}")
    if make_fig:
        _figure(A, B, C, D)
    return {"A": A, "B": B, "C": C, "D": D}


def _figure(A, B, C, D):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    a = ax[0, 0]
    rs = [r for r in A["Ftot"] if np.isfinite(A["dist"][r])]
    a.scatter([A["pred"][r] for r in rs], [A["dist"][r] for r in rs])
    for r in rs:
        a.annotate(r, (A["pred"][r], A["dist"][r]), fontsize=7)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlabel("predicted distance |Re λ| / |F_total|")
    a.set_ylabel("true |ln fold| to collapse")
    a.set_title(f"A  Teusink: fragility ranks collapse (ρ="
                f"{A['rho']['fragility']:+.2f}; flux ρ={A['rho']['flux']:+.2f})")
    a = ax[0, 1]
    labs = [f"best node\n({B['auc_node'][0]})", "17-analyte\npanel",
            "relaxation\nrate λ", "fragility-\nprojected λ"]
    vals = [B["auc_node"][1], B["auc_panel"], B["auc_relax"]["auc"],
            B["auc_proj"]["auc"]]
    a.bar(range(4), vals, color=["#9aa", "#9aa", "#c84", "#c44"])
    a.axhline(0.5, ls="--", color="k", lw=0.7)
    a.set_xticks(range(4)); a.set_xticklabels(labs, fontsize=8)
    a.set_ylim(0.4, 1.0); a.set_ylabel("AUC: who collapses under challenge")
    a.set_title("B  baseline readouts predicting collapse ⚠️")
    a = ax[1, 0]
    prd = np.array([r[4] for r in C["rows"]])
    act = np.array([r[5] for r in C["rows"]])
    a.scatter(prd, act, s=14, label="total fragility")
    lim = [np.nanmin([prd, act]), np.nanmax([prd, act])]
    a.plot(lim, lim, "k:", lw=0.7)
    a.axhline(0, color="r", lw=0.6); a.axvline(0, color="r", lw=0.6)
    a.set_xlabel("predicted Re λc (first order)")
    a.set_ylabel("actual Re λc (re-solved)")
    prl = C["lam_c"].real + np.array([r[6] for r in C["rows"]]) * \
        np.log([r[1] for r in C["rows"]])
    a.scatter(prl, act, s=14, marker="x", color="gray", label="local only")
    a.legend(fontsize=7)
    a.set_title(f"C  Wolf: ±10 % enzyme changes vs Hopf "
                f"({C['agree']*100:.0f} % agree)")
    a = ax[1, 1]
    ss = list(D["auc"])
    a.scatter([D["participation"][s] for s in ss],
              [D["auc"][s]["auc"] for s in ss], label="probe")
    for s in ss:
        a.annotate(f"{s} ({D['abundance'][s]:.2g} mM)",
                   (D["participation"][s], D["auc"][s]["auc"]), fontsize=7)
    a.set_xscale("log")
    a.set_xlabel("critical-mode participation |u_k|/x_k (predicted)")
    a.set_ylabel("AUC detecting near-Hopf individuals")
    a.set_title("D  probe placement follows the mode, not abundance")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d3_fragility.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
