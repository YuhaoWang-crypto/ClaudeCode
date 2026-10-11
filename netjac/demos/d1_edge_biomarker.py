"""
D1 — Edge biomarkers: same concentrations, different Jacobian.

System: Teusink 2000 yeast glycolysis, EXACT curated kinetics (BIOMD0000000064).
Teusink B et al. Eur J Biochem 267:5313 (2000).

Phenotype ("disease"): PFK loses its allosteric ATP inhibition (CiPFKATP 100->1;
an analogue of PFKFB3/F2,6BP-driven override in tumours, ⚠️ analogy only) and the
cell compensates by lowering PFK capacity until F6P is back at its healthy
value. Result: every metabolite concentration is identical, but one EDGE of the
network (PFK <- ATP) is gone.

Steps
  A  MCA engine validation: analytic dFlux/dEnzyme vs brute force; theorems.
  B  Healthy vs disease: concentrations, elasticities (edges), eigenvalues.
  C  Virtual cohort (inter-individual enzyme variation, a shared "basal rate"
     factor, 10 % assay CV). Node readouts vs a model-predicted EDGE readout:
     measure F6P and ATP before/after a standard ATP-demand challenge (ATPase
     x0.8) -> response coefficient R = dlnF6P / dlnATP. All data-driven
     selections use a discovery/validation split.
  D  Group-level discovery from snapshot covariance only (no perturbation):
     inverse Lyapunov (Nägele 2014 / COVRECON 2023) recovers the edges and
     points at PFK<-ATP.
"""
import os
import numpy as np
from scipy.optimize import brentq

from netjac import FIG
from netjac.core import KineticModel
from netjac import spectrum as sp, cohort as co, discriminate as di
from netjac import fluctuation as fl

CHALLENGE = ("vATP", 0.8)      # ATP-demand challenge: ATPase activity x0.8
# (x0.6 itself drives ~30 % of individuals into hexose-phosphate runaway: a
#  diagnostic perturbation must stay inside the basin -> keep it small)


def phenotypes(m, x0, ci=1.0):
    """ci = residual ATP-inhibition strength (healthy 100, total loss 1)."""
    base = m.snapshot()
    iF6P = m.state_ids.index("F6P")
    m.set_params({"CiPFKATP": ci})

    def g(f):
        m.set_params({"vPFK_VmPFK": base["vPFK_VmPFK"] * f})
        return np.log(m.steady_state(x0=x0)[iF6P] / x0[iF6P])
    f = brentq(g, 0.5, 1.0, xtol=1e-10)
    m.set_params(base)
    return base, {"CiPFKATP": ci, "vPFK_VmPFK": base["vPFK_VmPFK"] * f}, f


def part_a(m, x):
    c = m.control(x)
    cols, _ = m.control_fd(x)
    diff = max(np.abs(c["CJ_scaled"][:, m.rxn_ids.index(r)] - cols[r]).max()
               for r in cols)
    iglt = m.rxn_ids.index("vGLT")
    fcc = {r: c["CJ_scaled"][iglt, m.rxn_ids.index(r)] for r in m.rxn_ids}
    return {"fcc_flux": fcc, "fd_max_diff": diff,
            "summation": c["summation_residual"],
            "connectivity": c["connectivity_residual"],
            "row_sum": c["CJ_scaled"][iglt].sum()}


def part_b(m, base, dis, x0):
    m.set_params(base)
    xh = m.steady_state(x0=x0)
    eh = m.elasticities(xh)
    vh = m.v(xh)
    ch = sp.classify(m.reduced_jacobian(xh, eh))
    atp_h = m.rr["ATP"]
    m.set_params(dis)
    xd = m.steady_state(x0=x0)
    ed = m.elasticities(xd)
    vd = m.v(xd)
    cd = sp.classify(m.reduced_jacobian(xd, ed))
    m.set_params(base)
    sc = lambda e, v, x: e * x[None, :] / np.where(v == 0, 1, v)[:, None]
    dE = sc(ed, vd, xd) - sc(eh, vh, xh)
    order = [i for i in np.argsort(-np.abs(dE).ravel())[:6]
             if np.abs(dE).ravel()[i] > 1e-3 * np.abs(dE).max()]
    edges = [(m.rxn_ids[r], m.state_ids[k], sc(eh, vh, xh)[r, k],
              sc(ed, vd, xd)[r, k])
             for r, k in zip(*np.unravel_index(order, dE.shape))]
    return {"xh": xh, "xd": xd, "max_dlnx": np.abs(np.log(xd / xh)).max(),
            "lam_h": ch["lambda_lead"], "lam_d": cd["lambda_lead"],
            "kind_h": ch["kind"], "kind_d": cd["kind"], "edges": edges,
            "flux_h": vh[m.rxn_ids.index("vGLT")],
            "flux_d": vd[m.rxn_ids.index("vGLT")], "atp": atp_h}


def _simulate_group(m, base, pheno, x0, n, rng, cv_enzyme, cv_global):
    inds = co.make_individuals(m, n, rng, cv_enzyme, cv_global)
    out, fails = [], 0
    rid, fac = CHALLENGE
    for ind in inds:
        try:
            with co.applied(m, base, pheno, ind):
                xa = m.steady_state(x0=x0, max_rounds=3)
                pa = [m.rr[s] for s in ("ATP", "ADP", "AMP")]
                m.scale_enzyme(rid, fac)
                xb = m.steady_state(x0=xa, max_rounds=3)
                pb = [m.rr[s] for s in ("ATP", "ADP", "AMP")]
            out.append((np.r_[xa, pa], np.r_[xb, pb]))
        except RuntimeError:
            fails += 1
    pre = np.array([o[0] for o in out])
    post = np.array([o[1] for o in out])
    return pre, post, fails


def part_c(m, base, dis, x0, n=160, seed=7, cv_enzyme=0.10, cv_global=0.10,
           cv_assay=0.10):
    rng = np.random.default_rng(seed)
    names = m.state_ids + ["ATP", "ADP", "AMP"]
    Hpre, Hpost, fh = _simulate_group(m, base, {}, x0, n, rng, cv_enzyme,
                                      cv_global)
    Dpre, Dpost, fd = _simulate_group(m, base, dis, x0, n, rng, cv_enzyme,
                                      cv_global)
    # assays (independent noise on every measurement)
    A = lambda X: co.assay(X, rng, cv_assay)
    Hpre_m, Hpost_m, Dpre_m, Dpost_m = A(Hpre), A(Hpost), A(Dpre), A(Dpost)
    keep = [k for k, s in enumerate(names) if s != "BPG"]  # ~0 nM, unmeasurable
    # discovery / validation split
    sh = rng.permutation(len(Hpre)); sd = rng.permutation(len(Dpre))
    hd, hv = sh[: len(sh) // 2], sh[len(sh) // 2:]
    dd, dv = sd[: len(sd) // 2], sd[len(sd) // 2:]

    def pick_and_validate(featH, featD, labels):
        disc = [di.auc(featD[dd, j], featH[hd, j], n_boot=1)["auc"]
                for j in range(featH.shape[1])]
        j = int(np.argmax(disc))
        val = di.auc(featD[dv, j], featH[hv, j])
        return labels[j], disc[j], val

    # 1) single static nodes
    LH, LD = np.log(Hpre_m[:, keep]), np.log(Dpre_m[:, keep])
    node_auc = {names[k]: di.auc(LD[:, i], LH[:, i], n_boot=1)["auc"]
                for i, k in enumerate(keep)}
    best_node = pick_and_validate(LH, LD, [names[k] for k in keep])
    # 2) all static pairwise log-ratios
    pairs = [(i, j) for i in range(len(keep)) for j in range(i + 1, len(keep))]
    RH = np.array([LH[:, i] - LH[:, j] for i, j in pairs]).T
    RD = np.array([LD[:, i] - LD[:, j] for i, j in pairs]).T
    best_ratio = pick_and_validate(
        RH, RD, [f"{names[keep[i]]}/{names[keep[j]]}" for i, j in pairs])
    # 3) full static panel, cross-validated logistic regression
    X = np.vstack([LH, LD]); y = np.r_[np.zeros(len(LH)), np.ones(len(LD))]
    panel_auc = di.cv_auc(X, y)
    # 4) pre-specified edge readout from the model (NO data selection)
    iF, iA = names.index("F6P"), names.index("ATP")

    def R(pre, post):
        return (np.log(post[:, iF] / pre[:, iF])
                / np.log(post[:, iA] / pre[:, iA]))
    RH_e, RD_e = R(Hpre_m, Hpost_m), R(Dpre_m, Dpost_m)
    edge = di.auc(RD_e, RH_e)
    # 5) noise-free ceiling of the same readout
    edge_true = di.auc(R(Dpre, Dpost), R(Hpre, Hpost), n_boot=200)
    # 6) data-driven best challenge delta (discovery/validation)
    dH = np.log(Hpost_m[:, keep] / Hpre_m[:, keep])
    dD = np.log(Dpost_m[:, keep] / Dpre_m[:, keep])
    best_delta = pick_and_validate(dH, dD, [f"dln {names[k]}" for k in keep])
    # 7) fragility side-result: collapsed individuals per group
    return {"node_auc": node_auc, "best_node": best_node,
            "best_ratio": best_ratio, "panel_auc": panel_auc,
            "edge": edge, "edge_true": edge_true, "best_delta": best_delta,
            "R_h": RH_e, "R_d": RD_e, "fails": (fh, fd), "n": n,
            "F6P_h": LH[:, keep.index(iF)], "F6P_d": LD[:, keep.index(iF)]}


def severity_sweep(m, x0, cis=(60, 30, 10, 3), n=100, seed=11):
    """Partial loss of the feedback: AUC of the edge readout vs the static
    18-analyte panel as the phenotype gets milder (⚠️ synthetic)."""
    out = []
    for ci in cis:
        base, dis, _ = phenotypes(m, x0, ci)
        c = part_c(m, base, dis, x0, n=n, seed=seed)
        out.append((ci, c["edge"]["auc"], c["panel_auc"],
                    c["best_node"][2]["auc"]))
    return out


def part_d(m, base, dis, x0, Ns=(1000, 10000, 100000, 1000000), omega=1e5,
           seed=3):
    rng = np.random.default_rng(seed)
    res = {}
    for lbl, ph in (("healthy", base), ("disease", dis)):
        m.set_params(base); m.set_params(ph)
        x = m.steady_state(x0=x0)
        e = m.elasticities(x)
        mask = fl.sparsity(m, x, e)
        v = m.v(x)
        C_exact, _ = fl.covariance(m, x, omega, e)
        e_exact, _ = fl.inverse_elasticities(m, C_exact, v, mask, omega)
        errs, est = [], {}
        for N in Ns:
            S = fl.sample_lna(m, x, N, omega, rng, e)
            Ch = np.cov(S.T)
            eh, _ = fl.inverse_elasticities(m, Ch, v, mask, omega)
            errs.append(np.linalg.norm(eh - e) / np.linalg.norm(e))
            est[N] = eh
        # conditioning: error from a 1e-4 relative perturbation of exact C
        P = rng.standard_normal(C_exact.shape); P = (P + P.T) / 2
        Cp = C_exact + 1e-4 * np.linalg.norm(C_exact) * P / np.linalg.norm(P)
        ep, _ = fl.inverse_elasticities(m, Cp, v, mask, omega)
        lam = np.linalg.eigvals(m.reduced_jacobian(x, e)).real
        res[lbl] = {"x": x, "v": v, "eps": e, "mask": mask, "errs": errs,
                    "est": est, "exact_err": np.linalg.norm(e_exact - e)
                    / np.linalg.norm(e),
                    "err_1e4": np.linalg.norm(ep - e) / np.linalg.norm(e),
                    "stiffness": np.abs(lam).max() / np.abs(lam).min()}
    m.set_params(base)
    # differential Jacobian from the largest-N estimate: rank edges
    Nmax = Ns[-1]
    sc = lambda r, E: E * r["x"][None, :] / np.where(r["v"] == 0, 1,
                                                     r["v"])[:, None]
    dhat = sc(res["disease"], res["disease"]["est"][Nmax]) - \
        sc(res["healthy"], res["healthy"]["est"][Nmax])
    mask = res["healthy"]["mask"] | res["disease"]["mask"]
    dhat[~mask] = 0
    order = np.argsort(-np.abs(dhat).ravel())[:5]
    top = [(m.rxn_ids[r], m.state_ids[k], dhat[r, k])
           for r, k in zip(*np.unravel_index(order, dhat.shape))]
    return {"Ns": Ns, "res": res, "top_recovered": top}


def report(make_fig=True):
    m = KineticModel("BIOMD0000000064", "Teusink2000")
    x0 = m.steady_state()
    base, dis, fac = phenotypes(m, x0)
    A = part_a(m, x0)
    B = part_b(m, base, dis, x0)
    C = part_c(m, base, dis, x0)
    D = part_d(m, base, dis, x0)
    SW = severity_sweep(m, x0)

    print("=" * 68)
    print("D1  EDGE BIOMARKERS — Teusink 2000 glycolysis (BIOMD0000000064)")
    print("=" * 68)
    print(f"[A] dFlux/dEnzyme analytic vs brute-force max|diff| = "
          f"{A['fd_max_diff']:.1e}; summation residual {A['summation']:.1e}; "
          f"connectivity residual {A['connectivity']:.1e}  ✅")
    top = sorted(A["fcc_flux"].items(), key=lambda kv: -abs(kv[1]))[:4]
    print("    flux control on glucose uptake: " +
          ", ".join(f"{r} {c:+.3f}" for r, c in top) +
          f"  (sum = {A['row_sum']:.4f})")
    print(f"[B] disease = PFK loses ATP inhibition, VmPFK x{fac:.4f} "
          f"compensation")
    print(f"    max |dln concentration| over all 16 species = "
          f"{B['max_dlnx']:.1e}; flux {B['flux_h']:.2f} vs {B['flux_d']:.2f}  ✅")
    print(f"    leading eigenvalue {B['lam_h'].real:.3f} -> "
          f"{B['lam_d'].real:.3f} 1/min")
    for r, k, a, b in B["edges"][:4]:
        print(f"    edge {r:7s}<-{k:6s} scaled elasticity {a:+.3f} -> {b:+.3f}")
    print(f"[C] virtual cohort n={C['n']}/group (enzyme CV 10 %, shared CV "
          f"10 %, assay CV 10 %); collapsed (no steady state): "
          f"healthy {C['fails'][0]}, disease {C['fails'][1]}  ⚠️ synthetic")
    bn, bnd, bnv = C["best_node"]
    br, brd, brv = C["best_ratio"]
    bd, bdd, bdv = C["best_delta"]
    print(f"    best single node  ({bn}): validation AUC {bnv['auc']:.2f} "
          f"[{bnv['ci'][0]:.2f},{bnv['ci'][1]:.2f}]")
    print(f"    best static ratio ({br}): validation AUC {brv['auc']:.2f} "
          f"[{brv['ci'][0]:.2f},{brv['ci'][1]:.2f}]  (of "
          f"{len(m.state_ids) + 2} choose 2)")
    print(f"    full 18-analyte static panel, 5-fold CV logistic AUC "
          f"{C['panel_auc']:.2f}")
    e = C["edge"]
    print(f"    PRE-SPECIFIED edge readout R=dlnF6P/dlnATP after challenge: "
          f"AUC {e['auc']:.2f} [{e['ci'][0]:.2f},{e['ci'][1]:.2f}] "
          f"(noise-free {C['edge_true']['auc']:.2f}); median R healthy "
          f"{np.median(C['R_h']):+.2f} vs disease {np.median(C['R_d']):+.2f}")
    print(f"    data-driven best challenge delta ({bd}): validation AUC "
          f"{bdv['auc']:.2f}")
    print("    severity sweep (residual ATP-inhibition Ci; healthy=100): " +
          "; ".join(f"Ci={ci}: edge {e:.2f} / panel {p:.2f} / best node "
                    f"{b:.2f}" for ci, e, p, b in SW))
    print(f"[D] inverse Jacobian from snapshot covariance (LNA, Omega=1e5):")
    for lbl in ("healthy", "disease"):
        r = D["res"][lbl]
        print(f"    {lbl:8s} stiffness |λ|max/|λ|min = {r['stiffness']:.1e}; "
              f"1e-4 relative error in C -> {r['err_1e4']:.2f} error in ε")
        print(f"    {lbl:8s} exact-C error {r['exact_err']:.1e}; sample error "
              + ", ".join(f"N={N:.0e}:{err:.2f}" for N, err in
                          zip(D["Ns"], r["errs"])))
    print("    top recovered differential edges: " +
          "; ".join(f"{r}<-{k} {d:+.2f}" for r, k, d in D["top_recovered"]))
    if make_fig:
        _figure(m, A, B, C, D, SW)
    return {"A": A, "B": B, "C": C, "D": D, "sweep": SW, "compensation": fac}


def _figure(m, A, B, C, D, SW):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    a = ax[0, 0]
    idx = np.arange(len(m.state_ids))
    a.bar(idx - 0.2, np.log10(B["xh"]), 0.4, label="healthy")
    a.bar(idx + 0.2, np.log10(B["xd"]), 0.4, label="lost PFK-ATP feedback")
    a.set_xticks(idx); a.set_xticklabels(m.state_ids, rotation=70, fontsize=8)
    a.set_ylabel("log10 concentration (mM)")
    a.set_title(f"A  identical metabolome (max |Δln x| = {B['max_dlnx']:.0e})")
    a.legend(fontsize=8)
    a = ax[0, 1]
    labs = [f"{r}←{k}" for r, k, _, _ in B["edges"][:6]]
    hv = [h for _, _, h, _ in B["edges"][:6]]
    dv = [d for _, _, _, d in B["edges"][:6]]
    yy = np.arange(len(labs))
    a.barh(yy + 0.2, hv, 0.4, label="healthy")
    a.barh(yy - 0.2, dv, 0.4, label="disease")
    a.set_yticks(yy); a.set_yticklabels(labs, fontsize=8)
    a.axvline(0, color="k", lw=0.5)
    a.set_xlabel("scaled elasticity (edge strength)")
    a.set_title("B  …but the network edges differ"); a.legend(fontsize=8)
    a = ax[1, 0]
    names = ["best node\n(" + C["best_node"][0] + ")",
             "best static\nratio", "18-analyte\npanel (CV)",
             "edge readout\nR=ΔlnF6P/ΔlnATP"]
    vals = [C["best_node"][2]["auc"], C["best_ratio"][2]["auc"],
            C["panel_auc"], C["edge"]["auc"]]
    lo = [C["best_node"][2]["ci"][0], C["best_ratio"][2]["ci"][0], np.nan,
          C["edge"]["ci"][0]]
    hi = [C["best_node"][2]["ci"][1], C["best_ratio"][2]["ci"][1], np.nan,
          C["edge"]["ci"][1]]
    cols = ["#9aa", "#9aa", "#9aa", "#c44"]
    a.bar(range(4), vals, color=cols)
    for i in (0, 1, 3):
        a.plot([i, i], [lo[i], hi[i]], color="k")
    a.axhline(0.5, ls="--", color="k", lw=0.7)
    a.set_xticks(range(4)); a.set_xticklabels(names, fontsize=8)
    a.set_ylim(0.4, 1.0); a.set_ylabel("validation AUC (10 % assay CV)")
    a.set_title("C  node vs edge readouts (virtual cohort ⚠️)")
    a = ax[1, 1]
    for lbl, c in (("healthy", "C0"), ("disease", "C1")):
        a.loglog(D["Ns"], D["res"][lbl]["errs"], "o-", color=c, label=lbl)
    a.set_xlabel("snapshot samples N"); a.set_ylabel("relative error of ε̂")
    a.axhline(0.1, ls=":", color="gray")
    a.set_title("D  passive inverse Jacobian needs ~1e6 snapshots (stiff)")
    a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d1_edge_biomarker.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
