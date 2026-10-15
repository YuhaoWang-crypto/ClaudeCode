"""
D11  Defining a state by its perturbation-recovery behaviour, validated on
     known positives, and the three control questions.

Framework (from the uploaded planning report, sections 7, 8, 11): nested
models on the same subjects and folds
    M0 static readout
    M1 + one standardised stimulated endpoint
    M2 + a second stimulus strength (dose response)
    M3 + recovery kinetics after the stimulus is removed
and a model step counts only if it reduces misclassification: paired
bootstrap of the AUC / C-index increment and of the false-positive rate at
a fixed sensitivity. Then:
    Q1 which components are the core,
    Q2 the fewest knobs that control the network,
    Q3 whether adjusting them returns the recovery to the original
       equilibrium.

  A  Ground truth (exact model, known defect): Teusink 2000 yeast
     glycolysis (BIOMD64), virtual cells. Positives = loss of the PFK-ATP
     feedback, compensated so the static metabolome is unchanged
     (as in D1); decoy negatives = intact network with a large static shift
     (glucose transporter x1.15). Q1-Q3 on the same model, with the true fix
     (feedback strength CiPFKATP) as the positive control for Q2.
  B  Cells (real): Replogle K562 Perturb-seq. Known positives = the
     textbook transcriptional homeostats (proteasome, tRNA synthetases,
     UPR, sterol, actin) vs co-regulated machines (ribosome, mitoribosome,
     spliceosome, RNA Pol II ...): static co-fluctuation features vs
     perturbation-response features. Q1/Q2 on the 3,426-node genome-wide
     network (core regulators, structural driver nodes). Q3: Norman 2019
     double perturbations: do single-perturbation responses predict which
     second perturbation pulls the state back toward control?
  C  Human physiology (real): NHANES 2005-2016 oral glucose tolerance test
     (a standard stimulus, M1) vs fasting glucose + HbA1c (static, M0),
     known positive = death in the linked mortality file through 2019; and
     CGMacros meal challenges (M3 recovery kinetics) vs fasting glucose for
     type-2 diabetes.

Run:  python3 -m netjac.demos.d11_recovery_state [--part A,B,C]
"""
import argparse
import os
import pickle
import re
import time

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import recovery as rc

RES = os.path.join(DATA_CACHE, "d11_results.pkl")
MODELS = ["M0 static", "M1 +stimulated endpoint", "M2 +dose response", "M3 +recovery kinetics"]


def _load_res():
    return pickle.load(open(RES, "rb")) if os.path.exists(RES) else {}


def _save_res(R):
    pickle.dump(R, open(RES, "wb"))


# ======================================================================= A
OBS = ["F6P", "ATP", "G6P"]               # what a lab would track over time
T_REC = np.array([0.05, 0.1, 0.2, 0.4, 0.8, 1.6])   # s after washout


def _measure(m, x0, chal=("vATP", (0.9, 0.8))):
    """Static state, two challenge strengths, and the recovery after the
    stronger challenge is removed (parameters already applied)."""
    names = m.state_ids + ["ATP", "ADP", "AMP"]
    xa = m.steady_state(x0=x0, max_rounds=3)
    pre = np.r_[xa, [m.rr[s] for s in ("ATP", "ADP", "AMP")]]
    snap = m.snapshot()
    posts = []
    for f in chal[1]:
        m.scale_enzyme(chal[0], f)
        xb = m.steady_state(x0=xa, max_rounds=3)
        posts.append(np.r_[xb, [m.rr[s] for s in ("ATP", "ADP", "AMP")]])
        m.set_params(snap)
    # washout: challenge removed, start from the challenged state
    m.set_x(xb)
    r = m.rr.simulate(0, T_REC[-1], 161, selections=["time"] + OBS)
    tt = r[:, 0]
    rec = np.array([[np.interp(t, tt, r[:, 1 + j]) for t in T_REC]
                    for j in range(len(OBS))])
    m.set_x(xa)
    return names, pre, posts, rec


def _cohort(m, base, groups, x0, n, rng, cv_enzyme=0.10, cv_global=0.10):
    from netjac import cohort as co
    out = []
    for label, (pheno, extra) in groups.items():
        inds = co.make_individuals(m, n, rng, cv_enzyme, cv_global)
        for ind in inds:
            for k, v in extra.items():
                ind[k] = ind.get(k, 1.0) * v
            try:
                with co.applied(m, base, pheno, ind):
                    names, pre, posts, rec = _measure(m, x0)
                out.append((label, pre, posts, rec))
            except RuntimeError:
                pass
    return names, out


def _features(names, data, rng, cv_assay=0.10):
    from netjac import cohort as co
    keep = [k for k, s in enumerate(names) if s != "BPG"]
    A = lambda v: co.assay(v, rng, cv_assay)
    pre = np.log(np.array([A(d[1]) for d in data]))[:, keep]
    p08 = np.log(np.array([A(d[2][1]) for d in data]))[:, keep]
    p09 = np.log(np.array([A(d[2][0]) for d in data]))[:, keep]
    oi = [names.index(o) for o in OBS]
    rec = np.log(np.array([A(d[3]) for d in data]))           # n x obs x t
    base_obs = np.log(np.array([A(d[1][oi]) for d in data]))
    rec_rel = (rec - base_obs[:, :, None]).reshape(len(data), -1)
    F = {MODELS[0]: pre,
         MODELS[1]: np.c_[pre, p08 - pre],
         MODELS[2]: np.c_[pre, p08 - pre, p09 - pre],
         MODELS[3]: np.c_[pre, p08 - pre, p09 - pre, rec_rel]}
    return F


def _signature(m, x0):
    """Noise-free recovery signature used as the restoration target:
    static state + challenge responses + recovery curve (log units)."""
    names, pre, posts, rec = _measure(m, x0)
    keep = [k for k, s in enumerate(names) if s != "BPG"]
    oi = [names.index(o) for o in OBS]
    static = np.log(pre[keep])
    dyn = np.r_[np.log(posts[1][keep]) - static, np.log(posts[0][keep]) - static,
                (np.log(rec) - np.log(pre[oi])[:, None]).ravel()]
    return static, dyn, rec


def part_A(n=150, seed=3):
    from netjac.core import KineticModel
    from netjac import spectrum as sp
    from netjac.demos.d1_edge_biomarker import phenotypes
    rng = np.random.default_rng(seed)
    m = KineticModel("BIOMD0000000064")
    x0 = m.steady_state()
    base, dis, fcomp = phenotypes(m, x0, ci=1.0)
    groups = {"healthy": ({}, {}), "defect": (dis, {}),
              "decoy": ({}, {"vGLT": 1.15})}
    names, data = _cohort(m, base, groups, x0, n, rng)
    lab = np.array([d[0] for d in data])
    y = (lab == "defect").astype(int)
    F = _features(names, data, rng)
    inc = rc.incremental_classification(F, y, n_boot=1000, sens=0.8, reps=3)
    # decoy-specific false positives at 80 % sensitivity
    dec = {}
    for k, r in inc.items():
        thr = np.quantile(r["pred"][y == 1], 0.2)
        dec[k] = {"fp_decoy": float(np.mean(r["pred"][lab == "decoy"] >= thr)),
                  "fp_healthy": float(np.mean(r["pred"][lab == "healthy"] >= thr))}
        r.pop("pred")
    # static separation of the classes (how hidden is the defect?)
    out = {"n": {g: int((lab == g).sum()) for g in groups}, "inc": inc, "decoy": dec,
           "compensation": fcomp}
    # ---- Q1 core components (healthy reference)
    m.set_params(base)
    xh = m.steady_state(x0=x0)
    e = m.elasticities(xh)
    Jz = m.reduced_jacobian(xh, e)
    cls = sp.classify(Jz)
    out["eig_h"] = cls["eigenvalues"]
    m.set_params(dis)
    xd = m.steady_state(x0=x0)
    out["eig_d"] = sp.classify(m.reduced_jacobian(xd, m.elasticities(xd)))["eigenvalues"]
    m.set_params(base)
    # total fragility of the ringing (complex) pair and of the slowest mode
    knobs = list(m.enzyme_params())
    frag_slow, frag_ring = {}, {}
    ring = [l for l in cls["eigenvalues"] if abs(l.imag) > 1e-6]
    ring_t = ring[0] if ring else None
    for rid in knobs:
        try:
            frag_slow[rid] = float(np.real(sp.total_fragility(m, rid, xh, solver="integrate")))
            if ring_t is not None:
                frag_ring[rid] = float(np.real(sp.total_fragility(m, rid, xh, target=ring_t,
                                                                  solver="integrate")))
        except Exception:
            pass
        m.set_params(base)
    out["frag_slow"] = sorted(frag_slow.items(), key=lambda kv: -abs(kv[1]))
    out["frag_ring"] = sorted(frag_ring.items(), key=lambda kv: -abs(kv[1]))
    # edge fragility: which elasticity entries move the ringing pair most
    F_e, _, _ = sp.edge_fragility(m, xh, target=ring_t)
    E = np.abs(F_e)
    order = np.argsort(-E.ravel())
    out["edge_top"] = [(m.rxn_ids[i], m.state_ids[j], float(E[i, j]))
                       for i, j in zip(*np.unravel_index(order[:8], E.shape))]
    ip = m.rxn_ids.index("vPFK")
    pfk_rank = {m.state_ids[j]: int(np.where(order == ip * E.shape[1] + j)[0][0]) + 1
                for j in range(E.shape[1]) if E[ip, j] > 0}
    out["pfk_edge_ranks"] = pfk_rank
    out["n_edges"] = int((E > 0).sum())
    # ---- Q2/Q3 restoration on the defect phenotype (no individual noise)
    m.set_params(base)
    s_h, d_h, rec_h = _signature(m, x0)
    m.set_params(dis)
    s_d, d_d, rec_d = _signature(m, x0)
    snap_dis = m.snapshot()
    m.set_params(base)
    w_dyn = 1.0

    def make_obj(use_dyn):
        def obj(adj):
            m.set_params(snap_dis)
            for k, v in adj.items():
                if k == "CiPFKATP":
                    m.set_params({"CiPFKATP": snap_dis["CiPFKATP"] * np.exp(v)})
                else:
                    m.scale_enzyme(k, float(np.exp(v)))
            try:
                s, d, _ = _signature(m, x0)
            except Exception:
                m.set_params(base)
                return 1e3
            m.set_params(base)
            val = np.mean((s - s_h) ** 2)
            if use_dyn:
                val += w_dyn * np.mean((d - d_h) ** 2)
            return float(val)
        return obj
    # candidate knobs: the 6 enzymes with the largest ringing-pair fragility
    # + the true fix (feedback strength) as the positive control
    cand = [k for k, _ in out["frag_ring"][:6]]
    if "vPFK" not in cand:
        cand.append("vPFK")
    obj_dyn = make_obj(True)
    obj_sta = make_obj(False)
    out["restore_baseline"] = {"dyn": obj_dyn({}), "static": obj_sta({})}
    out["restore_vmax_only"] = rc.greedy_restoration_set(cand, obj_dyn, max_size=3,
                                                         starts=(0.0, -0.5, 0.5))
    out["restore_with_true_fix"] = rc.greedy_restoration_set(cand + ["CiPFKATP"], obj_dyn,
                                                             max_size=3, bounds=(-1, 6),
                                                             starts=(0.0, 2.0, 4.5))
    out["restore_static_objective"] = rc.greedy_restoration_set(cand + ["CiPFKATP"], obj_sta,
                                                                max_size=1, bounds=(-1, 6),
                                                                starts=(0.0, 4.5))
    # recovery curves: healthy, defect, best Vmax-only fix, true fix
    def curve(adj):
        m.set_params(snap_dis)
        for k, v in adj.items():
            if k == "CiPFKATP":
                m.set_params({"CiPFKATP": snap_dis["CiPFKATP"] * np.exp(v)})
            else:
                m.scale_enzyme(k, float(np.exp(v)))
        xa = m.steady_state(x0=x0)
        snap = m.snapshot()
        m.scale_enzyme("vATP", 0.8)
        xb = m.steady_state(x0=xa)
        m.set_params(snap)
        m.set_x(xb)
        r = m.rr.simulate(0, 1.0, 401, selections=["time", "F6P", "ATP"])
        m.set_params(base)
        return r[:, 0], r[:, 1] / xa[m.state_ids.index("F6P")], xa
    best_v = out["restore_vmax_only"][-1][1]
    best_t = out["restore_with_true_fix"][-1][1]
    out["curves"] = {"healthy": None, "defect": curve({}),
                     "best Vmax-only fix": curve(best_v), "true fix (feedback)": curve(best_t)}
    m.set_params(base)
    xa = m.steady_state(x0=x0)
    snap = m.snapshot()
    m.scale_enzyme("vATP", 0.8)
    xb = m.steady_state(x0=xa)
    m.set_params(snap)
    m.set_x(xb)
    r = m.rr.simulate(0, 1.0, 401, selections=["time", "F6P", "ATP"])
    out["curves"]["healthy"] = (r[:, 0], r[:, 1] / xa[m.state_ids.index("F6P")], xa)
    m.set_params(base)
    return out


# ======================================================================= C
NH = os.path.join(DATA_CACHE, "nhanes")
NH_CYCLES = {"D": 2005, "E": 2007, "F": 2009, "G": 2011, "H": 2013, "I": 2015}


def load_nhanes():
    """NHANES 2005-2016 OGTT subsample (fasting glucose, 2-h OGTT glucose,
    HbA1c, age, sex, BMI) + public-use linked mortality
    through 2019 (NCHS). Files cached in figures/_data_cache/nhanes/."""
    import subprocess
    import pandas as pd
    os.makedirs(NH, exist_ok=True)
    rows = []
    for c, y in NH_CYCLES.items():
        for f in ("OGTT", "GLU", "GHB", "DEMO", "BMX"):
            fn = os.path.join(NH, f"{f}_{c}.xpt")
            if not os.path.exists(fn):
                subprocess.run(["curl", "-sS", "-o", fn,
                                f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{y}/DataFiles/{f}_{c}.xpt"],
                               check=True)
        mf = os.path.join(NH, f"MORT_{c}.dat")
        if not os.path.exists(mf):
            subprocess.run(["curl", "-sS", "-o", mf,
                            "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/"
                            f"linked_mortality/NHANES_{y}_{y + 1}_MORT_2019_PUBLIC.dat"], check=True)
        r = lambda f: pd.read_sas(os.path.join(NH, f"{f}_{c}.xpt"))
        m = (r("OGTT")[["SEQN", "LBXGLT"]]
             .merge(r("GLU")[["SEQN", "LBXGLU"]], on="SEQN", how="left")
             .merge(r("GHB")[["SEQN", "LBXGH"]], on="SEQN", how="left")
             .merge(r("DEMO")[["SEQN", "RIDAGEYR", "RIAGENDR"]], on="SEQN")
             .merge(r("BMX")[["SEQN", "BMXBMI"]], on="SEQN", how="left"))
        mort = pd.read_fwf(mf, colspecs=[(0, 6), (14, 15), (15, 16), (16, 19), (19, 20),
                                         (20, 21), (42, 45), (45, 48)],
                           names=["SEQN", "eligstat", "mortstat", "ucod", "dm_mcod",
                                  "htn_mcod", "permth_int", "permth_exm"], na_values=["."])
        m = m.merge(mort, on="SEQN", how="left")
        m["cycle"] = c
        rows.append(m)
    A = pd.concat(rows, ignore_index=True)
    A = A[(A.RIDAGEYR >= 20) & (A.eligstat == 1)]
    return A.dropna(subset=["LBXGLT", "LBXGLU", "LBXGH", "BMXBMI", "mortstat", "permth_exm"])


def part_C():
    import pandas as pd
    out = {}
    A = load_nhanes()
    t = A.permth_exm.values.astype(float)
    ev = A.mortstat.values.astype(int)
    cov = A[["RIDAGEYR", "RIAGENDR", "BMXBMI"]].values.astype(float)
    fg, a1c, g2 = A.LBXGLU.values, A.LBXGH.values, A.LBXGLT.values
    models = {"M0 clinical": cov,
              "M0 static (+fasting glucose, HbA1c)": np.c_[cov, fg, a1c],
              "M1 +2-h OGTT glucose": np.c_[cov, fg, a1c, g2],
              "M1b +OGTT increment shape": np.c_[cov, fg, a1c, g2, (g2 - fg) / fg]}
    out["nhanes_n"] = len(A)
    out["nhanes_deaths"] = int(ev.sum())
    out["nhanes_follow"] = float(np.median(t))
    out["nhanes_all"] = rc.incremental_survival(models, t, ev, n_boot=300, horizon=96)
    # cause-specific known positives: deaths with diabetes as a multiple cause,
    # and heart-disease deaths (UCOD_LEADING 001)
    ev_dm = ((A.mortstat == 1) & (A.dm_mcod == 1)).values.astype(int)
    ev_cv = ((A.mortstat == 1) & (A.ucod == 1)).values.astype(int)
    out["nhanes_dm"] = rc.incremental_survival(models, t, ev_dm, n_boot=300, horizon=96)
    out["nhanes_cv"] = rc.incremental_survival(models, t, ev_cv, n_boot=300, horizon=96)
    out["nhanes_dm_n"], out["nhanes_cv_n"] = int(ev_dm.sum()), int(ev_cv.sum())
    # DECODE-style comparison (no HbA1c): fasting glucose vs + 2-h OGTT
    dec = {"M0 clinical + fasting glucose": np.c_[cov, fg],
           "M1 +2-h OGTT glucose": np.c_[cov, fg, g2]}
    out["decode_all"] = rc.incremental_survival(dec, t, ev, n_boot=300, horizon=96)
    out["decode_dm"] = rc.incremental_survival(dec, t, ev_dm, n_boot=300, horizon=96)
    for k in ("decode_all", "decode_dm"):
        for v in out[k].values():
            v.pop("lp", None)
    for k in ("nhanes_all", "nhanes_dm", "nhanes_cv"):
        for v in out[k].values():
            v.pop("lp", None)
    # CGMacros: meal challenges with known diagnoses
    from netjac.demos.d2_perturbation_relaxation import cgm_metrics
    R = cgm_metrics()
    R = R.dropna(subset=["k", "peak", "fasting_glucose", "Age", "BMI"]).copy()
    R["t_half"] = R["t_half"].fillna(R["t_half"].median())
    y = (R.group.astype(str) == "T2D").astype(int).values
    y2 = (R.group.astype(str) != "healthy").astype(int).values
    base = R[["Age", "BMI"]].values
    mods = {"M0 static (fasting glucose)": np.c_[base, R.fasting_glucose],
            "M1 +meal peak": np.c_[base, R.fasting_glucose, R.peak],
            "M3 +recovery kinetics": np.c_[base, R.fasting_glucose, R.peak, R.k, R.t_half]}
    out["cgm_n"] = {"total": len(R), "T2D": int(y.sum()), "dysglycaemia": int(y2.sum())}
    out["cgm_t2d"] = rc.incremental_classification(mods, y, n_splits=5, n_boot=2000, reps=20)
    out["cgm_dys"] = rc.incremental_classification(mods, y2, n_splits=5, n_boot=2000, reps=20)
    for k in ("cgm_t2d", "cgm_dys"):
        for v in out[k].values():
            v.pop("pred", None)
    return out


# ======================================================================= B
POS_FAM = {
    "proteasome": r"^PSM[ABCD]\d+$",
    "tRNA synthetase": r"^([A-Z]ARS|EPRS|FARS[AB])$",
    "sterol synthesis": r"^(HMGCS1|HMGCR|MVK|PMVK|MVD|FDPS|FDFT1|SQLE|LSS|CYP51A1|"
                        r"IDI1|DHCR7|SC5D|MSMO1|NSDHL|HSD17B7|EBP|TM7SF2|DHCR24|SCAP|SREBF2)$",
    "ER chaperone / UPR": r"^(HSPA5|HSP90B1|PDIA\d|DNAJB11|CALR|CANX|HYOU1|SEL1L|"
                          r"DERL\d|SYVN1|DDOST|SEC61A1|SEC63)$",
    "actin (MRTF-SRF)": r"^(ACTB|ACTG1|CFL1|PFN1|CAPZB|CAPZA1|DSTN)$",
}
# co-regulated machines: large stoichiometric complexes without a known
# transcriptional feedback on their own subunits (negatives)
NEG_FAM = {
    "cytosolic ribosome": r"^RP[SL]\d+[A-Z]?\d*$",
    "mito ribosome": r"^MRP[SL]\d+$",
    "spliceosome core": r"^(SF3[AB]\d|SNRP[ABDEFG]\d*|SNRNP\d+|PRPF\d+|LSM[2-8])$",
    "RNA polymerase II": r"^POLR2[A-L]$",
    "CCT chaperonin": r"^(CCT\d|CCT6A|TCP1)$",
    "translation initiation": r"^(EIF2B\d|EIF2S\d|EIF3[A-M])$",
    "RNA exosome": r"^EXOSC\d+$",
    "mediator": r"^MED\d+$",
    "NADH dehydrogenase": r"^NDUF[ABSV]\d+$",
    "MCM / replication": r"^(MCM[2-7]|GINS\d|POLA\d|POLE\d?)$",
}


def _family(g, fams):
    for k, pat in fams.items():
        if re.match(pat, g):
            return k
    return None


def part_B(K=50):
    import pandas as pd
    from netjac.data import perturbseq as ps
    out = {}
    # ---------------- B1: homeostat vs co-regulated, static vs response
    r = ps.reduce_singlecell()
    H, pert, gid = r["half_means"], r["pert"], r["gene_id"]
    pg = np.array([p.split("_")[1] for p in pert])
    pid = np.array([p.split("_")[-1] for p in pert])
    ix = {e: i for i, e in enumerate(gid)}
    ctrl = pg == "non-targeting"
    meas = np.array([i for i in range(len(pert)) if pid[i] in ix and not ctrl[i]])
    own = np.array([ix[pid[i]] for i in meas])
    g = pg[meas]
    keep = np.array([_family(x, POS_FAM) is not None or _family(x, NEG_FAM) is not None
                     for x in g])
    meas, own, g = meas[keep], own[keep], g[keep]
    C = r["ctrl_cells"].astype(np.float64)
    mu_raw = C.mean(0)
    C = (C - C.mean(0)) / C.std(0)
    n = len(meas)
    Ccol = (C.T @ C[:, own]) / len(C)
    # all perturbed genes' partner sets for the specificity rank
    allm = np.array([i for i in range(len(pert)) if pid[i] in ix and not ctrl[i]])
    allo = np.array([ix[pid[i]] for i in allm])
    Call = (C.T @ C[:, allo]) / len(C)
    del C
    Pall = np.zeros((len(allm), len(gid)), np.float32)
    for k in range(len(allm)):
        c = Call[:, k].copy()
        c[allo[k]] = -np.inf
        Pall[k, np.argsort(-c)[:K]] = 1.0 / K
    feats_s, feats_d = [], []
    D = H[meas].mean(1)
    pos_in_all = {m_: i for i, m_ in enumerate(allm)}
    for k in range(n):
        c = Ccol[:, k].copy()
        c[own[k]] = -np.inf
        top = np.argsort(-c)[:K]
        d = D[k].copy()
        d[own[k]] = 0
        # static (unperturbed ensemble only)
        feats_s.append([np.log(mu_raw[own[k]] + 1e-3) if mu_raw[own[k]] > -1 else 0.0,
                        float(np.mean(Ccol[top, k])), float(np.max(Ccol[top, k])),
                        float(np.sum(Ccol[:, k] > 0.05))])
        # response (perturbation) features
        M = Pall @ d                                   # partner-set responses
        q = float(np.mean(M < M[pos_in_all[meas[k]]]))
        feats_d.append([float(np.mean(d[top])), q, float(np.linalg.norm(d)),
                        float(D[k][own[k]])])
    Fs, Fd = np.array(feats_s), np.array(feats_d)
    y = np.array([_family(x, POS_FAM) is not None for x in g]).astype(int)
    fam = np.array([_family(x, POS_FAM) or _family(x, NEG_FAM) for x in g])
    models = {"M0 static co-fluctuation": Fs,
              "M1 +knockdown response": np.c_[Fs, Fd]}
    out["B1_n"] = {"pos": int(y.sum()), "neg": int((1 - y).sum()),
                   "families": pd.Series(fam).value_counts().to_dict()}
    out["B1_random_cv"] = rc.incremental_classification(models, y, n_splits=5, n_boot=2000, reps=10)
    # leave-family-out: train on other families, predict an unseen family
    ug = np.unique(fam)
    gnum = np.array([list(ug).index(f) for f in fam])
    out["B1_leave_family_out"] = rc.incremental_classification(models, y, groups=gnum,
                                                               n_splits=len(ug), n_boot=2000)
    out["B1_response_only"] = rc.incremental_classification(
        {"response only": Fd}, y, groups=gnum, n_splits=len(ug), n_boot=10)
    for kk in ("B1_random_cv", "B1_leave_family_out", "B1_response_only"):
        for v in out[kk].values():
            v.pop("pred", None)
    out["B1_genes"] = list(zip(g, fam, y, Fd[:, 1]))
    # ---------------- B2: core regulators and driver nodes (genome-wide network)
    from netjac.demos.d10_signs_genomewide_thresholds import _bulk_z
    Bg = ps.load_bulk("K562_gwps_bulk")
    X = Bg["X"].copy()
    X[~np.isfinite(X)] = 0
    Zs = _bulk_z(Bg)
    gi = {e: i for i, e in enumerate(Bg["gene_id"])}
    ctrl2 = Bg["control"]
    m2 = np.array([i for i in range(len(X)) if Bg["pert_id"][i] in gi and not ctrl2[i]])
    o2 = np.array([gi[Bg["pert_id"][i]] for i in m2])
    g2 = Bg["pert_gene"][m2]
    selfz = X[m2, o2]
    best = {}
    for k in range(len(m2)):
        if g2[k] not in best or selfz[k] < selfz[best[g2[k]]]:
            best[g2[k]] = k
    ks = np.array(sorted(best.values()))
    sel = ks[selfz[ks] < -0.5]
    Z = Zs[m2[sel]][:, o2[sel]].T
    gS = g2[sel]
    sig = (np.abs(Z) > 6.0)
    np.fill_diagonal(sig, False)
    drivers = rc.driver_nodes(sig)
    outdeg, indeg = sig.sum(0), sig.sum(1)
    out["B2"] = {"nodes": len(gS), "edges": int(sig.sum()), "n_drivers": len(drivers),
                 "frac_drivers": len(drivers) / len(gS),
                 "drivers_isolated": int(np.sum((outdeg[drivers] == 0) & (indeg[drivers] == 0))),
                 "drivers_with_in": int(np.sum(indeg[drivers] > 0)),
                 "top_out": [(gS[i], int(outdeg[i])) for i in np.argsort(-outdeg)[:15]]}
    # driver set on the connected part only (nodes with any edge)
    conn = (outdeg + indeg) > 0
    sub = sig[np.ix_(conn, conn)]
    dsub = rc.driver_nodes(sub)
    out["B2"]["connected_nodes"] = int(conn.sum())
    out["B2"]["connected_drivers"] = len(dsub)
    # random-graph null with the same in/out degree sequence (configuration)
    rng = np.random.default_rng(0)
    nd_null = []
    tgt, src = np.nonzero(sub)
    for _ in range(5):
        perm = rng.permutation(len(src))
        A = np.zeros_like(sub)
        A[tgt, src[perm]] = True
        np.fill_diagonal(A, False)
        nd_null.append(len(rc.driver_nodes(A)))
    out["B2"]["driver_null_degree_preserving"] = nd_null
    # hub regulators among drivers? (Liu 2011: drivers avoid hubs)
    sub_out = sub.sum(0)
    out["B2"]["driver_outdeg_median"] = float(np.median(sub_out[dsub]))
    out["B2"]["all_outdeg_median"] = float(np.median(sub_out))
    cg = gS[conn]
    out["B2"]["driver_examples_high_out"] = [(cg[i], int(sub_out[i]))
                                             for i in sorted(dsub, key=lambda i: -sub_out[i])[:12]]
    # ---------------- B3: restoration toward control from double perturbations
    out["B3"] = _norman_restoration()
    return out


def _norman_restoration(thr=0.75, seed=0):
    """Norman et al. 2019 (Science 365:786) CRISPRa singles + doubles in
    K562 (GEARS processed release, Harvard Dataverse doi:10.7910/DVN/Q2ZV3E).
    Observed restoration ratio rho = |dAB| / max(|dA|, |dB|) (distance to
    the control centroid); positives = doubles that end closer to control
    than their stronger single (rho < thr). Predictors from singles only:
    additive superposition (linear Jacobian response) vs static
    co-expression of the two target genes in control cells."""
    import anndata as ad
    from scipy.stats import spearmanr
    fn = os.path.join(DATA_CACHE, "norman_gears", "norman", "perturb_processed.h5ad")
    a = ad.read_h5ad(fn)
    cond = a.obs["condition"].astype(str).values
    X = a.X.tocsr() if hasattr(a.X, "tocsr") else a.X
    gnames = np.array(a.var["gene_name"].astype(str))
    rng = np.random.default_rng(seed)
    means, halves = {}, {}
    for c in np.unique(cond):
        idx = np.where(cond == c)[0]
        Xi = X[idx]
        means[c] = np.asarray(Xi.mean(0)).ravel()
        h = rng.permutation(len(idx))
        halves[c] = (np.asarray(X[idx[h[: len(h) // 2]]].mean(0)).ravel(),
                     np.asarray(X[idx[h[len(h) // 2:]]].mean(0)).ravel())
    ctrl = means["ctrl"]
    ctrl_cells = X[np.where(cond == "ctrl")[0]]
    d = {c: v - ctrl for c, v in means.items()}

    def single(gn):
        for c in (f"{gn}+ctrl", f"ctrl+{gn}"):
            if c in d:
                return c
        return None
    gidx = {gname: i for i, gname in enumerate(gnames)}
    Cc = ctrl_cells.toarray() if hasattr(ctrl_cells, "toarray") else np.asarray(ctrl_cells)
    rows = []
    for c in d:
        if c == "ctrl" or "ctrl" in c:
            continue
        A_, B_ = c.split("+")
        sa, sb = single(A_), single(B_)
        if sa is None or sb is None:
            continue
        dA, dB, dAB = d[sa], d[sb], d[c]
        nA, nB, nAB = map(np.linalg.norm, (dA, dB, dAB))
        rho = nAB / max(nA, nB)
        rho_add = np.linalg.norm(dA + dB) / max(nA, nB)
        cos = float(dA @ dB / (nA * nB))
        # split-half noise floor of the double's distance
        h1, h2 = halves[c]
        noise = np.linalg.norm(h1 - h2) / 2
        # static: co-expression of the two target genes among control cells
        if A_ in gidx and B_ in gidx:
            x1, x2 = Cc[:, gidx[A_]], Cc[:, gidx[B_]]
            st = float(np.corrcoef(x1, x2)[0, 1]) if x1.std() > 0 and x2.std() > 0 else 0.0
        else:
            st = 0.0
        rows.append({"double": c, "rho": rho, "rho_add": rho_add, "cos": cos,
                     "static_coexpr": st, "nA": nA, "nB": nB, "nAB": nAB,
                     "noise": noise})
    import pandas as pd
    T = pd.DataFrame(rows)
    T["restores"] = T.rho < thr
    y = T.restores.values.astype(int)
    res = {"n_doubles": len(T), "n_restore": int(y.sum()),
           "spearman_add": spearmanr(T.rho_add, T.rho),
           "spearman_static": spearmanr(T.static_coexpr, T.rho)}
    if 0 < y.sum() < len(y):
        res["auc_add"] = rc.auc(y, -T.rho_add.values)
        res["auc_cos"] = rc.auc(y, -T.cos.values)
        res["auc_static"] = rc.auc(y, -T.static_coexpr.values)
        b_add, b_st = [], []
        for _ in range(2000):
            b = rng.integers(0, len(y), len(y))
            if y[b].min() == y[b].max():
                continue
            b_add.append(rc.auc(y[b], -T.rho_add.values[b]))
            b_st.append(rc.auc(y[b], -T.static_coexpr.values[b]))
        res["auc_add_ci"] = np.percentile(b_add, [2.5, 97.5])
        res["auc_static_ci"] = np.percentile(b_st, [2.5, 97.5])
        res["d_auc"] = np.percentile(np.array(b_add) - np.array(b_st), [2.5, 50, 97.5])
    res["table"] = T.sort_values("rho")
    return res


# ================================================================== report
def _ci(t):
    return f"{t[0]:+.3f} [{t[1]:+.3f}, {t[2]:+.3f}]"


def report(R):
    if "A" in R:
        A = R["A"]
        print(f"\nA  Ground truth (Teusink glycolysis): {A['n']}")
        for k, v in A["inc"].items():
            s = f"   {k:28s} AUC {v['auc']:.3f}  FP@80%sens {v['fpr']:.1%}"
            if "d_prev" in v:
                s += f"  dAUC {_ci(v['d_prev'])} p={v['d_prev'][3]:.3f}  dFPR {_ci(v['d_fpr'])}"
            d = A["decoy"][k]
            s += f"  | decoys flagged {d['fp_decoy']:.0%}, healthy flagged {d['fp_healthy']:.0%}"
            print(s)
        ring = [l for l in A["eig_h"] if abs(l.imag) > 1e-6]
        print(f"   healthy has a ringing pair {ring[0]:.1f}; defect eigenvalues all real: "
              f"{all(abs(l.imag) < 1e-6 for l in A['eig_d'])}")
        print("   core (total fragility of the ringing pair): " + ", ".join(
            f"{k} {v:+.1f}" for k, v in A["frag_ring"][:6]))
        print("   core (slowest mode): " + ", ".join(f"{k} {v:+.1f}" for k, v in A["frag_slow"][:6]))
        print(f"   edge fragility top: " + ", ".join(f"{a}<-{b}" for a, b, _ in A["edge_top"][:6]) +
              f"; PFK edges rank {A['pfk_edge_ranks']} of {A['n_edges']}")
        b = A["restore_baseline"]
        print(f"   restoration objective before fix: dynamic {b['dyn']:.4f}, static {b['static']:.1e}")
        for nm in ("restore_vmax_only", "restore_with_true_fix", "restore_static_objective"):
            for st in A[nm][1:]:
                print(f"   {nm:26s} {st[0]} -> objective {st[2]:.5f} "
                      f"({1 - st[2] / A[nm][0][2]:.1%} reduction)" if A[nm][0][2] > 0 else
                      f"   {nm:26s} {st[0]} -> objective {st[2]:.1e} (static already matched)")
    if "B" in R:
        B = R["B"]
        print(f"\nB  Perturb-seq: {B['B1_n']['pos']} homeostat-family vs {B['B1_n']['neg']} "
              f"co-regulated-family knockdowns")
        for kk in ("B1_random_cv", "B1_leave_family_out", "B1_response_only"):
            for k, v in B[kk].items():
                s = f"   {kk:22s} {k:26s} AUC {v['auc']:.3f} FP@80% {v['fpr']:.1%}"
                if "d_prev" in v:
                    s += f"  dAUC {_ci(v['d_prev'])}  dFPR {_ci(v['d_fpr'])}"
                print(s)
        b2 = B["B2"]
        print(f"   network {b2['nodes']} nodes / {b2['edges']} edges; minimum driver nodes "
              f"{b2['n_drivers']} ({b2['frac_drivers']:.0%}); connected part {b2['connected_drivers']}/"
              f"{b2['connected_nodes']} vs degree-preserving random {b2['driver_null_degree_preserving']}")
        print("   top regulators (out-degree): " + ", ".join(f"{g}({d})" for g, d in b2["top_out"][:10]))
        b3 = B["B3"]
        print(f"   Norman 2019: {b3['n_doubles']} doubles, {b3['n_restore']} end closer to control "
              f"(rho<0.75); additive single responses: Spearman {b3['spearman_add'][0]:.2f}, "
              f"AUC {b3['auc_add']:.3f} {np.round(b3['auc_add_ci'], 3)}; static co-expression AUC "
              f"{b3['auc_static']:.3f} {np.round(b3['auc_static_ci'], 3)}; dAUC {np.round(b3['d_auc'], 3)}")
        T = b3["table"]
        print("   restoring doubles: " + "; ".join(
            f"{r.double} rho {r.rho:.2f} (additive {r.rho_add:.2f})" for r in T.head(8).itertuples()))
    if "C" in R:
        C = R["C"]
        print(f"\nC  NHANES 2005-2016 OGTT: {C['nhanes_n']} adults, {C['nhanes_deaths']} deaths "
              f"(median follow-up {C['nhanes_follow']:.0f} months); diabetes-associated {C['nhanes_dm_n']}, "
              f"heart disease {C['nhanes_cv_n']}")
        for kk in ("nhanes_all", "nhanes_dm", "nhanes_cv", "decode_all", "decode_dm"):
            for k, v in C[kk].items():
                s = f"   {kk:11s} {k:38s} C {v['c']:.4f}  8-y AUC {v['auc_h']:.4f} FP@80% {v['fp_h'][1]:.1%}"
                if "d_prev" in v:
                    s += f"  dC {_ci(v['d_prev'])} p={v['d_prev'][3]:.2f}  dFPR {_ci(v['d_fpr_h'])}"
                print(s)
        print(f"   CGMacros {C['cgm_n']}")
        for kk in ("cgm_t2d", "cgm_dys"):
            for k, v in C[kk].items():
                s = f"   {kk:8s} {k:30s} AUC {v['auc']:.3f} FP@80% {v['fpr']:.1%}"
                if "d_prev" in v:
                    s += f"  dAUC {_ci(v['d_prev'])}"
                print(s)


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 4, figsize=(21, 9.4))
    col = ["0.55", "C0", "C2", "C3"]
    if "A" in R:
        A = R["A"]
        a = ax[0, 0]
        ks = list(A["inc"])
        a.bar(range(4), [A["inc"][k]["auc"] for k in ks], color=col)
        for i, k in enumerate(ks):
            a.text(i, A["inc"][k]["auc"] + 0.01, f"FP {A['inc'][k]['fpr']:.0%}", ha="center", fontsize=8)
        a.set_ylim(0.4, 1.06)
        a.set_xticks(range(4), ["M0\nstatic", "M1\n+challenge", "M2\n+dose", "M3\n+recovery"])
        a.set_ylabel("AUC (defect vs healthy + decoys)")
        a.set_title("A  ground truth: hidden feedback defect\nstatic blind, one challenge sees it")
        a = ax[0, 1]
        for k, c in (("healthy", "k"), ("defect", "C3"), ("best Vmax-only fix", "C1"),
                     ("true fix (feedback)", "C2")):
            t, f6p, _ = A["curves"][k]
            a.plot(t * 1000, f6p, color=c, lw=2 if k != "healthy" else 3,
                   ls="--" if k == "true fix (feedback)" else "-", label=k)
        a.set_xlim(0, 400)
        a.set_xlabel("ms after the challenge is removed")
        a.set_ylabel("F6P / own baseline")
        a.set_title("A  restoring the recovery, not the set point")
        a.legend(fontsize=8)
        a = ax[0, 2]
        path_v = A["restore_vmax_only"]
        path_t = A["restore_with_true_fix"]
        a.plot(range(len(path_v)), [p[2] for p in path_v], "o-", color="C1", label="enzyme capacities only")
        a.plot(range(len(path_t)), [p[2] for p in path_t], "o-", color="C2", label="+ feedback strength allowed")
        for i, p in enumerate(path_t[1:], 1):
            a.annotate(p[0][-1], (i, p[2]), fontsize=8, xytext=(4, 4), textcoords="offset points")
        a.set_yscale("log")
        a.set_xlabel("number of knobs (greedy)")
        a.set_ylabel("distance to healthy recovery signature")
        a.set_title("A  minimal control set: the edge itself\n(static objective is already 0)")
        a.legend(fontsize=8)
    if "B" in R:
        B = R["B"]
        a = ax[0, 3]
        labs, vals, fps = [], [], []
        for kk, tag in (("B1_random_cv", "mixed folds"), ("B1_leave_family_out", "unseen family")):
            for k, v in B[kk].items():
                labs.append(f"{'static' if 'M0' in k else '+response'}\n{tag}")
                vals.append(v["auc"])
                fps.append(v["fpr"])
        a.bar(range(4), vals, color=["0.55", "C3", "0.55", "C3"])
        for i, (v, f) in enumerate(zip(vals, fps)):
            a.text(i, v + 0.01, f"FP {f:.0%}", ha="center", fontsize=8)
        a.axhline(0.5, color="k", ls=":", lw=1)
        a.set_xticks(range(4), labs, fontsize=8)
        a.set_ylabel("AUC: homeostat vs co-regulated family")
        a.set_title("B  cells: response beats co-fluctuation,\nand transfers to unseen components")
        a = ax[1, 0]
        T = B["B3"]["table"]
        a.scatter(T.rho_add, T.rho, c=np.where(T.restores, "C2", "0.6"), s=22)
        lim = [0.3, max(1.6, float(T.rho.max()) + 0.05)]
        a.plot(lim, lim, "k:", lw=1)
        a.axhline(0.75, color="C2", ls="--", lw=1)
        for j, r in enumerate(T.head(3).itertuples()):
            a.annotate(r.double, (r.rho_add, r.rho), fontsize=7,
                       xytext=(40, -14 - 11 * j), textcoords="offset points",
                       arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
        a.set_xlabel("predicted from singles (additive response)")
        a.set_ylabel("observed |double| / max |single|")
        a.set_title(f"B  which second perturbation restores?\nAUC {B['B3']['auc_add']:.2f} "
                    f"vs static co-expression {B['B3']['auc_static']:.2f}")
        a = ax[1, 1]
        b2 = B["B2"]
        a.bar([0, 1], [b2["connected_drivers"] / b2["connected_nodes"],
                       np.mean(b2["driver_null_degree_preserving"]) / b2["connected_nodes"]],
              color=["C3", "0.6"])
        a.set_xticks([0, 1], ["K562 response\nnetwork", "degree-preserving\nrandom"])
        a.set_ylabel("fraction of nodes that must be driven")
        a.set_title(f"B  full structural control is expensive\n({b2['connected_drivers']} of "
                    f"{b2['connected_nodes']} nodes)")
    if "C" in R:
        C = R["C"]
        a = ax[1, 2]
        ks = list(C["nhanes_dm"])
        x = np.arange(len(ks))
        for j, (kk, c, lab) in enumerate((("nhanes_all", "0.5", "all-cause"),
                                          ("nhanes_dm", "C3", "diabetes-associated"),
                                          ("nhanes_cv", "C0", "heart disease"))):
            a.plot(x, [C[kk][k]["c"] for k in ks], "o-", color=c, label=lab)
        a.set_xticks(x, ["clinical", "+fasting\n+HbA1c", "+2-h\nOGTT", "+OGTT\nshape"], fontsize=8)
        a.set_ylabel("out-of-fold C-index")
        a.set_title(f"C  NHANES OGTT → mortality (n={C['nhanes_n']})\nthe stimulated endpoint adds ~nothing")
        a.legend(fontsize=8)
        a = ax[1, 3]
        for kk, c, lab in (("cgm_t2d", "C3", "T2D"), ("cgm_dys", "C1", "dysglycaemia")):
            ks = list(C[kk])
            a.plot(range(len(ks)), [C[kk][k]["auc"] for k in ks], "o-", color=c, label=lab)
        a.set_xticks(range(3), ["fasting\nglucose", "+meal\npeak", "+recovery\nkinetics"])
        a.set_ylabel("cross-validated AUC")
        a.set_ylim(0.8, 1.0)
        a.set_title(f"C  CGMacros meals (n={C['cgm_n']['total']}):\nstatic already saturates")
        a.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(FIG, "d11_recovery_state.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return [p]


def main(parts="A,B,C"):
    R = _load_res()
    for p in parts.split(","):
        t0 = time.time()
        R[p] = {"A": part_A, "B": part_B, "C": part_C}[p]()
        print(f"[part {p} done in {time.time() - t0:.0f}s]", flush=True)
        _save_res(R)
    report(R)
    print(figures(R))
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="A,B,C")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.report:
        R = _load_res()
        report(R)
        print(figures(R))
    else:
        main(a.part)
