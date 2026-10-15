"""
D12  Real time-series validation of perturbation-recovery state definition:
     Myd88-deficient macrophages (GEO GSE75210; Tong et al., Cell Systems
     2016, PMID 27211859).

Design: bone-marrow-derived macrophages, wild type vs Myd88-/-, 12 innate
stimuli x 8 time points (0-4 h) x 2 replicate 96-well plates per genotype
(384 low-input RNA-seq libraries). Each (plate, stimulus) is one time
series sampled from neighbouring wells.

Known positives (textbook, independent of these data):
  MyD88-dependent inputs   Pam3Cys (TLR2), flagellin (TLR5), R848 (TLR7),
                           CpG (TLR9)
  MyD88-independent        Poly(I:C) (TLR3-TRIF), Sendai virus (RIG-I),
                           HT-DNA (cGAS-STING), MDP (NOD2), TDB (Mincle)
  mixed (not labelled)     LPS (TLR4: MyD88 + TRIF), MDP+LPS, zymosan
                           (Dectin-1 + TLR2)
  and the IRF3 / interferon arm downstream of TRIF stays intact in Myd88-/-.

Questions (same as D11):
  1  Which inputs route through MyD88? (the cell's functional state) -
     static vs single endpoints vs trajectory, transfer from replicate 1 to 2.
  2  Is this macrophage Myd88-deficient? (genotype state) - nested
     M0 static -> M3 trajectory, leave-replicate-plate-out.
  3  Recovery kinetics of immediate-early genes (peak and return).
  4  Core components: genes whose induction is lost in both replicates.
  5  Restoring the state: in Myd88-/- cells, which available input brings
     the response back to the wild-type MyD88-dependent response, and at
     what off-target cost (interferon arm)?

Modules are fixed from prior knowledge (no data-driven selection):
NF-kB primary-response genes and IRF3 / interferon-stimulated genes.

Run:  python3 -m netjac.demos.d12_timecourse_cell_state
"""
import os
import pickle
import re
import subprocess

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import recovery as rc

DIR = os.path.join(DATA_CACHE, "geo_gse75210")
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75210/suppl/GSE75210_myd88_exp.csv.gz"
RES = os.path.join(DATA_CACHE, "d12_results.pkl")

NFKB = ["Tnf", "Il1b", "Il1a", "Cxcl1", "Cxcl2", "Nfkbia", "Nfkbiz", "Nfkbie", "Tnfaip3",
        "Ier3", "Egr1", "Egr2", "Egr3", "Fos", "Junb", "Ccl3", "Ccl4", "Ptgs2", "Il6", "Il10",
        "Il12b", "Ccl20", "Zfp36", "Dusp1", "Dusp2", "Ier2", "Nr4a1", "Irf1", "Tnfaip2",
        "Icam1", "Sod2", "Bcl2a1a", "Socs3", "Ccl5", "Cxcl3", "Pim1", "Map3k8", "Gadd45b",
        "Rel", "Nfkb1", "Nfkb2", "Relb"]
IRF = ["Ifit1", "Ifit2", "Ifit3", "Isg15", "Rsad2", "Cxcl10", "Ifnb1", "Oasl1", "Mx1", "Irf7",
       "Ifi44", "Usp18", "Rtp4", "Gbp2", "Stat1", "Ccl12", "Mx2", "Oasl2", "Ifi205", "Zbp1"]
IEG = ["Egr1", "Egr2", "Egr3", "Fos", "Ier2", "Ier3", "Nr4a1", "Dusp1", "Zfp36", "Junb"]
DEP = ["Pam3Cys", "FLA", "R848", "CpG"]
INDEP = ["Poly.I.C.", "SeV", "HT.DNA", "MDP", "TDB"]
MIXED = ["LPS", "MDP.LPS", "Zymosan"]


def load():
    import pandas as pd
    fn = os.path.join(DIR, "myd88_exp.csv.gz")
    if not os.path.exists(fn):
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sS", "-o", fn, URL], check=True)
    X = pd.read_csv(fn, index_col=0)
    meta = pd.DataFrame([re.match(r"(WT|KO)(\d)\.(.*)\.T([\d.]+)$", c).groups()
                         for c in X.columns], columns=["g", "rep", "stim", "t"], index=X.columns)
    meta["t"] = meta.t.astype(float)
    lib = X.sum(axis=0)
    L = np.log2(X / lib * 1e6 + 1)
    expressed = (X >= 5).mean(axis=1) > 0.2
    L = L[expressed]
    T = sorted(meta.t.unique())
    traj = {}
    for key, idx in meta.groupby(["g", "rep", "stim"]).groups.items():
        m = meta.loc[idx].sort_values("t")
        traj[key] = L[m.index].T.values             # time x gene
    genes = np.array(L.index)
    return {"traj": traj, "genes": genes, "T": np.array(T), "meta": meta,
            "stims": list(meta.stim.unique()), "lib": lib}


def _mod(D, names):
    gi = {g: i for i, g in enumerate(D["genes"])}
    return np.array([gi[g] for g in names if g in gi])


def _induction(D, key, idx):
    """Module trajectory relative to the series' own t=0 well."""
    Z = D["traj"][key][:, idx].mean(1)
    return Z - Z[0]


# ------------------------------------------------------------ question 1
def q1_routing(D):
    """Per stimulus: lost induction in Myd88-/- = WT induction - KO
    induction of the NF-kB module. Readouts: static t=0 difference, each
    single time point, the trajectory area. Threshold learned on replicate
    1, applied to replicate 2."""
    nf = _mod(D, NFKB)
    T = D["T"]
    readouts = ["static t=0"] + [f"{t:g} h" for t in T[1:]] + ["trajectory area"]
    lab = {s: 1 for s in DEP} | {s: 0 for s in INDEP}

    def score(s, r, ro):
        w = _induction(D, ("WT", r, s), nf)
        k = _induction(D, ("KO", r, s), nf)
        if ro == "static t=0":
            W0 = D["traj"][("WT", r, s)][0, nf].mean()
            K0 = D["traj"][("KO", r, s)][0, nf].mean()
            return W0 - K0
        if ro == "trajectory area":
            return float(np.trapezoid(w - k, T) / T[-1])
        i = list(T).index(float(ro.split()[0]))
        return float(w[i] - k[i])
    out = {}
    for ro in readouts:
        sc = {r: {s: score(s, r, ro) for s in D["stims"]} for r in "12"}
        y = np.array([lab[s] for s in DEP + INDEP])
        a1 = rc.auc(y, [sc["1"][s] for s in DEP + INDEP])
        a2 = rc.auc(y, [sc["2"][s] for s in DEP + INDEP])
        # threshold: midpoint between classes on replicate 1 -> replicate 2
        thr = (min(sc["1"][s] for s in DEP) + max(sc["1"][s] for s in INDEP)) / 2
        mis = [s for s in DEP + INDEP if (sc["2"][s] > thr) != bool(lab[s])]
        rep_r = np.corrcoef([sc["1"][s] for s in D["stims"]], [sc["2"][s] for s in D["stims"]])[0, 1]
        out[ro] = {"auc_rep1": a1, "auc_rep2": a2, "threshold": thr, "misclassified_rep2": mis,
                   "rep_agreement": float(rep_r), "scores": sc,
                   "mixed": {s: (sc["1"][s] + sc["2"][s]) / 2 for s in MIXED}}
    return out


# ------------------------------------------------------------ question 2
def q2_genotype(D, n_top=300):
    """Genotype classification of the 48 (plate, stimulus) series, folds =
    replicate plate pair. Features: (a) the two prior modules, (b) the
    most time-variable genes of the replicate used for training only."""
    nf, ir = _mod(D, NFKB), _mod(D, IRF)
    units = [(g, r, s) for g in ("WT", "KO") for r in "12" for s in D["stims"]]
    y = np.array([u[0] == "KO" for u in units]).astype(int)
    grp = np.array([int(u[1]) for u in units])
    T = D["T"]
    it = {t: i for i, t in enumerate(T)}
    P = np.array([D["traj"][u] for u in units])               # unit x time x gene
    Dl = P - P[:, :1, :]
    mods = np.stack([P[:, :, nf].mean(2), P[:, :, ir].mean(2)], 2)      # unit x time x 2
    dmods = mods - mods[:, :1, :]
    res = {}
    feats_m = {"M0 static t=0": mods[:, 0, :],
               "M1 +1 h response": np.c_[mods[:, 0, :], dmods[:, it[1.0], :]],
               "M2 +2.5 h response": np.c_[mods[:, 0, :], dmods[:, it[1.0], :], dmods[:, it[2.5], :]],
               "M3 +full trajectory": np.c_[mods[:, 0, :], dmods[:, 1:, :].reshape(len(units), -1)]}
    res["modules"] = rc.incremental_classification(feats_m, y, groups=grp, n_splits=2,
                                                   n_boot=1000, C=1.0)
    # transcriptome: top time-variable genes chosen on WT series only
    var = np.array([D["traj"][u].var(0) for u in units if u[0] == "WT"]).mean(0)
    top = np.argsort(-var)[:n_top]
    feats_g = {"M0 static t=0": P[:, 0, top],
               "M1 +1 h response": np.c_[P[:, 0, top], Dl[:, it[1.0], top]],
               "M3 +full trajectory": np.c_[P[:, 0, top], Dl[:, 1:, top].reshape(len(units), -1)]}
    res["transcriptome"] = rc.incremental_classification(feats_g, y, groups=grp, n_splits=2,
                                                         n_boot=1000, C=0.1)
    # response-only, restricted to series with a MyD88-dependent stimulus
    dep = np.array([u[2] in DEP for u in units])
    feats_r = {"t=0 static": mods[dep, 0, :],
               "response trajectory only": dmods[dep, 1:, :].reshape(dep.sum(), -1)}
    res["dependent_only"] = rc.incremental_classification(feats_r, y[dep], groups=grp[dep],
                                                          n_splits=2, n_boot=1000, C=1.0)
    for k in res:
        for v in res[k].values():
            v.pop("pred", None)
    # what differs at baseline?
    gi = D["genes"]
    W0 = np.array([D["traj"][u][0] for u in units if u[0] == "WT"])
    K0 = np.array([D["traj"][u][0] for u in units if u[0] == "KO"])
    diff = W0.mean(0) - K0.mean(0)
    res["baseline_top"] = [(gi[i], float(diff[i])) for i in np.argsort(-np.abs(diff))[:15]]
    res["baseline_module"] = {"NF-kB": float(diff[nf].mean()), "IRF": float(diff[ir].mean())}
    return res


# ------------------------------------------------------------ question 3
def q3_recovery(D):
    """Immediate-early genes: time of peak and fraction recovered by the
    last time point, WT vs KO, for inputs where both genotypes respond."""
    ieg = _mod(D, IEG)
    T = D["T"]
    out = {}
    for s in ["LPS", "MDP.LPS", "Zymosan", "Poly.I.C.", "R848", "Pam3Cys"]:
        row = {}
        for g in ("WT", "KO"):
            curves = [_induction(D, (g, r, s), ieg) for r in "12"]
            c = np.mean(curves, 0)
            ip = int(np.argmax(c))
            peak = c[ip]
            row[g] = {"curve": c, "peak": float(peak), "t_peak": float(T[ip]),
                      "recovered": float(1 - c[-1] / peak) if peak > 0.3 else np.nan}
        out[s] = row
    return out


# ------------------------------------------------------------ question 4
def q4_core(D, thr=1.0):
    from netjac.data import gprofiler as gp
    T = D["T"]
    loss = {}
    for r in "12":
        L = []
        for s in DEP:
            w = D["traj"][("WT", r, s)]
            k = D["traj"][("KO", r, s)]
            L.append(np.trapezoid((w - w[0]) - (k - k[0]), T, axis=0) / T[-1])
        loss[r] = np.mean(L, 0)
    core = (loss["1"] > thr) & (loss["2"] > thr)
    genes = D["genes"]
    # genes induced by LPS in both genotypes (TRIF route kept)
    kept = np.ones(len(genes), bool)
    for g in ("WT", "KO"):
        for r in "12":
            z = D["traj"][(g, r, "LPS")]
            kept &= (z - z[0]).max(0) > thr
    out = {"n_core": int(core.sum()), "core": sorted(genes[core],
                                                     key=lambda x: -loss["1"][list(genes).index(x)])[:40],
           "rep_agreement": float(np.corrcoef(loss["1"], loss["2"])[0, 1]),
           "n_kept_lps": int(kept.sum()), "kept_lps": list(genes[kept])[:40]}
    bg = list(genes)
    out["enrich_core"] = [x for x in gp.enrich(list(genes[core]), bg, organism="mmusculus")
                          if x["term_size"] < 800][:8]
    out["enrich_kept"] = [x for x in gp.enrich(list(genes[kept]), bg, organism="mmusculus")
                          if x["term_size"] < 800][:8]
    out["nfkb_in_core"] = [g for g in NFKB if g in set(genes[core])]
    out["irf_in_kept"] = [g for g in IRF if g in set(genes[kept])]
    return out


# ------------------------------------------------------------ question 5
def q5_restore(D):
    """Target: the wild-type NF-kB-module response to MyD88-dependent inputs
    (mean of the four; genes x time, relative to t=0). For each input given
    to Myd88-/- cells:
      gain    projection on the target / |target|^2   (1 = full amplitude)
      cos     cosine with the target                  (shape)
      dist    1 - |KO_x - target| / |KO_dep - target| (Euclidean restoration;
              kept to show that it is noise-dominated, see `noise`)
    Positive control: each WT dependent input against the mean of the other
    three (what "fully restored" looks like at this noise level). Negative
    control: the knockout's own response to the dependent inputs.
    Off-target cost: IRF-module induction area of KO_x minus that of the
    WT target (log2 x h; > 0 = interferon programme the target does not
    have)."""
    nf, ir = _mod(D, NFKB), _mod(D, IRF)
    T = np.asarray(D["T"], float)
    out, noise = {}, []
    for r in "12":
        def resp(g, s, idx=nf):
            z = D["traj"][(g, r, s)][:, idx]
            return (z - z[0]).ravel()
        def irf_area(g, s):
            return float(np.trapezoid(_induction(D, (g, r, s), ir), T))
        tgt = np.mean([resp("WT", s) for s in DEP], 0)
        base = np.mean([np.linalg.norm(resp("KO", s) - tgt) for s in DEP])
        irf_t = np.mean([irf_area("WT", s) for s in DEP])
        noise.append(np.mean([np.linalg.norm(resp("KO", s)) for s in DEP]) / np.linalg.norm(tgt))

        def score(x, t, irf):
            nt = np.linalg.norm(t)
            return {"gain": x @ t / nt ** 2, "cos": x @ t / (np.linalg.norm(x) * nt + 1e-12),
                    "dist": 1 - np.linalg.norm(x - t) / base, "irf_excess": irf - irf_t}
        rows = {}
        for s in DEP:                               # positive control (leave-one-out WT)
            t = np.mean([resp("WT", q) for q in DEP if q != s], 0)
            rows[f"+ctrl WT {s}"] = score(resp("WT", s), t, irf_area("WT", s))
            rows[f"-ctrl KO {s}"] = score(resp("KO", s), tgt, irf_area("KO", s))
        for s in D["stims"]:
            if s not in DEP:
                rows[s] = score(resp("KO", s), tgt, irf_area("KO", s))
        rows["(no input) KO at rest"] = score(np.zeros_like(tgt), tgt, 0.0)
        for k, v in rows.items():
            for m, x in v.items():
                out.setdefault(k, {}).setdefault(m, []).append(float(x))
    res = {k: {m: (float(np.mean(x)), float(np.min(x)), float(np.max(x))) for m, x in v.items()}
           for k, v in out.items()}
    return {"rows": res, "noise": float(np.mean(noise))}


# ================================================================= report
def run():
    D = load()
    R = {"n_series": len(D["traj"]), "n_genes": len(D["genes"]),
         "q1": q1_routing(D), "q2": q2_genotype(D), "q3": q3_recovery(D),
         "q4": q4_core(D), "q5": q5_restore(D)}
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    print(f"\nD12 GSE75210 Myd88 macrophages: {R['n_series']} time series, {R['n_genes']} expressed genes")
    print("Q1 which inputs route through MyD88 (4 dependent vs 5 independent; threshold from replicate 1):")
    for ro, v in R["q1"].items():
        print(f"   {ro:16s} AUC rep1 {v['auc_rep1']:.2f} rep2 {v['auc_rep2']:.2f}  replicate r {v['rep_agreement']:+.2f}"
              f"  misclassified in rep2: {v['misclassified_rep2'] or 'none'}  | mixed inputs "
              + ", ".join(f"{s} {x:+.2f}" for s, x in v["mixed"].items()))
    print("Q2 is this macrophage Myd88-/- (48 series, leave-replicate-plate-out):")
    for k in ("modules", "transcriptome", "dependent_only"):
        for m, v in R["q2"][k].items():
            s = f"   {k:15s} {m:26s} AUC {v['auc']:.3f} FP@80% {v['fpr']:.0%}"
            if "d_prev" in v:
                s += f"  dAUC {v['d_prev'][0]:+.3f} [{v['d_prev'][1]:+.3f}, {v['d_prev'][2]:+.3f}]"
            print(s)
    print(f"   baseline WT-KO module difference: {R['q2']['baseline_module']}; top baseline genes: " +
          ", ".join(f"{g} {d:+.1f}" for g, d in R["q2"]["baseline_top"][:10]))
    print("Q3 immediate-early genes (peak time h, fraction recovered by 4 h):")
    for s, v in R["q3"].items():
        print(f"   {s:10s} WT peak {v['WT']['peak']:+.2f} at {v['WT']['t_peak']:g} h, recovered {v['WT']['recovered']:.0%} | "
              f"KO peak {v['KO']['peak']:+.2f} at {v['KO']['t_peak']:g} h, recovered {v['KO']['recovered']:.0%}")
    q = R["q4"]
    print(f"Q4 core: {q['n_core']} genes lose >1 log2 induction-area in BOTH replicates (replicate r {q['rep_agreement']:.2f});"
          f" NF-kB prior genes among them: {q['nfkb_in_core']}")
    print("   core enrichment: " + "; ".join(f"{x['name'][:45]} ({x['source']}, p={x['p']:.0e}, k={x['k']})"
                                       for x in q["enrich_core"][:5]))
    print(f"   LPS-induced in both genotypes: {q['n_kept_lps']} genes; IRF prior genes among them: {q['irf_in_kept']}")
    print("   kept enrichment: " + "; ".join(f"{x['name'][:45]} ({x['source']}, p={x['p']:.0e}, k={x['k']})"
                                       for x in q["enrich_kept"][:5]))
    q5 = R["q5"]
    print("Q5 restoring the MyD88-dependent NF-kB response in Myd88-/- cells "
          f"(noise floor: |KO response to dependent inputs| = {q5['noise']:.2f} x |target|)")
    print("   gain = amplitude along the WT target (1 = full), cos = shape, dist = Euclidean restoration, "
          "IRF excess = interferon area beyond the target (log2*h); ranges over replicates")
    for s, v in sorted(q5["rows"].items(), key=lambda kv: -kv[1]["gain"][0]):
        print(f"   {s:22s} gain {v['gain'][0]:+.2f} ({v['gain'][1]:+.2f}..{v['gain'][2]:+.2f})  cos {v['cos'][0]:+.2f}"
              f"  dist {v['dist'][0]:+.2f}  IRF excess {v['irf_excess'][0]:+5.1f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    q1 = R["q1"]
    a = ax[0, 0]
    tr = q1["trajectory area"]["scores"]
    order = DEP + MIXED + INDEP
    for j, r in enumerate("12"):
        a.barh(np.arange(len(order)) + (j - 0.5) * 0.38, [tr[r][s] for s in order], 0.38,
               color=["C3" if s in DEP else ("C1" if s in MIXED else "C0") for s in order],
               alpha=1 if r == "1" else 0.55, label=f"replicate {r}")
    a.axvline(q1["trajectory area"]["threshold"], color="k", ls="--", lw=1)
    a.set_yticks(range(len(order)), order)
    a.invert_yaxis()
    a.set_xlabel("NF-κB induction lost in Myd88-/- (trajectory area, log2)")
    a.set_title("Q1  which inputs route through MyD88\n(red = textbook dependent, blue = independent)")
    a.legend(fontsize=8)
    a = ax[0, 1]
    ros = list(q1)
    a.plot(range(len(ros)), [q1[r]["auc_rep1"] for r in ros], "o-", label="replicate 1")
    a.plot(range(len(ros)), [q1[r]["auc_rep2"] for r in ros], "s--", label="replicate 2")
    a.set_xticks(range(len(ros)), [r.replace("trajectory area", "traj.").replace("static t=0", "static")
                                   for r in ros], rotation=45, fontsize=8)
    a.set_ylabel("AUC dependent vs independent")
    a.set_ylim(0, 1.05)
    a.set_title("Q1  static t=0 carries no routing information;\nany stimulated readout separates")
    a.legend(fontsize=8)
    a = ax[0, 2]
    q2 = R["q2"]
    for k, c in (("modules", "C3"), ("transcriptome", "C0")):
        ks = list(q2[k])
        a.plot(range(len(ks)), [q2[k][m]["auc"] for m in ks], "o-", color=c, label=f"{k}")
    a.set_xticks(range(4), ["M0 static", "M1 +1 h", "M2 +2.5 h", "M3 traj."], fontsize=8)
    a.set_ylim(0.4, 1.02)
    a.set_ylabel("AUC Myd88-/- vs WT (leave-plate-pair-out)")
    a.set_title("Q2  genotype: baseline already differs\n(tonic signalling or plate batch)")
    a.legend(fontsize=8)
    a = ax[1, 0]
    T = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 4])
    for s, c in (("LPS", "C3"), ("Poly.I.C.", "C0"), ("R848", "C2")):
        for g, ls in (("WT", "-"), ("KO", "--")):
            a.plot(T, R["q3"][s][g]["curve"], ls, color=c, lw=2, label=f"{s} {g}")
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel("hours after stimulus")
    a.set_ylabel("immediate-early genes (log2 vs own t=0)")
    a.set_title("Q3  peak and return of immediate-early genes")
    a.legend(fontsize=7, ncol=2)
    a = ax[1, 1]
    rows = R["q5"]["rows"]
    pos = [k for k in rows if k.startswith("+ctrl")]
    neg = [k for k in rows if k.startswith("-ctrl")]
    byp = sorted([k for k in rows if k not in pos + neg], key=lambda k: -rows[k]["gain"][0])
    for grp, c in ((pos, "C2"), (byp, "C1"), (neg, "C3")):
        g = np.array([rows[k]["gain"] for k in grp])
        cs = [rows[k]["cos"][0] for k in grp]
        a.errorbar(g[:, 0], cs, xerr=[g[:, 0] - g[:, 1], g[:, 2] - g[:, 0]], fmt="o", color=c, ms=6, lw=1)
        for k, x, y in zip(grp, g[:, 0], cs):
            if grp is byp or k.endswith("R848"):
                a.annotate(f"{k.replace('+ctrl ', '').replace('-ctrl ', '')} (IFN {rows[k]['irf_excess'][0]:+.0f})"
                           if grp is byp else k, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    a.plot([], [], "o", color="C2", label="+ control: WT dependent input (LOO)")
    a.plot([], [], "o", color="C1", label="bypass input in Myd88-/-  (IFN excess)")
    a.plot([], [], "o", color="C3", label="- control: Myd88-/- dependent input")
    a.axvline(1, color="k", lw=0.8, ls=":")
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel("NF-κB gain along WT MyD88-dependent response (1 = full)")
    a.set_ylabel("shape (cosine)")
    a.set_title("Q5  bypass inputs restore NF-κB amplitude and shape\nbut add an interferon programme")
    a.legend(fontsize=7, loc="lower right")
    a = ax[1, 2]
    ec = R["q4"]["enrich_core"][:6]
    ek = R["q4"]["enrich_kept"][:6]
    labs = [f"lost: {x['name'][:38]}" for x in ec] + [f"kept: {x['name'][:38]}" for x in ek]
    vals = [-np.log10(x["p"]) for x in ec + ek]
    a.barh(range(len(labs)), vals, color=["C3"] * len(ec) + ["C0"] * len(ek))
    a.set_yticks(range(len(labs)), labs, fontsize=7)
    a.invert_yaxis()
    a.set_xlabel("-log10 p (g:Profiler, expressed-gene background)")
    a.set_title(f"Q4  core: {R['q4']['n_core']} genes lost in both replicates")
    fig.tight_layout()
    p = os.path.join(FIG, "d12_timecourse_cell_state.png")
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    R = pickle.load(open(RES, "rb")) if a.report and os.path.exists(RES) else run()
    report(R)
    print(figures(R))
