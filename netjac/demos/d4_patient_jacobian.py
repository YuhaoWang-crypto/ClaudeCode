"""
D4 — Patient-specific Jacobians from TCGA (clear-cell renal carcinoma).

Chassis: Holzhütter 2004 human glycolysis + pentose-phosphate + glutathione
kinetic model (BIOMD0000000070; erythrocyte enzymes, exact curated
kinetics). It is a Warburg-type chassis (glucose -> lactate, no OXPHOS).
Data: TCGA-KIRC STAR-Counts TPM (GDC open access), 541 primary tumours + 72
solid-tissue normals; clinical OS from GDC.

Instantiation (⚠️ assumptions, stated up front):
  * enzyme level of reaction r = sum of TPM of its isoform genes (GENE_MAP);
  * Vmax_r(sample) = Vmax_r(model) * level_r(sample) / median level_r(normals)
    i.e. the median kidney normal IS the curated model; tumours are
    expression-scaled perturbations of it;
  * mRNA ~ enzyme activity; isoform kinetics are the erythrocyte ones
    (HK1, PFK, PKLR ...) even when the tumour switches isoform (HK2, PKM2);
  * ATP demand (ATPase), oxidative load, glucose, Mg are held fixed.
Every per-sample number below is a ✅ deterministic computation on that
instantiated model; what it means for the patient is ⚠️.

Per sample: steady state existence, glycolytic flux, lactate export, PPP
share, leading eigenvalue (resilience; τ = 1/|Re λ|), full flux-control
matrix, and TOTAL fragility of every expression-mapped enzyme.

Analyses
  A  tumour vs paired normal: which network properties shift (Wilcoxon)
  B  patient-specific fragile nodes: top fragile enzyme per sample
  C  pre-specified survival test (OS, Cox): λ_lead, glycolytic flux, PPP
     share vs the boring baselines (glycolysis / PPP gene-expression scores,
     stage, age), incl. "adds beyond baseline" likelihood-ratio tests.
"""
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from netjac import FIG

MODEL = "BIOMD0000000070"
GENE_MAP = {
    "vGLT": ["SLC2A1", "SLC2A3"], "vHEX": ["HK1", "HK2", "HK3"],
    "vGPI": ["GPI"], "vPFK": ["PFKL", "PFKM", "PFKP"],
    "vALD": ["ALDOA", "ALDOB", "ALDOC"], "vTPI": ["TPI1"],
    "vGAPDH": ["GAPDH"], "vPGK": ["PGK1"], "vBPGM": ["BPGM"],
    "vBPGP": ["BPGM"], "vPGM": ["PGAM1", "PGAM2"],
    "vENO": ["ENO1", "ENO2", "ENO3"], "vPK": ["PKLR", "PKM"],
    "vLDHNADH": ["LDHA", "LDHB"], "vLDHNADPH": ["LDHA", "LDHB"],
    "vAK": ["AK1", "AK2"], "vG6PDH": ["G6PD"], "vPGLDH": ["PGD"],
    "vGSSGRD": ["GSR"], "vRibPepi": ["RPE"], "vRibPiso": ["RPIA"],
    "vTrKet1": ["TKT"], "vTrKet2": ["TKT"], "vTrAld": ["TALDO1"],
    "vPPRPPS": ["PRPS1", "PRPS2"], "vLacexch": ["SLC16A1", "SLC16A3"],
    "vPyrexch": ["SLC16A1", "SLC16A3"],
}
GLYC_GENES = ["SLC2A1", "HK1", "HK2", "GPI", "PFKL", "PFKM", "PFKP", "ALDOA",
              "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA"]
PPP_GENES = ["G6PD", "PGD", "RPE", "RPIA", "TKT", "TALDO1"]
FRAG_ENZ = ["vGLT", "vHEX", "vGPI", "vPFK", "vALD", "vTPI", "vGAPDH", "vPGK",
            "vBPGM", "vBPGP", "vPGM", "vENO", "vPK", "vLDHNADH", "vG6PDH",
            "vPGLDH", "vGSSGRD", "vTrKet1", "vTrAld", "vPPRPPS", "vLacexch"]
CLIP = (0.02, 50.0)
# species groups used to label eigenmodes by pathway: the slowest mode of the
# erythrocyte chassis is the 2,3-BPG/Mg pool (RBC-specific, ~absent in kidney),
# so "the" lead eigenvalue is chassis physiology, not tumour glycolysis
MODE_GROUPS = {
    "shunt": ["Gri23P2f", "MgGri23P2", "Mgf", "MgAMP", "MgATP", "MgADP"],
    "glycolysis": ["Glcin", "Glc6P", "Fru6P", "Fru16P2", "GraP", "DHAP",
                   "Gri13P2", "Gri3P", "Gri2P", "PEP", "Pyr", "Lac"],
    "ppp": ["GlcA6P", "Rul5P", "Xul5P", "Rib5P", "Sed7P", "E4P", "NADPHf",
            "NADPf", "P1NADP", "P1NADPH", "P2NADP", "P2NADPH", "GSH", "GSSG"],
}


def pathway_modes(m, x, Jz):
    """Slowest eigenvalue whose (relative) mode shape puts >50 % weight on each
    species group. Returns {group: eigenvalue or nan}."""
    vals, vr = np.linalg.eig(Jz)
    out = {g: np.nan for g in MODE_GROUPS}
    idx = {g: [m.state_ids.index(s) for s in ss if s in m.state_ids]
           for g, ss in MODE_GROUPS.items()}
    for i in np.argsort(-vals.real):
        u = np.abs(m.B @ vr[:, i]) / np.maximum(np.abs(x), 1e-12)
        w = u ** 2 / np.sum(u ** 2)
        for g, ii in idx.items():
            if np.isnan(out[g]) and w[ii].sum() > 0.5:
                out[g] = vals[i]
    return out

_W = {}


def _worker_model():
    if "m" not in _W:
        from netjac.core import KineticModel
        m = KineticModel(MODEL, "Holzhutter2004")
        _W["x0"] = m.steady_state()
        _W["base"] = m.snapshot()
        m.enzyme_params()
        _W["m"] = m
    return _W["m"], _W["x0"], _W["base"]


def analyse_sample(scales, with_fragility=True):
    """scales: {reaction: factor}. Returns dict of network properties, or
    {'ok': False} if the instantiated network has no steady state."""
    from netjac import spectrum as sp
    m, x0, base = _worker_model()
    m.set_params(base)
    for r, f in scales.items():
        m.scale_enzyme(r, f)
    try:
        x = m.steady_state(x0=x0, max_rounds=4)
    except RuntimeError:
        m.set_params(base)
        return {"ok": False}
    c = m.control(x)
    vals = np.linalg.eigvals(c["Jz"])
    lam = vals.real.max()
    modes = pathway_modes(m, x, c["Jz"])
    v = c["v"]
    ix = m.rxn_ids.index
    out = {"ok": True, "lam": lam,
           "lam_shunt": np.real(modes["shunt"]),
           "lam_glyc": np.real(modes["glycolysis"]),
           "lam_ppp": np.real(modes["ppp"]),
           "fcc_ppp": c["CJ_scaled"][ix("vG6PDH")].copy(),
           "glyc": v[ix("vHEX")],
           "lactate": -v[ix("vLacexch")],
           "ppp_share": v[ix("vG6PDH")] / v[ix("vHEX")],
           "atp": m.rr["MgATP"] + m.rr["ATPf"],
           "fcc": c["CJ_scaled"][ix("vHEX")].copy(),
           "fcc_theorem": c["summation_residual"]}
    if with_fragility:
        fr = {}
        for r in FRAG_ENZ:
            try:
                # fragility of the GLYCOLYTIC mode (not the RBC shunt mode)
                fr[r] = sp.total_fragility(m, r, x, target=modes["glycolysis"],
                                           solver="integrate")
            except RuntimeError:
                fr[r] = np.nan
        out["frag"] = fr
    m.set_params(base)
    return out


def enzyme_scales(expr):
    """Per-sample Vmax factors relative to the median normal."""
    normals = expr[expr.sample_type == "Solid Tissue Normal"]
    S = {}
    for r, genes in GENE_MAP.items():
        lvl = expr[genes].sum(axis=1)
        ref = normals[genes].sum(axis=1).median()
        S[r] = np.clip(lvl / ref, *CLIP)
    import pandas as pd
    return pd.DataFrame(S, index=expr.index)


def run_cohort(expr, workers=8):
    S = enzyme_scales(expr)
    tasks = [{r: float(v) for r, v in S.iloc[i].items()}
             for i in range(len(S))]
    res = []
    with ProcessPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(analyse_sample, tasks, chunksize=4)):
            res.append(r)
            if (i + 1) % 100 == 0:
                print(f"    ... {i + 1}/{len(tasks)} samples", flush=True)
    return S, res


def _cindex(time, event, risk):
    """Harrell's C (higher risk -> shorter survival)."""
    t, e, r = map(np.asarray, (time, event, risk))
    num = den = 0.0
    for i in np.where(e == 1)[0]:
        later = t > t[i]
        den += later.sum()
        num += (r[i] > r[later]).sum() + 0.5 * (r[i] == r[later]).sum()
    return num / den


def _cox(df, cols):
    from statsmodels.duration.hazard_regression import PHReg
    X = df[cols].values.astype(float)
    X = (X - X.mean(0)) / X.std(0)
    res = PHReg(df["time"].values, X, status=df["event"].values).fit()
    return res


def survival(tab):
    """tab: tumour rows with features + time/event/stage/age."""
    from scipy.stats import chi2
    d = tab.dropna(subset=["time", "event"]).copy()
    d = d[d.time > 0]
    stage_map = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4}
    d["stage_n"] = d["stage"].map(stage_map)
    d["age"] = d["age_days"] / 365.25
    d = d.dropna(subset=["stage_n", "age", "lam_glyc", "lam_ppp"])
    feats = {"lam_lead [pre]": "lam", "glycolytic flux [pre]": "glyc",
             "PPP share [pre]": "ppp_share",
             "lam_glyc [post-hoc]": "lam_glyc", "lam_ppp [post-hoc]": "lam_ppp",
             "glycolysis gene score": "glyc_score",
             "PPP gene score": "ppp_score"}
    out = {}
    for lab, col in feats.items():
        uni = _cox(d, [col])
        adj = _cox(d, [col, "stage_n", "age"])
        out[lab] = {"hr": float(np.exp(uni.params[0])),
                    "p": float(uni.pvalues[0]),
                    "hr_adj": float(np.exp(adj.params[0])),
                    "p_adj": float(adj.pvalues[0]),
                    "c": _cindex(d.time, d.event, d[col])}
    # does the network feature add beyond the expression score + stage + age?
    base_cols = ["glyc_score", "ppp_score", "stage_n", "age"]
    base = _cox(d, base_cols)
    add = {}
    for lab, col in [("lam_lead [pre]", "lam"), ("glycolytic flux [pre]", "glyc"),
                     ("PPP share [pre]", "ppp_share"),
                     ("lam_glyc [post-hoc]", "lam_glyc"),
                     ("lam_ppp [post-hoc]", "lam_ppp")]:
        full = _cox(d, base_cols + [col])
        lr = 2 * (full.llf - base.llf)
        add[lab] = {"LR": lr, "p": float(chi2.sf(lr, 1))}
    return {"n": len(d), "events": int(d.event.sum()), "uni": out,
            "adds": add, "c_stage": _cindex(d.time, d.event, d.stage_n)}


def report(make_fig=True, project="TCGA-KIRC", workers=4):
    import pandas as pd
    from scipy.stats import wilcoxon, mannwhitneyu
    from netjac.data import tcga
    from netjac import discriminate as di

    genes = sorted({g for v in GENE_MAP.values() for g in v} |
                   set(GLYC_GENES) | set(PPP_GENES))
    expr = tcga.expression(project, genes)
    expr = expr[expr.sample_type.isin(["Primary Tumor",
                                       "Solid Tissue Normal"])]
    n_dup = int(expr.index.duplicated().sum())
    expr = expr[~expr.index.duplicated()]      # one aliquot per sample
    S, res = run_cohort(expr, workers)
    from netjac.core import KineticModel
    m0 = KineticModel(MODEL)
    rxn = m0.rxn_ids

    ok = np.array([r["ok"] for r in res])
    tab = pd.DataFrame(index=expr.index)
    tab["sample_type"] = expr.sample_type
    tab["patient"] = expr.patient
    tab["ok"] = ok
    for k in ("lam", "lam_shunt", "lam_glyc", "lam_ppp", "glyc", "lactate",
              "ppp_share", "atp", "fcc_theorem"):
        tab[k] = [r.get(k, np.nan) for r in res]
    fr = pd.DataFrame([r.get("frag", {}) for r in res], index=expr.index)
    fcc = pd.DataFrame([r["fcc"] if r["ok"] else np.full(len(rxn), np.nan)
                        for r in res], index=expr.index, columns=rxn)
    z = lambda g: (np.log2(expr[g] + 1) - np.log2(expr[g] + 1).mean()) \
        / np.log2(expr[g] + 1).std()
    tab["glyc_score"] = pd.concat([z(g) for g in GLYC_GENES], axis=1).mean(axis=1)
    tab["ppp_score"] = pd.concat([z(g) for g in PPP_GENES], axis=1).mean(axis=1)
    T = tab[tab.sample_type == "Primary Tumor"]
    N = tab[tab.sample_type == "Solid Tissue Normal"]

    # ---------- A: paired tumour vs normal
    pairs = sorted(set(T.patient) & set(N.patient))
    tp = T.drop_duplicates("patient").set_index("patient").loc[pairs]
    npr = N.drop_duplicates("patient").set_index("patient").loc[pairs]
    okp = tp.ok.values & npr.ok.values
    paired = {}
    for k in ("lam", "lam_glyc", "lam_ppp", "glyc", "lactate", "ppp_share",
              "atp"):
        a, b = tp[k].values[okp], npr[k].values[okp]
        paired[k] = {"tumour": np.median(a), "normal": np.median(b),
                     "p": wilcoxon(a, b).pvalue,
                     "auc": di.auc(T[k][T.ok], N[k][N.ok], n_boot=300)}
    # expression baseline: best single gene AUC tumour vs normal
    gene_auc = {g: di.auc(expr.loc[T.index, g], expr.loc[N.index, g],
                          n_boot=1)["auc"] for g in genes}
    best_gene = max(gene_auc, key=gene_auc.get)
    # flux-control redistribution (median scaled C^J on glycolytic flux)
    fccT, fccN = fcc.loc[T.index[T.ok]].median(), fcc.loc[N.index[N.ok]].median()
    # PPP mechanism: is the model's PPP shift visible in G6PD expression alone?
    from scipy.stats import spearmanr
    g6pd_auc = di.auc(expr.loc[T.index, "G6PD"], expr.loc[N.index, "G6PD"])
    fppT = pd.DataFrame([r["fcc_ppp"] for r in res if r["ok"]],
                        index=[i for i, r in zip(expr.index, res) if r["ok"]],
                        columns=rxn)
    fppT_t = fppT.loc[[i for i in fppT.index if i in T.index]].median()
    fppT_n = fppT.loc[[i for i in fppT.index if i in N.index]].median()
    To = T[T.ok]
    drivers = {r: spearmanr(S.loc[To.index, r], To.ppp_share).statistic
               for r in GENE_MAP}
    # ---------- B: fragile nodes
    frT, frN = fr.loc[T.index[T.ok]], fr.loc[N.index[N.ok]]
    top_T = frT.abs().idxmax(1).value_counts(normalize=True)
    top_N = frN.abs().idxmax(1).value_counts(normalize=True)
    dfrag = {}
    trow = [T.index[T.patient == p][0] for p in pairs]
    nrow = [N.index[N.patient == p][0] for p in pairs]
    for r in FRAG_ENZ:
        a = fr.loc[trow, r].values[okp]
        b = fr.loc[nrow, r].values[okp]
        good = np.isfinite(a) & np.isfinite(b)
        dfrag[r] = {"tumour": np.median(a[good]), "normal": np.median(b[good]),
                    "p": wilcoxon(a[good], b[good]).pvalue}
    # ---------- C: survival
    cl = tcga.clinical(project)
    Tt = T[T.ok].drop_duplicates("patient").merge(cl, on="patient")
    surv = survival(Tt)

    print("=" * 68)
    print(f"D4  PATIENT-SPECIFIC JACOBIANS — {project} on Holzhütter 2004 "
          f"(BIOMD0000000070)")
    print("=" * 68)
    print(f"samples: {len(T)} tumour, {len(N)} normal ({n_dup} duplicate "
          f"aliquots dropped); steady state exists: "
          f"tumour {T.ok.mean()*100:.1f} %, normal {N.ok.mean()*100:.1f} %; "
          f"MCA summation residual max {np.nanmax(tab.fcc_theorem):.1e}  ✅")
    print(f"[A] {okp.sum()} matched tumour/normal pairs (paired Wilcoxon):")
    for k, d in paired.items():
        print(f"    {k:10s} tumour {d['tumour']:.4g} vs normal {d['normal']:.4g}"
              f"  p={d['p']:.1g}  unpaired AUC {d['auc']['auc']:.2f}")
    print(f"    baseline: best single gene ({best_gene}) tumour-vs-normal AUC "
          f"{gene_auc[best_gene]:.2f} (expression alone separates them)")
    top = (fccT - fccN).abs().sort_values(ascending=False).index[:5]
    print("    flux control on glycolysis, normal -> tumour: " + ", ".join(
        f"{r} {fccN[r]:+.3f}->{fccT[r]:+.3f}" for r in top))
    print(f"    PPP: model PPP share x{paired['ppp_share']['tumour']/paired['ppp_share']['normal']:.1f}"
          f" in tumours, but G6PD expression alone AUC {g6pd_auc['auc']:.2f} "
          f"(median TPM ratio {expr.loc[T.index,'G6PD'].median()/expr.loc[N.index,'G6PD'].median():.2f})")
    top_pc = (fppT_t.abs()).sort_values(ascending=False).index[:5]
    print("    flux control of PPP (vG6PDH) in tumours: " + ", ".join(
        f"{r} {fppT_t[r]:+.2f} (normal {fppT_n[r]:+.2f})" for r in top_pc))
    dr = sorted(drivers, key=lambda r: -abs(drivers[r]))[:5]
    print("    tumour PPP-share drivers (Spearman with enzyme scale): " +
          ", ".join(f"{r} {drivers[r]:+.2f}" for r in dr))
    print("    pathway modes (median, tumour vs normal): " + ", ".join(
        f"{k} {T[k].median():.3f} vs {N[k].median():.3f}"
        for k in ("lam_shunt", "lam_glyc", "lam_ppp")))
    print("[B] most fragile enzyme for the GLYCOLYTIC mode (share of samples):")
    print("    tumour: " + ", ".join(f"{r} {p*100:.0f}%" for r, p in
                                      top_T.head(4).items()))
    print("    normal: " + ", ".join(f"{r} {p*100:.0f}%" for r, p in
                                      top_N.head(4).items()))
    big = sorted(dfrag, key=lambda r: dfrag[r]["p"])[:5]
    print("    largest paired fragility shifts: " + "; ".join(
        f"{r} {dfrag[r]['normal']:+.3f}->{dfrag[r]['tumour']:+.3f} "
        f"(p={dfrag[r]['p']:.0e})" for r in big))
    print(f"[C] OS, n={surv['n']} tumours with stage+age, {surv['events']} "
          f"deaths (Cox, per SD):  C-index of stage alone {surv['c_stage']:.2f}")
    for lab, d in surv["uni"].items():
        print(f"    {lab:24s} HR {d['hr']:.2f} (p={d['p']:.1g}); adj. stage+age "
              f"HR {d['hr_adj']:.2f} (p={d['p_adj']:.1g}); C {d['c']:.2f}")
    for lab, d in surv["adds"].items():
        print(f"    adds beyond expression scores+stage+age: {lab:22s} "
              f"LR {d['LR']:.2f}, p={d['p']:.2g}")
    if make_fig:
        _figure(tab, T, N, fccT, fccN, frT, frN, surv, Tt)
    return {"tab": tab, "paired": paired, "gene_auc": gene_auc,
            "g6pd_auc": g6pd_auc, "fcc_ppp_t": fppT_t, "fcc_ppp_n": fppT_n,
            "drivers": drivers,
            "fccT": fccT, "fccN": fccN, "top_T": top_T, "top_N": top_N,
            "dfrag": dfrag, "surv": surv, "S": S}


def _figure(tab, T, N, fccT, fccN, frT, frN, surv, Tt):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    for k, (lab, d) in enumerate([("normal", N), ("tumour", T)]):
        v = d.lam_glyc[d.ok].dropna()
        a.scatter(np.full(len(v), k) + np.random.default_rng(k).uniform(
            -.15, .15, len(v)), -1 / v, s=6, alpha=.5)
        a.plot([k - .25, k + .25], [np.median(-1 / v)] * 2, color="k")
    a.set_xticks([0, 1]); a.set_xticklabels(["normal (72)", "tumour (541)"])
    a.set_yscale("log"); a.set_ylabel("recovery time τ = 1/|Re λ| (model h)")
    a.set_title("A  glycolytic-mode recovery time per patient")
    a = ax[0, 1]
    idx = (fccT - fccN).abs().sort_values(ascending=False).index[:8]
    y = np.arange(len(idx))
    a.barh(y + .2, fccN[idx], .4, label="normal")
    a.barh(y - .2, fccT[idx], .4, label="tumour")
    a.set_yticks(y); a.set_yticklabels(idx, fontsize=8)
    a.axvline(0, color="k", lw=.5)
    a.set_xlabel("median flux-control coefficient on glycolysis")
    a.set_title("A  who controls glycolytic flux"); a.legend(fontsize=8)
    a = ax[1, 0]
    cols = sorted(set(frT.columns), key=lambda r: -frT[r].abs().median())[:10]
    a.boxplot([frN[c].dropna() for c in cols], positions=np.arange(10) - .2,
              widths=.35, showfliers=False)
    b = a.boxplot([frT[c].dropna() for c in cols], positions=np.arange(10) + .2,
                  widths=.35, showfliers=False, patch_artist=True)
    for p in b["boxes"]:
        p.set_facecolor("#f4a582")
    a.set_xticks(range(10)); a.set_xticklabels(cols, rotation=60, fontsize=8)
    a.axhline(0, color="k", lw=.5)
    a.set_ylabel("total fragility dReλ/dlnE")
    a.set_title("B  fragility profiles (white normal, orange tumour)")
    a = ax[1, 1]
    labs = list(surv["uni"])
    hr = [surv["uni"][l]["hr_adj"] for l in labs]
    p = [surv["uni"][l]["p_adj"] for l in labs]
    a.barh(range(len(labs)), np.log2(hr),
           color=["#c44" if q < 0.05 else "#9aa" for q in p])
    a.set_yticks(range(len(labs))); a.set_yticklabels(labs, fontsize=8)
    a.axvline(0, color="k", lw=.5)
    a.set_xlabel("log2 HR per SD (adjusted for stage + age); red p<0.05")
    a.set_title(f"C  overall survival, n={surv['n']}, {surv['events']} deaths")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "d4_patient_jacobian.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    report()
