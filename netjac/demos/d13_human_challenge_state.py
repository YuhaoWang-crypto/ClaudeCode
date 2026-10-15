"""
D13  Human tissue-level real time series: controlled influenza challenge
     (blood transcriptome) as a perturbation-recovery state definition.

Data (GEO, open):
  GSE52428  Woods et al., PLoS One 2013 (PMID 23326326): 41 healthy adults
            inoculated with H3N2 (DEE2, n=17) or H1N1 (DEE3, n=24); blood
            sampled twice before inoculation (about -1 day and 0 h) and at
            14 times up to 108 h. Affymetrix U133A 2.0 (GPL571).
  GSE17156  Zaas et al., Cell Host Microbe 2009 (PMID 19664979): rhinovirus
            (n=20) and RSV (n=20) challenges, baseline + time of peak
            symptoms (matched time in asymptomatic subjects). Used only as
            an independent transfer test (its H3N2 arm overlaps GSE52428 and
            is dropped).

Known positive, independent of the transcriptome: the GEO 'clinical
phenotype' label (symptomatic 21 vs asymptomatic 20). The paper's stricter
'symptomatic infection' (symptoms + viral shedding) counts 18 of 41; the
results here use the GEO label.

Questions (same as D11 / D12):
  1  Does the pre-challenge (static) state define who becomes ill? Nested
     M0 baseline -> M1 + one post-challenge time -> M2 + second time ->
     M3 + full trajectory/recovery; random CV and leave-strain-out.
  2  Personal reference: is the change from one's OWN baseline better than
     the raw level against the population, time point by time point?
  3  Recovery: which modules return to the personal baseline by 108 h,
     relative to the baseline-to-baseline noise floor?
  4  Core components: genes whose induction separates symptomatic from
     asymptomatic in BOTH strains; enrichment; transfer to HRV / RSV.
  5  Minimal readout: fewest genes reaching the full-model accuracy across
     strains and viruses, against random gene sets.

Modules are fixed from textbook immunology (no data-driven selection).

Run:  python3 -m netjac.demos.d13_human_challenge_state
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import recovery as rc
from netjac.data import geo

RES = os.path.join(DATA_CACHE, "d13_results.pkl")

MODULES = {
    "ISG (type I IFN)": ["IFI44L", "IFI44", "RSAD2", "ISG15", "IFI27", "IFIT1", "IFIT2", "IFIT3",
                         "IFIT5", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFI6", "HERC5",
                         "HERC6", "EPSTI1", "CMPK2", "XAF1", "STAT1", "IRF7", "SIGLEC1", "LY6E",
                         "USP18", "DDX58", "IFIH1", "LAMP3", "SAMD9L", "PLSCR1", "RTP4", "BST2",
                         "IFITM3", "SP110", "PARP9", "DDX60", "EIF2AK2"],
    "neutrophil / inflammation": ["S100A8", "S100A9", "S100A12", "MMP9", "FCGR3B", "CXCR2",
                                  "CSF3R", "IL1R2", "ARG1", "MMP8", "LCN2", "CEACAM8"],
    "T / NK lymphocyte": ["CD3E", "CD3D", "CD3G", "CD2", "CD247", "CD8A", "LCK", "ZAP70", "IL7R",
                          "CCR7", "GZMK", "NKG7", "KLRB1", "CD28", "LEF1", "TCF7"],
    "B cell": ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "BLK", "BANK1", "FCRLA", "VPREB3"],
    "plasmablast": ["IGJ", "JCHAIN", "MZB1", "TXNDC5", "TNFRSF17", "IGLL1", "POU2AF1", "PRDM1"],
    "cell cycle": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "TYMS", "RRM2", "CCNA2", "MCM2"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1", "PPIA", "RPLP0",
                                        "GUSB", "B2M"],
}
ISG = "ISG (type I IFN)"


# ================================================================== data
def _time(v):
    if v == "Baseline":
        return -24.0
    if v == "Pre-Challenge Baseline":
        return 0.0
    return float(v.split()[0])


def load(n_genes=8000):
    X, M = geo.series_matrix("GSE52428")
    sym = geo.gpl_symbols("GPL571")
    L = np.log2(X.clip(lower=16))
    L = geo.quantile_normalise(L)
    G = geo.collapse_to_genes(L, sym)
    G = G.loc[G.mean(axis=1).sort_values(ascending=False).index[:n_genes]]
    M = M.assign(t=M["time point"].map(_time))
    grid = np.array(sorted(t for t in M.t.unique() if 0 < t <= 108))
    subj = sorted(M["subject id"].unique())
    info = M.drop_duplicates("subject id").set_index("subject id")
    y = np.array([info.loc[s, "clinical phenotype"] == "Symptomatic" for s in subj])
    virus = np.array([info.loc[s, "challenge virus"] for s in subj])
    score = np.array([int(info.loc[s, "total symptom score"]) for s in subj])
    V = G.values
    B1, B0, T = [], [], []                          # baseline -1 d, 0 h, trajectory
    for s in subj:
        m = M[M["subject id"] == s]
        get = lambda t: V[:, G.columns.get_loc(m.index[m.t == t][0])] if (m.t == t).any() else None
        b0 = get(0.0)
        b1 = get(-24.0)
        B0.append(b0)
        B1.append(b1 if b1 is not None else b0)
        ts = m[(m.t > 0) & (m.t <= 108)].sort_values("t")
        Z = V[:, [G.columns.get_loc(i) for i in ts.index]].T
        Ti = np.array([np.interp(grid, ts.t.values, Z[:, k]) for k in range(Z.shape[1])]).T
        T.append(Ti)
    return {"genes": np.array(G.index), "subj": subj, "y": y, "virus": virus, "score": score,
            "B0": np.array(B0), "B1": np.array(B1), "traj": np.array(T), "grid": grid}


def load_transfer(genes):
    X, M = geo.series_matrix("GSE17156")
    sym = geo.gpl_symbols("GPL571")
    G = geo.collapse_to_genes(geo.quantile_normalise(X), sym)
    G = G.reindex(genes)
    M = M.assign(v=M["subject id"].str.extract(r"^([A-Za-z]+)")[0])
    M = M[M.v.isin(["Rhino", "RSV"])]
    out = {"y": [], "virus": [], "base": [], "peak": []}
    for s, m in M.groupby("subject id"):
        if set(m.timepoint) != {"baseline", "T - peak symptoms"}:
            continue
        out["y"].append(m["symptom group"].iloc[0] == "Symptomatic")
        out["virus"].append(m.v.iloc[0])
        out["base"].append(G[m.index[m.timepoint == "baseline"][0]].values)
        out["peak"].append(G[m.index[m.timepoint == "T - peak symptoms"][0]].values)
    return {k: np.array(v) for k, v in out.items()}


def _idx(D, names):
    gi = {g: i for i, g in enumerate(D["genes"])}
    return np.array([gi[g] for g in names if g in gi])


def _base(D):
    return (D["B0"] + D["B1"]) / 2


def _mod_traj(D, name):
    """Module mean, change from own baseline: subjects x time."""
    i = _idx(D, MODULES[name])
    return D["traj"][:, :, i].mean(2) - _base(D)[:, i].mean(1)[:, None]


# ============================================================ question 1
def q1_nested(D):
    """Nested module models, per subject. Features are pre-specified module
    summaries (7 modules). M1 adds the change at 36 h (before most symptom
    peaks), M2 the change at 60 h, M3 the trajectory (area, peak, residual
    at 108 h)."""
    g = D["grid"]
    names = list(MODULES)
    base = np.column_stack([_base(D)[:, _idx(D, MODULES[n])].mean(1) for n in names])
    tr = {n: _mod_traj(D, n) for n in names}
    at = lambda t: np.column_stack([tr[n][:, np.argmin(abs(g - t))] for n in names])
    traj = np.column_stack([np.column_stack([np.trapezoid(tr[n], g) / g[-1], tr[n].max(1),
                                             tr[n][:, -1]]) for n in names])
    models = {"M0 pre-challenge baseline": base,
              "M1 + change at 36 h": np.hstack([base, at(36)]),
              "M2 + change at 60 h": np.hstack([base, at(36), at(60)]),
              "M3 + trajectory / recovery": np.hstack([base, at(36), at(60), traj])}
    out = {"random CV": rc.incremental_classification(models, D["y"], n_splits=5, reps=5,
                                                      n_boot=1000, C=0.3),
           "leave-strain-out": rc.incremental_classification(models, D["y"], groups=D["virus"],
                                                             n_splits=2, n_boot=1000, C=0.3)}
    # transcriptome-wide versions (top-variance genes, label-free selection)
    v = np.var(_base(D), 0)
    top = np.argsort(-v)[:500]
    Bt = _base(D)[:, top]
    d60 = D["traj"][:, np.argmin(abs(g - 60)), top] - Bt
    out["transcriptome"] = rc.incremental_classification(
        {"M0 baseline (500 genes)": Bt, "M1 + change at 60 h": np.hstack([Bt, d60])},
        D["y"], n_splits=5, reps=5, n_boot=1000, C=0.01)
    out["transcriptome leave-strain-out"] = rc.incremental_classification(
        {"M0 baseline (500 genes)": Bt, "M1 + change at 60 h": np.hstack([Bt, d60])},
        D["y"], groups=D["virus"], n_splits=2, n_boot=1000, C=0.01)
    return out


# ============================================================ question 2
def q2_personal(D):
    """Per time point AUC of symptomatic vs asymptomatic: ISG module as raw
    level (population reference) vs change from own baseline; and a
    transcriptome logistic (leave-subject-out) on raw vs baseline-subtracted
    profiles. Also the noise floor: baseline(-1 d) vs baseline(0 h)."""
    g = D["grid"]
    i = _idx(D, MODULES[ISG])
    raw = D["traj"][:, :, i].mean(2)
    delta = _mod_traj(D, ISG)
    v = np.var(_base(D), 0)
    top = np.argsort(-v)[:500]
    res = {"t": g, "isg_raw": [], "isg_delta": [], "tx_raw": [], "tx_delta": []}
    y = D["y"]
    for k in range(len(g)):
        res["isg_raw"].append(rc.auc(y, raw[:, k]))
        res["isg_delta"].append(rc.auc(y, delta[:, k]))
        Xr = D["traj"][:, k, top]
        Xd = Xr - _base(D)[:, top]
        res["tx_raw"].append(rc.auc(y, rc.oof_predictions(Xr, y, n_splits=len(y) // 2, C=0.01)))
        res["tx_delta"].append(rc.auc(y, rc.oof_predictions(Xd, y, n_splits=len(y) // 2, C=0.01)))
    res["isg_static"] = rc.auc(y, _base(D)[:, i].mean(1))
    res["isg_static_b0"] = rc.auc(y, D["B0"][:, i].mean(1))
    res["tx_static"] = rc.auc(y, rc.oof_predictions(_base(D)[:, top], y, n_splits=len(y) // 2, C=0.01))
    # baseline-to-baseline noise of the ISG module
    nb = D["B0"][:, i].mean(1) - D["B1"][:, i].mean(1)
    res["noise_sd"] = float(np.std(nb))
    res["first_t"] = {k: float(next((t for t, a in zip(g, res[k]) if a >= 0.85), np.nan))
                      for k in ("isg_raw", "isg_delta", "tx_raw", "tx_delta")}
    rng = np.random.default_rng(0)
    boots = []
    for _ in range(2000):
        b = rng.integers(0, len(y), len(y))
        if y[b].min() == y[b].max():
            continue
        boots.append([rc.auc(y[b], delta[b, k]) - rc.auc(y[b], raw[b, k]) for k in range(len(g))])
    boots = np.array(boots)
    res["isg_d_personal"] = [(float(boots[:, k].mean()), *np.percentile(boots[:, k], [2.5, 97.5]))
                             for k in range(len(g))]
    return res


# ============================================================ question 3
def q3_recovery(D, noise):
    """Per module, symptomatic and asymptomatic means of the change from own
    baseline; peak; residual at 108 h in units of the module's own
    baseline-to-baseline SD; fraction of symptomatic subjects back within
    2 SD."""
    g = D["grid"]
    out = {}
    for n in MODULES:
        i = _idx(D, MODULES[n])
        tr = _mod_traj(D, n)
        sd = float(np.std(D["B0"][:, i].mean(1) - D["B1"][:, i].mean(1)) / np.sqrt(2))
        s, a = tr[D["y"]], tr[~D["y"]]
        ms = s.mean(0)
        kp = int(np.argmax(abs(ms)))
        out[n] = {"n_genes": len(i), "sym": ms, "asym": a.mean(0),
                  "sym_se": s.std(0) / np.sqrt(len(s)), "asym_se": a.std(0) / np.sqrt(len(a)),
                  "peak": float(ms[kp]), "t_peak": float(g[kp]), "sd": sd,
                  "end_z": float(ms[-1] / sd),
                  "frac_back": float(np.mean(abs(s[:, -1]) < 2 * sd)),
                  "recovered": float(1 - ms[-1] / ms[kp]) if abs(ms[kp]) > 2 * sd else np.nan}
    # quantitative relevance: ISG area vs total symptom score in symptomatic
    from scipy.stats import spearmanr
    area = np.trapezoid(_mod_traj(D, ISG), g)
    r = spearmanr(area[D["y"]], D["score"][D["y"]])
    out["_score_rho"] = (float(r.statistic), float(r.pvalue))
    out["_score_rho_all"] = tuple(map(float, spearmanr(area, D["score"])))
    return out


# ============================================================ question 4
def q4_core(D, T2):
    """Gene-level induction area (change from own baseline, 0-108 h);
    Welch t symptomatic vs asymptomatic in each strain separately; core =
    same sign and p < 0.01 in both. Enrichment against expressed genes.
    Transfer: core-up score (peak - baseline) in rhinovirus / RSV."""
    from scipy.stats import ttest_ind
    from netjac.data import gprofiler
    g = D["grid"]
    A = np.trapezoid(D["traj"] - _base(D)[:, None, :], g, axis=1) / g[-1]
    st = {}
    for v in ("H3N2", "H1N1"):
        m = D["virus"] == v
        t, p = ttest_ind(A[m & D["y"]], A[m & ~D["y"]], equal_var=False)
        st[v] = (t, p)
    t3, p3 = st["H3N2"]
    t1, p1 = st["H1N1"]
    core_up = (t3 > 0) & (t1 > 0) & (p3 < 0.01) & (p1 < 0.01)
    core_dn = (t3 < 0) & (t1 < 0) & (p3 < 0.01) & (p1 < 0.01)
    # null: how many genes would replicate by chance (permute labels within strain)
    rng = np.random.default_rng(0)
    null = []
    for _ in range(50):
        yy = D["y"].copy()
        for v in ("H3N2", "H1N1"):
            m = np.where(D["virus"] == v)[0]
            yy[m] = rng.permutation(yy[m])
        n = np.ones(A.shape[1], bool)
        sg = None
        for v in ("H3N2", "H1N1"):
            m = D["virus"] == v
            t, p = ttest_ind(A[m & yy], A[m & ~yy], equal_var=False)
            s_ = np.sign(t)
            n &= p < 0.01
            sg = s_ if sg is None else sg
            n &= s_ == sg
        null.append(int(n.sum()))
    genes = D["genes"]
    up, dn = list(genes[core_up]), list(genes[core_dn])
    enr_up = gprofiler.enrich(up, background=list(genes), sources=("GO:BP", "REAC", "KEGG")) if up else []
    enr_dn = gprofiler.enrich(dn, background=list(genes), sources=("GO:BP", "REAC", "KEGG")) if len(dn) >= 5 else []
    isg = set(MODULES[ISG])
    # transfer to rhinovirus / RSV (independent subjects and viruses)
    tr = {}
    iu = np.where(core_up)[0]
    ii = _idx(D, MODULES[ISG])
    ok = ~np.isnan(T2["peak"]).any(0)
    iu = iu[ok[iu]]
    ii = ii[ok[ii]]
    for v in ("Rhino", "RSV", "both"):
        m = np.ones(len(T2["y"]), bool) if v == "both" else T2["virus"] == v
        yv = T2["y"][m]
        tr[v] = {"n": (int(yv.sum()), int((~yv).sum())),
                 "core static": rc.auc(yv, T2["base"][m][:, iu].mean(1)),
                 "core peak - baseline": rc.auc(yv, (T2["peak"][m] - T2["base"][m])[:, iu].mean(1)),
                 "core peak raw": rc.auc(yv, T2["peak"][m][:, iu].mean(1)),
                 "ISG static": rc.auc(yv, T2["base"][m][:, ii].mean(1)),
                 "ISG peak - baseline": rc.auc(yv, (T2["peak"][m] - T2["base"][m])[:, ii].mean(1)),
                 "ISG peak raw": rc.auc(yv, T2["peak"][m][:, ii].mean(1))}
    return {"i_up": np.where(core_up)[0], "i_dn": np.where(core_dn)[0], "i_isg": _idx(D, MODULES[ISG]),
            "t3": t3, "t1": t1, "n_up": len(up), "n_dn": len(dn), "up": up, "dn": dn,
            "null_median": float(np.median(null)), "null_max": int(max(null)),
            "isg_in_up": [x for x in up if x in isg], "enr_up": enr_up[:8], "enr_dn": enr_dn[:8],
            "rep_r": float(np.corrcoef(t3, t1)[0, 1]), "transfer": tr}


# ============================================================ question 5
def q5_minimal(D, T2, t_eval=60, ks=(1, 2, 3, 5, 8, 13, 20, 34), n_rand=500):
    """Fewest genes that carry the state. Readout = mean change from own
    baseline at t_eval of k genes, each signed by its H3N2 direction (no
    fitted weights). Genes are ranked on H3N2 only (standardised effect);
    scored on H1N1 (unseen strain) and on rhinovirus / RSV peak - baseline
    (unseen viruses). Prior alternative: single textbook ISGs and the whole
    prior ISG module. Null: random gene sets of the same size, signed the
    same way."""
    g = D["grid"]
    k = np.argmin(abs(g - t_eval))
    Dx = D["traj"][:, k, :] - _base(D)
    m3, m1 = D["virus"] == "H3N2", D["virus"] == "H1N1"
    ok = ~np.isnan(T2["peak"]).any(0)
    T2d = T2["peak"] - T2["base"]
    y1, y2 = D["y"][m1], T2["y"]
    sd = Dx[m3].std(0) + 1e-9
    eff = (Dx[m3 & D["y"]].mean(0) - Dx[m3 & ~D["y"]].mean(0)) / sd
    pool = np.where(ok)[0]
    rank = pool[np.argsort(-abs(eff[pool]))]

    def ev(sel):
        sel = list(sel)
        w = np.sign(eff[sel])
        return (rc.auc(y1, (Dx[m1][:, sel] * w).mean(1)), rc.auc(y2, (T2d[:, sel] * w).mean(1)))
    rng = np.random.default_rng(1)
    out = {"t_eval": float(g[k]), "ks": list(ks), "top": [], "null": [], "top_genes": list(D["genes"][rank[:8]])}
    for n in ks:
        out["top"].append(ev(rank[:n]))
        r = np.array([ev(rng.choice(pool, n, replace=False)) for _ in range(n_rand)])
        out["null"].append((np.percentile(r, 50, 0).tolist(), np.percentile(r, 95, 0).tolist()))
    prior = {}
    gi = {x: i for i, x in enumerate(D["genes"])}
    for gname in ("RSAD2", "IFI44L", "ISG15", "IFI27", "SIGLEC1"):
        if gname in gi and ok[gi[gname]]:
            i = gi[gname]
            prior[gname] = (rc.auc(y1, Dx[m1][:, i]), rc.auc(y2, T2d[:, i]))
    ii = [i for i in _idx(D, MODULES[ISG]) if ok[i]]
    prior["ISG module (34)"] = (rc.auc(y1, Dx[m1][:, ii].mean(1)), rc.auc(y2, T2d[:, ii].mean(1)))
    out["prior"] = prior
    # collective mode: label-free PCA of the change from own baseline (all 41)
    Z = Dx[:, pool] - Dx[:, pool].mean(0)
    U, S, Vt = np.linalg.svd(Z, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    isg = Dx[:, ii].mean(1)
    sgn = np.sign(np.corrcoef(pc1, isg)[0, 1])
    out["pc1"] = {"var_frac": float(S[0] ** 2 / (S ** 2).sum()),
                  "r_isg": float(abs(np.corrcoef(pc1, isg)[0, 1])),
                  "auc_all": rc.auc(D["y"], sgn * pc1),
                  "frac_genes_r>0.5": float(np.mean(abs(np.corrcoef(np.column_stack([pc1, Z]).T)[0, 1:]) > 0.5))}
    return out


# ================================================================= report
def run():
    D = load()
    T2 = load_transfer(D["genes"])
    q2 = q2_personal(D)
    R = {"n": len(D["subj"]), "n_sym": int(D["y"].sum()), "n_genes": len(D["genes"]),
         "grid": D["grid"], "q1": q1_nested(D), "q2": q2, "q3": q3_recovery(D, q2["noise_sd"]),
         "q4": q4_core(D, T2), "q5": q5_minimal(D, T2), "n_transfer": len(T2["y"])}
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    print(f"\nD13 influenza challenge (GSE52428): {R['n']} adults ({R['n_sym']} symptomatic), "
          f"{R['n_genes']} genes, {len(R['grid'])} post-inoculation times; transfer cohort {R['n_transfer']} (HRV/RSV)")
    print("Q1 nested models (subjects; AUC, FP at 80% sensitivity, dAUC vs previous [95% CI]):")
    for cv, res in R["q1"].items():
        for m, v in res.items():
            s = f"   {cv:16s} {m:28s} AUC {v['auc']:.3f} FP {v['fpr']:4.0%}"
            if "d_prev" in v:
                s += f"  dAUC {v['d_prev'][0]:+.3f} [{v['d_prev'][1]:+.3f}, {v['d_prev'][2]:+.3f}]  dFP {v['d_fpr'][0]:+.2f}"
            print(s)
    q = R["q2"]
    print(f"Q2 time-resolved AUC (ISG module; transcriptome LOSO). static: ISG {q['isg_static']:.2f} "
          f"(0 h only {q['isg_static_b0']:.2f}), transcriptome {q['tx_static']:.2f}; ISG baseline-baseline SD {q['noise_sd']:.3f}")
    for k, t in enumerate(q["t"]):
        print(f"   {t:6.1f} h  ISG raw {q['isg_raw'][k]:.2f}  ISG vs own baseline {q['isg_delta'][k]:.2f}  "
              f"| transcriptome raw {q['tx_raw'][k]:.2f}  vs own baseline {q['tx_delta'][k]:.2f}")
    print("   first time AUC >= 0.85:", {k: v for k, v in q["first_t"].items()})
    print("   ISG own-baseline minus raw dAUC [95% CI]: " + ", ".join(
        f"{t:g} h {d[0]:+.2f} [{d[1]:+.2f}, {d[2]:+.2f}]" for t, d in zip(q["t"], q["isg_d_personal"]) if 21 < t < 70))
    print("Q3 recovery in symptomatic subjects (change from own baseline, log2):")
    for n, v in R["q3"].items():
        if n.startswith("_"):
            continue
        print(f"   {n:32s} genes {v['n_genes']:2d}  peak {v['peak']:+.2f} at {v['t_peak']:5.1f} h  "
              f"108 h residual {v['sym'][-1]:+.2f} ({v['end_z']:+.1f} SD)  back within 2 SD {v['frac_back']:.0%}  "
              f"asym peak {v['asym'][np.argmax(abs(v['asym']))]:+.2f}")
    print(f"   ISG area vs symptom score: symptomatic only rho {R['q3']['_score_rho'][0]:+.2f} (p={R['q3']['_score_rho'][1]:.2g}); "
          f"all rho {R['q3']['_score_rho_all'][0]:+.2f} (p={R['q3']['_score_rho_all'][1]:.1g})")
    q = R["q4"]
    print(f"Q4 core: {q['n_up']} up / {q['n_dn']} down replicate in both strains (p<0.01 each; permutation null median "
          f"{q['null_median']:.0f}, max {q['null_max']}); strain t-stat r {q['rep_r']:.2f}")
    print(f"   ISG among core-up: {len(q['isg_in_up'])}/{len(MODULES[ISG])}: {q['isg_in_up']}")
    print("   up enrichment: " + "; ".join(f"{x['name'][:45]} ({x['source']}, p={x['p']:.0e}, k={x['k']})" for x in q["enr_up"][:5]))
    print("   down enrichment: " + ("; ".join(f"{x['name'][:45]} ({x['source']}, p={x['p']:.0e}, k={x['k']})" for x in q["enr_dn"][:5]) or "n.s."))
    print(f"   down genes: {q['dn'][:20]}")
    for v, d in q["transfer"].items():
        print(f"   transfer {v:6s} n sym/asym {d['n']}: " + ", ".join(f"{k} {x:.2f}" for k, x in d.items() if k != "n"))
    q = R["q5"]
    print(f"Q5 minimal readout (change from own baseline at {q['t_eval']:g} h; genes ranked on H3N2 only; top genes {q['top_genes']}):")
    for n, (a1, a2), (med, p95) in zip(q["ks"], q["top"], q["null"]):
        print(f"   top {n:2d} genes  H1N1 {a1:.2f}  HRV+RSV {a2:.2f}  | random sets median {med[0]:.2f}/{med[1]:.2f}, 95th pct {p95[0]:.2f}/{p95[1]:.2f}")
    c = q["pc1"]
    print(f"   collective mode (label-free PC1 of change at {q['t_eval']:g} h): {c['var_frac']:.0%} of variance, |r| with ISG {c['r_isg']:.2f}, "
          f"AUC {c['auc_all']:.2f}, genes with |r|>0.5 to PC1 {c['frac_genes_r>0.5']:.0%}")
    print("   prior single genes / module (H1N1, HRV+RSV): " + ", ".join(f"{g} {a:.2f}/{b:.2f}" for g, (a, b) in q["prior"].items()))


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    g = R["grid"]
    q3 = R["q3"]
    a = ax[0, 0]
    v = q3[ISG]
    for key, c, lab in (("sym", "C3", "symptomatic"), ("asym", "C0", "asymptomatic")):
        a.plot(g, v[key], "o-", color=c, label=f"{lab} (n={R['n_sym'] if key == 'sym' else R['n'] - R['n_sym']})")
        a.fill_between(g, v[key] - v[key + "_se"], v[key] + v[key + "_se"], color=c, alpha=0.2)
    a.axhspan(-2 * v["sd"], 2 * v["sd"], color="0.85", label="±2 SD baseline-to-baseline")
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel("hours after inoculation")
    a.set_ylabel("ISG module, change from own baseline (log2)")
    a.set_title("Perturbation: influenza inoculation (H3N2 + H1N1)\nresponse and return of the interferon module")
    a.legend(fontsize=8)
    a = ax[0, 1]
    q = R["q2"]
    a.plot(g, q["isg_delta"], "o-", color="C3", label="ISG vs own baseline")
    a.plot(g, q["isg_raw"], "s--", color="C3", alpha=0.5, label="ISG raw level")
    a.plot(g, q["tx_delta"], "o-", color="C0", label="transcriptome vs own baseline (LOSO)")
    a.plot(g, q["tx_raw"], "s--", color="C0", alpha=0.5, label="transcriptome raw (LOSO)")
    a.scatter([-8, -8], [q["isg_static"], q["tx_static"]], c=["C3", "C0"], marker="D", s=50, zorder=5)
    a.annotate("pre-challenge\n(static)", (-8, 0.36), fontsize=8, ha="center")
    a.axhline(0.5, color="k", lw=0.6, ls=":")
    a.set_ylim(0.3, 1.02)
    a.set_xlabel("hours after inoculation")
    a.set_ylabel("AUC symptomatic vs asymptomatic")
    a.set_title("Static state does not predict illness;\nthe response does from ~45 h (own baseline vs raw: n.s.)")
    a.legend(fontsize=7, loc="lower right")
    a = ax[0, 2]
    for cv, c in (("random CV", "C3"), ("leave-strain-out", "C1")):
        res = R["q1"][cv]
        a.plot(range(4), [res[m]["auc"] for m in res], "o-", color=c, label=f"{cv}: AUC")
        a.plot(range(4), [res[m]["fpr"] for m in res], "s--", color=c, alpha=0.6, label=f"{cv}: FP at 80% sens.")
    a.set_xticks(range(4), ["M0 baseline", "M1 +36 h", "M2 +60 h", "M3 traj."], fontsize=8)
    a.set_ylim(0, 1.02)
    a.set_title("Nested models M0 -> M3 (7 prior modules)")
    a.legend(fontsize=7)
    a = ax[1, 0]
    for k, n in enumerate(MODULES):
        if n == ISG:
            continue
        a.plot(g, q3[n]["sym"], "-", lw=2, color=f"C{k}", label=f"{n} ({q3[n]['end_z']:+.1f} SD at 108 h)")
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel("hours after inoculation")
    a.set_ylabel("symptomatic mean, change from own baseline")
    a.set_title("Other modules in symptomatic subjects:\nwho has not returned by 108 h")
    a.legend(fontsize=7)
    a = ax[1, 1]
    q4 = R["q4"]
    a.scatter(q4["t3"], q4["t1"], s=3, c="0.7", alpha=0.4)
    a.scatter(q4["t3"][q4["i_up"]], q4["t1"][q4["i_up"]], s=8, c="C3", label="core up (both strains)")
    a.scatter(q4["t3"][q4["i_dn"]], q4["t1"][q4["i_dn"]], s=8, c="C0", label="core down")
    a.scatter(q4["t3"][q4["i_isg"]], q4["t1"][q4["i_isg"]], s=30, facecolors="none", edgecolors="k", lw=0.8,
              label="prior ISG module")
    a.legend(fontsize=7, loc="lower right")
    a.axhline(0, color="k", lw=0.5)
    a.axvline(0, color="k", lw=0.5)
    a.set_xlabel("H3N2: t symptomatic vs asymptomatic (induction area)")
    a.set_ylabel("H1N1: t")
    lab = q4["enr_up"][0]["name"][:40] if q4["enr_up"] else "n.s."
    a.set_title(f"Core genes: {q4['n_up']} up / {q4['n_dn']} down in both strains\n"
                f"(null median {q4['null_median']:.0f}); top term: {lab}", fontsize=10)
    a = ax[1, 2]
    q5 = R["q5"]
    ks = q5["ks"]
    for j, (lab, c) in enumerate((("H1N1 (unseen strain)", "C3"), ("HRV + RSV (unseen viruses)", "C2"))):
        a.plot(ks, [t[j] for t in q5["top"]], "o-", color=c, label=f"top-k genes ranked on H3N2 -> {lab}")
        a.fill_between(ks, [n[0][j] for n in q5["null"]], [n[1][j] for n in q5["null"]], color=c, alpha=0.12)
    a.axhline(q5["prior"]["ISG module (34)"][0], color="C3", ls=":", lw=1)
    a.axhline(q5["prior"]["ISG module (34)"][1], color="C2", ls=":", lw=1)
    a.set_xscale("log")
    a.set_xticks(ks, [str(k) for k in ks])
    a.set_xlabel("genes in readout")
    a.set_ylabel("AUC")
    a.set_ylim(0.4, 1.02)
    a.set_title("Minimal readout (band: random sets median-95th pct;\ndotted: prior ISG module)")
    a.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    p = os.path.join(FIG, "d13_human_challenge_state.png")
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
