"""
D14  Does treatment return a human tissue to its balance state?
     Real longitudinal treatment data with clinical response as the known
     positive and, in one study, a placebo arm as the negative control.

Data (GEO, open):
  A  GSE85034  psoriasis skin, Correa da Rosa et al., J Invest Dermatol 2017
               (PMID 27667537): 29 patients, non-lesional (NL) + lesional
               (LS) at week 0, LS at weeks 1, 2, 4, 16 on adalimumab (ADA) or
               methotrexate (MTX); PASI at every visit. Illumina HT-12.
  B  GSE130588 atopic dermatitis skin, Guttman-Yassky et al., J Allergy Clin
               Immunol 2019 (PMID 30194992): dupilumab vs PLACEBO, LS at
               weeks 0, 4, 16, NL at weeks 0, 16, plus skin of 20 healthy
               volunteers; EASI at every visit. Affymetrix U133 Plus 2.
  C  GSE16879  inflammatory bowel disease mucosa, Arijs et al., PLoS One 2009
               (PMID 19956723): UC, colonic and ileal Crohn's, before and
               4-6 weeks after first infliximab, response by endoscopy /
               histology; 12 control biopsies.

Balance state = the patient's own non-lesional skin (A, B) or control
mucosa (C). Restoration is measured by PROJECTION on the patient's own
disease axis d = LS0 - reference (the D12 lesson: distance is noise
dominated):
    R(t) = 1 - <x_t - ref, d> / |d|^2      (0 = untreated, 1 = back to balance)
    off-axis = |(x_t - ref) - (1 - R) d| / |d|   (new abnormality)
Disease genes are chosen from LS0 vs reference only (label-free with
respect to response).

Questions:
  1  Does the static pre-treatment state predict who responds, or does the
     early restoration (weeks 1-4) define it? (nested M0 -> M3, against the
     clinical early-improvement baseline)
  2  Does treatment restore the balance state? Responders vs non-responders,
     drug vs drug, drug vs placebo; molecular restoration vs clinical score.
  3  Which modules restore first and which stay (residual disease)?
  4  Which module's early restoration is the lever that predicts response?
  5  Transfer to another tissue / drug (IBD, infliximab) and the choice of
     reference (NL vs healthy skin).

Run:  python3 -m netjac.demos.d14_treatment_restoration
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import recovery as rc
from netjac.data import geo

RES = os.path.join(DATA_CACHE, "d14_results.pkl")

SKIN = {
    "IL-23/IL-17 cytokines": ["IL17A", "IL17F", "IL22", "IL23A", "IL12B", "IL19", "IL36G", "IL36A"],
    "IL-17 keratinocyte response": ["DEFB4A", "DEFB4B", "S100A7A", "S100A7", "S100A8", "S100A9",
                                    "S100A12", "LCN2", "PI3", "SERPINB4", "SERPINB3", "CXCL1",
                                    "CCL20", "KRT16", "KRT6A", "KRT17"],
    "Th2": ["CCL17", "CCL18", "CCL22", "CCL26", "CCL13", "IL13", "IL4R", "IL31", "TNFRSF4"],
    "interferon": ["MX1", "ISG15", "IFI27", "IFI44L", "RSAD2", "OAS1", "OAS2", "STAT1", "IFIT1",
                   "IFIT3", "CXCL10", "CXCL9"],
    "T cell / DC infiltrate": ["CD3E", "CD3D", "CD2", "CD8A", "ITGAX", "LAMP3", "CCR7", "GZMB",
                               "CD274", "IL2RA"],
    "proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "TYMS", "RRM2", "CCNA2"],
    "barrier / lipid (down in lesions)": ["CCL27", "BTC", "WIF1", "THRSP", "ELOVL3", "KRT77",
                                          "KRT79", "FLG", "LOR", "CLDN23", "FA2H", "AQP9"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1", "PPIA", "RPLP0",
                                        "GUSB"],
}
GUT = {
    "inflammation": ["CXCL1", "CXCL8", "IL1B", "S100A8", "LCN2", "DUOX2", "MMP3", "MMP1",
                     "TNFAIP6", "OSM", "TREM1", "REG1A", "REG3A", "CHI3L1", "SAA1"],
    "epithelial metabolism (down in disease)": ["AQP8", "SLC26A2", "HMGCS2", "CA1", "ABCG2",
                                                "PCK1", "CLDN8", "GUCA2B", "SLC30A10", "UGT2A3"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1", "PPIA", "RPLP0",
                                        "GUSB"],
}


# ================================================================== helpers
def _genes(X, gpl, log=False):
    if log:
        X = np.log2(X.clip(lower=16))
    return geo.collapse_to_genes(geo.quantile_normalise(X), geo.gpl_symbols(gpl))


def disease_genes(LS, REF, thr=1.0, q=0.05):
    """Paired LS - reference: |mean| > thr log2 and BH q < `q`."""
    from scipy.stats import ttest_1samp
    D = LS - REF
    t, p = ttest_1samp(D, 0, axis=0)
    p = np.nan_to_num(p, nan=1.0)                    # constant genes (floor values)
    o = np.argsort(p)
    bh = np.empty_like(p)
    bh[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    return np.where((abs(D.mean(0)) > thr) & (bh < q))[0]


def loo_axes(diffs):
    """Leave-one-out cohort disease axis for each patient: the mean
    disease deviation of all OTHER patients. Using the patient's own
    deviation as the axis would turn measurement noise in the pre-treatment
    sample into apparent restoration (regression to the mean)."""
    diffs = np.asarray(diffs)
    tot = diffs.sum(0)
    return (tot[None] - diffs) / (len(diffs) - 1)


def restoration(x, ref, ls0, idx, u):
    """R = 1 - <x - ref, u> / <ls0 - ref, u> on genes `idx` (u: cohort axis).
    Returns R and the off-axis deviation at t relative to the on-axis
    baseline deviation (new abnormality not along the disease direction)."""
    u = u[idx] / np.linalg.norm(u[idx])
    e0, et = (ls0 - ref)[idx], (x - ref)[idx]
    a0, at = e0 @ u, et @ u
    return 1 - at / a0, float(np.linalg.norm(et - at * u) / abs(a0))


def baseline_projection(ls0, ref, idx, u):
    """On-axis size of the pre-treatment deviation. A biopsy whose
    projection is < 25% of the cohort median carries no detectable disease
    signature, so its restoration ratio is undefined (pre-specified
    exclusion)."""
    u = u[idx] / np.linalg.norm(u[idx])
    return float((ls0 - ref)[idx] @ u)


def auc_ci(y, s, n_boot=2000, seed=0):
    y, s = np.asarray(y), np.asarray(s)
    m = ~np.isnan(s)
    y, s = y[m], s[m]
    rng = np.random.default_rng(seed)
    b = []
    for _ in range(n_boot):
        i = rng.integers(0, len(y), len(y))
        if y[i].min() != y[i].max():
            b.append(rc.auc(y[i], s[i]))
    return (rc.auc(y, s), *np.percentile(b, [2.5, 97.5]))


def module_index(gi, names, u, thr=0.585, kmin=3):
    """Module genes that are disrupted in the cohort (|mean deviation| >
    `thr` log2); None if fewer than `kmin` (module not part of the disease
    state, so 'restoration' is undefined)."""
    idx = [gi[g] for g in names if g in gi and abs(u[gi[g]]) > thr]
    return np.array(idx) if len(idx) >= kmin else None


def module_restoration(x, ref, ls0, midx, u):
    return restoration(x, ref, ls0, midx, u)[0] if midx is not None else np.nan


# ========================================================== A  psoriasis
def load_psoriasis():
    X, M = geo.series_matrix("GSE85034")
    G = _genes(X, "GPL10558")
    M = M.assign(subj=M.title.str.extract(r"Subject (\d+)")[0],
                 tp=M.timepoint.replace({"WK2": "WK 2", "WK 1LS": "WK 1"}))
    M = M[M.tp != "WK 1NL"]
    pts = []
    for s, m in M.groupby("subj"):
        tp = dict(zip(m.tp, m.index))
        if not {"NL", "LS"} <= set(tp):
            continue
        pasi = {k: float(v) for k, v in zip(m.tp, m.pasi) if v not in ("NA", "")}
        pts.append({"id": s, "drug": m.treatment.iloc[0],
                    "x": {k: G[v].values for k, v in tp.items()}, "pasi": pasi})
    return {"genes": np.array(G.index), "pts": pts}


def analyse_psoriasis(P):
    gi = {g: i for i, g in enumerate(P["genes"])}
    pts = [p for p in P["pts"] if "WK 16" in p["pasi"] and "LS" in p["pasi"]]
    LS = np.array([p["x"]["LS"] for p in P["pts"]])
    NL = np.array([p["x"]["NL"] for p in P["pts"]])
    dis = disease_genes(LS, NL)
    dset = set(dis)
    up = int(((LS - NL).mean(0)[dis] > 0).sum())
    U = dict(zip([p["id"] for p in P["pts"]], loo_axes(LS - NL)))
    ucoh = (LS - NL).mean(0)
    midx = {mn: module_index(gi, gl, ucoh) for mn, gl in SKIN.items()}
    T = ["WK 1", "WK 2", "WK 4", "WK 16"]
    rows = []
    for p in pts:
        r = {"id": p["id"], "drug": p["drug"], "pasi0": p["pasi"]["LS"],
             "off WK 0": restoration(p["x"]["LS"], p["x"]["NL"], p["x"]["LS"], dis, U[p["id"]])[1],
             "a0": baseline_projection(p["x"]["LS"], p["x"]["NL"], dis, U[p["id"]]),
             "resp": p["pasi"]["WK 16"] <= 0.25 * p["pasi"]["LS"],
             "pasi_impr16": 1 - p["pasi"]["WK 16"] / p["pasi"]["LS"],
             "pasi_impr4": 1 - p["pasi"].get("WK 4", np.nan) / p["pasi"]["LS"],
             "mag": float(np.linalg.norm((p["x"]["LS"] - p["x"]["NL"])[dis]))}
        for t in T:
            if t in p["x"]:
                u = U[p["id"]]
                r[f"R {t}"], r[f"off {t}"] = restoration(p["x"][t], p["x"]["NL"], p["x"]["LS"], dis, u)
                for mname in SKIN:
                    r[f"{mname} | {t}"] = module_restoration(p["x"][t], p["x"]["NL"], p["x"]["LS"],
                                                             midx[mname], u)
            else:
                r[f"R {t}"] = r[f"off {t}"] = np.nan
        rows.append(r)
    a0 = np.array([r["a0"] for r in rows])
    excluded = [r["id"] for r in rows if r["a0"] < 0.25 * np.median(a0)]
    keep = [r["id"] not in excluded for r in rows]
    rows = [r for r, k in zip(rows, keep) if k]
    pts = [p for p, k in zip(pts, keep) if k]
    y = np.array([r["resp"] for r in rows])
    get = lambda k: np.array([r.get(k, np.nan) for r in rows], float)
    ok = ~np.isnan(get("R WK 4")) & ~np.isnan(get("pasi_impr4")) & ~np.isnan(get("R WK 1")) & ~np.isnan(get("R WK 2"))
    drug = np.array([r["drug"] == "ADA" for r in rows], float)
    base = np.column_stack([get("pasi0"), get("mag"), drug])
    models = {"M0 static: PASI0, lesion magnitude, drug": base[ok],
              "M1 + molecular restoration wk 1": np.column_stack([base, get("R WK 1")])[ok],
              "M2 + wk 2": np.column_stack([base, get("R WK 1"), get("R WK 2")])[ok],
              "M3 + wk 4": np.column_stack([base, get("R WK 1"), get("R WK 2"), get("R WK 4")])[ok]}
    nested = rc.incremental_classification(models, y[ok], n_splits=5, reps=10, n_boot=1000, C=0.5)
    clin = rc.incremental_classification(
        {"M0 static": base[ok], "+ clinical PASI improvement wk 4": np.column_stack([base, get("pasi_impr4")])[ok],
         "+ molecular restoration wk 4": np.column_stack([base, get("pasi_impr4"), get("R WK 4")])[ok]},
        y[ok], n_splits=5, reps=10, n_boot=1000, C=0.5)
    # single pre-specified readouts (no fitting)
    single = {"static PASI0": auc_ci(y, get("pasi0")), "static lesion magnitude": auc_ci(y, get("mag")),
              "clinical PASI improvement wk 4": auc_ci(y, get("pasi_impr4"))}
    for t in T[:3]:
        single[f"molecular restoration {t}"] = auc_ci(y, get(f"R {t}"))
    for t in T[:3]:
        single[f"IL-17 keratinocyte module {t}"] = auc_ci(y, get(f"IL-17 keratinocyte response | {t}"))
    # restoration validity: molecular vs clinical at week 16
    from scipy.stats import spearmanr, mannwhitneyu
    m16 = ~np.isnan(get("R WK 16"))
    rho = spearmanr(get("R WK 16")[m16], get("pasi_impr16")[m16])
    curves = {}
    for grp, sel in (("responders", y), ("non-responders", ~y)):
        curves[grp] = [float(np.nanmedian(get(f"R {t}")[sel])) for t in T]
        curves[grp + " q"] = [np.nanpercentile(get(f"R {t}")[sel], [25, 75]).tolist() for t in T]
    by_drug = {d: {"n": int((get("pasi0") > 0)[[r["drug"] == d for r in rows]].sum()),
                   "responders": int(y[[r["drug"] == d for r in rows]].sum()),
                   "R16": float(np.nanmedian(get("R WK 16")[[r["drug"] == d for r in rows]]))}
               for d in ("ADA", "MTX")}
    rr = get("R WK 16")
    resid = {"R16 responders median": float(np.nanmedian(rr[y])),
             "R16 non-responders median": float(np.nanmedian(rr[~y])),
             "p": float(mannwhitneyu(rr[y & m16], rr[~y & m16]).pvalue),
             "frac responders R16 >= 0.75": float(np.nanmean(rr[y] >= 0.75)),
             "frac responders R16 >= 0.9": float(np.nanmean(rr[y] >= 0.9)),
             "off-axis wk0": float(np.nanmedian(get("off WK 0"))),
             "off-axis wk16 responders": float(np.nanmedian(get("off WK 16")[y])),
             "off-axis wk16 non-responders": float(np.nanmedian(get("off WK 16")[~y]))}
    # module restoration curves in responders; module early AUC
    mods = {}
    for mname in SKIN:
        mods[mname] = {"resp": [float(np.nanmedian(get(f"{mname} | {t}")[y])) for t in T],
                       "nonresp": [float(np.nanmedian(get(f"{mname} | {t}")[~y])) for t in T],
                       "n_dis": 0 if midx[mname] is None else len(midx[mname]),
                       "genes": [] if midx[mname] is None else list(P["genes"][midx[mname]])}
        for t in ("WK 1", "WK 2", "WK 4"):
            v = get(f"{mname} | {t}")
            m = ~np.isnan(v)
            mods[mname][f"auc {t}"] = rc.auc(y[m], v[m]) if m.sum() > 5 and y[m].min() != y[m].max() else np.nan
    # gene-level residual disease in responders at wk 16
    resp_pts = [p for p, r in zip(pts, rows) if r["resp"] and "WK 16" in p["x"]]
    D0 = np.array([(p["x"]["LS"] - p["x"]["NL"])[dis] for p in resp_pts])
    D16 = np.array([(p["x"]["WK 16"] - p["x"]["NL"])[dis] for p in resp_pts])
    gene_rest = 1 - D16.mean(0) / D0.mean(0)
    from scipy.stats import ttest_1samp
    _, pres = ttest_1samp(D16, 0, axis=0)
    residual = dis[(gene_rest < 0.5) & (pres < 0.05)]
    from netjac.data import gprofiler
    rg = list(P["genes"][residual])
    enr = gprofiler.enrich(rg, background=list(P["genes"]), sources=("GO:BP", "REAC", "KEGG")) if len(rg) >= 5 else []
    return {"n": len(rows), "n_resp": int(y.sum()), "n_dis": len(dis), "n_up": up, "excluded": excluded,
            "nested": nested, "clin": clin, "single": single, "rho16": (float(rho.statistic), float(rho.pvalue)),
            "curves": curves, "T": T, "by_drug": by_drug, "resid": resid, "mods": mods,
            "residual_genes": rg, "residual_frac": len(rg) / len(dis), "residual_enr": enr[:6],
            "gene_rest_median": float(np.median(gene_rest)),
            "rows": rows}


# ===================================================== B  atopic dermatitis
def load_ad():
    X, M = geo.series_matrix("GSE130588")
    G = _genes(X, "GPL570")
    return G, M


def analyse_ad(G, M):
    from scipy.stats import mannwhitneyu, spearmanr
    genes = np.array(G.index)
    gi = {g: i for i, g in enumerate(genes)}
    V = lambda i: G[i].values
    normal = np.array([V(i) for i in M.index[M.tissue == "Normal"]])
    nl_all = np.array([V(i) for i in M.index[(M.tissue == "NL") & (M.time == "Week0")]])
    pts = []
    for s, m in M[M.tissue != "Normal"].groupby("subject id"):
        key = {(t, ti): i for t, ti, i in zip(m.time, m.tissue, m.index)}
        if ("Week0", "LS") not in key:
            continue
        e = {(t, ti): float(v) for t, ti, v in zip(m.time, m.tissue, m.easi) if v not in ("NA", "")}
        pts.append({"id": s, "arm": m.treatment.iloc[0], "key": key, "easi": e})
    LS0 = np.array([V(p["key"][("Week0", "LS")]) for p in pts if ("Week0", "NL") in p["key"]])
    NL0 = np.array([V(p["key"][("Week0", "NL")]) for p in pts if ("Week0", "NL") in p["key"]])
    dis = disease_genes(LS0, NL0)
    dset = set(dis)
    nl_mean = nl_all.mean(0)
    paired = [p["id"] for p in pts if ("Week0", "NL") in p["key"]]
    U = dict(zip(paired, loo_axes(LS0 - NL0)))
    ucoh = (LS0 - NL0).mean(0)
    midx = {mn: module_index(gi, gl, ucoh) for mn, gl in SKIN.items()}
    rows = []
    for p in pts:
        k = p["key"]
        ref = V(k[("Week0", "NL")]) if ("Week0", "NL") in k else nl_mean
        ls0 = V(k[("Week0", "LS")])
        u0 = U.get(p["id"], ucoh)
        r = {"id": p["id"], "arm": p["arm"], "easi0": p["easi"].get(("Week0", "LS"), np.nan),
             "a0": baseline_projection(ls0, ref, dis, u0), "off0": restoration(ls0, ref, ls0, dis, u0)[1]}
        for t in ("Week4", "Week16"):
            if (t, "LS") in k:
                u = U.get(p["id"], ucoh)
                r[f"R {t}"], r[f"off {t}"] = restoration(V(k[(t, "LS")]), ref, ls0, dis, u)
                for mname in SKIN:
                    r[f"{mname} | {t}"] = module_restoration(V(k[(t, "LS")]), ref, ls0, midx[mname], u)
                r[f"easi {t}"] = p["easi"].get((t, "LS"), np.nan)
            else:
                r[f"R {t}"] = r[f"off {t}"] = r[f"easi {t}"] = np.nan
        # non-lesional skin under treatment: change toward healthy skin
        if ("Week0", "NL") in k and ("Week16", "NL") in k:
            a = (normal.mean(0) - V(k[("Week0", "NL")]))[dis]
            b = (V(k[("Week16", "NL")]) - V(k[("Week0", "NL")]))[dis]
            r["NL toward healthy"] = float(b @ a / (a @ a))
        else:
            r["NL toward healthy"] = np.nan
        rows.append(r)
    a0 = np.array([r["a0"] for r in rows])
    excluded = [r["id"] for r in rows if r["a0"] < 0.25 * np.median(a0)]
    rows = [r for r in rows if r["id"] not in excluded]
    arm = np.array([r["arm"] for r in rows])
    get = lambda key: np.array([r.get(key, np.nan) for r in rows], float)
    out = {"excluded": excluded, "off0": float(np.median(get("off0"))),
           "off16": {a: float(np.nanmedian(get("off Week16")[arm == a])) for a in ("Dupilumab", "Placebo")},
           "n_dis": len(dis), "n": {a: int((arm == a).sum()) for a in ("Dupilumab", "Placebo")},
           "mod_genes": {mn: (0 if v is None else len(v)) for mn, v in midx.items()}}
    for t in ("Week4", "Week16"):
        R = get(f"R {t}")
        d, pl = R[(arm == "Dupilumab") & ~np.isnan(R)], R[(arm == "Placebo") & ~np.isnan(R)]
        out[t] = {"dupi": (float(np.median(d)), len(d)), "placebo": (float(np.median(pl)), len(pl)),
                  "p": float(mannwhitneyu(d, pl).pvalue), "auc_arm": rc.auc(np.r_[np.ones(len(d)), np.zeros(len(pl))], np.r_[d, pl]),
                  "mods": {mn: (float(np.nanmedian(get(f"{mn} | {t}")[arm == "Dupilumab"])),
                                float(np.nanmedian(get(f"{mn} | {t}")[arm == "Placebo"])))
                           for mn in SKIN}}
    # molecular vs clinical, all patients with week-16 LS
    impr = 1 - get("easi Week16") / get("easi0")
    m = ~np.isnan(impr) & ~np.isnan(get("R Week16"))
    rho = spearmanr(get("R Week16")[m], impr[m])
    out["rho16"] = (float(rho.statistic), float(rho.pvalue), int(m.sum()))
    # early (week 4) restoration predicts week-16 EASI-75 within dupilumab
    mm = (arm == "Dupilumab") & ~np.isnan(impr) & ~np.isnan(get("R Week4"))
    y75 = impr[mm] >= 0.75
    out["dupi_easi75"] = (int(y75.sum()), int((~y75).sum()))
    if y75.min() != y75.max():
        out["auc_R4_easi75"] = rc.auc(y75, get("R Week4")[mm])
        out["auc_easi0_easi75"] = rc.auc(y75, get("easi0")[mm])
        e4 = 1 - get("easi Week4")[mm] / get("easi0")[mm]
        out["auc_clin4_easi75"] = rc.auc(y75, e4)
    # reference choice: NL vs healthy skin
    out["nl_vs_normal_n_genes"] = int(len(disease_genes(nl_all[:min(len(nl_all), len(normal))],
                                                        normal[:min(len(nl_all), len(normal))], thr=1.0)))
    nt = get("NL toward healthy")
    out["NL toward healthy"] = {a: (float(np.nanmedian(nt[arm == a])), int(np.sum(~np.isnan(nt[arm == a]))))
                                for a in ("Dupilumab", "Placebo")}
    # where does treated lesional skin end up: fraction of the LS0 -> healthy gap closed vs LS0 -> NL gap
    gap = []
    for p in pts:
        k = p["key"]
        if p["arm"] == "Dupilumab" and ("Week16", "LS") in k and ("Week0", "NL") in k:
            ls0, nl, x = V(k[("Week0", "LS")]), V(k[("Week0", "NL")]), V(k[("Week16", "LS")])
            u = U[p["id"]]
            gap.append((restoration(x, nl, ls0, dis, u)[0], restoration(x, normal.mean(0), ls0, dis, u)[0]))
    out["toward NL vs healthy"] = np.median(np.array(gap), 0).tolist() if gap else []
    out["rows"] = rows
    return out


# ================================================================= C  IBD
def analyse_ibd():
    from scipy.stats import mannwhitneyu
    X, M = geo.series_matrix("GSE16879")
    G = _genes(X, "GPL570", log=True)
    genes = np.array(G.index)
    gi = {g: i for i, g in enumerate(genes)}
    M = M.assign(pid=M.title.str.replace(r"_(before|after)T$", "", regex=True),
                 when=M.title.str.extract(r"_(before|after)T$")[0])
    out = {}
    for tissue, ctrl_prefix, groups in (("Colon", "CO", ("UC", "CDc")), ("Ileum", "IL", ("CDi",))):
        ctrl = np.array([G[i].values for i in M.index[M.title.str.match(rf"^{ctrl_prefix}\d+$")]])
        ref = ctrl.mean(0)
        P = M[M.pid.str.match("^(" + "|".join(groups) + ")")]
        pairs = []
        for pid, m in P.groupby("pid"):
            if set(m.when) == {"before", "after"}:
                pairs.append({"pid": pid, "resp": m["response to infliximab"].iloc[0] == "Yes",
                              "before": G[m.index[m.when == "before"][0]].values,
                              "after": G[m.index[m.when == "after"][0]].values})
        B = np.array([p["before"] for p in pairs])
        for thr in (1.0, 0.585, 0.3):                # ileal Crohn's is heterogeneous
            dis = disease_genes(B, np.repeat(ref[None], len(B), 0), thr=thr)
            if len(dis) >= 50:
                break
        y = np.array([p["resp"] for p in pairs])
        UU = loo_axes(B - ref)
        ucoh = (B - ref).mean(0)
        midx = {mn: module_index(gi, gl, ucoh) for mn, gl in GUT.items()}
        a0 = np.array([baseline_projection(p["before"], ref, dis, u) for p, u in zip(pairs, UU)])
        keep = a0 >= 0.25 * np.median(a0)
        pairs = [p for p, k in zip(pairs, keep) if k]
        UU = UU[keep]
        B = B[keep]
        y = y[keep]
        off0 = np.array([restoration(p["before"], ref, p["before"], dis, u)[1] for p, u in zip(pairs, UU)])
        R = np.array([restoration(p["after"], ref, p["before"], dis, u) for p, u in zip(pairs, UU)])
        mag = np.array([np.linalg.norm((p["before"] - ref)[dis]) for p in pairs])
        modR = {mn: np.array([module_restoration(p["after"], ref, p["before"], midx[mn], u)
                              for p, u in zip(pairs, UU)]) for mn in GUT}
        infl0 = np.array([np.mean([p["before"][gi[g]] for g in GUT["inflammation"] if g in gi]) for p in pairs])
        # static transcriptome prediction (top-variance genes among patients, LOO)
        top = np.argsort(-B.var(0))[:500]
        ptx = rc.oof_predictions(B[:, top], y, n_splits=min(len(y) // 2, 10), C=0.01)
        out[tissue] = {"n": (int(y.sum()), int((~y).sum())), "n_dis": len(dis), "thr": thr,
                       "n_excluded": int((~keep).sum()), "off0": float(np.median(off0)),
                       "auc_R_ci": auc_ci(y, R[:, 0]), "auc_mag_ci": auc_ci(y, -np.array([np.linalg.norm((p["before"] - ref)[dis]) for p in pairs])),
                       "R resp": float(np.median(R[y, 0])), "R nonresp": float(np.median(R[~y, 0])),
                       "auc_R": rc.auc(y, R[:, 0]), "p": float(mannwhitneyu(R[y, 0], R[~y, 0]).pvalue),
                       "off resp": float(np.median(R[y, 1])), "off nonresp": float(np.median(R[~y, 1])),
                       "frac resp R>=0.75": float(np.mean(R[y, 0] >= 0.75)),
                       "auc static magnitude": rc.auc(y, -mag), "auc static inflammation": rc.auc(y, -infl0),
                       "auc static transcriptome": rc.auc(y, ptx),
                       "mods": {mn: (float(np.nanmedian(v[y])), float(np.nanmedian(v[~y]))) for mn, v in modR.items()},
                       "R": R[:, 0], "y": y}
    return out


# ================================================================= report
def run():
    R = {"A": analyse_psoriasis(load_psoriasis())}
    R["B"] = analyse_ad(*load_ad())
    R["C"] = analyse_ibd()
    pickle.dump(R, open(RES, "wb"))
    return R


def _fmt(v):
    return f"{v:+.2f}" if isinstance(v, float) else str(v)


def report(R):
    A = R["A"]
    print(f"\nD14 A psoriasis (GSE85034): {A['n']} patients with week-16 PASI, {A['n_resp']} PASI75 responders; "
          f"{A['n_dis']} disease genes (LS vs own NL; {A['n_up']} up); excluded (no lesional signature): {A['excluded']}")
    print("   by drug:", A["by_drug"])
    print("   Q1 single pre-specified readouts, AUC for PASI75 at week 16:")
    for k, v in A["single"].items():
        print(f"      {k:38s} {v[0]:.2f} [{v[1]:.2f}, {v[2]:.2f}]")
    for name, res in (("nested molecular", A["nested"]), ("vs clinical early improvement", A["clin"])):
        print(f"   Q1 {name}:")
        for m, v in res.items():
            s = f"      {m:44s} AUC {v['auc']:.3f} FP {v['fpr']:4.0%}"
            if "d_prev" in v:
                s += f"  dAUC {v['d_prev'][0]:+.3f} [{v['d_prev'][1]:+.3f}, {v['d_prev'][2]:+.3f}]"
            print(s)
    print(f"   Q2 restoration R(t) (0 untreated, 1 = own NL): responders {np.round(A['curves']['responders'], 2)} "
          f"non-responders {np.round(A['curves']['non-responders'], 2)} at {A['T']}")
    print(f"      week 16: {A['resid']}; molecular R16 vs PASI improvement rho {A['rho16'][0]:+.2f} (p={A['rho16'][1]:.1g})")
    print(f"   Q3 residual disease in responders: {len(A['residual_genes'])} of {A['n_dis']} disease genes restored < 50% "
          f"(gene median restoration {A['gene_rest_median']:.2f}); examples {A['residual_genes'][:25]}")
    print("      enrichment: " + ("; ".join(f"{x['name'][:45]} ({x['source']}, p={x['p']:.0e}, k={x['k']})" for x in A["residual_enr"]) or "n.s."))
    print("   Q3/Q4 modules (median restoration in responders over weeks 1,2,4,16 | non-responders | AUC of early restoration for PASI75 wk1/2/4):")
    for mn, v in A["mods"].items():
        print(f"      {mn:36s} n_dis {v['n_dis']:2d}  resp {np.round(v['resp'], 2)}  nonresp {np.round(v['nonresp'], 2)}  "
              f"AUC {v['auc WK 1']:.2f}/{v['auc WK 2']:.2f}/{v['auc WK 4']:.2f}")
    B = R["B"]
    print(f"\nD14 B atopic dermatitis (GSE130588): {B['n']} ; {B['n_dis']} disease genes; module genes {B['mod_genes']}; "
          f"excluded {B['excluded']}; off-axis wk0 {B['off0']:.2f}, wk16 {B['off16']}")
    for t in ("Week4", "Week16"):
        v = B[t]
        print(f"   {t}: restoration dupilumab {v['dupi'][0]:+.2f} (n={v['dupi'][1]}) vs placebo {v['placebo'][0]:+.2f} "
              f"(n={v['placebo'][1]}), p={v['p']:.1g}, AUC {v['auc_arm']:.2f}")
        print("      modules (dupilumab, placebo): " + "; ".join(f"{mn.split(' (')[0]} {a:+.2f}/{b:+.2f}" for mn, (a, b) in v["mods"].items()))
    print(f"   molecular R16 vs EASI improvement: rho {B['rho16'][0]:+.2f} (p={B['rho16'][1]:.1g}, n={B['rho16'][2]})")
    print(f"   within dupilumab, EASI-75 at wk16 (yes/no {B['dupi_easi75']}): AUC R(wk4) {B.get('auc_R4_easi75', np.nan):.2f}, "
          f"static EASI0 {B.get('auc_easi0_easi75', np.nan):.2f}, clinical EASI improvement wk4 {B.get('auc_clin4_easi75', np.nan):.2f}")
    print(f"   reference choice: NL vs healthy skin differ in {B['nl_vs_normal_n_genes']} genes (|log2FC|>1, q<0.05); "
          f"NL moves toward healthy (fraction of gap, wk16): {B['NL toward healthy']}; dupilumab LS wk16 restoration "
          f"toward own NL vs toward healthy skin: {np.round(B['toward NL vs healthy'], 2)}")
    print("\nD14 C IBD infliximab (GSE16879):")
    for tissue, v in R["C"].items():
        print(f"   {tissue}: responders/non {v['n']}, {v['n_dis']} disease genes (|log2FC|>{v['thr']}); restoration median resp {v['R resp']:+.2f} vs "
              f"non {v['R nonresp']:+.2f} (AUC {v['auc_R']:.2f}, p={v['p']:.1g}); responders R>=0.75 {v['frac resp R>=0.75']:.0%}; "
              f"off-axis after {v['off resp']:.2f}/{v['off nonresp']:.2f} (before {v['off0']:.2f}); excluded {v['n_excluded']}")
        print(f"      AUC restoration {np.round(v['auc_R_ci'], 2)}, static magnitude {np.round(v['auc_mag_ci'], 2)}")
        print(f"      static before treatment, AUC for response: disease magnitude {v['auc static magnitude']:.2f}, "
              f"inflammation module {v['auc static inflammation']:.2f}, transcriptome (CV) {v['auc static transcriptome']:.2f}")
        print("      modules (resp, non-resp): " + "; ".join(f"{mn.split(' (')[0]} {a:+.2f}/{b:+.2f}" for mn, (a, b) in v["mods"].items()))


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    A = R["A"]
    wk = [1, 2, 4, 16]
    a = ax[0, 0]
    for grp, c in (("responders", "C2"), ("non-responders", "C3")):
        m, q = np.array(A["curves"][grp]), np.array(A["curves"][grp + " q"])
        a.plot([0] + wk, [0] + list(m), "o-", color=c, label=f"PASI75 {grp} (median, IQR)")
        a.fill_between(wk, q[:, 0], q[:, 1], color=c, alpha=0.2)
    a.axhline(1, color="k", ls=":", lw=1)
    a.text(16, 1.02, "own non-lesional skin", ha="right", fontsize=8)
    a.set_xlabel("week of treatment (adalimumab or methotrexate)")
    a.set_ylabel("restoration along own disease axis")
    a.set_title(f"A psoriasis: does treatment return skin to balance?\n"
                f"(n={A['n']}, {A['n_resp']} responders; {A['n_dis']} disease genes)")
    a.legend(fontsize=8)
    a = ax[0, 1]
    keys = list(A["single"])
    v = np.array([A["single"][k] for k in keys])
    a.barh(range(len(keys)), v[:, 0], xerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]],
           color=["0.6", "0.6", "C1"] + ["C2"] * 3 + ["C4"] * 3, error_kw={"lw": 1})
    a.axvline(0.5, color="k", lw=0.8, ls=":")
    a.set_yticks(range(len(keys)), keys, fontsize=8)
    a.invert_yaxis()
    a.set_xlim(0.2, 1.05)
    a.set_xlabel("AUC for PASI75 at week 16")
    a.set_title("A static state vs early restoration\n(pre-specified single readouts)")
    a = ax[0, 2]
    mods = [m for m in A["mods"] if A["mods"][m]["n_dis"] > 0]
    for k, mn in enumerate(mods):
        a.plot([0] + wk, [0] + A["mods"][mn]["resp"], "o-", color=f"C{k}", label=mn)
    a.axhline(1, color="k", ls=":", lw=1)
    a.set_xlabel("week")
    a.set_ylabel("module restoration (responders, median)")
    a.set_title("A which modules restore first, which stay")
    a.legend(fontsize=7)
    a = ax[1, 0]
    B = R["B"]
    for j, t in enumerate(("Week4", "Week16")):
        rows = B["rows"]
        for k, (armn, c) in enumerate((("Dupilumab", "C2"), ("Placebo", "0.5"))):
            v = [r[f"R {t}"] for r in rows if r["arm"] == armn and not np.isnan(r[f"R {t}"])]
            x = j * 3 + k + np.random.default_rng(k).uniform(-0.15, 0.15, len(v))
            a.scatter(x, v, s=14, color=c, alpha=0.8)
            a.plot([j * 3 + k - 0.3, j * 3 + k + 0.3], [np.median(v)] * 2, color="k", lw=2)
    a.set_xticks([0, 1, 3, 4], ["dupilumab\nwk 4", "placebo\nwk 4", "dupilumab\nwk 16", "placebo\nwk 16"])
    a.axhline(0, color="k", lw=0.6)
    a.axhline(1, color="k", ls=":", lw=1)
    a.set_ylabel("restoration toward own non-lesional skin")
    a.set_title(f"B atopic dermatitis: drug vs placebo (negative control)\n"
                f"wk16 p={B['Week16']['p']:.0g}; R16 vs EASI improvement rho {B['rho16'][0]:.2f}")
    a = ax[1, 1]
    mn = [m for m in SKIN if B["mod_genes"][m] > 0]
    d = [B["Week4"]["mods"][m][0] for m in mn]
    p = [B["Week4"]["mods"][m][1] for m in mn]
    yy = np.arange(len(mn))
    a.barh(yy - 0.2, d, 0.4, color="C2", label="dupilumab")
    a.barh(yy + 0.2, p, 0.4, color="0.6", label="placebo")
    a.set_yticks(yy, mn, fontsize=8)
    a.invert_yaxis()
    a.axvline(1, color="k", ls=":", lw=1)
    a.axvline(0, color="k", lw=0.6)
    a.set_xlabel("module restoration at week 4 (median; placebo n=17)")
    a.set_title("B which modules dupilumab restores (specificity)")
    a.legend(fontsize=8)
    a = ax[1, 2]
    C = R["C"]
    for j, (tissue, v) in enumerate(C.items()):
        for k, (sel, c, lab) in enumerate(((v["y"], "C2", "responder"), (~v["y"], "C3", "non-responder"))):
            vals = v["R"][sel]
            x = j * 3 + k + np.random.default_rng(j * 2 + k).uniform(-0.15, 0.15, len(vals))
            a.scatter(x, vals, s=14, color=c, label=lab if j == 0 else None)
            a.plot([j * 3 + k - 0.3, j * 3 + k + 0.3], [np.median(vals)] * 2, color="k", lw=2)
    a.set_xticks([0.5, 3.5], [f"colon (UC + CD)\nAUC {C['Colon']['auc_R']:.2f}", f"ileum (CD)\nAUC {C['Ileum']['auc_R']:.2f}"])
    a.axhline(1, color="k", ls=":", lw=1)
    a.axhline(0, color="k", lw=0.6)
    a.set_ylabel("restoration toward control mucosa")
    a.set_title("C transfer: IBD mucosa after infliximab")
    a.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(FIG, "d14_treatment_restoration.png")
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
