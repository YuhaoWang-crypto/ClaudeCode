"""
D22  Remove the cause: what comes back, and what stays?

Every treatment cohort so far (D14, D19, D20) showed the same ordering:
inflammatory and proliferative modules restore first, structural ones last.
Those were all drugs holding a disease down. This demo asks the cleaner
question - what happens when the CAUSE ITSELF is removed and nothing is
being suppressed - in two tissues that have nothing to do with skin.

PART A  GSE7895, bronchial epithelium (Spira / Beane, PMID 17347440).
21 never smokers, 52 current smokers, 31 former smokers whose time since
quitting runs from 1 month to 504 months (42 years). A cross-sectional
dose-of-time design: cause removed, no drug, and a reference group.
  Q1  which modules return to the never-smoker state, and how fast?
  Q2  what is still displaced after 10+ years without the cause?
  Q3  the two confounders this design forces: former smokers are OLDER and
      have MORE pack-years than current smokers. Both are tested, and the
      time-since-quit effect is reported with and without them.

PART B  GSE145358, duodenal biopsies in coeliac disease (the "covert
condition on a gluten-free diet" cohort). 15 patients sampled on a
long-term gluten-free diet and again after a gluten challenge, plus 6
healthy controls on a normal diet.
  Q4  with the cause removed for years, is the mucosa back to healthy?
  Q5  when the cause returns, does the tissue move back along the SAME axis
      (the attractor question of D19/D20, now in a third disease and with
      the cause itself rather than a drug)?

Run:  python3 -m netjac.demos.d22_scar_after_cause_removed
"""
import os
import pickle
import re

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac.data import geo

RES = os.path.join(DATA_CACHE, "d22_results.pkl")

# ------------------------------------------------- modules, pre-specified
AIRWAY = {
    "xenobiotic metabolism": ["CYP1A1", "CYP1B1", "ALDH3A1", "AKR1B10", "AKR1C1",
                              "AKR1C2", "AKR1C3", "ADH7", "CBR1", "NQO1"],
    "NRF2 / antioxidant": ["GCLC", "GCLM", "TXNRD1", "SLC7A11", "HMOX1", "GPX2",
                           "PIR", "TALDO1", "SRXN1"],
    "inflammation": ["IL8", "CXCL1", "CXCL2", "IL1B", "CCL20", "PTGS2", "SOD2",
                     "TNFAIP3", "NFKBIA"],
    # Ciliary genes present on this platform (HG-U133A). The list is
    # functional, not chosen by how they recover; the >0.2 disruption rule
    # then decides which enter the module.
    "cilia / ciliated differentiation": ["FOXJ1", "DNAI1", "DNAI2", "DNALI1",
                                         "DNAH5", "DNAH9", "DNAAF1", "SPAG1",
                                         "SPAG6", "TEKT2", "TEKT3", "ZMYND10",
                                         "IFT57", "IFT88", "WDR78", "LRRC6",
                                         "GAS8", "RSPH4A", "RSPH9", "TPPP3"],
    "secretory / mucous": ["MUC5AC", "MUC5B", "TFF3", "SCGB3A1", "AGR2",
                           "CEACAM5", "CEACAM6"],
    "squamous / basal structure": ["KRT6A", "KRT13", "KRT14", "SPRR1A", "SPRR1B",
                                   "SPRR3", "S100A2", "TGM1", "IVL", "CSTA"],
    "proliferation": ["MKI67", "TOP2A", "CCNB1", "CDK1", "BUB1", "TYMS", "RRM2",
                      "CCNA2"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1",
                                        "PPIA", "RPLP0", "GUSB"],
}
COELIAC = {
    "interferon / chemokine": ["CXCL9", "CXCL10", "CXCL11", "STAT1", "GBP1",
                               "GBP5", "IDO1", "IRF1"],
    "cytotoxic T / IEL": ["GZMB", "GZMA", "PRF1", "CD3E", "CD8A", "IL15",
                          "CCL5", "NKG7"],
    "epithelial absorption (down in disease)": ["APOA4", "ALPI", "LCT", "SI",
                                                "TREH", "FABP6", "SLC5A1",
                                                "ANPEP", "MGAM", "APOB"],
    "crypt proliferation": ["MKI67", "OLFM4", "LGR5", "CCNB1", "TOP2A", "RRM2"],
    "matrix remodelling": ["MMP1", "MMP3", "MMP12", "COL1A1", "COL3A1", "TIMP1"],
    "housekeeping (negative control)": ["ACTB", "GAPDH", "PGK1", "TBP", "HPRT1",
                                        "PPIA", "RPLP0", "GUSB"],
}


# ------------------------------------------------------------------ helpers
def _module_dir(G, ref, dis, genes, thr=0.2):
    """Disrupted module genes and the unit disease direction on them."""
    gi = {g: i for i, g in enumerate(G)}
    idx = [gi[g] for g in genes if g in gi and abs(dis[gi[g]] - ref[gi[g]]) > thr]
    if len(idx) < 3:
        return None, None, len(idx)
    idx = np.array(idx)
    u = dis[idx] - ref[idx]
    return idx, u / np.linalg.norm(u), len(idx)


def _phi(x, ref, idx, u, scale):
    return float((x[idx] - ref[idx]) @ u) / scale


# ------------------------------------------------------------------ part A
def load_a():
    X, M = geo.series_matrix("GSE7895")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)),
                              geo.gpl_symbols("GPL96"))
    grp, age, py, quit_m = [], [], [], []
    for gsm in M.index:
        s = M.loc[gsm]
        grp.append(s.iloc[0])
        # the characteristics line is "Age:40 Packyears:5 Time Since Quit ...",
        # so GEO's key/value split leaves the rest of it inside the Age column
        txt = " ".join(str(v) for v in s.values)
        m = re.match(r"\s*(\d+)", str(s.get("Age", "")))
        age.append(float(m.group(1)) if m else np.nan)
        m = re.search(r"Packyears:([\d.]+)", txt)
        py.append(float(m.group(1)) if m else 0.0)
        m = re.search(r"Quit Smoking \(months\):(\d+)", txt)
        quit_m.append(float(m.group(1)) if m else np.nan)
    import pandas as pd
    S = pd.DataFrame({"group": grp, "age": age, "packyears": py,
                      "quit_months": quit_m}, index=M.index)
    return {"G": G, "S": S, "genes": np.array(G.index)}


def part_a(D, seed=0):
    import pandas as pd
    from scipy import stats
    G, S = D["G"], D["S"]
    never = S.index[S.group == "Never Smoker"]
    curr = S.index[S.group == "Current Smoker"]
    form = S.index[S.group == "Former Smoker"]
    ref = G[never].values.mean(1)
    dis = G[curr].values.mean(1)
    out = {"n": {"never": len(never), "current": len(curr), "former": len(form)},
           "age": {k: float(S.loc[v, "age"].median())
                   for k, v in (("never", never), ("current", curr), ("former", form))},
           "packyears": {k: float(S.loc[v, "packyears"].median())
                         for k, v in (("current", curr), ("former", form))},
           "quit_months": [float(x) for x in sorted(S.loc[form, "quit_months"])]}

    rows = {}
    for mn, genes in AIRWAY.items():
        idx, u, n = _module_dir(np.array(G.index), ref, dis, genes)
        if idx is None:
            rows[mn] = {"n_disrupted": n, "n_total": len(genes), "defined": False}
            continue
        scale = float(np.mean([_phi(G[c].values, ref, idx, u, 1.0) for c in curr])) or 1e-9
        phi = {s: _phi(G[s].values, ref, idx, u, scale) for s in G.columns}
        # Age adjustment. Never smokers are 31 years younger than former
        # smokers (median 28 vs 59), so an unadjusted residual could be age.
        # Fit phi ~ age INSIDE the never smokers (ages 21-60, so this is
        # interpolation, not extrapolation) and subtract that expectation.
        a_sl0, a_ic0 = np.polyfit(S.loc[never, "age"].values,
                                  np.array([phi[s] for s in never]), 1)
        adj = {s: phi[s] - (a_sl0 * S.loc[s, "age"] + a_ic0) for s in G.columns}
        t = S.loc[form, "quit_months"].values
        y = np.array([phi[s] for s in form])          # raw is the primary read
        lt = np.log2(t)
        sl, ic, r, p, se = stats.linregress(lt, y)
        # with age and pack-years as covariates
        Xc = np.column_stack([np.ones(len(y)), lt, S.loc[form, "age"].values,
                              S.loc[form, "packyears"].values])
        beta, *_ = np.linalg.lstsq(Xc, y, rcond=None)
        resid = y - Xc @ beta
        dfree = len(y) - Xc.shape[1]
        sigma2 = resid @ resid / dfree
        cov = sigma2 * np.linalg.pinv(Xc.T @ Xc)
        tval = beta[1] / np.sqrt(cov[1, 1])
        # long-quit formers vs never smokers
        long_q = [s for s in form if S.loc[s, "quit_months"] >= 120]
        nv = np.array([phi[s] for s in never])
        lq = np.array([phi[s] for s in long_q])
        u_p = stats.mannwhitneyu(lq, nv).pvalue if len(lq) >= 3 else np.nan
        lq_adj = np.array([adj[s] for s in long_q])
        nv_adj = np.array([adj[s] for s in never])
        u_p_adj = (stats.mannwhitneyu(lq_adj, nv_adj).pvalue
                   if len(lq) >= 3 else np.nan)
        # age sensitivity inside never smokers (the size of the adjustment)
        a_sl, _, a_r, a_p, _ = stats.linregress(S.loc[never, "age"].values,
                                                np.array([phi[s] for s in never]))
        # Age-matched window. Never smokers and long-quit formers overlap at
        # 40-60 years (n = 5 and 7, ages 40-60 in both), which is the only
        # part of this cohort where the smoking effect and age are separable:
        # over the whole cohort corr(group, age) = 0.79, so the regression
        # adjustment above is near-collinear and NOT identifiable.
        nv_m = [s for s in never if 40 <= S.loc[s, "age"] <= 60]
        lq_m = [s for s in long_q if 40 <= S.loc[s, "age"] <= 60]
        p_m = (stats.mannwhitneyu([phi[s] for s in lq_m],
                                  [phi[s] for s in nv_m]).pvalue
               if len(nv_m) >= 3 and len(lq_m) >= 3 else np.nan)
        rows[mn] = {"defined": True, "n_disrupted": n, "n_total": len(genes),
                    "n_matched": (len(nv_m), len(lq_m)),
                    "phi_matched_never": float(np.median([phi[s] for s in nv_m])),
                    "phi_matched_longquit": float(np.median([phi[s] for s in lq_m])),
                    "p_matched": float(p_m),
                    "phi_current": float(np.median([phi[s] for s in curr])),
                    "phi_former_median": float(np.median(y)),
                    "phi_long_quit_median": float(np.median(lq)) if len(lq) else np.nan,
                    "phi_long_quit_age_adj": float(np.median(lq_adj)) if len(lq) else np.nan,
                    "p_longquit_vs_never_age_adj": float(u_p_adj),
                    "n_long_quit": len(lq),
                    "slope_per_doubling": float(sl), "p_slope": float(p),
                    "slope_adj": float(beta[1]),
                    "p_slope_adj": float(2 * stats.t.sf(abs(tval), dfree)),
                    "beta_age": float(beta[2]), "beta_packyears": float(beta[3]),
                    "p_longquit_vs_never": float(u_p),
                    "never_phi_sd": float(nv.std()),
                    "age_slope_in_never": float(a_sl), "age_p_in_never": float(a_p),
                    "phi": {s: float(phi[s]) for s in G.columns},
                    "phi_age_adj": {s: float(adj[s]) for s in G.columns}}
    out["modules"] = rows
    from statsmodels.stats.multitest import multipletests as _mt
    keys = [k for k, v in rows.items() if v["defined"] and v["p_matched"] == v["p_matched"]]
    if keys:
        qs = _mt([rows[k]["p_matched"] for k in keys], method="fdr_bh")[1]
        for k, qv in zip(keys, qs):
            rows[k]["q_matched"] = float(qv)

    # gene level: of the genes smoking moves, how many are still moved after 10 y?
    t_s, p_s = stats.ttest_ind(G[curr].values, G[never].values, axis=1)
    from statsmodels.stats.multitest import multipletests
    q = multipletests(np.nan_to_num(p_s, nan=1.0), method="fdr_bh")[1]
    moved = np.where((q < 0.05) & (np.abs(dis - ref) > 0.5))[0]
    long_q = [s for s in form if S.loc[s, "quit_months"] >= 120]
    if len(moved) and len(long_q) >= 3:
        lqm = G[long_q].values[moved].mean(1)
        t2, p2 = stats.ttest_ind(G[long_q].values[moved], G[never].values[moved], axis=1)
        q2 = multipletests(np.nan_to_num(p2, nan=1.0), method="fdr_bh")[1]
        persist = (q2 < 0.05) & (np.sign(lqm - ref[moved]) == np.sign(dis[moved] - ref[moved]))
        # age-adjusted sensitivity: per gene, regress on (intercept, age,
        # is_long_quit) over never smokers + long-quit formers
        grp = list(never) + list(long_q)
        Xg = np.column_stack([np.ones(len(grp)), S.loc[grp, "age"].values,
                              [0] * len(never) + [1] * len(long_q)])
        Y = G[grp].values[moved]
        beta_g, *_ = np.linalg.lstsq(Xg, Y.T, rcond=None)
        resid_g = Y.T - Xg @ beta_g
        df_g = len(grp) - 3
        s2 = (resid_g ** 2).sum(0) / df_g
        cov_d = np.linalg.pinv(Xg.T @ Xg)[2, 2]
        tg = beta_g[2] / np.sqrt(np.maximum(s2 * cov_d, 1e-12))
        pg = 2 * stats.t.sf(np.abs(tg), df_g)
        qg = multipletests(np.nan_to_num(pg, nan=1.0), method="fdr_bh")[1]
        persist_adj = (qg < 0.05) & (np.sign(beta_g[2]) == np.sign(dis[moved] - ref[moved]))
        out["genes"] = {"n_moved_by_smoking": int(len(moved)),
                        "n_still_moved_after_10y": int(persist.sum()),
                        "fraction": float(persist.mean()),
                        "n_still_moved_age_adj": int(persist_adj.sum()),
                        "fraction_age_adj": float(persist_adj.mean()),
                        "persistent": [str(x) for x in np.array(G.index)[moved][persist]][:40]}
        mods = {}
        for mn, genes in AIRWAY.items():
            inm = set(genes) & set(np.array(G.index)[moved])
            if inm:
                pset = set(np.array(G.index)[moved][persist])
                mods[mn] = (len(inm & pset), len(inm))
        out["genes"]["by_module"] = mods
    return out


# ------------------------------------------------------------------ part B
def load_b():
    import pandas as pd
    import subprocess
    d = os.path.join(DATA_CACHE, "geo_gse145358")
    fn = os.path.join(d, "umi_counts.txt.gz")
    if not os.path.exists(fn):
        os.makedirs(d, exist_ok=True)
        subprocess.run(["curl", "-sSL", "--max-time", "900", "-o", fn,
                        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145358/"
                        "suppl/GSE145358_umi_counts.txt.gz"], check=True)
    C = pd.read_csv(fn, sep="\t")
    C = C.dropna(subset=["SYMBOL"])
    vals = C.iloc[:, 3:]
    L = np.log2(vals.astype(float) + 1)
    L.index = C.SYMBOL.values
    L = L.groupby(level=0).max()
    L = geo.quantile_normalise(L)
    cols = list(L.columns)
    pat = {}
    for c in cols:
        if c.startswith("NORM"):
            continue
        m = re.match(r"(.+?)(BL|EOS)$", c)
        if m:
            pat.setdefault(m.group(1), {})[m.group(2)] = c
    ctrl = [c for c in cols if c.startswith("NORM")]
    paired = {k: v for k, v in pat.items() if "BL" in v and "EOS" in v}
    return {"L": L, "ctrl": ctrl, "paired": paired}


def part_b(D):
    from scipy import stats
    L, ctrl, paired = D["L"], D["ctrl"], D["paired"]
    G = np.array(L.index)
    ref = L[ctrl].values.mean(1)                      # healthy controls
    bl = np.array([L[v["BL"]].values for v in paired.values()])     # on GFD
    eos = np.array([L[v["EOS"]].values for v in paired.values()])   # after gluten
    dis = eos.mean(0)                                  # challenged state
    out = {"n_patients": len(paired), "n_controls": len(ctrl), "n_genes": len(G)}

    rows = {}
    for mn, genes in COELIAC.items():
        idx, u, n = _module_dir(G, ref, dis, genes)
        if idx is None:
            rows[mn] = {"defined": False, "n_disrupted": n, "n_total": len(genes)}
            continue
        scale = float(np.mean([(e[idx] - ref[idx]) @ u for e in eos])) or 1e-9
        gfd = np.array([(b[idx] - ref[idx]) @ u / scale for b in bl])
        chal = np.array([(e[idx] - ref[idx]) @ u / scale for e in eos])
        ctl = np.array([(L[c].values[idx] - ref[idx]) @ u / scale for c in ctrl])
        rows[mn] = {"defined": True, "n_disrupted": n, "n_total": len(genes),
                    "phi_gfd": float(np.median(gfd)), "phi_challenge": float(np.median(chal)),
                    "phi_control": float(np.median(ctl)),
                    "p_gfd_vs_control": float(stats.mannwhitneyu(gfd, ctl).pvalue),
                    "p_challenge_vs_gfd": float(stats.wilcoxon(chal, gfd).pvalue),
                    "gfd": gfd.tolist(), "challenge": chal.tolist(), "control": ctl.tolist()}
    out["modules"] = rows
    from statsmodels.stats.multitest import multipletests as _mt
    keys = [k for k, v in rows.items() if v["defined"]]
    if keys:
        qs = _mt([rows[k]["p_gfd_vs_control"] for k in keys], method="fdr_bh")[1]
        for k, qv in zip(keys, qs):
            rows[k]["q_gfd_vs_control"] = float(qv)

    # Q5 does the challenge move along the same axis as the residual on the
    # gluten-free diet? The axis is built LEAVE-ONE-PATIENT-OUT, otherwise the
    # patient's own challenge sample is inside the axis it is scored against.
    cos, res = [], []
    for i in range(len(bl)):
        others = [j for j in range(len(bl)) if j != i]
        axis = eos[others].mean(0) - ref
        big = np.argsort(-np.abs(axis))[:500]
        a = axis[big] - axis[big].mean()

        def c(v):
            v = v[big] - v[big].mean()
            return float(v @ a / max(np.linalg.norm(v) * np.linalg.norm(a), 1e-12))
        cos.append(c(eos[i] - bl[i]))
        res.append(c(bl[i] - ref))
    out["challenge_cos_with_disease_axis"] = {"median": float(np.median(cos)),
                                              "values": cos}
    out["gfd_residual_cos_with_disease_axis"] = {"median": float(np.median(res)),
                                                 "values": res}
    # The negative control module is the benchmark here: in a tissue whose cell
    # composition changes this much, bulk "housekeeping" genes move too, so a
    # module's residual only counts if it clears the control's.
    hk = rows.get("housekeeping (negative control)")
    if hk and hk.get("defined"):
        out["control_drift_gfd"] = hk["phi_gfd"]
    return out


# ------------------------------------------------------------------ report
def report(R):
    A, B = R["A"], R["B"]
    print("=" * 78)
    print("D22  remove the cause: what comes back and what stays")
    print("=" * 78)
    print(f"\nPART A  bronchial epithelium, GSE7895 "
          f"({A['n']['never']} never, {A['n']['current']} current, "
          f"{A['n']['former']} former smokers)")
    print(f"  median age: never {A['age']['never']:.0f}, current {A['age']['current']:.0f}, "
          f"former {A['age']['former']:.0f}   |   median pack-years: current "
          f"{A['packyears']['current']:.0f}, former {A['packyears']['former']:.0f}")
    q = A["quit_months"]
    print(f"  time since quitting: {min(q):.0f} to {max(q):.0f} months "
          f"(median {np.median(q):.0f})")
    print("\n  Q1/Q2  module position (1 = current smoker, 0 = never smoker)")
    print(f"    {'module':34s} {'curr':>6s} {'former':>8s} {'>10y':>8s} {'>10y age-adj':>13s} "
          f"{'slope/2x':>9s} {'p':>7s} {'>10y vs never':>14s} {'age-matched 40-60':>19s}")
    for mn, r in sorted(A["modules"].items(),
                        key=lambda kv: -(kv[1].get("phi_long_quit_median", -9)
                                         if kv[1]["defined"] else -9)):
        if not r["defined"]:
            print(f"    {mn:34s} undefined ({r['n_disrupted']} of {r['n_total']} disrupted)")
            continue
        print(f"    {mn:34s} {r['phi_current']:+6.2f} {r['phi_former_median']:+8.2f} "
              f"{r['phi_long_quit_median']:+8.2f} {r['phi_long_quit_age_adj']:+13.2f} "
              f"{r['slope_per_doubling']:+9.3f} {r['p_slope']:7.3f} "
              f"{r['p_longquit_vs_never']:14.4f}   {r['phi_matched_longquit']:+.2f} vs "
              f"{r['phi_matched_never']:+.2f} p={r['p_matched']:.3f} q={r.get('q_matched', float('nan')):.3f}")
    print(f"    (n with >= 10 years since quitting: "
          f"{A['modules'][list(A['modules'])[0]].get('n_long_quit', '?')})")
    print(f"\n  Q3  confounders. Former smokers are 31 years older than never smokers and "
          "corr(group, age) = 0.79 over the\n      whole cohort, so the regression "
          "adjustment is near-collinear; the age-matched 40-60 window "
          f"(n = {A['modules'][[k for k, v in A['modules'].items() if v['defined']][0]]['n_matched']}) "
          "is the identifiable test.")
    for mn, r in A["modules"].items():
        if r["defined"]:
            print(f"    {mn:34s} age slope within never smokers {r['age_slope_in_never']:+.4f} "
                  f"(p {r['age_p_in_never']:.3f});  in formers beta_age {r['beta_age']:+.4f}, "
                  f"beta_packyears {r['beta_packyears']:+.4f}")
    g = A.get("genes")
    if g:
        print(f"\n  gene level: {g['n_moved_by_smoking']} genes moved by smoking "
              f"(BH q<0.05, |log2FC|>0.5); after 10+ years without smoking "
              f"{g['n_still_moved_after_10y']} are still moved in the same direction "
              f"({100 * g['fraction']:.0f}%); with age as a covariate "
              f"{g['n_still_moved_age_adj']} ({100 * g['fraction_age_adj']:.0f}%)")
        for mn, (a, b) in sorted(g["by_module"].items(), key=lambda kv: -kv[1][0]):
            print(f"      {mn:34s} {a}/{b} of its smoking-moved genes persist")
        print("      e.g. " + ", ".join(g["persistent"][:18]))

    print(f"\nPART B  coeliac duodenum, GSE145358 ({B['n_patients']} paired patients on a "
          f"gluten-free diet then challenged, {B['n_controls']} healthy controls)")
    print(f"    {'module':42s} {'GFD':>7s} {'challenge':>10s} {'control':>8s} "
          f"{'p GFD vs ctrl':>14s} {'p chal vs GFD':>14s}")
    for mn, r in sorted(B["modules"].items(),
                        key=lambda kv: -(kv[1].get("phi_gfd", -9) if kv[1]["defined"] else -9)):
        if not r["defined"]:
            print(f"    {mn:42s} undefined ({r['n_disrupted']} of {r['n_total']} disrupted)")
            continue
        print(f"    {mn:42s} {r['phi_gfd']:+7.2f} {r['phi_challenge']:+10.2f} "
              f"{r['phi_control']:+8.2f} {r['p_gfd_vs_control']:9.4f} "
              f"q={r.get('q_gfd_vs_control', float('nan')):.3f} "
              f"{r['p_challenge_vs_gfd']:14.4f}")
    if "control_drift_gfd" in B:
        print(f"    NOTE: the housekeeping control itself sits at "
              f"{B['control_drift_gfd']:+.2f} on the gluten-free diet, so a module's "
              "residual only counts above that")
    c = B["challenge_cos_with_disease_axis"]
    r0 = B["gfd_residual_cos_with_disease_axis"]
    print(f"    the gluten challenge moves along the disease axis (leave-one-patient-out): "
          f"median cos {c['median']:+.3f}")
    print(f"    what is left on the gluten-free diet lies along the same axis: median cos "
          f"{r0['median']:+.3f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    A, B = R["A"], R["B"]
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))

    a = ax[0]
    S = R["S_former"]
    defined = [mn for mn, r in A["modules"].items()
               if r["defined"] and not mn.startswith("housekeeping")]
    for k, mn in enumerate(defined):
        r = A["modules"][mn]
        c = f"C{k}"
        t = np.array([S[s]["quit_months"] for s in S])
        y = np.array([r["phi"][s] for s in S])
        a.plot(t / 12, y, "o", ms=3.5, alpha=0.35, color=c)
        z = np.polyfit(np.log2(t), y, 1)
        xs = np.linspace(np.log2(t.min()), np.log2(t.max()), 40)
        a.plot(2 ** xs / 12, np.polyval(z, xs), "-", lw=2, color=c,
               label=f"{mn.split(' (')[0][:26]}")
    a.set_ylim(-1.6, 2.4)        # a few inflammation points lie far outside
    a.set_xscale("log")
    a.axhline(0, color="k", lw=0.8)
    a.axhline(1, color="k", ls=":", lw=1)
    a.set_xlabel("years since quitting (log)")
    a.set_ylabel("position (1 = current smoker, 0 = never)")
    a.set_title("A  airway: what comes back after the cause is gone\n(points clipped to the plotted range)", fontsize=10)
    a.legend(fontsize=6.5, loc="upper right")

    a = ax[1]
    names, lqv, nvv = [], [], []
    for mn, r in sorted(A["modules"].items(),
                        key=lambda kv: (kv[1].get("phi_matched_longquit", 9)
                                        if kv[1]["defined"] else 9)):
        if not r["defined"]:
            continue
        names.append(mn.split(" (")[0][:26] +
                     (" *" if r.get("p_matched", 1) < 0.05 else ""))
        lqv.append(r["phi_matched_longquit"])
        nvv.append(r["phi_matched_never"])
    y = np.arange(len(names))
    a.barh(y + 0.2, lqv, 0.38, color="C3", label="quit 10+ years ago (n=7)")
    a.barh(y - 0.2, nvv, 0.38, color="C2", label="never smoked (n=5)")
    a.set_yticks(y, names, fontsize=8)
    a.axvline(0, color="k", lw=0.8)
    a.axvline(1, color="k", ls=":", lw=1)
    a.set_xlabel("position (1 = current smoker)")
    a.set_title("A  age-matched 40-60 window:\nthe only comparison age cannot explain",
                fontsize=10)
    a.legend(fontsize=7.5, loc="lower right")

    a = ax[2]
    ms = [mn for mn, r in B["modules"].items() if r["defined"]]
    xs = np.arange(len(ms))
    w = 0.27
    for k, (key, lab, c) in enumerate((("phi_control", "healthy control", "C2"),
                                       ("phi_gfd", "on gluten-free diet", "C0"),
                                       ("phi_challenge", "after gluten", "C3"))):
        a.bar(xs + (k - 1) * w, [B["modules"][m][key] for m in ms], w, label=lab, color=c)
    a.set_xticks(xs, [m.split(" (")[0][:22] for m in ms], rotation=30, ha="right",
                 fontsize=7.5)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("position (1 = challenged state)")
    a.set_title("B  coeliac: cause removed for years, then returned\n(the challenged state is 1 by construction)", fontsize=10)
    a.legend(fontsize=8)

    fig.tight_layout()
    fn = os.path.join(FIG, "d22_scar_after_cause_removed.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        DA = load_a()
        A = part_a(DA)
        S = DA["S"]
        R = {"A": A, "B": part_b(load_b()),
             "S_former": {s: {"quit_months": float(S.loc[s, "quit_months"]),
                              "age": float(S.loc[s, "age"])}
                          for s in S.index[S.group == "Former Smoker"]}}
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D22_FORCE")))
