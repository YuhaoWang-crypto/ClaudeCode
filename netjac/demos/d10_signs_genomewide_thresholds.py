"""
D10  The three D9 follow-ups.

  A  Kinase -> substrate SIGNS from CPTAC phosphoproteomes (8 indications,
     same cases as the D9 proteomes).
     - kinase activity = mean z of its known substrate sites (OmniPath
       enz_sub), leave-one-out;
     - phosphosites normalised to their own protein (removes abundance);
     - top-3 collective modes subtracted (the D9 "mean-field" lesson);
     - does the SIGNOR sign (up- vs down-regulates) of kinase -> kinase
       edges show up in activity correlations, when it did NOT show up in
       protein abundance?
     - 2-node kinase loops: are negative-feedback loops hidden here too?
  B  Directed Jacobian genome-wide: Replogle 2022 K562 genome-wide
     Perturb-seq (11,258 knockdowns x 8,248 genes).
     - noise calibrated on 585 non-targeting pseudo-bulk rows;
     - total-response network among strongly knocked-down genes;
     - replication in the independent essential-gene screen;
     - signed loops;
     - genome-wide FDT / anti-FDT homeostat atlas (control-cell covariance
       from the essential screen);
     - macro response modes.
  C  Beta-cell threshold distribution from measured data.
     - Per-cell Ca2+ recruitment (Jonkers & Henquin 2001, mouse): log-normal
       threshold fits for uncoupled single cells and coupled clusters.
     - Per-cell human exocytosis Patch-seq (Camunas-Soler 2020, MIT licence):
       amplitude / glucose-amplification distributions, ND vs T2D, and
       donor-level micro -> macro tests.
     - Both plugged into the D9 beta-population Topp model.

Run:  python3 -m netjac.demos.d10_signs_genomewide_thresholds [--part A,B,C]
"""
import argparse
import collections
import csv
import os
import pickle
import re
import subprocess
import time

import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, norm, spearmanr, wilcoxon

from netjac import DATA_CACHE, FIG
from netjac import statmech as sm

RES = os.path.join(DATA_CACHE, "d10_results.pkl")
PATCHSEQ = os.path.join(DATA_CACHE, "patchseq")


def _load_res():
    return pickle.load(open(RES, "rb")) if os.path.exists(RES) else {}


def _save_res(R):
    pickle.dump(R, open(RES, "wb"))


# ======================================================================= A
def _enz_sub():
    from netjac.data import omnipath as op
    fn = os.path.join(op.DIR, "enz_sub.tsv")
    if not os.path.exists(fn):
        subprocess.run(["curl", "-sS", "-o", fn,
                        "https://omnipathdb.org/enz_sub?genesymbols=1&organisms=9606"],
                       check=True)
    ES = collections.defaultdict(set)
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["modification"] == "phosphorylation":
            ES[r["enzyme_genesymbol"]].add(
                f"{r['substrate_genesymbol']}_{r['residue_type']}{r['residue_offset']}")
    return ES


def _signor_phospho():
    """SIGNOR protein->protein (de)phosphorylation edges with a sign:
    +1 if the upstream enzyme ACTIVATES the substrate."""
    from netjac.data import omnipath as op
    op._file("signor_human.tsv")
    out = []
    for ln in open(os.path.join(op.DIR, "signor_human.tsv")):
        f = ln.rstrip("\n").split("\t")
        if (len(f) > 10 and f[1] == "protein" and f[5] == "protein"
                and f[9] in ("phosphorylation", "dephosphorylation")
                and ("up-regulates" in f[8] or "down-regulates" in f[8])):
            out.append((f[0], f[4], 1 if "up-regulates" in f[8] else -1, f[9], f[10]))
    return out


def _prep_phospho(name):
    from netjac.data import cptac
    ph = cptac.load_phospho(name, min_present=0.5)
    pr = cptac.load(name)
    pi = {(c, bool(t)): i for i, (c, t) in enumerate(zip(pr["cases"], pr["tumour"]))}
    keys = list(zip(ph["cases"], map(bool, ph["tumour"])))
    if np.mean([k in pi for k in keys]) < 0.5:
        # HBV-HCC: proteome aliquots are 'T471'/'P472', phospho aliquots
        # '471'/'472' (same aliquot numbering); match on the number
        pi = {(a.lstrip("TP"), bool(t)): i for i, (a, t) in
              enumerate(zip(pr["aliquots"], pr["tumour"]))}
        keys = list(zip(ph["aliquots"], map(bool, ph["tumour"])))
    rows = [i for i, k in enumerate(keys) if k in pi]
    P = cptac.impute(pr["X"])[[pi[keys[i]] for i in rows]]
    X = cptac.impute(ph["X"])[rows]
    gi = {g: j for j, g in enumerate(pr["genes"])}
    has = np.array([g in gi for g in ph["site_gene"]])
    Xn = X[:, has] - P[:, [gi[g] for g in ph["site_gene"][has]]]
    return {"X": X, "site": ph["site"], "Xn": Xn, "site_n": ph["site"][has],
            "t": ph["tumour"][rows], "P": P, "pgenes": pr["genes"]}


def _activity_tests(Z, site, ES, SIG, rng, Pz=None, pgenes=None):
    ix = {s: i for i, s in enumerate(site)}
    sub = {k: [ix[s] for s in v if s in ix and not s.startswith(k + "_")]
           for k, v in ES.items()}
    sub = {k: v for k, v in sub.items() if len(v) >= 5}
    ks = list(sub)
    obs, nul = [], []
    for k in ks:
        for j in sub[k]:
            rest = [x for x in sub[k] if x != j]
            obs.append(np.corrcoef(Z[:, rest].mean(1), Z[:, j])[0, 1])
            k2 = ks[rng.integers(len(ks))]
            r2 = [x for x in sub[k2] if x != j]
            nul.append(np.corrcoef(Z[:, r2].mean(1), Z[:, j])[0, 1])
    pgi = {g: i for i, g in enumerate(pgenes)} if pgenes is not None else {}
    up, dn, up_ab, dn_ab, pairs = [], [], [], [], {}
    for a, b, s, mech, _ in SIG:
        if a in sub and b in sub and a != b:
            sa = [x for x in sub[a] if x not in sub[b]]
            sb = [x for x in sub[b] if x not in sub[a]]
            if len(sa) >= 3 and len(sb) >= 3:
                sgn = s * (1 if mech == "phosphorylation" else -1)
                r = np.corrcoef(Z[:, sa].mean(1), Z[:, sb].mean(1))[0, 1]
                (up if sgn > 0 else dn).append(r)
                pairs[(a, b)] = (sgn, r)
                if Pz is not None and a in pgi and b in pgi:
                    rab = np.corrcoef(Pz[:, pgi[a]], Pz[:, pgi[b]])[0, 1]
                    (up_ab if sgn > 0 else dn_ab).append(rab)
    # 2-node kinase loops among the testable pairs
    loops = {"positive": [], "negative": []}
    for (a, b), (s1, r) in pairs.items():
        if a < b and (b, a) in pairs:
            s2 = pairs[(b, a)][0]
            loops["positive" if s1 * s2 > 0 else "negative"].append(abs(r))
    return {"n_kinase": len(sub), "loo": np.array(obs), "loo_null": np.array(nul),
            "up": np.array(up), "dn": np.array(dn), "up_ab": np.array(up_ab),
            "dn_ab": np.array(dn_ab),
            "loops": {k: np.array(v) for k, v in loops.items()},
            "activity": {k: Z[:, v].mean(1) for k, v in sub.items()}}


def part_A():
    from netjac.data import cptac
    ES, SIG = _enz_sub(), _signor_phospho()
    rng = np.random.default_rng(0)
    out = {"n_enz": len(ES), "n_signor": len(SIG), "ind": {}}
    for name in cptac.PHOSPHO:
        t0 = time.time()
        D = _prep_phospho(name)
        res = {"n_sites": D["X"].shape[1], "n_sites_norm": D["Xn"].shape[1],
               "nT": int(D["t"].sum()), "nN": int((~D["t"]).sum())}
        for lab, m in (("T", D["t"]), ("N", ~D["t"])):
            Pz = sm.remove_modes(D["P"][m], 5)
            for mode in ("raw", "protnorm", "protnorm_k3"):
                if mode == "raw":
                    Z, s = sm.zscore(D["X"][m]), D["site"]
                elif mode == "protnorm":
                    Z, s = sm.zscore(D["Xn"][m]), D["site_n"]
                else:
                    Z, s = sm.remove_modes(D["Xn"][m], 3), D["site_n"]
                r = _activity_tests(Z, s, ES, SIG, rng, Pz, D["pgenes"])
                if not (lab == "T" and mode == "protnorm_k3"):
                    r.pop("activity")
                res[(lab, mode)] = r
        # tumour vs normal rewiring of kinase-activity couplings (k3 layer)
        aT = res[("T", "protnorm_k3")]["activity"]
        res[("T", "protnorm_k3")].pop("activity")
        out["ind"][name] = res
        print(f"  A {name}: {time.time() - t0:.0f}s", flush=True)
    return out


# ======================================================================= B
def _bulk_z(B):
    obs = B["obs"]
    n = obs["num_cells_filtered"].values.astype(float)
    n = np.where(np.isnan(n), obs["num_cells_unfiltered"].values.astype(float), n)
    Z = B["X"] * np.sqrt(n)[:, None]
    Z[~np.isfinite(Z)] = 0.0
    return Z


def part_B(thr=6.0, K=50):
    from netjac.data import perturbseq as ps, gprofiler as gp
    from netjac.demos.d9_micro_to_macro import FAMILIES
    out = {}
    Bg = ps.load_bulk("K562_gwps_bulk")
    X = Bg["X"].copy()
    out["n_nonfinite"] = int((~np.isfinite(X)).sum())
    X[~np.isfinite(X)] = 0.0
    ctrl = Bg["control"]
    Zs = _bulk_z(Bg)
    nt = Zs[ctrl]
    out["shape"] = X.shape
    out["n_ctrl_rows"] = int(ctrl.sum())
    out["nt_rate"] = {t: float(np.mean(np.abs(nt) > t)) for t in (4, 5, 6, 8)}
    out["gauss_rate"] = {t: float(2 * norm.sf(t)) for t in (4, 5, 6, 8)}
    gi = {e: i for i, e in enumerate(Bg["gene_id"])}
    meas = np.array([i for i in range(len(X)) if Bg["pert_id"][i] in gi and not ctrl[i]])
    own = np.array([gi[Bg["pert_id"][i]] for i in meas])
    g = Bg["pert_gene"][meas]
    selfz = X[meas, own]
    out["n_perturbed_measured"] = len(meas)
    best = {}
    for k in range(len(meas)):
        if g[k] not in best or selfz[k] < selfz[best[g[k]]]:
            best[g[k]] = k
    ks = np.array(sorted(best.values()))
    sel = ks[selfz[ks] < -0.5]
    S = len(sel)
    oS, gS = own[sel], g[sel]
    Z = Zs[meas[sel]][:, oS].T                         # readout i x knockdown j
    off = ~np.eye(S, dtype=bool)
    rate = out["nt_rate"][thr] if thr in out["nt_rate"] else float(np.mean(np.abs(nt) > thr))
    sig = (np.abs(Z) > thr) & off
    eff = -np.sign(Z)                                  # + : j activates i
    out["net"] = {"nodes": S, "edges": int(sig.sum()), "tests": int(off.sum()),
                  "expected_false": rate * off.sum(), "fdr": rate * off.sum() / sig.sum(),
                  "frac_activating": float(np.mean(eff[sig] > 0)),
                  "out_degree": sorted(zip(sig.sum(0), gS), reverse=True)[:15],
                  "in_degree": sorted(zip(sig.sum(1), gS), reverse=True)[:15]}
    mut = np.argwhere(np.triu(sig & sig.T, 1))
    kinds = {"mutual activation": [], "mutual inhibition": [], "negative": []}
    for i, j in mut:
        s = eff[i, j] * eff[j, i]
        if s < 0:
            act, rep = (gS[j], gS[i]) if eff[i, j] > 0 else (gS[i], gS[j])
            kinds["negative"].append((act, rep, float(min(abs(Z[i, j]), abs(Z[j, i])))))
        else:
            kinds["mutual activation" if eff[i, j] > 0 else "mutual inhibition"].append(
                (gS[i], gS[j], float(min(abs(Z[i, j]), abs(Z[j, i])))))
    out["loops"] = {k: sorted(v, key=lambda x: -x[2]) for k, v in kinds.items()}
    out["loops_chance"] = float(sig[off].mean() ** 2 * off.sum() / 2)
    # replication in the independent essential-gene screen
    E = ps.load_bulk("K562_essential_bulk")
    ZE = _bulk_z(E)
    egi = {e: i for i, e in enumerate(E["gene_id"])}
    erow = {}
    for i, pg in enumerate(E["pert_gene"]):
        if pg != "non-targeting":
            erow.setdefault(pg, []).append(i)
    ii, jj = np.where(sig)
    agree, resig = [], []
    for i, j in zip(ii, jj):
        gid_i = Bg["gene_id"][oS[i]]
        if gS[j] in erow and gid_i in egi:
            v = ZE[erow[gS[j]], egi[gid_i]].mean()
            agree.append(np.sign(v) == np.sign(Z[i, j]))
            resig.append(abs(v) > 3 and np.sign(v) == np.sign(Z[i, j]))
    out["replication"] = {"tested": len(agree), "sign": float(np.mean(agree)),
                          "resig": float(np.mean(resig))}
    # replication of loops (both arms re-tested)
    rep_loop = {"tested": 0, "both_same_sign": 0}
    for i, j in mut:
        a_ = [(i, j), (j, i)]
        ok = True
        vals = []
        for (p, q) in a_:
            gid_p = Bg["gene_id"][oS[p]]
            if gS[q] in erow and gid_p in egi:
                vals.append(np.sign(ZE[erow[gS[q]], egi[gid_p]].mean()) == np.sign(Z[p, q]))
            else:
                ok = False
        if ok:
            rep_loop["tested"] += 1
            rep_loop["both_same_sign"] += int(all(vals))
    out["loop_replication"] = rep_loop
    # macro response modes (genome-wide) vs essential screen
    strong = (Bg["obs"]["energy_test_p_value"].values < 0.01) & ~ctrl
    Rr = X[strong]
    s = np.linalg.svd(Rr - Rr.mean(0), compute_uv=False) ** 2
    Rn = sm.shuffle_columns(Rr[:3000], np.random.default_rng(0))
    s2 = np.linalg.svd(Rn - Rn.mean(0), compute_uv=False) ** 2 * (len(Rr) / 3000)
    out["modes"] = {"n_strong": int(strong.sum()), "PR": sm.participation_ratio(s),
                    "k_null": int((s > s2[0]).sum()), "top10": float(s[:10].sum() / s.sum()),
                    "lam": s[:40] / s.sum()}
    # genome-wide FDT / anti-FDT atlas (control cells from the essential screen)
    r = ps.reduce_singlecell()
    cgid = r["gene_id"]
    common = np.intersect1d(Bg["gene_id"], cgid)
    bi = {e: i for i, e in enumerate(Bg["gene_id"])}
    ci = {e: i for i, e in enumerate(cgid)}
    C = r["ctrl_cells"][:, [ci[e] for e in common]].astype(np.float32)
    C = (C - C.mean(0)) / C.std(0)
    Xc = X[:, [bi[e] for e in common]]
    ix = {e: i for i, e in enumerate(common)}
    m2 = np.array([i for i in range(len(Xc)) if Bg["pert_id"][i] in ix and not ctrl[i]])
    o2 = np.array([ix[Bg["pert_id"][i]] for i in m2])
    keep = Xc[m2, o2] < -0.3
    m2, o2 = m2[keep], o2[keep]
    n = len(m2)
    Ccol = (C.T @ C[:, o2]) / len(C)
    del C
    P = np.zeros((n, len(common)), np.float32)
    for k in range(n):
        c = Ccol[:, k].copy()
        c[o2[k]] = -np.inf
        P[k, np.argsort(-c)[:K]] = 1.0 / K
    Dm = Xc[m2].copy()
    Dm[np.arange(n), o2] = 0
    M = Dm @ P.T
    q = np.array([(M[k] < M[k, k]).mean() for k in range(n)])
    sh = np.random.default_rng(1).permutation(n)
    qn = np.array([(M[k] < M[k, sh[k]]).mean() for k in range(n)])
    g2 = Bg["pert_gene"][m2]
    anti, fdt = q > 0.98, q < 0.02
    fam = {}
    for k, pat in FAMILIES.items():
        mk = np.array([bool(re.match(pat, x)) for x in g2])
        a, f, b = int((anti & mk).sum()), int((fdt & mk).sum()), int(mk.sum())
        fam[k] = {"tested": b, "anti": a, "fdt": f,
                  "p_anti": fisher_exact([[a, b - a], [anti.sum() - a,
                                          (~mk).sum() - anti.sum() + a]], alternative="greater")[1],
                  "p_fdt": fisher_exact([[f, b - f], [fdt.sum() - f,
                                         (~mk).sum() - fdt.sum() + f]], alternative="greater")[1]}
    d9 = pickle.load(open(os.path.join(DATA_CACHE, "d9_results.pkl"), "rb"))["D"]
    qA, qB, ge = d9["fdt_q"]
    eq = dict(zip(ge, (qA + qB) / 2))
    pr = np.array([(q[k], eq[g2[k]]) for k in range(n) if g2[k] in eq])
    ess = set(ge)
    anti_genes = sorted(set(g2[anti]))
    bg = sorted(set(g2))
    out["fdt"] = {"n": n, "counts": {t: {"fdt": int(np.sum(q < t)), "fdt_null": int(np.sum(qn < t)),
                                         "anti": int(np.sum(q > 1 - t)),
                                         "anti_null": int(np.sum(qn > 1 - t))}
                                     for t in (0.02, 0.01)},
                  "families": fam, "rep_n": len(pr),
                  "rep_rho": spearmanr(pr[:, 0], pr[:, 1]),
                  "anti_genes": anti_genes,
                  "anti_new": [x for x in anti_genes if x not in ess],
                  "enrich_anti": [x for x in gp.enrich(anti_genes, bg) if x["term_size"] < 500][:10],
                  "enrich_fdt": [x for x in gp.enrich(sorted(set(g2[fdt])), bg)
                                 if x["term_size"] < 500][:6],
                  "q": q, "qn": qn, "g": g2}
    return out


# ======================================================================= C
# Jonkers & Henquin 2001 Diabetes 50:540 (doi:10.2337/diabetes.50.3.540),
# mouse islet cells, cytosolic Ca2+ responses, values as stated in the
# abstract: single cells: thresholds between 6 and 10 mM, max ~65 %
# responsive at 10 mM, 13 % unresponsive even to glucose + tolbutamide.
# Clusters: responding clusters ~20 % at 6 mM -> ~95 % at 10 mM; within
# responsive clusters 75 % of cells active at 6 mM, 95-100 % at 8-10 mM.
JH2001 = {
    "single cells": [(6.0, 0.05), (10.0, 0.65 / 0.87)],   # among excitable cells
    "clusters": [(6.0, 0.20), (10.0, 0.95)],
    "cells within clusters": [(6.0, 0.75), (8.0, 0.975)],
}


def lognormal_from_two(points):
    """Exact log-normal threshold CDF through two (glucose, fraction) points."""
    (g1, f1), (g2, f2) = points
    z1, z2 = norm.ppf(f1), norm.ppf(f2)
    w = (np.log(g2) - np.log(g1)) / (z2 - z1)
    mu = np.log(g1) - z1 * w
    return np.exp(mu), w


def _patchseq():
    fn = os.path.join(PATCHSEQ, "patchclamp_wcryo_human.metadata.tab")
    if not os.path.exists(fn):
        os.makedirs(PATCHSEQ, exist_ok=True)
        subprocess.run(["curl", "-sSL", "-o", fn,
                        "https://raw.githubusercontent.com/jcamunas/patchseq/master/data/"
                        "patchclamp_wcryo_human.metadata.tab"], check=True)
    import pandas as pd
    m = pd.read_csv(fn, sep="\t", index_col=0)
    b = m[(m.cell_type == "beta") & m["NormalizedTotalCapacitance_fF/pF"].notna()].copy()
    b["x"] = b["NormalizedTotalCapacitance_fF/pF"]
    grp = {"healthy": "ND", "T2D": "T2D", "T2Dreversed": "T2D",
           "prediabetic": "pre", "elevated HbA1c": "pre"}
    b["g"] = b.DiabetesStatus.map(grp)
    return b[b.g.notna()]


def part_C(n_boot=2000):
    from netjac.demos.d9_micro_to_macro import (BetaPopulation, TOPP, fast_glucose,
                                                 fate, fit_population, topp_state)
    out = {"fits": {}}
    for k, pts in JH2001.items():
        med, w = lognormal_from_two(pts)
        out["fits"][k] = {"median_mM": med, "w": w, "points": pts}
    # model consequences: same macro curve shape, different micro widths
    m0, w0, *_ = fit_population()
    beta0 = topp_state(BetaPopulation(m0, w0))["beta_h"]
    rows = {}
    for lab, w in [("Topp-implied (pure recruitment)", w0)] + [
            (k, v["w"]) for k, v in out["fits"].items()]:
        P = BetaPopulation(m0, w)
        st = topp_state(P)
        from scipy.optimize import brentq
        crit = brentq(lambda x: fast_glucose(BetaPopulation(x * m0, w), beta0) - st["G_s"],
                      1.0, 80.0)
        rows[lab] = {"w": w, "zeta": st["zeta"], "tau": st["tau_crit_min"],
                     "gain": st["gain"], "crit_shift": crit}
    out["width_consequences"] = rows
    # recruitment-only macro curve vs Topp curve: how graded can it be?
    Gs = np.linspace(40, 600, 200)
    target = TOPP["sigma"] * Gs ** 2 / (TOPP["alpha"] + Gs ** 2)
    curves = {}
    for k, v in out["fits"].items():
        P = BetaPopulation(v["median_mM"] * 18.0, v["w"])
        curves[k] = P.sigma(Gs)
    out["curves"] = {"G": Gs, "topp": target, **curves}
    # 10-90 % dynamic range of recruitment vs of the Topp macro curve
    def rng1090(F):
        F = F / F.max()
        return Gs[np.argmax(F >= 0.9)] / Gs[np.argmax(F >= 0.1)]
    out["dynamic_range"] = {k: rng1090(c) for k, c in curves.items()}
    out["dynamic_range"]["Topp macro curve"] = rng1090(target)
    # Patch-seq human per-cell exocytosis
    b = _patchseq()
    ps_ = {}
    rng = np.random.default_rng(0)
    for g in ("ND", "T2D"):
        s = b[b.g == g]
        med = {G: float(s[s.Glucose_mM == G].x.median()) for G in (1.0, 5.0, 10.0)}
        n = {G: int((s.Glucose_mM == G).sum()) for G in (1.0, 5.0, 10.0)}
        donors = s.DonorID.unique()
        boots = []
        for _ in range(n_boot):
            dd = rng.choice(donors, len(donors), replace=True)
            ss = __import__("pandas").concat([s[s.DonorID == d] for d in dd])
            hi = ss[ss.Glucose_mM == 10.0].x
            lo = ss[ss.Glucose_mM == 1.0].x
            if len(hi) > 3 and len(lo) > 3:
                boots.append(np.log2((hi.median() + 1) / (max(lo.median(), 0) + 1)))
        x10 = s[s.Glucose_mM == 10.0].x.values
        q25, q75 = np.percentile(x10, [25, 75])
        ps_[g] = {"median": med, "n": n, "donors": len(donors),
                  "amp_log2": float(np.log2((med[10.0] + 1) / (med[1.0] + 1))),
                  "amp_ci": np.percentile(boots, [2.5, 97.5]),
                  "iqr10": (q25, q75), "x10": x10,
                  "frac_neg10": float(np.mean(x10 <= 0))}
    ps_["MWU_10mM_ND_vs_T2D"] = mannwhitneyu(ps_["ND"]["x10"], ps_["T2D"]["x10"]).pvalue
    ps_["MWU_1mM_ND_vs_T2D"] = mannwhitneyu(b[(b.g == "ND") & (b.Glucose_mM == 1)].x,
                                            b[(b.g == "T2D") & (b.Glucose_mM == 1)].x).pvalue
    # recruitment proxy: fraction above the ND 1-mM 90th percentile
    cut = np.percentile(b[(b.g == "ND") & (b.Glucose_mM == 1)].x, 90)
    ps_["active_cut"] = float(cut)
    ps_["active_frac"] = {g: {G: float((b[(b.g == g) & (b.Glucose_mM == G)].x > cut).mean())
                              for G in (1.0, 5.0, 10.0)} for g in ("ND", "T2D")}
    # donor-level micro -> macro
    def stim(v):
        if not isinstance(v, str):
            return np.nan
        mm = re.match(r"\s*([\d.]+)", v)
        return float(mm.group(1)) if mm else np.nan
    rows = []
    for d, s in b.groupby("DonorID"):
        hi = s[s.Glucose_mM == 10.0].x
        if len(hi) >= 3:
            rows.append({"donor": d, "g": s.g.iloc[0], "n10": len(hi),
                         "med10": hi.median(), "frac10": float((hi > cut).mean()),
                         "hba1c": float(__import__("pandas").to_numeric(
                             s["Glycated hemoglobin (HbA1c)"], errors="coerce").iloc[0]),
                         "SI16": stim(s["16.7 mM to 1 mM stimulation index"].iloc[0]),
                         "SI10": stim(s["10 mM to 1 mM stimulation index"].iloc[0])})
    import pandas as pd
    Dd = pd.DataFrame(rows)
    tests = {}
    for x in ("med10", "frac10"):
        for y in ("SI16", "SI10", "hba1c"):
            ok = Dd[[x, y]].dropna()
            r_ = spearmanr(ok[x], ok[y])
            tests[(x, y)] = (float(r_[0]), float(r_[1]), len(ok))
    ps_["donor_tests"] = tests
    ps_["donors"] = Dd
    out["patchseq"] = ps_
    # plug measured amplification into the model: T2D amplification loss
    # rescales sigma_max (independent of thresholds) -> what does Topp do?
    amp = 2 ** (ps_["T2D"]["amp_log2"] - ps_["ND"]["amp_log2"])
    P = BetaPopulation(m0, w0, sigma_max=TOPP["sigma"] * amp)
    st = topp_state(P)
    out["amp_model"] = {"rel_amp": amp, "G_fast": fast_glucose(P, beta0),
                        "beta_needed": st["beta_h"] / beta0, "margin": st["margin"],
                        "fate": fate(P, beta0)}
    return out


# ================================================================== report
def report(R):
    if "A" in R:
        A = R["A"]
        print(f"\nA  Phosphoproteome kinase signs ({A['n_enz']} enzymes with substrates, "
              f"{A['n_signor']} signed SIGNOR phospho edges)")
        pooled = {"up": [], "dn": [], "up_ab": [], "dn_ab": [], "pos": [], "neg": []}
        for nm, r in A["ind"].items():
            s = [f"{nm:6s} sites {r['n_sites']} (protein-matched {r['n_sites_norm']}), T {r['nT']} N {r['nN']}"]
            for mode in ("raw", "protnorm", "protnorm_k3"):
                o = r[("T", mode)]
                p = mannwhitneyu(o["up"], o["dn"], alternative="greater").pvalue
                s.append(f"{mode}: LOO {np.median(o['loo']):.2f}/null {np.median(o['loo_null']):.2f}; "
                         f"K->K up {o['up'].mean():+.3f} dn {o['dn'].mean():+.3f} p={p:.0e}")
            o = r[("T", "protnorm_k3")]
            pab = mannwhitneyu(o["up_ab"], o["dn_ab"], alternative="greater").pvalue
            s.append(f"abundance layer (res5): up {o['up_ab'].mean():+.3f} dn {o['dn_ab'].mean():+.3f} p={pab:.0e}")
            lp = o["loops"]
            s.append(f"kinase 2-loops |r|: pos {np.median(lp['positive']) if len(lp['positive']) else np.nan:.3f} "
                     f"(n={len(lp['positive'])}) neg {np.median(lp['negative']) if len(lp['negative']) else np.nan:.3f} "
                     f"(n={len(lp['negative'])})")
            print("   " + " | ".join(s))
            pooled["up"].append(o["up"].mean() - o["dn"].mean())
            pooled["up_ab"].append(o["up_ab"].mean() - o["dn_ab"].mean())
            if len(lp["positive"]) and len(lp["negative"]):
                pooled["pos"].append(np.median(lp["positive"]) - np.median(lp["negative"]))
        d = np.array(pooled["up"])
        da = np.array(pooled["up_ab"])
        print(f"   activity layer: up-dn > 0 in {np.sum(d > 0)}/8 (Wilcoxon p={wilcoxon(d, alternative='greater').pvalue:.3f}); "
              f"abundance layer: {np.sum(da > 0)}/8 (p={wilcoxon(da, alternative='greater').pvalue:.3f})")
        if pooled["pos"]:
            dp = np.array(pooled["pos"])
            print(f"   kinase loops: |r| pos > neg in {np.sum(dp > 0)}/{len(dp)}")
    if "B" in R:
        B = R["B"]
        print(f"\nB  Genome-wide Perturb-seq {B['shape'][0]} knockdowns x {B['shape'][1]} genes; "
              f"{B['n_ctrl_rows']} non-targeting rows")
        print("   empirical vs Gaussian per-entry rate: " + ", ".join(
            f"|z|>{t}: {B['nt_rate'][t]:.1e} vs {B['gauss_rate'][t]:.1e}" for t in (4, 5, 6, 8)))
        n_ = B["net"]
        print(f"   network: {n_['nodes']} nodes, {n_['edges']} edges (expected false {n_['expected_false']:.0f}, "
              f"FDR {n_['fdr']:.3f}); activating {n_['frac_activating']:.2f}")
        print("   top out-degree (regulators): " + ", ".join(f"{g}({d})" for d, g in n_["out_degree"][:12]))
        print("   top in-degree (sensors): " + ", ".join(f"{g}({d})" for d, g in n_["in_degree"][:12]))
        rp = B["replication"]
        print(f"   replication in essential screen: {rp['tested']} edges, sign {rp['sign']:.3f}, "
              f"re-significant {rp['resig']:.3f}")
        L = B["loops"]
        print(f"   loops: " + ", ".join(f"{k} {len(v)}" for k, v in L.items()) +
              f" (chance {B['loops_chance']:.0f}); both arms replicated in "
              f"{B['loop_replication']['both_same_sign']}/{B['loop_replication']['tested']}")
        print("   strongest negative loops (activator -> repressor ⊣ activator): " + "; ".join(
            f"{a}->{b}" for a, b, _ in L["negative"][:15]))
        m = B["modes"]
        print(f"   response modes: {m['n_strong']} strong knockdowns, PR {m['PR']:.1f}, "
              f"{m['k_null']} modes > null, top-10 {m['top10']:.0%}")
        f = B["fdt"]
        for t, c in f["counts"].items():
            print(f"   FDT q<{t}: FDT {c['fdt']} (null {c['fdt_null']}), anti {c['anti']} (null {c['anti_null']})")
        for k, v in f["families"].items():
            print(f"   {k:20s} tested {v['tested']:3d}: anti {v['anti']:2d} (p={v['p_anti']:.0e}), "
                  f"FDT {v['fdt']:2d} (p={v['p_fdt']:.0e})")
        print(f"   replication vs essential-screen ranks: n={f['rep_n']}, Spearman "
              f"{f['rep_rho'][0]:.2f} (p={f['rep_rho'][1]:.0e})")
        print(f"   anti-FDT genes {len(f['anti_genes'])}, new {len(f['anti_new'])}: " + ", ".join(f["anti_new"]))
        print("   anti-FDT enrichment: " + "; ".join(
            f"{x['name'][:40]} ({x['source']}, p={x['p']:.0e}, k={x['k']})" for x in f["enrich_anti"][:6]))
    if "C" in R:
        C = R["C"]
        print("\nC  Beta-cell threshold distributions from data")
        for k, v in C["fits"].items():
            print(f"   {k:22s} log-normal median {v['median_mM']:.1f} mM ({v['median_mM'] * 18:.0f} mg/dl), "
                  f"width w = {v['w']:.3f}")
        for k, v in C["width_consequences"].items():
            print(f"   {k:32s} w {v['w']:.2f}: gain {v['gain']:.0f}, zeta {v['zeta']:.2f}, "
                  f"tau* {v['tau']:.0f} min, tipping threshold shift x{v['crit_shift']:.1f}")
        print("   10-90 % dynamic range (G90/G10): " + ", ".join(
            f"{k} {v:.2f}" for k, v in C["dynamic_range"].items()))
        p = C["patchseq"]
        for g in ("ND", "T2D"):
            e = p[g]
            print(f"   Patch-seq {g}: {e['donors']} donors; median exocytosis (fF/pF) 1/5/10 mM "
                  f"{e['median'][1.0]:.1f}/{e['median'][5.0]:.1f}/{e['median'][10.0]:.1f} "
                  f"(n {e['n'][1.0]}/{e['n'][5.0]}/{e['n'][10.0]}); amplification log2 {e['amp_log2']:.2f} "
                  f"[donor-bootstrap {e['amp_ci'][0]:.2f}, {e['amp_ci'][1]:.2f}]; IQR at 10 mM "
                  f"{e['iqr10'][0]:.1f}-{e['iqr10'][1]:.1f}; active (> ND 1 mM P90) " +
                  "/".join(f"{p['active_frac'][g][G]:.2f}" for G in (1.0, 5.0, 10.0)))
        print(f"   ND vs T2D at 10 mM MWU p={p['MWU_10mM_ND_vs_T2D']:.2f}; at 1 mM p={p['MWU_1mM_ND_vs_T2D']:.2f}")
        for (x, y), (r_, pv, n_) in p["donor_tests"].items():
            print(f"   donor {x} vs {y}: rho {r_:+.2f} p={pv:.2f} n={n_}")
        a = C["amp_model"]
        print(f"   T2D/ND amplification {a['rel_amp']:.2f} -> Topp: acute G {a['G_fast']:.0f}, beta needed "
              f"x{a['beta_needed']:.2f}, fate {a['fate']['kind']} (G {a['fate']['G']:.0f})")


# ================================================================= figures
def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(2, 4, figsize=(20, 9))
    if "A" in R:
        A = R["A"]
        names = list(A["ind"])
        x = np.arange(len(names))
        a = ax[0, 0]
        for j, (mode, c) in enumerate((("raw", "0.6"), ("protnorm", "C7"), ("protnorm_k3", "C3"))):
            v = [np.median(A["ind"][n][("T", mode)]["loo"]) - np.median(A["ind"][n][("T", mode)]["loo_null"])
                 for n in names]
            a.bar(x + (j - 1) * 0.27, v, 0.27, color=c, label=mode)
        a.set_xticks(x, names, rotation=30)
        a.set_ylabel("LOO r(kinase activity, its site) − null")
        a.set_title("A  kinase → substrate visibility\n(signal above a random kinase)")
        a.legend(fontsize=8)
        a = ax[0, 1]
        for j, (lay, ku, kd, c) in enumerate((("activity (k3)", "up", "dn", "C3"),
                                               ("abundance (res5)", "up_ab", "dn_ab", "C0"))):
            v = [A["ind"][n][("T", "protnorm_k3")][ku].mean() - A["ind"][n][("T", "protnorm_k3")][kd].mean()
                 for n in names]
            a.bar(x + (j - 0.5) * 0.4, v, 0.4, color=c, label=lay)
        a.axhline(0, color="k", lw=0.8)
        a.set_xticks(x, names, rotation=30)
        a.set_ylabel("mean r(up-regulating) − mean r(down-regulating)")
        a.set_title("A  SIGNOR kinase→kinase sign:\nvisible in activity, not in abundance")
        a.legend(fontsize=8)
    if "B" in R:
        B = R["B"]
        a = ax[0, 2]
        f = B["fdt"]["families"]
        ks = list(f)
        y = np.arange(len(ks))
        a.barh(y - 0.2, [f[k]["anti"] / max(f[k]["tested"], 1) for k in ks], 0.4, color="C3",
               label="anti-FDT (homeostat)")
        a.barh(y + 0.2, [f[k]["fdt"] / max(f[k]["tested"], 1) for k in ks], 0.4, color="C0",
               label="FDT-consistent")
        a.set_yticks(y, ks, fontsize=8)
        a.set_xlabel("fraction of family knockdowns")
        a.set_title(f"B  genome-wide homeostat atlas\n({B['fdt']['n']} knockdowns)")
        a.legend(fontsize=8)
        a = ax[0, 3]
        L = B["loops"]
        kinds = list(L)
        a.bar(range(len(kinds)), [len(L[k]) for k in kinds], color=["C2", "C1", "C0"])
        a.axhline(B["loops_chance"] / 3, color="k", ls="--", lw=1, label="chance (per class)")
        a.set_xticks(range(len(kinds)), kinds, fontsize=8)
        a.set_ylabel("mutual pairs")
        n_ = B["net"]
        a.set_title(f"B  {n_['nodes']} nodes, {n_['edges']} edges (FDR {n_['fdr']:.1%})\n"
                    f"replicated sign {B['replication']['sign']:.0%}")
        a.legend(fontsize=8)
    if "C" in R:
        C = R["C"]
        a = ax[1, 0]
        cv = C["curves"]
        a.plot(cv["G"] / 18, cv["topp"] / cv["topp"].max(), "k", lw=3, label="Topp macro curve")
        for k, c in zip(C["fits"], ("C0", "C1", "C2")):
            a.plot(cv["G"] / 18, cv[k] / cv[k].max(), color=c, lw=2,
                   label=f"recruitment, {k} (w={C['fits'][k]['w']:.2f})")
            for gx, fx in C["fits"][k]["points"]:
                a.plot(gx, fx, "o", color=c)
        a.set_xlim(0, 25)
        a.set_xlabel("glucose (mM)")
        a.set_ylabel("fraction of maximum")
        a.set_title("C  measured recruitment is a narrow switch;\nthe macro curve is graded")
        a.legend(fontsize=7)
        a = ax[1, 1]
        wc = C["width_consequences"]
        ks = list(wc)
        a.bar(range(len(ks)), [min(v["tau"], 2000) for v in wc.values()], color="C3")
        a.set_yscale("log")
        a.set_xticks(range(len(ks)), [k.replace(" (pure recruitment)", "\n(pure recr.)") for k in ks],
                  fontsize=7, rotation=15)
        a.set_ylabel("critical insulin-action delay (min)")
        a.set_title("C  if recruitment alone set the gain:\nmeasured widths → oscillation-prone")
        p = C["patchseq"]
        a = ax[1, 2]
        for g, c in (("ND", "C0"), ("T2D", "C3")):
            a.plot([1, 5, 10], [p[g]["median"][G] for G in (1.0, 5.0, 10.0)], "o-", color=c,
                   lw=2, label=f"{g} ({p[g]['donors']} donors)")
        a.set_xlabel("glucose in dish (mM)")
        a.set_ylabel("median exocytosis (fF/pF)")
        a.set_title("C  human Patch-seq: per-cell amplification\n(voltage-clamped: not a trigger threshold)")
        a.legend(fontsize=8)
        a = ax[1, 3]
        Dd = p["donors"]
        for g, c in (("ND", "C0"), ("T2D", "C3"), ("pre", "C1")):
            s = Dd[Dd.g == g]
            a.scatter(s["med10"], s["SI16"], color=c, label=g)
        r_, pv, n_ = p["donor_tests"][("med10", "SI16")]
        a.set_xlabel("donor median β-cell exocytosis at 10 mM (micro)")
        a.set_ylabel("islet perifusion stimulation index 16.7/1 mM (macro)")
        a.set_title(f"C  micro → macro across donors\nρ = {r_:+.2f}, p = {pv:.2f}, n = {n_}")
        a.legend(fontsize=8)
    fig.tight_layout()
    p = os.path.join(FIG, "d10_signs_genomewide_thresholds.png")
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
