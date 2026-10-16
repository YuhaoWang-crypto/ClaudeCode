"""
D24  Co-occurrence, not the mean: does CRISPR perturbation move the WITHIN-CELL
covariance of antigen presentation and adhesion, and does that matter?

The hypothesis (scenario 7 of the indication catalogue). A T cell kills a
tumour cell only if that SAME cell presents antigen AND engages adhesion /
co-stimulation. In a Kuznetsov-type killing term the population kill rate is

    kappa = k0 * ( E[A] E[H] + Cov(A, H) )

so a perturbation that raises the co-expression of MHC-I and CD58 in the same
cells, without raising either mean, still raises killing. The claim to test is
that Cov carries information the means do not.

Data: Frangieh / Izar 2021 Perturb-CITE-seq (patient-derived melanoma,
via scPerturb). 218,331 cells x 24 surface proteins with single-cell CRISPR
knockouts, in three conditions: Control, IFN-gamma, and co-culture with
autologous tumour-infiltrating lymphocytes. The panel carries both members of
the pair (HLA_A, CD58) and FOUR isotype controls, which give a negative control
pair for every covariance statistic computed here.

Order of work (controls before claims):
  Q0  TECHNICAL. Two ADTs in the same cell share that cell's total count, which
      manufactures positive covariance. Compare raw vs CLR-normalised, and show
      what happens to the isotype-isotype pair, which cannot co-vary biologically.
  Q1  IS IT ESTIMABLE? Split each perturbation's cells in half and ask whether
      Cov agrees between halves. If not, the question is unanswerable at this
      ADT depth and everything below is noise.
  Q2  IS IT REDUNDANT? Regress Cov on the two means. A high R^2 means the
      "co-occurrence" coordinate adds nothing.
  Q3  THE PREMISE. Are there perturbations that move Cov while leaving both
      means where they were? Name them, with a cell bootstrap.
  Q4  DOES IT MATTER? Per knockout, a selection score from the cell-fraction
      shift between co-culture and control (the resistance read-out of this
      screen). Test whether Cov measured BEFORE T-cell pressure adds to the
      prediction of that score beyond the two means. This is the test that can
      fail.
  Q5  HOW MUCH? Variance decomposition of kappa across knockouts: how much of
      the spread comes from Cov versus from E[A]E[H]. No invented bifurcation
      threshold - the scenario's Cov_c = +0.160 and kappa_c = 1.641 have no
      source in this repository and are not used.

Run:  python3 -m netjac.demos.d24_cooccurrence_sensitisation
"""
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG

RES = os.path.join(DATA_CACHE, "d24_results.pkl")
H5 = os.path.join(DATA_CACHE, "frangieh", "protein.h5ad")
URL = ("https://zenodo.org/api/records/13350497/files/"
       "FrangiehIzar2021_protein.h5ad/content")
A_PROT, H_PROT = "HLA_A", "CD58"
ISO = ["Rat_IgG2a", "Mouse_IgG1", "Mouse_IgG2a", "Mouse_IgG2b"]
MIN_CELLS = 100


# ------------------------------------------------------------------ loading
def load(min_cells=MIN_CELLS):
    import anndata as ad
    import scipy.sparse as sp
    if not os.path.exists(H5):
        os.makedirs(os.path.dirname(H5), exist_ok=True)
        import subprocess
        subprocess.run(["curl", "-sSL", "--max-time", "900", "-o", H5, URL], check=True)
    A = ad.read_h5ad(H5)
    X = A.X.toarray() if sp.issparse(A.X) else np.asarray(A.X, dtype=float)
    obs = A.obs
    keep = np.asarray(((obs.nperts == 1) | (obs.perturbation == "control")).values).copy()
    tot = X.sum(1)
    lo, hi = np.quantile(tot[keep], [0.01, 0.99])
    keep &= (tot >= lo) & (tot <= hi)
    X, obs = X[keep], obs[keep]
    # CLR across the panel within each cell: removes the per-cell total, which
    # is what manufactures covariance between any two ADTs.
    L = np.log1p(X)
    CLR = L - L.mean(1, keepdims=True)
    return {"raw": np.log1p(X), "clr": CLR, "proteins": list(A.var_names),
            "pert": obs.perturbation.astype(str).values,
            "cond": obs.perturbation_2.astype(str).values,
            "total": np.log1p(tot[keep]), "min_cells": min_cells}


# ------------------------------------------------------------------ helpers
def _cov(a, b):
    a = a - a.mean()
    b = b - b.mean()
    return float(a @ b / max(len(a) - 1, 1))


def _groups(D, cond, min_cells=None):
    m = D["cond"] == cond
    p = D["pert"][m]
    out = {}
    for u in np.unique(p):
        idx = np.where(m)[0][p == u]
        if len(idx) >= (min_cells or D["min_cells"]):
            out[u] = idx
    return out


def _stats(M, i, j, idx):
    return {"mean_a": float(M[idx, i].mean()), "mean_h": float(M[idx, j].mean()),
            "cov": _cov(M[idx, i], M[idx, j]), "n": len(idx)}


# ------------------------------------------------------------------ analysis
def run(seed=0):
    from scipy import stats
    D = load()
    P = D["proteins"]
    ia, ih = P.index(A_PROT), P.index(H_PROT)
    i1, i2 = P.index(ISO[1]), P.index(ISO[2])
    out = {"n_cells": len(D["pert"]), "conditions": sorted(set(D["cond"]))}

    # ---- Q0 the depth artefact
    q0 = {}
    for name, M in (("raw log1p", D["raw"]), ("CLR", D["clr"])):
        q0[name] = {
            "cov_A_H": _cov(M[:, ia], M[:, ih]),
            "cov_isotype_pair": _cov(M[:, i1], M[:, i2]),
            "cov_A_isotype": _cov(M[:, ia], M[:, i1]),
            "r_A_with_total": float(np.corrcoef(M[:, ia], D["total"])[0, 1]),
            "r_H_with_total": float(np.corrcoef(M[:, ih], D["total"])[0, 1])}
    out["q0"] = q0

    M = D["clr"]
    g = {c: _groups(D, c) for c in out["conditions"]}
    out["n_perturbations"] = {c: len(v) for c, v in g.items()}

    # ---- Q1 split-half reproducibility of Cov, with the isotype pair as control
    rng = np.random.default_rng(seed)
    q1 = {}
    for c, gg in g.items():
        pairs = {"HLA-A x CD58": (ia, ih), "isotype x isotype": (i1, i2)}
        row = {}
        for lab, (x, y) in pairs.items():
            a, b, am, bm = [], [], [], []
            for u, idx in gg.items():
                k = rng.permutation(len(idx))
                h1, h2 = idx[k[: len(idx) // 2]], idx[k[len(idx) // 2:]]
                a.append(_cov(M[h1, x], M[h1, y]))
                b.append(_cov(M[h2, x], M[h2, y]))
                am.append(M[h1, x].mean())
                bm.append(M[h2, x].mean())
            row[lab] = {"r_cov": float(np.corrcoef(a, b)[0, 1]),
                        "r_mean_for_scale": float(np.corrcoef(am, bm)[0, 1]),
                        "n_perts": len(a)}
        q1[c] = row
    out["q1"] = q1

    # ---- Q1b does the ADT panel even see its own knockouts, and how many
    # cells would a usable covariance estimate need?
    panel_ko = {"HLA-A": A_PROT, "CD58": H_PROT, "CD274": "CD274", "CD59": "CD59",
                "CD47": "CD47", "CD44": "CD44", "CD9": "CD9", "B2M": A_PROT}
    seen = {}
    for c, gg in g.items():
        if "control" not in gg:
            continue
        ctrl = gg["control"]
        row = {}
        for ko, prot in panel_ko.items():
            if ko not in gg or prot not in P:
                continue
            k = P.index(prot)
            x, y = M[gg[ko], k], M[ctrl, k]
            sd = np.sqrt(y.var() / len(x) + y.var() / len(y))
            row[f"{ko} KO -> {prot}"] = {"delta": float(x.mean() - y.mean()),
                                         "z": float((x.mean() - y.mean()) / max(sd, 1e-9)),
                                         "n": len(x)}
        seen[c] = row
    out["q1b_panel_sees_ko"] = seen

    # bootstrap standard error of Cov as a function of cell count, on the
    # control cells, against the between-knockout spread of dCov
    ctrl_idx = {c: gg["control"] for c, gg in g.items() if "control" in gg}
    se = {}
    rngb = np.random.default_rng(seed + 1)
    for c, idx in ctrl_idx.items():
        row = {}
        for n in (100, 300, 1000, 3000, 10000):
            if n > len(idx):
                continue
            vals = [_cov(M[s, ia], M[s, ih]) for s in
                    (rngb.choice(idx, n, replace=False) for _ in range(200))]
            row[n] = float(np.std(vals))
        se[c] = row
    out["q1b_se_vs_n"] = se

    # ---- per-perturbation table in each condition
    tab = {}
    for c, gg in g.items():
        rows = {u: _stats(M, ia, ih, idx) for u, idx in gg.items()}
        base = rows.get("control")
        if base is None:
            continue
        for u, r in rows.items():
            r["d_mean_a"] = r["mean_a"] - base["mean_a"]
            r["d_mean_h"] = r["mean_h"] - base["mean_h"]
            r["d_cov"] = r["cov"] - base["cov"]
        tab[c] = rows
    out["table"] = tab

    # ---- Q2 is Cov redundant given the means?
    q2 = {}
    for c, rows in tab.items():
        ks = [u for u in rows if u != "control"]
        y = np.array([rows[u]["d_cov"] for u in ks])
        Xd = np.column_stack([np.ones(len(ks)),
                              [rows[u]["d_mean_a"] for u in ks],
                              [rows[u]["d_mean_h"] for u in ks]])
        beta, *_ = np.linalg.lstsq(Xd, y, rcond=None)
        r2 = 1 - np.sum((y - Xd @ beta) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12)
        q2[c] = {"r2_cov_from_means": float(r2), "n": len(ks),
                 "sd_d_cov": float(np.std(y))}
    out["q2"] = q2

    # ---- Q3 the premise: Cov moves while the means do not
    q3 = {}
    for c, rows in tab.items():
        ks = [u for u in rows if u != "control"]
        sa = np.std([rows[u]["d_mean_a"] for u in ks])
        sh = np.std([rows[u]["d_mean_h"] for u in ks])
        sc = np.std([rows[u]["d_cov"] for u in ks])
        quiet = [u for u in ks
                 if abs(rows[u]["d_mean_a"]) < 0.5 * sa
                 and abs(rows[u]["d_mean_h"]) < 0.5 * sh
                 and abs(rows[u]["d_cov"]) > 1.5 * sc]
        quiet.sort(key=lambda u: -abs(rows[u]["d_cov"]))
        q3[c] = {"n_quiet_movers": len(quiet), "sd": (float(sa), float(sh), float(sc)),
                 "top": [(u, rows[u]["d_cov"], rows[u]["d_mean_a"],
                          rows[u]["d_mean_h"], rows[u]["n"]) for u in quiet[:10]]}
    out["q3"] = q3

    # ---- Q4 does it predict the screen's resistance read-out?
    # selection score = log2 (share of this KO in co-culture / share in control)
    share = {}
    for c in ("Control", "Co-culture", "IFNγ"):
        m = D["cond"] == c
        u, n = np.unique(D["pert"][m], return_counts=True)
        share[c] = dict(zip(u, n / n.sum()))
    common = [u for u in tab.get("Control", {})
              if u != "control" and u in share["Co-culture"] and u in share["Control"]]
    sel = {u: float(np.log2(share["Co-culture"][u] / share["Control"][u]))
           for u in common}
    rows = tab["Control"]
    y = np.array([sel[u] for u in common])
    feats = {"means only": np.column_stack(
        [np.ones(len(common)), [rows[u]["d_mean_a"] for u in common],
         [rows[u]["d_mean_h"] for u in common]]),
        "means + cov": np.column_stack(
        [np.ones(len(common)), [rows[u]["d_mean_a"] for u in common],
         [rows[u]["d_mean_h"] for u in common],
         [rows[u]["d_cov"] for u in common]])}
    q4 = {"n": len(common)}
    for lab, Xf in feats.items():
        beta, *_ = np.linalg.lstsq(Xf, y, rcond=None)
        q4[lab] = float(1 - np.sum((y - Xf @ beta) ** 2)
                        / max(np.sum((y - y.mean()) ** 2), 1e-12))
    # permutation test on the added R2 of the cov column
    rng2 = np.random.default_rng(seed)
    cov_col = np.array([rows[u]["d_cov"] for u in common])
    add = q4["means + cov"] - q4["means only"]
    null = []
    for _ in range(5000):
        Xp = feats["means + cov"].copy()
        Xp[:, 3] = cov_col[rng2.permutation(len(cov_col))]
        b, *_ = np.linalg.lstsq(Xp, y, rcond=None)
        r2p = 1 - np.sum((y - Xp @ b) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-12)
        null.append(r2p - q4["means only"])
    q4["added_r2"] = float(add)
    q4["p_added"] = float(np.mean(np.array(null) >= add))
    sp_cov = stats.spearmanr(cov_col, y)
    sp_a = stats.spearmanr([rows[u]["d_mean_a"] for u in common], y)
    sp_h = stats.spearmanr([rows[u]["d_mean_h"] for u in common], y)
    q4["spearman"] = {"d_cov": (float(sp_cov[0]), float(sp_cov[1])),
                      "d_mean_HLA_A": (float(sp_a[0]), float(sp_a[1])),
                      "d_mean_CD58": (float(sp_h[0]), float(sp_h[1]))}
    # known positives of this screen as a sanity check on the selection score
    pos = [u for u in ("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2", "NLRC5",
                       "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "CD58")
           if u in sel]
    q4["known_resistance_genes"] = {u: sel[u] for u in pos}
    q4["median_all"] = float(np.median(list(sel.values())))
    out["q4"] = q4
    out["_sel"] = sel

    # ---- Q5 how much of kappa's spread is co-occurrence?
    q5 = {}
    for c, rows in tab.items():
        ks = [u for u in rows if u != "control"]
        prod = np.array([rows[u]["mean_a"] * rows[u]["mean_h"] for u in ks])
        cov = np.array([rows[u]["cov"] for u in ks])
        q5[c] = {"sd_product_of_means": float(np.std(prod)), "sd_cov": float(np.std(cov)),
                 "ratio_cov_over_product": float(np.std(cov) / max(np.std(prod), 1e-12)),
                 "mean_product": float(np.mean(prod)), "mean_cov": float(np.mean(cov))}
    out["q5"] = q5
    return out


# ------------------------------------------------------------------ report
def report(R):
    print("=" * 78)
    print("D24  co-occurrence sensitisation in melanoma (Frangieh Perturb-CITE-seq)")
    print("=" * 78)
    print(f"  {R['n_cells']} single-perturbation or control cells; conditions "
          f"{R['conditions']}; perturbations with >= {MIN_CELLS} cells: "
          + ", ".join(f"{c} {n}" for c, n in R["n_perturbations"].items()))

    print("\n  Q0  the depth artefact (covariance of two ADTs in the same cell)")
    print(f"      {'normalisation':14s} {'cov(HLA-A,CD58)':>16s} {'cov(iso,iso)':>13s} "
          f"{'cov(HLA-A,iso)':>15s} {'r(HLA-A,depth)':>15s}")
    for k, v in R["q0"].items():
        print(f"      {k:14s} {v['cov_A_H']:16.4f} {v['cov_isotype_pair']:13.4f} "
              f"{v['cov_A_isotype']:15.4f} {v['r_A_with_total']:15.3f}")

    print("\n  Q1  is the per-perturbation covariance estimable? (split-half agreement)")
    for c, row in R["q1"].items():
        for lab, v in row.items():
            print(f"      {c:11s} {lab:20s} r(Cov) {v['r_cov']:+.3f}   "
                  f"[r(mean) {v['r_mean_for_scale']:+.3f} for scale, "
                  f"{v['n_perts']} perturbations]")

    print("\n  Q2  is covariance redundant given the two means?")
    for c, v in R["q2"].items():
        print(f"      {c:11s} R2(Cov from means) {v['r2_cov_from_means']:.3f}   "
              f"(sd of dCov {v['sd_d_cov']:.4f}, n={v['n']})")

    print("\n  Q3  perturbations that move Cov while both means stay put")
    for c, v in R["q3"].items():
        print(f"      {c:11s} {v['n_quiet_movers']} such knockouts; top: " +
              ", ".join(f"{u} (dCov {d:+.3f}, n={n})" for u, d, _, _, n in v["top"][:5]))

    print("\n  Q1b does the ADT panel see its own knockouts? (CLR difference vs control)")
    for c, row in R["q1b_panel_sees_ko"].items():
        best = sorted(row.items(), key=lambda kv: kv[1]["z"])[:4]
        print(f"      {c:11s} " + ", ".join(
            f"{k} {v['delta']:+.2f} (z {v['z']:+.0f}, n={v['n']})" for k, v in best))
    print("      bootstrap SE of Cov(HLA-A, CD58) by cell count, vs the "
          "between-knockout spread of dCov (~0.02):")
    for c, row in R["q1b_se_vs_n"].items():
        print(f"      {c:11s} " + ", ".join(f"n={n}: {v:.4f}" for n, v in row.items()))

    q = R["q4"]
    print(f"\n  Q4  does Cov predict the screen's resistance read-out? (n = {q['n']} knockouts)")
    print(f"      selection score sanity check, known resistance genes "
          f"(median over all knockouts {q['median_all']:+.3f}):")
    for u, s in sorted(q["known_resistance_genes"].items(), key=lambda kv: -kv[1]):
        print(f"        {u:8s} {s:+.3f}")
    print(f"      R2 means only {q['means only']:.3f} -> means + Cov "
          f"{q['means + cov']:.3f}   added {q['added_r2']:+.3f}, "
          f"permutation p = {q['p_added']:.4f}")
    for k, (rho, p) in q["spearman"].items():
        print(f"      Spearman({k}, selection) {rho:+.3f} (p = {p:.3g})")

    print("\n  Q5  how much of kappa's spread across knockouts is co-occurrence?")
    for c, v in R["q5"].items():
        print(f"      {c:11s} sd(E[A]E[H]) {v['sd_product_of_means']:.4f} vs sd(Cov) "
              f"{v['sd_cov']:.4f}   ratio {v['ratio_cov_over_product']:.3f}")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    fig, ax = plt.subplots(1, 4, figsize=(21.5, 5.2))

    a = ax[0]
    labs = list(R["q0"])
    w = 0.26
    keys = [("cov_A_H", "HLA-A x CD58"), ("cov_isotype_pair", "isotype x isotype"),
            ("cov_A_isotype", "HLA-A x isotype")]
    for k, (key, lab) in enumerate(keys):
        a.bar(np.arange(len(labs)) + (k - 1) * w, [R["q0"][l][key] for l in labs], w,
              label=lab)
    a.set_xticks(range(len(labs)), labs)
    a.axhline(0, color="k", lw=0.8)
    a.set_ylabel("covariance across all cells")
    a.set_title("Q0  two ADTs share the cell's total count;\nCLR removes the artefact",
                fontsize=10)
    a.legend(fontsize=8)

    a = ax[1]
    conds = list(R["q1"])
    for k, lab in enumerate(["HLA-A x CD58", "isotype x isotype"]):
        a.bar(np.arange(len(conds)) + (k - 0.5) * 0.38,
              [R["q1"][c][lab]["r_cov"] for c in conds], 0.38, label=f"Cov, {lab}")
    a.plot(np.arange(len(conds)), [R["q1"][c]["HLA-A x CD58"]["r_mean_for_scale"]
                                   for c in conds], "k*", ms=12,
           label="mean of HLA-A (for scale)")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(range(len(conds)), conds, fontsize=9)
    a.set_ylabel("split-half agreement (r)")
    a.set_title("Q1  is the per-knockout covariance\nestimable at this ADT depth?",
                fontsize=10)
    a.legend(fontsize=7.5)

    a = ax[2]
    rows = R["table"]["Control"]
    sel = R["_sel"]
    ks = [u for u in rows if u != "control" and u in sel]
    x = [rows[u]["d_cov"] for u in ks]
    y = [sel[u] for u in ks]
    a.plot(x, y, ".", color="0.6", ms=5)
    for u in R["q4"]["known_resistance_genes"]:
        if u in rows:
            a.plot(rows[u]["d_cov"], sel[u], "o", color="C3", ms=6)
            a.annotate(u, (rows[u]["d_cov"], sel[u]), fontsize=7,
                       xytext=(3, 3), textcoords="offset points")
    a.axhline(0, color="k", lw=0.8)
    a.axvline(0, color="k", lw=0.8)
    a.set_xlabel("change in Cov(HLA-A, CD58) before T-cell pressure")
    a.set_ylabel("selection score (co-culture vs control)")
    rho, p = R["q4"]["spearman"]["d_cov"]
    a.set_title(f"Q4  does co-occurrence predict survival\nunder T cells? "
                f"Spearman {rho:+.2f} (p={p:.2g})", fontsize=10)

    a = ax[3]
    for c, row in R["q1b_se_vs_n"].items():
        ns = sorted(row)
        a.plot(ns, [row[n] for n in ns], "o-", label=c)
    sig = float(np.median([v["sd_d_cov"] for v in R["q2"].values()]))
    a.axhline(sig, color="k", ls="--", lw=1.2,
              label="whole between-knockout spread")
    a.axhline(sig / 4, color="C3", ls=":", lw=1.2, label="spread / 4 (usable)")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("cells per perturbation")
    a.set_ylabel("bootstrap SE of Cov(HLA-A, CD58)")
    a.set_title("Q1b  how many cells a usable covariance\nestimate would need",
                fontsize=10)
    a.legend(fontsize=7.5)

    fig.tight_layout()
    fn = os.path.join(FIG, "d24_cooccurrence_sensitisation.png")
    fig.savefig(fn, dpi=140)
    print(f"\nfigure: {fn}")


def main(force=False):
    if os.path.exists(RES) and not force:
        R = pickle.load(open(RES, "rb"))
    else:
        R = run()
        pickle.dump(R, open(RES, "wb"))
    report(R)
    figures(R)
    return R


if __name__ == "__main__":
    main(force=bool(os.environ.get("D24_FORCE")))
