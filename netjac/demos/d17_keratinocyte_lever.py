"""
D17  The psoriasis lever question, with a keratinocyte perturbation library.

D15 asked which module to adjust in psoriasis and got the wrong answer
(barrier, not IL-17), because its J came from one observational treatment.
D16 showed the fix works when the response library is measured, but only in
K562, where a genome-wide screen exists. This demo closes the loop in the
relevant cell type.

Perturbation library (GSE305982, open): cultured primary KERATINOCYTES and
fibroblasts from one donor, stimulated with six T-cell cytokines - IL-9,
IL-13, IL-17A, IL-22, IL-26, IFN-gamma - each on top of ImmunoCult (IC), with
IC alone as the matched control, at 3 h and 24 h, 3 replicates. Each cytokine
is a different lever, so the columns of R are measured response directions.

Disease displacements (already used in D14/D15):
  psoriasis         GSE85034  lesional - own non-lesional
  atopic dermatitis GSE130588 lesional - own non-lesional

CROSSED GROUND TRUTH, from the clinic and independent of every dataset here:
  psoriasis         IL-17A blockade is highly effective (secukinumab,
                    ixekizumab); IL-13 blockade is not.
  atopic dermatitis IL-13 blockade is effective (dupilumab via IL-4Ralpha,
                    tralokinumab, lebrikizumab); IL-17 blockade is not.
So the prediction is falsifiable in both directions at once: the top lever
must be IL-17A for psoriasis and IL-13 for atopic dermatitis.

The sign is part of the prediction: a driver cytokine must be aligned WITH
the lesional displacement, so restoring the state requires a negative
coefficient - removing that input, which is what a blocking antibody does.

Run:  python3 -m netjac.demos.d17_keratinocyte_lever
"""
import gzip
import os
import pickle

import numpy as np

from netjac import DATA_CACHE, FIG
from netjac import statepipe as sp
from netjac.data import geo
from netjac.demos import d14_treatment_restoration as d14

DIR = os.path.join(DATA_CACHE, "geo_gse305982")
RES = os.path.join(DATA_CACHE, "d17_results.pkl")
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305982/suppl/GSE305982_geneCounts.csv.gz"
CYT = {"IL-9": "IL-9", "IL-13": "IL-13", "IL-17A": "IL-17a", "IL-22": "IL-22",
       "IL-26": "IL-26", "IFN-gamma": "gamma IFN"}
TRUTH = {"psoriasis": "IL-17A", "atopic dermatitis": "IL-13"}


def _ens2sym():
    """Ensembl gene id -> symbol, from the NCBI gene_info already cached."""
    fn = os.path.join(DATA_CACHE, "islet_xin", "Homo_sapiens.gene_info.gz")
    out = {}
    with gzip.open(fn, "rt") as f:
        f.readline()
        for line in f:
            p = line.split("\t")
            for x in p[5].split("|"):
                if x.startswith("Ensembl:"):
                    out[x.split(":")[1]] = p[2]
    return out


def load_library():
    """Response library: per cell type, time and cytokine, the mean log2 CPM
    difference against the matched ImmunoCult-only control."""
    import pandas as pd
    fn = os.path.join(DIR, "geneCounts.csv.gz")
    if not os.path.exists(fn):
        import subprocess
        os.makedirs(DIR, exist_ok=True)
        subprocess.run(["curl", "-sSL", "-o", fn, URL], check=True)
    C = pd.read_csv(fn, index_col=0)
    _, M = geo.series_matrix("GSE305982")
    if "description" not in M.columns:
        raise RuntimeError("sample description with library name not found")
    libname = M["description"].str.extract(r"(lib\d+)")[0]   # counts columns
    keep = libname.notna() & libname.isin(C.columns)
    M = M[keep]
    C = C[libname[keep].values]
    cpm = C / C.sum(0) * 1e6
    L = np.log2(cpm + 1)
    L.columns = M.index
    sym = _ens2sym()
    L = L.loc[[g for g in L.index if g in sym]]
    L.index = [sym[g] for g in L.index]
    L = L.groupby(level=0).max()                     # one row per symbol
    stim = M["stimulation"].values
    cell = M["cell type"].values
    time = M["time"].values
    lib_out, ctrl_out = {}, {}
    for ct in ("Keratinocytes", "Fibroblasts"):
        for t in np.unique(time):
            m0 = (cell == ct) & (time == t) & np.array(["Immunocult" in s and "+" not in s for s in stim])
            if m0.sum() == 0:
                continue
            base = L.loc[:, M.index[m0]].mean(1)
            cols = {}
            for name, tag in CYT.items():
                m = (cell == ct) & (time == t) & np.array([tag in s for s in stim])
                if m.sum():
                    cols[name] = (L.loc[:, M.index[m]].mean(1) - base).values
            lib_out[(ct, t)] = {"genes": np.array(L.index), "R": np.column_stack([cols[c] for c in CYT if c in cols]),
                                "names": [c for c in CYT if c in cols], "n_ctrl": int(m0.sum())}
            ctrl_out[(ct, t)] = base.values
    # genes detectable in keratinocytes (mean CPM >= 1 across KC samples)
    kc = cell == "Keratinocytes"
    kc_expr = set(np.array(L.index)[(2 ** L.loc[:, M.index[kc]].values - 1).mean(1) >= 1])
    return {"lib": lib_out, "kc_expressed": kc_expr}


def disease_vectors():
    """Cohort-mean lesional - non-lesional displacement per indication."""
    out = {}
    P = d14.load_psoriasis()
    LS = np.array([p["x"]["LS"] for p in P["pts"]])
    NL = np.array([p["x"]["NL"] for p in P["pts"]])
    out["psoriasis"] = {"genes": P["genes"], "d": (LS - NL).mean(0), "n": len(LS)}
    G, M = d14.load_ad()
    pairs = []
    for s, m in M[M.tissue != "Normal"].groupby("subject id"):
        k = {(t, ti): i for t, ti, i in zip(m.time, m.tissue, m.index)}
        if ("Week0", "LS") in k and ("Week0", "NL") in k:
            pairs.append(G[k[("Week0", "LS")]].values - G[k[("Week0", "NL")]].values)
    out["atopic dermatitis"] = {"genes": np.array(G.index), "d": np.mean(pairs, 0), "n": len(pairs)}
    return out


def rank(d, genes_d, R, genes_R, names, restrict=None, n_perm=2000, seed=0):
    """Project a disease displacement on each measured cytokine response."""
    gi = {g: i for i, g in enumerate(genes_R)}
    shared = [g for g in genes_d if g in gi and (restrict is None or g in restrict)]
    idx_d = [i for i, g in enumerate(genes_d) if g in gi and (restrict is None or g in restrict)]
    idx_r = [gi[g] for g in shared]
    y = np.asarray(d, float)[idx_d]
    A = np.asarray(R, float)[idx_r]
    y = y - y.mean()
    A = A - A.mean(0)
    proj = (A.T @ y) / np.maximum((A ** 2).sum(0), 1e-12)        # signed coefficient
    energy = (A.T @ y) ** 2 / np.maximum((A ** 2).sum(0), 1e-12)
    cos = (A.T @ y) / np.maximum(np.linalg.norm(A, axis=0) * np.linalg.norm(y), 1e-12)
    rng = np.random.default_rng(seed)
    null = np.array([((A.T @ y[rng.permutation(len(y))]) ** 2
                      / np.maximum((A ** 2).sum(0), 1e-12)) for _ in range(n_perm)])
    p = [(float(np.mean(null[:, j] >= energy[j])) if energy[j] > 0 else 1.0) for j in range(A.shape[1])]
    o = np.argsort(-energy)
    return {"names": list(names), "energy": energy, "cos": cos, "coef": proj,
            "p": np.array(p), "order": [names[j] for j in o], "n_genes": len(shared),
            "top": names[int(o[0])]}


def crossed(LIB, DIS, ct="Keratinocytes", t="24 hours", n_boot=2000, seed=0):
    """Specificity, which does not depend on who is top-1: is the IL-17A
    lever more aligned with psoriasis than with atopic dermatitis, AND the
    IL-13 lever more aligned with atopic dermatitis than with psoriasis?
    Bootstrap over the shared genes."""
    lib = LIB["lib"][(ct, t)]
    gi = {g: i for i, g in enumerate(lib["genes"])}
    keep = LIB["kc_expressed"]
    shared = sorted(set(gi) & keep & set(DIS["psoriasis"]["genes"]) & set(DIS["atopic dermatitis"]["genes"]))
    A = np.asarray(lib["R"], float)[[gi[g] for g in shared]]
    A = A - A.mean(0)
    Y = {}
    for ind, dv in DIS.items():
        di = {g: i for i, g in enumerate(dv["genes"])}
        y = np.asarray(dv["d"], float)[[di[g] for g in shared]]
        Y[ind] = y - y.mean()
    names = lib["names"]

    def cosv(a, y):
        return float(a @ y / max(np.linalg.norm(a) * np.linalg.norm(y), 1e-12))
    rng = np.random.default_rng(seed)
    out = {"n_genes": len(shared), "cell": ct, "time": t, "contrasts": {}}
    for lever, pro in (("IL-17A", "psoriasis"), ("IL-13", "atopic dermatitis")):
        j = names.index(lever)
        other = "atopic dermatitis" if pro == "psoriasis" else "psoriasis"
        obs = cosv(A[:, j], Y[pro]) - cosv(A[:, j], Y[other])
        bs = []
        for _ in range(n_boot):
            b = rng.integers(0, len(shared), len(shared))
            bs.append(cosv(A[b, j], Y[pro][b]) - cosv(A[b, j], Y[other][b]))
        out["contrasts"][lever] = {"for": pro, "delta": obs,
                                   "ci": tuple(np.percentile(bs, [2.5, 97.5])),
                                   "p": float(np.mean(np.array(bs) <= 0))}
    return out


def robustness(LIB, DIS, ct="Keratinocytes", t="24 hours", drop=200):
    """Is the psoriasis -> IL-17A result carried by a handful of huge genes?
    Rank correlation, and dropping the most extreme disease genes."""
    from scipy.stats import spearmanr
    lib = LIB["lib"][(ct, t)]
    gi = {g: i for i, g in enumerate(lib["genes"])}
    keep = LIB["kc_expressed"]
    out = {}
    for ind in DIS:
        dv = DIS[ind]
        di = {g: i for i, g in enumerate(dv["genes"])}
        shared = sorted(set(gi) & keep & set(di))
        A = np.asarray(lib["R"], float)[[gi[g] for g in shared]]
        y = np.asarray(dv["d"], float)[[di[g] for g in shared]]
        sp_ = [spearmanr(A[:, j], y).statistic for j in range(A.shape[1])]
        o = np.argsort(-np.abs(y))[drop:]            # drop the strongest disease genes
        A2, y2 = A[o] - A[o].mean(0), y[o] - y[o].mean()
        cos2 = [float(A2[:, j] @ y2 / max(np.linalg.norm(A2[:, j]) * np.linalg.norm(y2), 1e-12))
                for j in range(A.shape[1])]
        out[ind] = {"names": lib["names"], "spearman_top": lib["names"][int(np.argmax(sp_))],
                    "spearman": sp_, "drop_top": lib["names"][int(np.argmax(cos2))], "drop_cos": cos2}
    return out


def run():
    LIB = load_library()
    DIS = disease_vectors()
    R = {"n_genes_disease": {k: len(v["genes"]) for k, v in DIS.items()},
         "n_subjects": {k: v["n"] for k, v in DIS.items()}, "res": {}, "kc_n": len(LIB["kc_expressed"])}
    for (ct, t), lib in LIB["lib"].items():
        for ind, dv in DIS.items():
            for restrict, tag in ((None, "all shared"), (LIB["kc_expressed"], "keratinocyte-expressed")):
                key = (ind, ct, t, tag)
                R["res"][key] = rank(dv["d"], dv["genes"], lib["R"], lib["genes"], lib["names"],
                                     restrict=restrict)
    R["truth"] = TRUTH
    R["crossed"] = {t: crossed(LIB, DIS, t=t) for t in ("03 hours", "24 hours")}
    R["robust"] = robustness(LIB, DIS)
    pickle.dump(R, open(RES, "wb"))
    return R


def report(R):
    print(f"\nD17 keratinocyte cytokine response library (GSE305982) vs lesional displacement")
    print(f"   disease vectors: " + ", ".join(f"{k} n={v}" for k, v in R["n_subjects"].items()))
    print(f"   ground truth (clinic, independent of all data here): " +
          ", ".join(f"{k} -> {v}" for k, v in R["truth"].items()))
    main = [(k, v) for k, v in R["res"].items() if k[1] == "Keratinocytes" and k[3] == "keratinocyte-expressed"]
    print("\n   MAIN TEST  keratinocytes, genes expressed in keratinocytes:")
    for (ind, ct, t, tag), v in sorted(main):
        hit = "HIT" if v["top"] == R["truth"][ind] else "miss"
        print(f"      {ind:19s} {t:9s} top lever {v['top']:10s} [{hit}]  ({v['n_genes']} genes)")
        for j in np.argsort(-v["energy"]):
            n = v["names"][j]
            star = " <- truth" if n == R["truth"][ind] else ""
            print(f"         {n:10s} cos {v['cos'][j]:+.3f}  coef {v['coef'][j]:+.3f}  "
                  f"p={v['p'][j]:.4f}{star}")
    print("\n   CROSSED SPECIFICITY (does not depend on who is top-1):")
    for t, c in sorted(R["crossed"].items()):
        print(f"      {t}, keratinocytes, {c['n_genes']} genes:")
        for lever, v in c["contrasts"].items():
            ok = "OK" if v["ci"][0] > 0 else "n.s."
            print(f"         cos({lever}) is higher for {v['for']:18s} than for the other: "
                  f"{v['delta']:+.3f} [{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}] p={v['p']:.4f} [{ok}]")
    print("\n   ROBUSTNESS of the ranking (keratinocytes, 24 h):")
    for ind, v in R["robust"].items():
        print(f"      {ind:19s} rank-based top {v['spearman_top']:10s} | after dropping the 200 "
              f"strongest disease genes: {v['drop_top']:10s} (truth {R['truth'][ind]})")
    print("\n   TOP LEVER PER SETTING:")
    print(f"      {'indication':19s} {'cell':14s} {'time':9s} {'genes':22s} top")
    for (ind, ct, t, tag), v in sorted(R["res"].items()):
        hit = "HIT" if v["top"] == R["truth"][ind] else "miss"
        print(f"      {ind:19s} {ct:14s} {t:9s} {tag:22s} {v['top']:10s} [{hit}]")


def figures(R):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    inds = list(R["truth"])
    times = sorted({k[2] for k in R["res"]})
    fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.2))
    for a, t in zip(ax[:2], times):
        keys = [(i, "Keratinocytes", t, "keratinocyte-expressed") for i in inds]
        names = R["res"][keys[0]]["names"]
        w = 0.38
        for k, (ind, key) in enumerate(zip(inds, keys)):
            v = R["res"][key]
            col = ["C3" if n == R["truth"][ind] else "0.7" for n in v["names"]]
            a.bar(np.arange(len(names)) + (k - 0.5) * w, v["cos"], w,
                  color=col, edgecolor="k" if k else None, linewidth=0.6,
                  label=f"{ind} (truth: {R['truth'][ind]})", hatch="" if k == 0 else "//")
        a.axhline(0, color="k", lw=0.8)
        a.set_xticks(range(len(names)), names, rotation=30, ha="right", fontsize=8)
        a.set_ylabel("cosine with lesional displacement")
        a.set_title(f"Keratinocyte cytokine levers, {t}\n(red = clinically validated lever)", fontsize=10)
        a.legend(fontsize=7)
    a = ax[2]
    labs, vals, los, his, cols = [], [], [], [], []
    for t in sorted(R["crossed"]):
        for lever, v in R["crossed"][t]["contrasts"].items():
            labs.append(f"{lever} for\n{v['for'][:18]}\n({t.split()[0]} h)")
            vals.append(v["delta"])
            los.append(v["delta"] - v["ci"][0])
            his.append(v["ci"][1] - v["delta"])
            cols.append("C2" if v["ci"][0] > 0 else "C1")
    a.bar(range(len(vals)), vals, yerr=[los, his], color=cols, capsize=4)
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(range(len(labs)), labs, fontsize=7)
    a.set_ylabel("cos(lever, its own indication)\n- cos(lever, the other indication)")
    a.set_title("Crossed specificity: is each lever more aligned\nwith the indication it treats? (green = CI excludes 0)",
                fontsize=9)
    fig.tight_layout()
    p = os.path.join(FIG, "d17_keratinocyte_lever.png")
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
