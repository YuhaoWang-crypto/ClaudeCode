"""
Default keratinocyte cytokine LEVER LIBRARY for the state-restoration skill.

A lever library is a matrix whose columns are MEASURED responses to perturbing
one input each. `statepipe.rank_levers` and `minimal_adjustment_from_responses`
need exactly this: with the MRA Jacobian J = -R^-1, the effect direction of
lever j is J^-1 e_j = -R[:, j], the measured response itself.

Two open keratinocyte panels are combined, because neither alone covers the
levers the clinic actually uses:

  GSE36287   primary keratinocytes, 3 donors, GPL570 (Chiricozzi / Krueger):
             IFN-alpha, IFN-gamma, IL-4, IL-13, IL-17A, TNF vs untreated.
             This is the only source with IL-4 and TNF, and it carries the
             psoriasis and atopic-dermatitis levers that are approved drug
             targets. It is the DEFAULT.
  GSE305982  primary keratinocytes (and fibroblasts), 1 donor, RNA-seq:
             IL-9, IL-13, IL-17A, IL-22, IL-26, IFN-gamma over ImmunoCult,
             with ImmunoCult alone as the matched control, 3 h and 24 h.
             Adds IL-9, IL-22 and IL-26, which GSE36287 lacks.

Why IL-4 matters (demos/d18_il4_and_combinations.py): with IL-13 as the only
type-2 lever, the atopic-dermatitis question cannot be answered - IL-13's
crossed specificity is +0.004 (n.s.) while IL-4's is +0.097 [0.070, 0.125].
A library missing IL-4 is not a weaker library for atopic dermatitis, it is
the wrong one.

Cross-source note: each column is compared with the disease displacement on
its own (cosine / projection energy), so columns from different studies are
directionally comparable. They are NOT magnitude-comparable, so a multi-lever
least-squares solve (`minimal_adjustment_from_responses`) should be run within
one source; `lever_library` labels every column with its source for this.

Validated targets (from approved drugs, independent of any dataset here) are
in `VALIDATED`, so the crossed specificity check can be run automatically.
"""
import numpy as np

from netjac.data import geo

# indication -> levers that are approved drug targets for it
VALIDATED = {
    "psoriasis": {"IL-17A", "TNF"},                 # secukinumab/ixekizumab; etanercept/adalimumab
    "atopic dermatitis": {"IL-4", "IL-13"},         # dupilumab (IL-4Ralpha), tralokinumab, lebrikizumab
}

# GSE36287: treatment label -> lever name
_G36287 = {"IFNa-treated": "IFN-alpha", "IFNg-treated": "IFN-gamma", "IL4-treated": "IL-4",
           "IL13-treated": "IL-13", "IL17A-treated": "IL-17A", "TNF-treated": "TNF"}
# GSE305982: lever name -> substring of the stimulation label
_G305982 = {"IL-9": "IL-9", "IL-13": "IL-13", "IL-17A": "IL-17a", "IL-22": "IL-22",
            "IL-26": "IL-26", "IFN-gamma": "gamma IFN"}


def _gse36287():
    X, M = geo.series_matrix("GSE36287")
    G = geo.collapse_to_genes(geo.quantile_normalise(X.astype(float)), geo.gpl_symbols("GPL570"))
    ctrl = G[M.index[M.treatment == "untreated control"]].values.mean(1)
    cols, names = [], []
    for tag, name in _G36287.items():
        idx = M.index[M.treatment == tag]
        if len(idx):
            cols.append(G[idx].values.mean(1) - ctrl)
            names.append(name)
    return np.array(G.index), np.column_stack(cols), names, int((M.treatment == "untreated control").sum())


def _gse305982(cell="Keratinocytes", time="24 hours"):
    import pandas as pd
    import os
    import subprocess
    from netjac import DATA_CACHE
    d = os.path.join(DATA_CACHE, "geo_gse305982")
    fn = os.path.join(d, "geneCounts.csv.gz")
    if not os.path.exists(fn):
        os.makedirs(d, exist_ok=True)
        subprocess.run(["curl", "-sSL", "-o", fn,
                        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305982/suppl/"
                        "GSE305982_geneCounts.csv.gz"], check=True)
    C = pd.read_csv(fn, index_col=0)
    _, M = geo.series_matrix("GSE305982")
    lib = M["description"].str.extract(r"(lib\d+)")[0]
    keep = lib.notna() & lib.isin(C.columns)
    M, C = M[keep], C[lib[keep].values]
    L = np.log2(C / C.sum(0) * 1e6 + 1)
    L.columns = M.index
    sym = _ens2sym()
    L = L.loc[[g for g in L.index if g in sym]]
    L.index = [sym[g] for g in L.index]
    L = L.groupby(level=0).max()
    sel = (M["cell type"] == cell) & (M["time"] == time)
    base = sel & M["stimulation"].str.contains("Immunocult") & ~M["stimulation"].str.contains(r"\+")
    ctrl = L[M.index[base]].values.mean(1)
    cols, names = [], []
    for name, tag in _G305982.items():
        idx = M.index[sel & M["stimulation"].str.contains(tag, regex=False)]
        if len(idx):
            cols.append(L[idx].values.mean(1) - ctrl)
            names.append(name)
    return np.array(L.index), np.column_stack(cols), names, int(base.sum())


def _ens2sym():
    import gzip
    import os
    from netjac import DATA_CACHE
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


def lever_library(sources=("GSE36287", "GSE305982"), time="24 hours", prefer="GSE36287"):
    """Keratinocyte cytokine lever library.

    Returns dict with
      genes   gene symbols (intersection over the sources used)
      R       genes x levers, each column a measured response (log2)
      names   lever names
      source  per lever, which study it came from
      n_ctrl  per source, how many control samples the baseline averaged
    A lever present in several sources is taken from `prefer` (3 donors and
    the same platform as the atopic-dermatitis cohort), and the duplicate is
    reported in `also_in`.
    """
    loaded = {}
    if "GSE36287" in sources:
        loaded["GSE36287"] = _gse36287()
    if "GSE305982" in sources:
        loaded["GSE305982"] = _gse305982(time=time)
    if not loaded:
        raise ValueError("no sources requested")
    shared = None
    for g, _, _, _ in loaded.values():
        shared = set(g) if shared is None else (shared & set(g))
    shared = sorted(shared)
    order = [prefer] + [s for s in loaded if s != prefer]
    cols, names, src, also = [], [], [], {}
    for s in order:
        if s not in loaded:
            continue
        g, R, nm, _ = loaded[s]
        gi = {x: i for i, x in enumerate(g)}
        idx = [gi[x] for x in shared]
        for k, n in enumerate(nm):
            if n in names:
                also.setdefault(n, []).append(s)
                continue
            cols.append(R[idx, k])
            names.append(n)
            src.append(s)
    return {"genes": np.array(shared), "R": np.column_stack(cols), "names": names,
            "source": src, "also_in": also,
            "n_ctrl": {s: loaded[s][3] for s in loaded}, "time": time}


def crossed_specificity(lib, disease, indication, lever, n_boot=2000, seed=0):
    """Is `lever` more aligned with `indication`'s lesional displacement than
    with the other indication's? `disease` maps indication -> dict(genes, d).
    Returns (delta, lo, hi, p) by bootstrap over the shared genes."""
    others = [k for k in disease if k != indication]
    if not others:
        raise ValueError("need two indications for a crossed contrast")
    other = others[0]
    gi = {g: i for i, g in enumerate(lib["genes"])}
    shared = sorted(set(gi) & set(disease[indication]["genes"]) & set(disease[other]["genes"]))
    a = lib["R"][[gi[g] for g in shared], lib["names"].index(lever)]
    a = a - a.mean()
    Y = {}
    for k in (indication, other):
        di = {g: i for i, g in enumerate(disease[k]["genes"])}
        y = np.asarray(disease[k]["d"], float)[[di[g] for g in shared]]
        Y[k] = y - y.mean()

    def cos(u, v):
        return float(u @ v / max(np.linalg.norm(u) * np.linalg.norm(v), 1e-12))
    rng = np.random.default_rng(seed)
    obs = cos(a, Y[indication]) - cos(a, Y[other])
    bs = [cos(a[b], Y[indication][b]) - cos(a[b], Y[other][b])
          for b in (rng.integers(0, len(shared), len(shared)) for _ in range(n_boot))]
    return (obs, *np.percentile(bs, [2.5, 97.5]), float(np.mean(np.array(bs) <= 0)))


def self_check(verbose=True):
    """Build the default library and re-run the published checks on it, so a
    broken download or a changed accession is caught immediately.
    Expected (demos/d18_il4_and_combinations.py): IL-17A is the top lever for
    psoriasis, and IL-17A and IL-4 pass their crossed specificity contrast."""
    from netjac.demos import d17_keratinocyte_lever as d17
    lib = lever_library()
    dis = d17.disease_vectors()
    out = {"levers": list(zip(lib["names"], lib["source"])), "n_genes": len(lib["genes"]), "checks": {}}
    r = d17.rank(dis["psoriasis"]["d"], dis["psoriasis"]["genes"], lib["R"], lib["genes"], lib["names"])
    out["checks"]["psoriasis top lever is IL-17A"] = (r["top"] == "IL-17A", r["top"])
    for ind, lev in (("psoriasis", "IL-17A"), ("atopic dermatitis", "IL-4")):
        d, lo, hi, p = crossed_specificity(lib, dis, ind, lev)
        out["checks"][f"{lev} crossed specificity for {ind}"] = (lo > 0, f"{d:+.3f} [{lo:+.3f}, {hi:+.3f}]")
    out["ok"] = all(v[0] for v in out["checks"].values())
    if verbose:
        print(f"keratinocyte lever library: {len(lib['names'])} levers, {out['n_genes']} genes")
        for n, s in out["levers"]:
            print(f"   {n:10s} {s}")
        print("checks:")
        for k, (ok, val) in out["checks"].items():
            print(f"   [{'PASS' if ok else 'FAIL'}] {k}: {val}")
    return out


if __name__ == "__main__":
    import sys
    sys.exit(0 if self_check()["ok"] else 1)
