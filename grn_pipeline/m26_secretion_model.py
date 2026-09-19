"""M26 - a better secretion predictor, built by auditing the previous one.

The question this answers: given a T cell's state, how well can we predict how
much cytokine it secretes? Two earlier attempts disagreed, and the disagreement
is informative:

  * M24 (this repo) used the cell's OWN cytokine transcript and found the link
    weak: per-cell R^2 = 0.017 (IL-2), 0.119 (TNF), 0.252 (IFN-g), with
    secretion better predicted by the same cell's earlier secretion.
  * A parallel analysis used transcript + four pathway axes and reported
    high-secretor AUCs of 0.92 (IFN-g), 0.87 (TNF), 0.69 (IL-2) under random
    5-fold cross-validation on RAW capture counts.

Both can be right, because they predict different things. This module builds
one protocol that separates them, and then tries to beat both.

THE AUDIT. Raw antibody-capture counts scale with how deeply a cell was
captured and sequenced. So does almost every transcriptomic feature. A model
given raw counts as its target can score well by learning cell size. M26
therefore always fits a DEPTH-ONLY baseline first, and reports every other
model as a gain over it, on two targets:
    raw        top-25% by raw capture count      (the earlier definition)
    normalised top-25% by capture count / total surface-antibody counts

THE FEATURE LADDER, each tier added on top of the last:
    D  depth      log RNA UMI, log antibody total, log HTO total, genes detected
    G  own gene   + the cell's own cytokine transcript
    A  axes       + 13 T-cell activation axes, UCell rank scores
    P  protein    + the 11-antibody surface panel (CLR)
    X  cross      + the other two cytokine transcripts
    S  persistence + the same cell's EARLIER secretion window (Tracing only)

THE SPLIT. Random cell-level folds leak: cells in a hashtag share condition and
batch. M26 reports random 5-fold (for comparability with the earlier number)
AND leave-one-hashtag-out, which is the honest one.

Run:  python3 -m grn_pipeline.m26_secretion_model
"""

from __future__ import annotations

import gzip
import os

import numpy as np
from scipy import sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from grn_pipeline.m24_rna_secretion_coupling import (CACHE, END, TRACE, ensure_data,
                                                     _lines, _barcodes)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "figures")

# 13 T-cell activation axes. Curated from canonical TCR-signalling and
# activation biology; the same axis set the parallel analysis used, kept
# identical on purpose so the comparison is about protocol, not gene lists.
AXES = {
    "TCR_SIGNALING": ["CD3D", "CD3E", "CD3G", "CD247", "ZAP70", "LCK", "FYN", "LAT",
                      "LCP2", "PLCG1", "GRAP2", "VAV1", "PTPRC", "ITK", "TXK", "SH2D1A"],
    "NFAT": ["NFATC1", "NFATC2", "NFATC3", "IL2", "IFNG", "TNF", "CSF2", "CD69",
             "EGR2", "EGR3", "NR4A1", "NR4A2", "NR4A3", "DUSP2", "DUSP4", "DUSP5", "PPP3CA"],
    "NFKB": ["NFKB1", "NFKB2", "RELA", "RELB", "REL", "NFKBIA", "NFKBIZ", "TNFAIP3",
             "BIRC3", "BCL2A1", "TRAF1", "BCL10", "MALT1", "CARD11"],
    "AP1": ["JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3", "BATF", "BATF3"],
    "IL2_STAT5": ["IL2", "IL2RA", "IL2RB", "IL2RG", "JAK1", "JAK3", "STAT5A", "STAT5B",
                  "SOCS3", "CISH", "PIM1", "BCL2"],
    "CD28_COSTIM": ["CD28", "ICOS", "PIK3CD", "PIK3R1", "AKT1", "AKT2", "MTOR", "RHEB", "RPS6KB1"],
    "CELL_CYCLE": ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6",
                   "MCM7", "CDK1", "CDK2", "CDK4", "CCNB1", "CCND1", "CCNE1", "E2F1"],
    "GLYCOLYSIS": ["HK2", "GPI", "PFKL", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA",
                   "SLC2A1", "MYC", "HIF1A"],
    "AIM_MARKERS": ["CD69", "IL2RA", "CD40LG", "TNFRSF4", "TNFRSF9", "CD274", "ENTPD1",
                    "CXCR5", "PDCD1"],
    "EXHAUSTION": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "EOMES", "CTLA4", "LAYN"],
    "ANERGY": ["RNF128", "DGKZ", "DGKA", "CBLB", "ITCH", "NDRG1", "JDP2"],
    "CYTOKINE_EFF": ["IL2", "IFNG", "TNF", "CSF2", "IL3", "IL21", "XCL1", "CCL3", "CCL4",
                     "CCL5", "LTA", "LTB"],
    "EARLY_RESPONSE": ["EGR1", "EGR2", "EGR3", "FOS", "JUN", "NR4A1", "DUSP1", "DUSP2",
                       "IER2", "IER3", "MYC", "JUNB"],
}
CYT = [("IL2", "IL-2"), ("IFNG", "IFN-g"), ("TNF", "TNF-a")]
CYT_TR = [("IL2", "IL-2"), ("IFNG", "IFN-y"), ("TNF", "TNF-a")]


# ------------------------------------------------------------------ loading
def _load_sparse(prefix, tag):
    feats = _lines(f"{CACHE}/{prefix}_{tag}_features.tsv.gz")
    bcs = _barcodes(f"{CACHE}/{prefix}_{tag}_barcodes.tsv.gz")
    rows, cols, vals = [], [], []
    with gzip.open(f"{CACHE}/{prefix}_{tag}_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                dims = line.split()
                break
        n_f, n_c = int(dims[0]), int(dims[1])
        for line in fh:
            r, c, v = line.split()
            rows.append(int(r) - 1); cols.append(int(c) - 1); vals.append(float(v))
    M = sparse.csr_matrix((vals, (rows, cols)), shape=(n_f, n_c))
    return feats, bcs, M


def _load_rna(prefix):
    feats = _lines(f"{CACHE}/{prefix}_RNA_features.tsv.gz")
    syms = [l.split("\t")[1] if "\t" in l else l for l in feats]
    bcs = _barcodes(f"{CACHE}/{prefix}_RNA_barcodes.tsv.gz")
    rows, cols, vals = [], [], []
    with gzip.open(f"{CACHE}/{prefix}_RNA_matrix.mtx.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                dims = line.split()
                break
        n_g, n_c = int(dims[0]), int(dims[1])
        for line in fh:
            r, c, v = line.split()
            rows.append(int(r) - 1); cols.append(int(c) - 1); vals.append(float(v))
    M = sparse.csc_matrix((vals, (rows, cols)), shape=(n_g, n_c))
    return syms, bcs, M


# -------------------------------------------------------------------- UCell
def ucell(M_csc, syms, axes, max_rank=1500):
    """UCell (Andreatta & Carmona 2021): per cell, rank genes by expression and
    score each signature by its rank sum, so the score is depth-robust by
    construction. Genes outside the top max_rank get rank max_rank+1."""
    n_genes, n_cells = M_csc.shape
    idx_of = {g: i for i, g in enumerate(syms)}
    sig_idx = {}
    for name, genes in axes.items():
        ii = sorted({idx_of[g] for g in genes if g in idx_of})
        if len(ii) >= 3:
            sig_idx[name] = np.array(ii)
    names = sorted(sig_idx)
    out = np.zeros((n_cells, len(names)), dtype=np.float32)
    M_csr = M_csc.tocsr().T.tocsr()          # cells x genes
    for c in range(n_cells):
        s, e = M_csr.indptr[c], M_csr.indptr[c + 1]
        gi, gv = M_csr.indices[s:e], M_csr.data[s:e]
        if len(gv) == 0:
            continue
        r = stats.rankdata(-gv, method="average")
        keep = r <= max_rank
        rank_of = dict(zip(gi[keep], r[keep]))
        for j, nm in enumerate(names):
            ii = sig_idx[nm]
            n = len(ii)
            U = 0.0
            for g in ii:
                U += rank_of.get(g, max_rank + 1)
            out[c, j] = 1.0 - (U - n * (n + 1) / 2.0) / (n * max_rank)
    return names, out


def clr(M):
    """Centred log-ratio across features, the standard antibody-count transform."""
    L = np.log1p(M)
    return L - L.mean(axis=0, keepdims=True)


# ------------------------------------------------------------------ assembly
def build(sample="END"):
    P = END if sample == "END" else TRACE
    syms, rbc, RNA = _load_rna(P["rna"])
    cf, cbc, C = _load_sparse(P["cyt"], "Cytokine")
    sf, sbc, S = _load_sparse(P["srf"], "SurfaceMarker")
    hf, hbc, H = _load_sparse(P["hto"], "HTO")

    shared = set(cbc) & set(sbc) & set(hbc)
    common = [b for b in rbc if b in shared]
    pos = {b: i for i, b in enumerate(rbc)}
    ri = np.array([pos[b] for b in common])
    ci = np.array([{b: i for i, b in enumerate(cbc)}[b] for b in common])
    si = np.array([{b: i for i, b in enumerate(sbc)}[b] for b in common])
    hi = np.array([{b: i for i, b in enumerate(hbc)}[b] for b in common])

    R = RNA[:, ri]
    umi = np.asarray(R.sum(axis=0)).ravel()
    ngene = np.asarray((R > 0).sum(axis=0)).ravel()
    Cm = np.asarray(C[:, ci].todense())
    Sm = np.asarray(S[:, si].todense())
    Hm = np.asarray(H[:, hi].todense())
    ab_tot = Sm[:-1].sum(axis=0)                 # 11 antibodies, drop 'unmapped'
    hto_tot = Hm[:-1].sum(axis=0)
    ok = (umi > 0) & (ab_tot > 0) & (hto_tot > 0)

    R = R[:, ok]; Cm = Cm[:, ok]; Sm = Sm[:, ok]; Hm = Hm[:, ok]
    umi, ngene, ab_tot, hto_tot = umi[ok], ngene[ok], ab_tot[ok], hto_tot[ok]

    # log-normalised expression for the own-gene / cross-gene features
    Rn = R.multiply(sparse.csr_matrix(1e4 / umi))
    Rn = Rn.tocsc()
    Rn.data = np.log1p(Rn.data)

    ax_names, AX = ucell(R.tocsc(), syms, AXES)
    hto = Hm[:-1].argmax(axis=0)

    return dict(syms=syms, Rn=Rn, C=Cm, cf=cf, S=Sm, sf=sf, hto=hto,
                umi=umi, ngene=ngene, ab_tot=ab_tot, hto_tot=hto_tot,
                ax_names=ax_names, AX=AX, n=int(ok.sum()))


def gene_vec(d, sym):
    idx = [i for i, g in enumerate(d["syms"]) if g == sym]
    if not idx:
        return np.zeros(d["n"])
    return np.asarray(d["Rn"][idx[0], :].todense()).ravel()


# ------------------------------------------------------------------ evaluation
def _fit_eval(X, y, folds, model="lr"):
    aucs, aps = [], []
    for tr, te in folds:
        if y[tr].sum() < 5 or y[te].sum() < 5:
            continue
        sc = StandardScaler().fit(X[tr])
        Xtr, Xte = sc.transform(X[tr]), sc.transform(X[te])
        clf = (LogisticRegression(max_iter=2000, C=1.0) if model == "lr"
               else HistGradientBoostingClassifier(max_iter=200, learning_rate=0.06,
                                                   max_depth=4, random_state=0))
        clf.fit(Xtr, y[tr])
        p = clf.predict_proba(Xte)[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        aps.append(average_precision_score(y[te], p))
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(aps))


def _hvg(Rn, n=2000):
    """Indices of the n most variable log-normalised genes."""
    X = Rn.tocsr()
    mean = np.asarray(X.mean(axis=1)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=1)).ravel()
    var = sq - mean ** 2
    return np.argsort(-var)[:n]


def _fit_eval_pca(Rn, hv, extra, y, folds, n_comp=50):
    """Same protocol, but the transcriptome enters as PCs fitted on the training
    cells only. Leaking the PCA across folds would flatter this tier."""
    from sklearn.decomposition import TruncatedSVD
    X = np.asarray(Rn[hv, :].todense()).T
    aucs = []
    for tr, te in folds:
        if y[tr].sum() < 5 or y[te].sum() < 5:
            continue
        svd = TruncatedSVD(n_components=n_comp, random_state=0).fit(X[tr])
        Ztr = np.column_stack([extra[tr], svd.transform(X[tr])])
        Zte = np.column_stack([extra[te], svd.transform(X[te])])
        sc = StandardScaler().fit(Ztr)
        clf = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.06,
                                             max_depth=4, random_state=0)
        clf.fit(sc.transform(Ztr), y[tr])
        aucs.append(roc_auc_score(y[te], clf.predict_proba(sc.transform(Zte))[:, 1]))
    return float(np.mean(aucs)), float(np.std(aucs)), float("nan")


def _depth_residual(y, depth):
    """Residual of y after linear regression on the depth covariates. A model
    given only depth cannot beat chance on the resulting target, which makes
    every AUC above 0.5 attributable to something other than capture depth."""
    Z = np.column_stack([depth, np.ones(len(y))])
    beta = np.linalg.lstsq(Z, y, rcond=None)[0]
    return y - Z @ beta


def folds_random(y, k=5, seed=0):
    return list(StratifiedKFold(k, shuffle=True, random_state=seed).split(np.zeros(len(y)), y))


def folds_grouped(groups, y):
    out = []
    for g in np.unique(groups):
        te = np.where(groups == g)[0]
        tr = np.where(groups != g)[0]
        if y[te].sum() >= 5 and y[tr].sum() >= 5:
            out.append((tr, te))
    return out


def _figure(results, best, better, persist):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.6))
    labels = [ab for _, ab in CYT]
    x = np.arange(len(labels))
    w = 0.26

    ax = axes[0]
    combos = [("raw", "random5", "raw counts, random folds", "#c44e52"),
              ("raw", "by-hashtag", "raw counts, grouped folds", "#dd8452"),
              ("depth-residual", "by-hashtag", "depth-residual, grouped", "#4c72b0")]
    for i, (tgt, split, lab, col) in enumerate(combos):
        vals = [results[(sym, tgt, split, "D depth")][0] for sym, _ in CYT]
        ax.bar(x + (i - 1) * w, vals, w, label=lab, color=col)
    ax.axhline(0.5, ls="--", lw=1, color="grey")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("AUC from capture depth ALONE")
    ax.set_title("the audit: what a model scores\nknowing only how deeply the cell was read")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3, axis="y")

    ax = axes[1]
    tiers = ["D depth", "G +gene", "A +axes", "P +prot"]
    for i, t in enumerate(tiers):
        vals = [results[(sym, "depth-residual", "by-hashtag", t)][0] for sym, _ in CYT]
        ax.plot(x, vals, "o-", label=t)
    ax.plot(x, [better[sym]["pca_model"] for sym, _ in CYT], "s--",
            color="black", label="50 transcriptome PCs")
    ax.axhline(0.5, ls="--", lw=1, color="grey")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("AUC, depth-residual target, grouped folds")
    ax.set_title("the honest ladder: what each\nfeature tier actually adds")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3)

    ax = axes[2]
    keys = ["state", "earlier", "both"]
    cols = ["#4c72b0", "#dd8452", "#55a868"]
    for i, k in enumerate(keys):
        vals = [persist[sym][k] for sym, _ in CYT_TR]
        ax.bar(x + (i - 1) * w, vals, w,
               label={"state": "cell state", "earlier": "its earlier secretion",
                      "both": "both"}[k], color=cols[i])
    ax.axhline(0.5, ls="--", lw=1, color="grey")
    ax.set_xticks(x); ax.set_xticklabels([ab for _, ab in CYT_TR])
    ax.set_ylabel("AUC for later-window secretion")
    ax.set_title("the decisive test: cell state vs\n'this cell was already secreting'")
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3, axis="y")

    fig.suptitle("M26  secretion prediction, audited: depth confound removed, "
                 "grouped folds, measured ceiling", fontsize=12)
    fig.tight_layout()
    os.makedirs(FIGDIR, exist_ok=True)
    path = os.path.join(FIGDIR, "m26_secretion_model.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"\n  figure -> {path}")


def report() -> dict:
    print("=" * 80)
    print("M26  secretion prediction: audit the protocol, then try to beat it")
    print("     TRAPS-seq GSE200690")
    print("=" * 80)
    ensure_data()
    d = build("END")
    print(f"\nEndTimePoint: {d['n']} cells, {len(d['ax_names'])} activation axes scored")
    print(f"  median UMI {np.median(d['umi']):.0f}, median antibody counts "
          f"{np.median(d['ab_tot']):.0f}, hashtag groups {len(np.unique(d['hto']))}")

    cidx = {f.split("-")[0] + "-" + f.split("-")[1]: i for i, f in enumerate(d["cf"][:-1])}
    depth = np.column_stack([np.log(d["umi"]), np.log(d["ab_tot"]),
                             np.log(np.maximum(d["hto_tot"], 1)), d["ngene"]])
    prot = clr(d["S"][:-1]).T
    axes = d["AX"]

    results = {}
    print("\n" + "-" * 80)
    print("the audit: how much of the signal is capture depth?")
    print("-" * 80)
    print(f"{'cytokine':9s} {'target':11s} {'split':10s} "
          f"{'D depth':>9s} {'G +gene':>9s} {'A +axes':>9s} {'P +prot':>9s} {'X +cross':>9s}")

    for sym, ab in CYT:
        raw = d["C"][cidx[ab]]
        norm = raw / d["ab_tot"]
        own = gene_vec(d, sym)
        cross = np.column_stack([gene_vec(d, s) for s, _ in CYT if s != sym])
        tiers = {
            "D depth": depth,
            "G +gene": np.column_stack([depth, own]),
            "A +axes": np.column_stack([depth, own, axes]),
            "P +prot": np.column_stack([depth, own, axes, prot]),
            "X +cross": np.column_stack([depth, own, axes, prot, cross]),
        }
        # A ratio target is confounded the other way: dividing by the antibody
        # total makes log(antibody total) predict it mechanically. The clean
        # target is the residual of log secretion after depth is regressed out,
        # where a depth-only model must score ~0.5 by construction.
        resid = _depth_residual(np.log1p(raw), depth)
        for tname, tvals in (("raw", raw), ("normalised", norm), ("depth-residual", resid)):
            y = (tvals >= np.percentile(tvals, 75)).astype(int)
            for sname, folds in (("random5", folds_random(y)),
                                 ("by-hashtag", folds_grouped(d["hto"], y))):
                row = []
                for k, X in tiers.items():
                    auc, sd, ap = _fit_eval(X, y, folds)
                    results[(sym, tname, sname, k)] = (auc, sd, ap)
                    row.append(f"{auc:9.3f}")
                print(f"{ab:9s} {tname:11s} {sname:10s} " + " ".join(row))

    print("\nread the first column: that is what a model knows from capture depth alone.")

    # ---- the best honest model, and whether a nonlinear learner adds anything
    print("\n" + "-" * 80)
    print("best model on the honest target (depth-residual, leave-one-hashtag-out)")
    print("-" * 80)
    best = {}
    for sym, ab in CYT:
        raw = d["C"][cidx[ab]]
        resid = _depth_residual(np.log1p(raw), depth)
        y = (resid >= np.percentile(resid, 75)).astype(int)
        folds = folds_grouped(d["hto"], y)
        own = gene_vec(d, sym)
        cross = np.column_stack([gene_vec(d, s) for s, _ in CYT if s != sym])
        X = np.column_stack([depth, own, axes, prot, cross])
        a_lr, s_lr, ap_lr = _fit_eval(X, y, folds, "lr")
        a_gb, s_gb, ap_gb = _fit_eval(X, y, folds, "gb")
        a_d, _, ap_d = results[(sym, "depth-residual", "by-hashtag", "D depth")]
        a_g, _, _ = results[(sym, "depth-residual", "by-hashtag", "G +gene")]
        best[sym] = dict(depth=a_d, own_gene=a_g, linear=a_lr, boosted=a_gb,
                         ap_linear=ap_lr, ap_depth=ap_d, prevalence=float(y.mean()))
        print(f"  {ab:6s} depth-only {a_d:.3f} | own gene {a_g:.3f} | "
              f"full linear {a_lr:.3f}±{s_lr:.3f} | boosted {a_gb:.3f}±{s_gb:.3f} | "
              f"AP {ap_lr:.3f} (prevalence {y.mean():.2f})")
    results_best = best

    # ---- can we do better than 13 hand-curated axes?
    print("\n" + "-" * 80)
    print("can a data-driven representation beat the hand-curated axes?")
    print("  PCA is refit inside each training fold, so this is not leaked.")
    print("-" * 80)
    hv = _hvg(d["Rn"], n=2000)
    better = {}
    for sym, ab in CYT:
        raw = d["C"][cidx[ab]]
        resid = _depth_residual(np.log1p(raw), depth)
        y = (resid >= np.percentile(resid, 75)).astype(int)
        folds = folds_grouped(d["hto"], y)
        own = gene_vec(d, sym)
        cross = np.column_stack([gene_vec(d, s) for s, _ in CYT if s != sym])
        base = np.column_stack([depth, own, axes, prot, cross])
        a_base, _, _ = _fit_eval(base, y, folds, "gb")
        a_pca, _, _ = _fit_eval_pca(d["Rn"], hv, np.column_stack([depth, own, prot]),
                                    y, folds, n_comp=50)
        # diagnostic only: the other two cytokines this cell actually secreted.
        # Not usable inside a peptide -> cytokine chain, but it bounds how much
        # of secretion is a general secretory state rather than cytokine-specific.
        other = np.column_stack([_depth_residual(np.log1p(d["C"][cidx[o]]), depth)
                                 for _, o in CYT if o != ab])
        a_other, _, _ = _fit_eval(np.column_stack([depth, other]), y, folds, "gb")
        better[sym] = dict(axes_model=a_base, pca_model=a_pca, other_secretion=a_other)
        print(f"  {ab:6s} curated axes {a_base:.3f} | 50 transcriptome PCs {a_pca:.3f} | "
              f"[diagnostic] other cytokines this cell secreted {a_other:.3f}")

    print("\n" + "-" * 80)
    print("which features carry it (full linear model, standardised coefficients)")
    print("-" * 80)
    fnames = (["logUMI", "logAb", "logHTO", "nGene", "ownGene"] + list(d["ax_names"])
              + [f.split("-")[0] for f in d["sf"][:-1]]
              + [f"cross_{s}" for s, _ in CYT])
    for sym, ab in CYT:
        raw = d["C"][cidx[ab]]
        resid = _depth_residual(np.log1p(raw), depth)
        y = (resid >= np.percentile(resid, 75)).astype(int)
        own = gene_vec(d, sym)
        cross = np.column_stack([gene_vec(d, s) for s, _ in CYT if s != sym])
        X = np.column_stack([depth, own, axes, prot, cross])
        names = fnames[:X.shape[1]]
        sc = StandardScaler().fit(X)
        clf = LogisticRegression(max_iter=2000).fit(sc.transform(X), y)
        co = clf.coef_[0]
        top = np.argsort(-np.abs(co))[:6]
        print(f"  {ab:6s} " + "  ".join(f"{names[i]}:{co[i]:+.2f}" for i in top))

    # ---- the decisive test, on the sample that has two capture windows
    print("\n" + "-" * 80)
    print("the decisive test: does cell state beat 'this cell was already secreting'?")
    print("-" * 80)
    t = build("TRACE")
    tcidx = {f.rsplit("-", 1)[0]: i for i, f in enumerate(t["cf"][:-1])}
    tdepth = np.column_stack([np.log(t["umi"]), np.log(t["ab_tot"]),
                              np.log(np.maximum(t["hto_tot"], 1)), t["ngene"]])
    tprot = clr(t["S"][:-1]).T
    stim = t["hto"] == 1                      # Hashtag2_10hr
    print(f"  Tracing: {t['n']} cells, {int(stim.sum())} stimulated")
    persist = {}
    for sym, ab in CYT_TR:
        later_raw = t["C"][tcidx[f"{ab}_2nd"]][stim]
        earlier = t["C"][tcidx[f"{ab}_1st"]][stim] / t["ab_tot"][stim]
        later = _depth_residual(np.log1p(later_raw), tdepth[stim])
        y = (later >= np.percentile(later, 75)).astype(int)
        own = gene_vec(t, sym)[stim]
        cross = np.column_stack([gene_vec(t, s)[stim] for s, _ in CYT_TR if s != sym])
        Xstate = np.column_stack([tdepth[stim], own, t["AX"][stim], tprot[stim], cross])
        Xearly = np.column_stack([tdepth[stim], np.log1p(earlier * 1e3)])
        Xboth = np.column_stack([Xstate, np.log1p(earlier * 1e3)])
        f5 = folds_random(y)
        a_state, _, _ = _fit_eval(Xstate, y, f5)
        a_early, _, _ = _fit_eval(Xearly, y, f5)
        a_both, _, _ = _fit_eval(Xboth, y, f5)
        persist[sym] = dict(state=a_state, earlier=a_early, both=a_both)
        print(f"  {ab:6s} cell state {a_state:.3f} | earlier window {a_early:.3f} | "
              f"both {a_both:.3f}  (state adds {a_both-a_early:+.3f} over the "
              f"earlier window)")

    _figure(results, results_best, better, persist)

    print("\n" + "-" * 80)
    print("VERDICT")
    print("  The depth-only column is the number to compare against, not 0.5.")
    print("  Cell state does carry real signal above it, and the axes are where it")
    print("  lives. But a cell's own earlier secretion remains the single best")
    print("  predictor of its later secretion, which is what M24 found and what")
    print("  any peptide -> cytokine chain has to live with.")
    print("-" * 80)
    return dict(ladder={f"{k[0]}|{k[1]}|{k[2]}|{k[3]}": v for k, v in results.items()},
                best=results_best, persistence=persist, better=better,
                n_end=d["n"], n_trace=t["n"])


main = report

if __name__ == "__main__":
    report()
