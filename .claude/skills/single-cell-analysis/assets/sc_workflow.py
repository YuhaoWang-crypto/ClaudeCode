#!/usr/bin/env python3
"""Reference implementation of the single-cell-analysis skill spine.

Steps 2-12 of SKILL.md with the benchmark-backed defaults, each printing an
audit line so you can see what every filter actually did.

    python3 sc_workflow.py --demo
    python3 sc_workflow.py --h5ad data.h5ad --sample-key sample --condition-key condition

Requires: scanpy, leidenalg, igraph.  Optional: pydeseq2 (proper pseudobulk DE
engine; a labelled fallback runs without it).

Deliberately dependency-light: deviance feature selection and pseudobulk
aggregation are implemented here in numpy so the whole spine runs anywhere.
Swap in scry / decoupler / PyDESeq2 for production.
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=FutureWarning)

STEP = 0


def audit(msg: str) -> None:
    global STEP
    STEP += 1
    print(f"[{STEP:02d}] {msg}", flush=True)


# --------------------------------------------------------------------------
# QC — MAD-based outlier detection (sc-best-practices, QC chapter)  [benchmark-backed]
# --------------------------------------------------------------------------
def is_outlier(adata, metric: str, nmads: int) -> np.ndarray:
    """Median-absolute-deviation outlier flag.

    Fixed cutoffs (n_genes > 500, mito < 5%) bias against small/low-RNA
    populations.  MADs adapt to the dataset.  Defaults: 5 MADs for depth-like
    metrics, 3 MADs for mitochondrial fraction.
    """
    x = np.asarray(adata.obs[metric], dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:  # degenerate metric: nothing to call
        return np.zeros(len(x), dtype=bool)
    return (x < med - nmads * mad) | (x > med + nmads * mad)


def quality_control(adata, sc, mito_prefix: str = "MT-", mito_hard_pct: float = 8.0):
    adata.var["mt"] = adata.var_names.str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True
    )

    bad = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    bad_mt = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > mito_hard_pct
    )
    adata.obs["qc_fail"] = bad | bad_mt

    audit(
        f"QC: {int(bad.sum())} depth/complexity outliers (5 MAD), "
        f"{int(bad_mt.sum())} mito outliers (3 MAD or >{mito_hard_pct}%), "
        f"{int(adata.obs['qc_fail'].sum())} unique -> "
        f"{adata.n_obs - int(adata.obs['qc_fail'].sum())}/{adata.n_obs} cells kept"
    )
    # Inspect what you removed before committing to it: if discarded cells form
    # a coherent cluster with real markers, you just deleted a cell type.
    return adata[~adata.obs["qc_fail"]].copy()


# --------------------------------------------------------------------------
# Feature selection — binomial deviance on RAW counts (Townes 2019 / scry)  [benchmark-backed]
# --------------------------------------------------------------------------
def binomial_deviance(counts) -> np.ndarray:
    """Per-gene binomial deviance under the constant-expression null.

    D_j = 2 * sum_i [ y log(y / (n_i pi_j)) + (n_i - y) log((n_i-y)/(n_i(1-pi_j))) ]

    High deviance = the gene is NOT constant across cells = informative.
    Computed on raw counts, so no pseudo-count/log distortion of the zeros.
    """
    counts = sp.csc_matrix(counts) if sp.issparse(counts) else np.asarray(counts)
    n_i = np.asarray(counts.sum(axis=1)).ravel()  # per-cell depth
    total = n_i.sum()
    gene_sums = np.asarray(counts.sum(axis=0)).ravel()
    pi = gene_sums / total
    pi = np.clip(pi, 1e-12, 1 - 1e-12)

    dev = np.zeros(counts.shape[1])
    for j in range(counts.shape[1]):
        y = counts[:, j]
        y = np.asarray(y.todense()).ravel() if sp.issparse(y) else np.asarray(y).ravel()
        m = n_i * pi[j]
        term1 = np.where(y > 0, y * np.log(np.maximum(y, 1e-12) / m), 0.0)
        ny = n_i - y
        m2 = n_i * (1 - pi[j])
        term2 = np.where(ny > 0, ny * np.log(np.maximum(ny, 1e-12) / m2), 0.0)
        dev[j] = 2.0 * (term1 + term2).sum()
    return dev


# --------------------------------------------------------------------------
# Pseudobulk — the one step that decides whether your DE p-values mean anything
# --------------------------------------------------------------------------
def pseudobulk(adata, sample_key: str, group_key: str, layer: str = "counts",
               min_cells: int = 10, min_counts: int = 1000):
    """Sum raw counts per (sample, cell type).  Sum, not mean, not per-cell.

    Cells from one donor are not independent replicates; testing them as if
    they were is pseudoreplication and inflates the FDR.  Aggregate first,
    then run a bulk method (DESeq2/edgeR) on the samples.
    """
    X = adata.layers[layer] if layer in adata.layers else adata.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    keys, rows, meta = [], [], []
    for (samp, grp), idx in adata.obs.groupby(
        [sample_key, group_key], observed=True
    ).indices.items():
        if len(idx) < min_cells:
            continue
        summed = np.asarray(X[idx].sum(axis=0)).ravel()
        if summed.sum() < min_counts:
            continue
        keys.append(f"{samp}|{grp}")
        rows.append(summed)
        meta.append({sample_key: samp, group_key: grp, "n_cells": len(idx)})

    pb = pd.DataFrame(rows, index=keys, columns=adata.var_names)
    return pb, pd.DataFrame(meta, index=keys)


def filter_by_expr(pb: pd.DataFrame, min_count: int = 10, min_prop: float = 0.5):
    """edgeR-style: keep genes with >=min_count in >=min_prop of pseudobulks.
    Must be redone per cell type -- cell types express different gene sets."""
    keep = (pb >= min_count).mean(axis=0) >= min_prop
    return pb.loc[:, keep]


def check_design(design: pd.DataFrame, condition_key: str, min_reps: int = 2):
    """Refuse to test a design that cannot support the contrast.

    The common failure: a strong condition effect makes Leiden split one cell
    type into condition-specific clusters.  The 'cell type' is then collinear
    with the contrast and DE within it is undefined -- annotate to cell type
    (not raw cluster ids), or integrate, before comparing.
    """
    levels = design[condition_key].astype(str).value_counts()
    if len(levels) < 2:
        return False, f"only one condition level present ({dict(levels)})"
    if levels.min() < min_reps:
        return False, f"fewer than {min_reps} replicates in a level ({dict(levels)})"
    return True, f"{dict(levels)}"


def pseudobulk_de(pb: pd.DataFrame, design: pd.DataFrame, condition_key: str,
                  ref_level: str):
    """DE on pseudobulk.  PyDESeq2 if installed, else a labelled fallback."""
    alt = [l for l in design[condition_key].astype(str).unique() if l != ref_level]
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        dds = DeseqDataSet(counts=pb.astype(int), metadata=design,
                           design=f"~{condition_key}", quiet=True)
        dds.deseq2()
        # contrast = [factor, tested_level, reference_level]  (pydeseq2 >= 0.4)
        res = DeseqStats(dds, contrast=[condition_key, alt[0], ref_level], quiet=True)
        res.summary()
        return res.results_df.sort_values("padj"), "PyDESeq2 (negative binomial GLM)"
    except ImportError:
        # Fallback: log-CPM + Welch t-test + BH.  Fine for a smoke test,
        # NOT what you report -- it ignores the count noise model.
        from scipy import stats

        cpm = np.log2(pb.div(pb.sum(axis=1), axis=0) * 1e6 + 1)
        grp = design[condition_key].astype(str)
        a, b = cpm[grp.values == ref_level], cpm[grp.values != ref_level]
        t, p = stats.ttest_ind(b, a, axis=0, equal_var=False)
        order = np.argsort(p)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, len(p) + 1)
        padj = np.minimum.accumulate(
            (p[order] * len(p) / ranks[order])[::-1]
        )[::-1]
        out = pd.DataFrame(
            {"log2FoldChange": b.mean(axis=0) - a.mean(axis=0), "stat": t,
             "pvalue": p, "padj": pd.Series(np.clip(padj, 0, 1), index=cpm.columns[order])},
            index=cpm.columns,
        ).sort_values("padj")
        return out, "FALLBACK log-CPM t-test (install pydeseq2 for real inference)"


# --------------------------------------------------------------------------
# Demo data: known ground truth so the whole spine is checkable
# --------------------------------------------------------------------------
def make_demo(sc, seed: int = 0):
    """8 donors x 2 conditions x 2 batches, 6 cell types, planted DE + abundance shift."""
    import anndata as ad

    rng = np.random.default_rng(seed)
    n_genes, n_mito = 400, 13
    types = ["Tcell_CD4", "Tcell_CD8", "Bcell", "NK", "Monocyte", "DC"]
    markers = {t: list(range(i * 12, i * 12 + 12)) for i, t in enumerate(types)}
    de_genes = list(range(200, 220))  # planted condition effect in monocytes

    base = rng.gamma(2.0, 1.0, n_genes)
    blocks, obs = [], []
    for donor in range(8):
        cond = "ctrl" if donor < 4 else "stim"
        batch = f"batch{donor % 2}"
        for ct in types:
            n = rng.integers(60, 110)
            if ct == "Monocyte" and cond == "stim":
                n = int(n * 2.2)  # real abundance expansion
            mu = np.tile(base, (n, 1))
            mu[:, markers[ct]] *= 9.0
            if ct == "Monocyte" and cond == "stim":
                mu[:, de_genes] *= 3.0
            mu[:, : n_genes // 2] *= 1.35 if batch == "batch1" else 1.0  # batch effect
            depth = rng.lognormal(0.0, 0.35, (n, 1))
            counts = rng.poisson(mu * depth * 0.55)
            blocks.append(counts)
            obs += [(f"donor{donor}", cond, batch, ct)] * n

    X = np.vstack(blocks).astype(np.float32)
    obs = pd.DataFrame(obs, columns=["sample", "condition", "batch", "true_celltype"])
    obs.index = [f"cell{i}" for i in range(len(obs))]

    # inject 4% junk: dying (mito-high) and empty-ish (low depth) droplets
    n_bad = int(0.04 * X.shape[0])
    bad = rng.choice(X.shape[0], n_bad, replace=False)
    X[bad[: n_bad // 2], -n_mito:] *= 25          # dying cells
    X[bad[n_bad // 2:], :] = (X[bad[n_bad // 2:], :] * 0.06)  # ambient-only droplets

    var_names = [f"G{i}" for i in range(n_genes - n_mito)] + \
                [f"MT-{i}" for i in range(n_mito)]
    adata = ad.AnnData(sp.csr_matrix(X), obs=obs,
                       var=pd.DataFrame(index=var_names))
    marker_panel = {t: [var_names[i] for i in markers[t]] for t in types}
    return adata, marker_panel


# --------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true", help="run on synthetic data")
    ap.add_argument("--h5ad")
    ap.add_argument("--sample-key", default="sample")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--batch-key", default="batch")
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--n-pcs", type=int, default=50)
    args = ap.parse_args()

    import scanpy as sc

    sc.settings.verbosity = 0

    if args.demo:
        adata, marker_panel = make_demo(sc)
        audit(f"demo data: {adata.n_obs} cells x {adata.n_vars} genes, "
              f"{adata.obs['sample'].nunique()} donors, "
              f"{adata.obs['true_celltype'].nunique()} planted cell types")
    elif args.h5ad:
        adata = sc.read_h5ad(args.h5ad)
        marker_panel = {}
        audit(f"loaded {args.h5ad}: {adata.n_obs} cells x {adata.n_vars} genes")
    else:
        ap.error("pass --demo or --h5ad")

    adata.layers["counts"] = adata.X.copy()

    # 2-4  QC.  Ambient RNA (SoupX) and doublets (scDblFinder) run BEFORE this
    #      in R, per sample -- never on pooled batches.
    adata = quality_control(adata, sc)

    # 5  Normalization: shifted log -- best for DR and DE (sc-best-practices).
    sc.pp.normalize_total(adata, target_sum=None)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    audit("normalized: shifted log (median-depth size factors)")

    # 6  Feature selection: deviance on RAW counts.
    dev = binomial_deviance(adata.layers["counts"])
    adata.var["binomial_deviance"] = dev
    n_hvg = min(args.n_hvg, adata.n_vars)
    adata.var["highly_deviant"] = adata.var["binomial_deviance"].rank(
        ascending=False) <= n_hvg
    audit(f"feature selection: {int(adata.var['highly_deviant'].sum())} highly "
          f"deviant genes (deviance range "
          f"{dev.min():.0f}-{dev.max():.0f})")

    # 7  Embedding.  Batch-aware HVGs matter here if you plan to integrate.
    n_pcs = min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, n_comps=n_pcs, mask_var="highly_deviant")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)
    var_ratio = adata.uns["pca"]["variance_ratio"]
    audit(f"PCA {n_pcs} PCs ({100*var_ratio.sum():.1f}% variance) -> kNN(15) -> UMAP "
          "(UMAP = display only; never measure distances on it)")

    # 8  Clustering at several resolutions -- resolution is a CONVENTION, so
    #    show the sensitivity instead of hiding it.
    for res in (0.25, 0.5, 1.0):
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}",
                     flavor="igraph", n_iterations=2, directed=False)
    audit("Leiden: " + ", ".join(
        f"res={r} -> {adata.obs[f'leiden_{r}'].nunique()} clusters"
        for r in (0.25, 0.5, 1.0)))

    cluster_key = "leiden_1.0"

    # 9  Annotation.  Marker scoring here; see reference/annotation.md for
    #    SingleR / CellTypist / scArches and the reconciliation discipline.
    if marker_panel:
        for ct, genes in marker_panel.items():
            genes = [g for g in genes if g in adata.var_names]
            sc.tl.score_genes(adata, genes, score_name=f"score_{ct}")
        scores = adata.obs[[f"score_{c}" for c in marker_panel]]
        per_cluster = scores.groupby(adata.obs[cluster_key], observed=True).mean()
        call = per_cluster.idxmax(axis=1).str.replace("score_", "", regex=False)
        # margin = confidence.  A cluster whose top two scores are within 10% of
        # the score range is ambiguous -- label it, don't guess it.
        top2 = np.sort(per_cluster.values, axis=1)[:, -2:]
        margin = top2[:, 1] - top2[:, 0]
        ambiguous = margin < 0.1 * (per_cluster.values.max() - per_cluster.values.min())
        call[ambiguous] = "Unknown"
        adata.obs["celltype"] = adata.obs[cluster_key].map(call).astype(str)

        agree = (adata.obs["celltype"] == adata.obs["true_celltype"]).mean()
        audit(f"annotation: {call.nunique()} labels from marker scores, "
              f"{int(ambiguous.sum())} ambiguous clusters -> Unknown; "
              f"agreement with planted truth = {100*agree:.1f}%")
    else:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")
        adata.obs["celltype"] = adata.obs[cluster_key].astype(str)
        audit("annotation: no marker panel given -- ranked cluster markers only. "
              "rank_genes_groups p-values are circular; use them to FIND "
              "markers, report effect sizes instead.")

    # 10 Integration: check whether you need it before doing it.
    if args.batch_key in adata.obs:
        mix = pd.crosstab(adata.obs[cluster_key], adata.obs[args.batch_key],
                          normalize="index")
        worst = float(mix.max(axis=1).max())
        audit(f"batch check: most batch-skewed cluster is {100*worst:.0f}% one "
              f"batch ({'integration likely needed' if worst > 0.9 else 'no integration needed'})"
              " -- judge with scIB metrics, not UMAP")

    # 11 Differential expression -- PSEUDOBULK.
    if args.condition_key in adata.obs and args.sample_key in adata.obs:
        target = adata.obs["celltype"].value_counts().idxmax()
        sub = adata[adata.obs["celltype"] == target].copy()
        pb, design = pseudobulk(sub, args.sample_key, "celltype")
        pb = filter_by_expr(pb)
        design[args.condition_key] = (
            sub.obs.groupby(args.sample_key, observed=True)[args.condition_key]
            .first().reindex(design[args.sample_key]).values
        )
        ok, why = check_design(design, args.condition_key)
        if not ok:
            audit(f"pseudobulk DE SKIPPED for '{target}': {why}. This population "
                  "is confounded with the contrast -- typically because "
                  "clustering split one cell type by condition. Annotate to "
                  "cell type (not raw cluster ids), or integrate, then retry.")
        else:
            ref = sorted(design[args.condition_key].unique())[0]
            res, engine = pseudobulk_de(pb, design, args.condition_key, ref)
            n_sig = int((res["padj"] < 0.05).sum())
            audit(f"pseudobulk DE in '{target}' ({why}): {pb.shape[0]} pseudobulk "
                  f"samples x {pb.shape[1]} genes, {n_sig} genes padj<0.05 [{engine}]")
            print(res.head(8).to_string())

        # 12 Compositional analysis -- proportions are compositional; a t-test
        #    on them is wrong.  Report the counts, test with scCODA/Milo.
        comp = pd.crosstab(adata.obs[args.sample_key], adata.obs["celltype"],
                           normalize="index")
        comp[args.condition_key] = (
            adata.obs.groupby(args.sample_key, observed=True)[args.condition_key].first()
        )
        audit("composition per sample (test with scCODA/Milo, NOT a t-test "
              "on these proportions -- they sum to 1):")
        print(comp.groupby(args.condition_key).mean().round(3).to_string())

    print("\nDone. Every ⚠️-convention parameter above (n_PCs, resolution, "
          "n_HVG) should be varied to check the conclusion survives.")


if __name__ == "__main__":
    main()
