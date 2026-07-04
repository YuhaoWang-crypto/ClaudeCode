#!/usr/bin/env python3
"""
scRNA-seq reanalysis of Sade-Feldman et al., Cell 2018 (GEO: GSE120575)
Melanoma tumor biopsies profiled before / on checkpoint immunotherapy.

Goal:
  1. Load the CD45+ single-cell TPM matrix and per-cell clinical metadata.
  2. Cluster and annotate immune cell types.
  3. Quantify which populations expand / contract with (a) response and
     (b) treatment timepoint (Pre vs Post).
  4. Score two CD8 T-cell programs (memory-like "CD8_G" vs exhausted "CD8_B")
     and derive a per-lesion signature that stratifies responders.
  5. Export figures, tables and a results summary.

Data referenced by name (not committed): the raw TPM matrix and metadata are
downloaded from GEO GSE120575; only derived tables/figures are kept.
"""
import os
import gzip
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sc.settings.verbosity = 1
np.random.seed(0)

DATA = "/home/user/ClaudeCode/data"
OUT = "/home/user/ClaudeCode/results"
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
sc.settings.figdir = FIG

TPM_FILE = os.path.join(DATA, "GSE120575_TPM.txt.gz")
META_FILE = os.path.join(DATA, "GSE120575_patient_ID_single_cells.txt.gz")

# ---------------------------------------------------------------------------
# Signature gene sets — curated from the marker genes reported by
# Sade-Feldman et al. 2018 for the two CD8 T-cell states.
#   CD8_G  = memory-like / "good" state, enriched in responding lesions
#   CD8_B  = exhausted / "bad" state, enriched in non-responding lesions
# ---------------------------------------------------------------------------
CD8_G_MEMORY = [
    "TCF7", "IL7R", "SELL", "CCR7", "CD28", "LEF1", "CD27",
    "GZMK", "S1PR1", "FOXP1", "BCL2", "NELL2", "CD55",
]
CD8_B_EXHAUSTION = [
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1", "CD38",
    "TOX", "TNFRSF9", "LAYN", "CXCL13", "GZMB", "PRF1", "HLA-DRA",
]

# canonical lineage markers for annotation
LINEAGE = {
    "CD8 T": ["CD8A", "CD8B"],
    "CD4 T": ["CD4", "IL7R", "CD3D"],
    "Treg": ["FOXP3", "IL2RA"],
    "B": ["MS4A1", "CD79A", "CD79B"],
    "Plasma": ["MZB1", "IGHG1", "SDC1"],
    "NK": ["NKG7", "KLRD1", "GNLY"],
    "Myeloid": ["LYZ", "CD68", "CD14", "ITGAX"],
    "pDC": ["LILRA4", "IL3RA"],
}


def load_metadata():
    """Parse the GEO sample table -> per-cell clinical annotation."""
    rows = []
    with gzip.open(META_FILE, "rt", encoding="latin-1") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            if not parts[0].startswith("Sample "):
                continue
            barcode, patient_id, response, therapy = parts[1], parts[4], parts[5], parts[6]
            if not barcode or response not in ("Responder", "Non-responder"):
                continue
            timepoint = patient_id.split("_")[0]  # Pre / Post
            patient = patient_id.split("_", 1)[1] if "_" in patient_id else patient_id
            rows.append((barcode, patient, timepoint, response, therapy))
    meta = pd.DataFrame(
        rows, columns=["cell", "patient", "timepoint", "response", "therapy"]
    ).set_index("cell")
    return meta


def load_expression():
    """Load the genes x cells TPM matrix into an AnnData (cells x genes)."""
    with gzip.open(TPM_FILE, "rt") as fh:
        cells = fh.readline().rstrip("\n").split("\t")[1:]
    print(f"[load] {len(cells)} cell columns in TPM header")
    # skip header row (cells) + patient-label row; gene names in col 0
    # keep col 0 (gene names) as object, data columns as float32 (memory-light)
    dtypes = {0: "object"}
    dtypes.update({i: "float32" for i in range(1, len(cells) + 1)})
    dtypes[len(cells) + 1] = "float32"  # trailing empty column from a stray tab
    df = pd.read_csv(
        TPM_FILE, sep="\t", skiprows=2, header=None, index_col=0, dtype=dtypes
    )
    if df.shape[1] == len(cells) + 1:
        df = df.iloc[:, : len(cells)]  # drop trailing all-NaN column
    df.columns = cells
    df.index.name = "gene"
    df = df[~df.index.duplicated(keep="first")]
    print(f"[load] matrix genes x cells = {df.shape}")
    # memory-careful construction: dense -> sparse, freeing copies as we go
    from scipy import sparse
    genes = df.index.astype(str).to_numpy()
    cell_names = list(df.columns)
    arr = df.to_numpy(dtype="float32")  # genes x cells
    del df
    X = sparse.csr_matrix(arr.T)  # cells x genes (TPM is mostly zeros)
    del arr
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame(index=cell_names),
        var=pd.DataFrame(index=genes),
    )
    return adata


CACHE = os.path.join(DATA, "adata_processed.h5ad")


def build_processed():
    meta = load_metadata()
    adata = load_expression()

    # attach metadata, keep only annotated cells
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    adata.obs = adata.obs.join(meta.loc[common])
    print(f"[merge] cells with clinical annotation: {adata.n_obs}")
    adata.obs["group"] = adata.obs["response"] + " / " + adata.obs["timepoint"]

    # ---- QC + gene filtering (drop genes seen in <3 cells to save memory) ----
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(1)).ravel()
    sc.pp.filter_genes(adata, min_cells=3)
    print("[qc] median genes/cell:", int(np.median(adata.obs["n_genes"])),
          "| genes after filter:", adata.n_vars)

    # ---- normalize (input is already TPM) -> log1p in place ----
    sc.pp.log1p(adata)

    # ---- HVG / scale / PCA / neighbors / clustering / UMAP ----
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata_hv = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata_hv, max_value=10)
    sc.tl.pca(adata_hv, n_comps=30)
    sc.pp.neighbors(adata_hv, n_neighbors=15, n_pcs=30)
    sc.tl.leiden(adata_hv, resolution=1.0, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(adata_hv)
    # carry results back
    adata.obs["leiden"] = adata_hv.obs["leiden"]
    adata.obsm["X_umap"] = adata_hv.obsm["X_umap"]
    adata.obsm["X_pca"] = adata_hv.obsm["X_pca"]
    print("[cluster] n clusters:", adata.obs["leiden"].nunique())

    # ---- annotate clusters by lineage marker mean expression ----
    present = {ct: [g for g in gs if g in adata.var_names] for ct, gs in LINEAGE.items()}
    score_df = pd.DataFrame(index=adata.obs_names)
    for ct, gs in present.items():
        if gs:
            score_df[ct] = np.asarray(adata[:, gs].X.mean(1)).ravel()
    cl_scores = score_df.groupby(adata.obs["leiden"]).mean()
    # z-score across clusters so lineages are comparable
    cl_z = (cl_scores - cl_scores.mean()) / (cl_scores.std() + 1e-9)
    cluster_label = cl_z.idxmax(1).to_dict()
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_label).astype(str)
    cl_scores.to_csv(os.path.join(OUT, "cluster_lineage_scores.csv"))
    print("[annotate] cluster -> cell_type:", cluster_label)

    # ---- signature scoring ----
    g_genes = [g for g in CD8_G_MEMORY if g in adata.var_names]
    b_genes = [g for g in CD8_B_EXHAUSTION if g in adata.var_names]
    sc.tl.score_genes(adata, g_genes, score_name="CD8_G_memory")
    sc.tl.score_genes(adata, b_genes, score_name="CD8_B_exhaustion")
    print(f"[score] CD8_G genes used ({len(g_genes)}): {g_genes}")
    print(f"[score] CD8_B genes used ({len(b_genes)}): {b_genes}")
    adata.write(CACHE)  # cache the heavy compute for fast re-runs
    print("[cache] wrote", CACHE)
    return adata


def main():
    if os.path.exists(CACHE):
        print("[cache] loading", CACHE)
        adata = sc.read_h5ad(CACHE)
    else:
        adata = build_processed()

    # =====================================================================
    # 1) Population expand/contract:  fractions by response and timepoint
    # =====================================================================
    from scipy.stats import mannwhitneyu

    def fractions(by):
        tab = (
            adata.obs.groupby([by, "cell_type"], observed=True)
            .size().unstack(fill_value=0)
        )
        return tab.div(tab.sum(1), axis=0)

    frac_resp = fractions("response")
    frac_time = fractions("timepoint")
    frac_resp.to_csv(os.path.join(OUT, "fraction_by_response.csv"))
    frac_time.to_csv(os.path.join(OUT, "fraction_by_timepoint.csv"))

    # per-lesion composition (one row per patient x timepoint)
    per_lesion = (
        adata.obs.groupby(["patient", "timepoint", "response", "cell_type"],
                          observed=True)
        .size().unstack(fill_value=0)
    )
    per_lesion = per_lesion.div(per_lesion.sum(1), axis=0).reset_index()
    cell_types = sorted(adata.obs["cell_type"].unique())

    def enrichment(df, col, a_label, b_label, tag):
        """Per-lesion fraction contrast a_label vs b_label for each cell type."""
        rows = []
        for ct in cell_types:
            a = df.loc[df[col] == a_label, ct].dropna()
            b = df.loc[df[col] == b_label, ct].dropna()
            lfc = np.log2((a.mean() + 1e-3) / (b.mean() + 1e-3))
            try:
                p = mannwhitneyu(a, b, alternative="two-sided",
                                 nan_policy="omit").pvalue
            except ValueError:
                p = np.nan
            rows.append((ct, a.mean(), b.mean(), lfc, p))
        out = pd.DataFrame(
            rows, columns=["cell_type", f"mean_frac_{a_label}",
                           f"mean_frac_{b_label}", f"log2FC_{tag}",
                           "pval_mannwhitney"]
        ).sort_values(f"log2FC_{tag}", ascending=False)
        return out

    # (a) Responder vs Non-responder  (who responds)
    enrich_df = enrichment(per_lesion, "response", "Responder",
                           "Non-responder", "R_vs_NR")
    enrich_df.to_csv(os.path.join(OUT, "population_enrichment_R_vs_NR.csv"),
                     index=False)
    # (b) Post vs Pre  (what changes WITH treatment)
    enrich_time = enrichment(per_lesion, "timepoint", "Post", "Pre",
                             "Post_vs_Pre")
    enrich_time.to_csv(os.path.join(OUT, "population_enrichment_Post_vs_Pre.csv"),
                       index=False)
    # (c) Post vs Pre restricted to responders (treatment effect in responders)
    enrich_time_R = enrichment(per_lesion[per_lesion.response == "Responder"],
                               "timepoint", "Post", "Pre", "Post_vs_Pre")
    enrich_time_R.to_csv(
        os.path.join(OUT, "population_enrichment_Post_vs_Pre_responders.csv"),
        index=False)

    print("[fractions] R vs NR:\n", enrich_df.to_string(index=False))
    print("[fractions] Post vs Pre (all):\n", enrich_time.to_string(index=False))

    # =====================================================================
    # 2) Marker genes per cluster (Wilcoxon)
    # =====================================================================
    sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon", n_genes=30)
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(os.path.join(OUT, "cluster_marker_genes.csv"), index=False)

    # =====================================================================
    # 3) CD8 signature -> per-lesion memory/exhaustion ratio -> stratifier
    # =====================================================================
    cd8 = adata[adata.obs["cell_type"] == "CD8 T"].copy()
    print(f"[cd8] CD8 T cells: {cd8.n_obs}")
    cd8.obs["mem_minus_exh"] = cd8.obs["CD8_G_memory"] - cd8.obs["CD8_B_exhaustion"]

    lesion = (
        cd8.obs.groupby(["patient", "timepoint", "response"], observed=True)
        .agg(n_cd8=("CD8_G_memory", "size"),
             CD8_G=("CD8_G_memory", "mean"),
             CD8_B=("CD8_B_exhaustion", "mean"),
             score=("mem_minus_exh", "mean"))
        .reset_index()
    )
    lesion = lesion[lesion["n_cd8"] >= 10]  # require enough CD8 cells
    lesion["GB_ratio"] = lesion["CD8_G"] - lesion["CD8_B"]
    lesion.to_csv(os.path.join(OUT, "per_lesion_cd8_signature.csv"), index=False)

    # AUC of the signature separating R vs NR lesions
    from sklearn.metrics import roc_auc_score
    y = (lesion["response"] == "Responder").astype(int)
    auc = roc_auc_score(y, lesion["score"]) if y.nunique() == 2 else np.nan
    p_sig = mannwhitneyu(
        lesion.loc[y == 1, "score"], lesion.loc[y == 0, "score"],
        alternative="two-sided",
    ).pvalue
    print(f"[stratify] per-lesion signature AUC (R vs NR) = {auc:.3f}, "
          f"Mann-Whitney p = {p_sig:.2e}, n_lesions = {len(lesion)}")

    # single-marker TCF7 fraction per lesion as a portable proxy
    if "TCF7" in cd8.var_names:
        # log1p preserves zeros, so X>0 == TPM>0 (cell expresses TCF7)
        tcf7_x = cd8[:, "TCF7"].X
        tcf7_x = tcf7_x.toarray().ravel() if hasattr(tcf7_x, "toarray") else np.asarray(tcf7_x).ravel()
        cd8.obs["TCF7_pos"] = tcf7_x > 0
        tcf7 = cd8.obs.groupby(["patient", "timepoint", "response"],
                               observed=True)["TCF7_pos"].mean()
        tcf7 = tcf7.reset_index().rename(columns={"TCF7_pos": "TCF7_frac"})
        tcf7 = tcf7.merge(lesion[["patient", "timepoint", "n_cd8"]],
                          on=["patient", "timepoint"], how="inner")
        tcf7 = tcf7.dropna(subset=["TCF7_frac"])
        tcf7.to_csv(os.path.join(OUT, "per_lesion_TCF7_fraction.csv"), index=False)
        y2 = (tcf7["response"] == "Responder").astype(int)
        auc_tcf7 = roc_auc_score(y2, tcf7["TCF7_frac"]) if y2.nunique() == 2 else np.nan
        print(f"[stratify] per-lesion TCF7+ fraction AUC (R vs NR) = {auc_tcf7:.3f}")
    else:
        auc_tcf7 = np.nan

    # =====================================================================
    # Figures
    # =====================================================================
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", size=8,
               title="Immune cell types (GSE120575)", show=False,
               save="_celltypes.png")
    sc.pl.umap(adata, color=["CD8_G_memory", "CD8_B_exhaustion"], size=8,
               show=False, save="_cd8_signatures.png")
    sc.pl.umap(adata, color="response", size=8, show=False, save="_response.png")

    # stacked-bar of composition by group
    comp = (adata.obs.groupby(["group", "cell_type"], observed=True).size()
            .unstack(fill_value=0))
    comp = comp.div(comp.sum(1), axis=0)
    ax = comp.plot(kind="bar", stacked=True, figsize=(9, 5), colormap="tab20")
    ax.set_ylabel("fraction of cells")
    ax.set_title("Immune composition by response / timepoint")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "composition_by_group.png"), dpi=130)
    plt.close()

    # per-lesion signature boxplot R vs NR
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    data = [lesion.loc[y == 1, "score"], lesion.loc[y == 0, "score"]]
    ax.boxplot(data)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Responder", "Non-resp."])
    for i, d in enumerate(data, 1):
        ax.scatter(np.random.normal(i, 0.05, len(d)), d, s=14, alpha=0.6)
    ax.set_ylabel("CD8 memory − exhaustion (per lesion)")
    ax.set_title(f"Responder signature\nAUC={auc:.2f}  p={p_sig:.1e}")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, "signature_boxplot.png"), dpi=130)
    plt.close()

    # ---- write summary ----
    with open(os.path.join(OUT, "SUMMARY.txt"), "w") as fh:
        fh.write("GSE120575 reanalysis summary\n")
        fh.write("=" * 60 + "\n")
        fh.write(f"cells analyzed: {adata.n_obs}\n")
        fh.write(f"clusters: {adata.obs['leiden'].nunique()}  "
                 f"cell types: {sorted(adata.obs['cell_type'].unique())}\n")
        fh.write(f"CD8 T cells: {cd8.n_obs}\n")
        fh.write(f"per-lesion signature AUC (R vs NR): {auc:.3f}  p={p_sig:.2e}  "
                 f"n_lesions={len(lesion)}\n")
        fh.write(f"per-lesion TCF7+ fraction AUC: {auc_tcf7:.3f}\n\n")
        fh.write("Population enrichment (R vs NR), log2FC>0 = expanded in responders:\n")
        fh.write(enrich_df.to_string(index=False) + "\n\n")
        fh.write("Population change WITH treatment (Post vs Pre, all lesions),\n"
                 "log2FC>0 = expanded after therapy:\n")
        fh.write(enrich_time.to_string(index=False) + "\n\n")
        fh.write("Population change WITH treatment in responders only (Post vs Pre):\n")
        fh.write(enrich_time_R.to_string(index=False) + "\n")

    print("[done] outputs in", OUT)


if __name__ == "__main__":
    main()
