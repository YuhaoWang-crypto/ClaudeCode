#!/usr/bin/env python3
"""
Part A: fine sub-clustering of CD8 T cells (GSE120575) into transcriptional states.
Part B: cross-dataset validation of the CD8_G/CD8_B gene sets on Yost et al.
        (GSE123813, basal/squamous cell carcinoma, anti-PD-1) single cells.
Part C: TCR / clonal dynamics in Yost — clonal expansion by CD8 state and
        pre vs post anti-PD-1 clonal replacement.

Data referenced by name (not committed): GSE120575 (cache), GSE123813 counts/
metadata/tcr. Yost signature genes + per-cell library sizes were pre-extracted
into data/yost_sig_genes.tsv and data/yost_libsize.tsv.
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
np.random.seed(0)
DATA = "/home/user/ClaudeCode/data"
OUT = "/home/user/ClaudeCode/results"
FIG = os.path.join(OUT, "figures")
sc.settings.figdir = FIG

MEM = ["TCF7", "IL7R", "SELL", "CCR7", "CD28", "LEF1", "CD27", "GZMK"]
EXH = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "ENTPD1", "CD38", "TOX"]

STATE_MARKERS = {
    "Naive/memory": ["TCF7", "CCR7", "SELL", "IL7R", "LEF1"],
    "Effector-memory": ["GZMK", "CD27", "CXCR3"],
    "Cytotoxic/effector": ["GZMB", "GNLY", "NKG7", "FGFBP2", "PRF1"],
    "Exhausted": ["PDCD1", "HAVCR2", "LAG3", "TOX", "ENTPD1", "TIGIT"],
    "Proliferating": ["MKI67", "TOP2A", "STMN1"],
}


# ===========================================================================
# Part A — CD8 sub-states in GSE120575
# ===========================================================================
def cd8_substates():
    adata = sc.read_h5ad(os.path.join(DATA, "adata_processed.h5ad"))
    cd8 = adata[adata.obs["cell_type"] == "CD8 T"].copy()
    print(f"[A] CD8 T cells: {cd8.n_obs}")

    sc.pp.highly_variable_genes(cd8, n_top_genes=1500)
    sub = cd8[:, cd8.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=20)
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=20)
    sc.tl.leiden(sub, resolution=0.6, key_added="cd8_leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(sub)
    cd8.obs["cd8_leiden"] = sub.obs["cd8_leiden"]
    cd8.obsm["X_umap"] = sub.obsm["X_umap"]

    # annotate each cluster by the best-matching state module (z across clusters)
    mod = pd.DataFrame(index=cd8.obs_names)
    for st, gs in STATE_MARKERS.items():
        gg = [g for g in gs if g in cd8.var_names]
        mod[st] = np.asarray(cd8[:, gg].X.mean(1)).ravel() if gg else 0.0
    clmean = mod.groupby(cd8.obs["cd8_leiden"], observed=True).mean()
    clz = (clmean - clmean.mean()) / (clmean.std() + 1e-9)
    lab = clz.idxmax(1).to_dict()
    cd8.obs["cd8_state"] = cd8.obs["cd8_leiden"].map(lab).astype(str)
    print("[A] cluster->state:", lab)

    # per-state signature scores + response enrichment
    summ = cd8.obs.groupby("cd8_state", observed=True).agg(
        n=("CD8_G_memory", "size"),
        CD8_G=("CD8_G_memory", "mean"),
        CD8_B=("CD8_B_exhaustion", "mean"),
    )
    # fraction of each state within responder vs non-responder CD8 pools
    frac = (cd8.obs.groupby(["response", "cd8_state"], observed=True).size()
            .unstack(fill_value=0))
    frac = frac.div(frac.sum(1), axis=0).T
    summ = summ.join(frac)
    summ["log2FC_R_vs_NR"] = np.log2((summ["Responder"] + 1e-3) /
                                     (summ["Non-responder"] + 1e-3))
    summ.to_csv(os.path.join(OUT, "cd8_substate_summary.csv"))
    print("[A] per-state summary:\n", summ.round(3).to_string())

    # figures
    sc.pl.umap(cd8, color="cd8_state", size=12, title="CD8 T-cell states (GSE120575)",
               show=False, save="_cd8_states.png")
    sc.pl.umap(cd8, color=["CD8_G_memory", "CD8_B_exhaustion"], size=12,
               show=False, save="_cd8_state_scores.png")

    fig, ax = plt.subplots(figsize=(6.5, 4))
    summ_sorted = summ.sort_values("log2FC_R_vs_NR", ascending=False)
    colors = ["#4C78A8" if v > 0 else "#E45756" for v in summ_sorted["log2FC_R_vs_NR"]]
    ax.barh(summ_sorted.index, summ_sorted["log2FC_R_vs_NR"], color=colors)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("log2 fold-change of state fraction (Responder / Non-responder)")
    ax.set_title("CD8 states enriched in responders (blue) vs non-responders (red)")
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "cd8_state_enrichment.png"), dpi=130)
    plt.close()
    return summ


# ===========================================================================
# Part B — validate CD8_G / CD8_B gene sets on Yost CD8 states
# ===========================================================================
def yost_geneset_validation():
    expr = pd.read_csv(os.path.join(DATA, "yost_sig_genes.tsv"), sep="\t", index_col=0)
    lib = np.loadtxt(os.path.join(DATA, "yost_libsize.tsv"))
    assert expr.shape[1] == len(lib), f"align mismatch {expr.shape[1]} vs {len(lib)}"
    # CP10k normalize + log1p, then z-score each gene across cells
    norm = np.log1p(expr.div(lib, axis=1) * 1e4)
    z = norm.sub(norm.mean(1), axis=0).div(norm.std(1) + 1e-9, axis=0)

    meta = pd.read_csv(os.path.join(DATA, "GSE123813_bcc_tcell_metadata.txt.gz"),
                       sep="\t", index_col=0)
    mem_g = [g for g in MEM if g in z.index]
    exh_g = [g for g in EXH if g in z.index]
    cd8g = z.loc[mem_g].mean(0)
    cd8b = z.loc[exh_g].mean(0)
    df = meta.join(pd.DataFrame({"CD8_G": cd8g, "CD8_B": cd8b}))
    cd8_states = ["Naive", "CD8_mem", "CD8_eff", "CD8_act", "CD8_ex_act", "CD8_ex"]
    df = df[df["cluster"].isin(cd8_states)]
    summ = df.groupby("cluster", observed=True)[["CD8_G", "CD8_B"]].mean()
    summ = summ.reindex([c for c in cd8_states if c in summ.index])
    summ.to_csv(os.path.join(OUT, "yost_geneset_by_state.csv"))
    print("[B] Yost CD8 states, mean gene-set z-score:\n", summ.round(3).to_string())

    fig, ax = plt.subplots(figsize=(6.5, 4))
    summ.plot(kind="bar", ax=ax, color={"CD8_G": "#4C78A8", "CD8_B": "#E45756"})
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel("mean signature z-score")
    ax.set_title("Gene sets generalize to Yost anti-PD-1 CD8 states (GSE123813)")
    ax.legend(["CD8_G memory", "CD8_B exhaustion"])
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "yost_geneset_validation.png"), dpi=130)
    plt.close()
    return summ


# ===========================================================================
# Part C — TCR / clonal dynamics in Yost
# ===========================================================================
def yost_tcr_dynamics():
    tcr = pd.read_csv(os.path.join(DATA, "GSE123813_bcc_tcr.txt.gz"),
                      sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(DATA, "GSE123813_bcc_tcell_metadata.txt.gz"),
                       sep="\t", index_col=0)
    d = meta.join(tcr["cdr3s_aa"].rename("clonotype")).dropna(subset=["clonotype"])
    print(f"[C] T cells with a TCR: {len(d)}")

    # clone = (patient, clonotype); clone size across that patient's cells
    d["clone_id"] = d["patient"].astype(str) + "|" + d["clonotype"].astype(str)
    size = d.groupby("clone_id").size()
    d["clone_size"] = d["clone_id"].map(size)
    d["expanded"] = d["clone_size"] >= 2

    cd8_states = ["Naive", "CD8_mem", "CD8_eff", "CD8_act", "CD8_ex_act", "CD8_ex"]
    dd = d[d["cluster"].isin(cd8_states)]

    # (1) clonal expansion by state
    by_state = dd.groupby("cluster", observed=True)["expanded"].mean().reindex(
        [c for c in cd8_states if c in dd["cluster"].unique()])
    # (2) expansion pre vs post by state
    by_state_time = (dd.groupby(["cluster", "treatment"], observed=True)["expanded"]
                     .mean().unstack())
    tab = pd.concat([by_state.rename("expanded_frac"), by_state_time], axis=1)
    tab.to_csv(os.path.join(OUT, "yost_clonal_expansion_by_state.csv"))
    print("[C] fraction of cells in expanded clones (size>=2):\n", tab.round(3).to_string())

    # (3) clonal replacement: among POST expanded clonotypes, fraction NOT seen pre
    rep_rows = []
    for pat, g in d.groupby("patient"):
        pre = set(g.loc[g.treatment == "pre", "clonotype"])
        post = g.loc[g.treatment == "post"]
        post_exp = set(post.loc[post.clone_size >= 2, "clonotype"])
        if not post_exp:
            continue
        novel = len(post_exp - pre) / len(post_exp)
        rep_rows.append((pat, len(post_exp), round(novel, 3)))
    rep = pd.DataFrame(rep_rows, columns=["patient", "n_post_expanded_clones",
                                          "frac_novel_post"])
    rep.to_csv(os.path.join(OUT, "yost_clonal_replacement.csv"), index=False)
    print("[C] clonal replacement (novel fraction of post expanded clones):\n",
          rep.to_string(index=False))
    print(f"[C] mean novel fraction across patients = {rep['frac_novel_post'].mean():.2f}")

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    order = [c for c in cd8_states if c in by_state.index]
    axes[0].bar(order, by_state.reindex(order).values,
                color=["#54A24B" if "ex" in c else "#4C78A8" for c in order])
    axes[0].set_ylabel("fraction of cells in expanded clones")
    axes[0].set_title("Clonal expansion by CD8 state")
    axes[0].tick_params(axis="x", rotation=45)

    bt = by_state_time.reindex(order)
    x = np.arange(len(order)); w = 0.38
    axes[1].bar(x - w/2, bt["pre"].values, w, label="pre", color="#B0B0B0")
    axes[1].bar(x + w/2, bt["post"].values, w, label="post", color="#E45756")
    axes[1].set_xticks(x); axes[1].set_xticklabels(order, rotation=45)
    axes[1].set_ylabel("fraction in expanded clones")
    axes[1].set_title(f"Pre vs post anti-PD-1 (mean {rep['frac_novel_post'].mean():.0%} "
                      "of post clones are novel)")
    axes[1].legend()
    plt.tight_layout(); plt.savefig(os.path.join(FIG, "yost_tcr_dynamics.png"), dpi=130)
    plt.close()
    return tab, rep


if __name__ == "__main__":
    print("==== Part A: GSE120575 CD8 sub-states ====")
    cd8_substates()
    print("\n==== Part B: Yost gene-set validation ====")
    yost_geneset_validation()
    print("\n==== Part C: Yost TCR clonal dynamics ====")
    yost_tcr_dynamics()
    print("\n[done]")
