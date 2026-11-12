"""AlphaGenome v2 demo: cell embeddings and expression conditioning.

`OutputType.CELL_EMBEDDINGS` cannot be requested from the predict_* calls (the
server rejects it); embeddings come from `model.encode_cells(adata)`, which
runs only the single-cell autoencoder on raw UMI counts.

This script:
  1. downloads the public 10x PBMC 3k dataset (raw counts, Ensembl gene IDs),
  2. encodes all 2,700 cells into 256-d embeddings,
  3. checks that the embeddings separate marker-defined cell types,
  4. uses groups of real cells as `cell_conditioning=` to predict the effect of
     chr22:36201698 A>C on denoised single-cell gene expression.

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_cells.py [data_dir]
"""

import io
import os
import pathlib
import sys
import tarfile
import urllib.request

import anndata
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

OUT_DIR = pathlib.Path(__file__).parent
DATA_DIR = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DIR / "data"
PBMC_URL = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
MTX_DIR = DATA_DIR / "filtered_gene_bc_matrices" / "hg19"

# ---------------------------------------------------------------- load data
if not MTX_DIR.exists():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(PBMC_URL) as resp:
        tarfile.open(fileobj=io.BytesIO(resp.read())).extractall(DATA_DIR, filter="data")

counts = scipy.sparse.csr_matrix(scipy.io.mmread(MTX_DIR / "matrix.mtx").T)
genes = pd.read_csv(MTX_DIR / "genes.tsv", sep="\t", header=None, names=["gene_id", "symbol"])
adata = anndata.AnnData(counts, var=genes.set_index("gene_id"))
print("PBMC 3k:", adata.shape, "(cells x genes)")

# Approximate cell-type labels from canonical markers (no clustering).
MARKERS = {
    "T cell": ["CD3D", "CD3E", "IL7R"],
    "B cell": ["MS4A1", "CD79A", "CD79B"],
    "Monocyte": ["LYZ", "CD14", "S100A8", "FCGR3A"],
    "NK cell": ["GNLY", "NKG7", "KLRD1"],
    "Platelet": ["PPBP", "PF4"],
}
lib = np.asarray(counts.sum(axis=1)).ravel()
symbol_to_col = {s: i for i, s in enumerate(genes.symbol)}
scores = {}
for name, syms in MARKERS.items():
    cols = [symbol_to_col[s] for s in syms]
    x = np.log1p(counts[:, cols].toarray() / lib[:, None] * 1e4).mean(axis=1)
    scores[name] = (x - x.mean()) / x.std()
scores = pd.DataFrame(scores)
labels = scores.idxmax(axis=1).where(scores.max(axis=1) > 0.5, "Unassigned").to_numpy()
print("Marker-based labels:", pd.Series(labels).value_counts().to_dict())

# ------------------------------------------------------------- embeddings
model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
emb = model.encode_cells(adata, progress_bar=False)
print("Embeddings:", emb.shape, emb.dtype)
np.save(OUT_DIR / "v2_pbmc_embeddings.npy", emb)


def knn_purity(x, labels, k=15):
    """Fraction of each cell's k nearest neighbours sharing its label."""
    x = x / np.linalg.norm(x, axis=1, keepdims=True)
    sim = x @ x.T
    np.fill_diagonal(sim, -np.inf)
    nn = np.argpartition(-sim, k, axis=1)[:, :k]
    return (labels[nn] == labels[:, None]).mean(axis=1)


def pca(x, n):
    x = x - x.mean(axis=0)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    return u[:, :n] * s[:n]


keep = labels != "Unassigned"
# Baseline: 50 PCs of log-normalised counts over the 2,000 most variable genes.
lognorm = np.log1p(counts.toarray() / lib[:, None] * 1e4)
hv = np.argsort(-lognorm.var(axis=0))[:2000]
baseline = pca(lognorm[:, hv], 50)
rng = np.random.default_rng(0)
purity = pd.DataFrame(
    {
        "AlphaGenome embedding": knn_purity(emb[keep], labels[keep]),
        "log-norm PCA (50 PCs)": knn_purity(baseline[keep], labels[keep]),
        "shuffled labels": knn_purity(emb[keep], rng.permutation(labels[keep])),
        "label": labels[keep],
    }
)
print("\n15-NN label purity (cosine), mean per marker-defined type:")
by_type = purity.groupby("label").mean()
by_type.loc["ALL CELLS"] = purity.drop(columns="label").mean()
print(by_type.round(3).to_string())

# ------------------------------------------- expression-conditioned variant
interval = genome.Interval(chromosome="chr22", start=35677410, end=36725986)
variant = genome.Variant(
    chromosome="chr22", position=36201698, reference_bases="A", alternate_bases="C"
)
rows = []
for name in ["Monocyte", "T cell", "B cell", "NK cell"]:
    idx = np.flatnonzero(labels == name)
    idx = rng.choice(idx, size=min(len(idx), dna_client.MAX_CONDITIONING_CELLS), replace=False)
    out = model.predict_variant(
        interval,
        variant,
        requested_outputs=[dna_client.OutputType.SC_GENE_EXPRESSION_DENOISED],
        cell_conditioning=adata[np.sort(idx)],
        cell_condition_name=name,
    )
    ref, alt = out.reference.sc_gene_expression_denoised, out.alternate.sc_gene_expression_denoised
    if not rows:
        print("\nSC_GENE_EXPRESSION_DENOISED:", ref.values.shape, "(genes x tracks)")
        print(ref.track_metadata.to_string())
    # Gene-level CPM: sum the per-strand channels.
    r, a = ref.values.sum(axis=1), alt.values.sum(axis=1)
    for g, gene in enumerate(ref.gene_metadata["gene_name"]):
        rows.append({"cells": name, "n_cells": len(idx), "gene": gene, "ref_cpm": float(r[g]), "alt_cpm": float(a[g])})
den = pd.DataFrame(rows)
den["log2fc"] = np.log2((den.alt_cpm + 0.1) / (den.ref_cpm + 0.1))
den.to_csv(OUT_DIR / "v2_pbmc_conditioned_variant.csv", index=False)
show = ["APOL4", "APOL3", "APOL1", "APOL2", "MYH9"]
print("\nVariant effect conditioned on real PBMC cells (denoised CPM):")
print(den[den.gene.isin(show)].pivot(index="gene", columns="cells", values="ref_cpm").loc[show].round(2).to_string())
print("log2(ALT/REF):")
print(den[den.gene.isin(show)].pivot(index="gene", columns="cells", values="log2fc").loc[show].round(2).to_string())

# ------------------------------------------------------------------ figure
pcs = pca(emb, 2)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
palette = dict(zip(list(MARKERS) + ["Unassigned"], ["tab:blue", "tab:orange", "tab:red", "tab:green", "tab:purple", "lightgrey"]))
for name, color in palette.items():
    m = labels == name
    axes[0].scatter(pcs[m, 0], pcs[m, 1], s=4, color=color, label=f"{name} ({m.sum()})")
axes[0].legend(markerscale=3, fontsize=8)
axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")
axes[0].set_title("AlphaGenome v2 cell embeddings (PCA), PBMC 3k")
plot = by_type.drop(index="ALL CELLS")
plot.plot.bar(ax=axes[1], color=["tab:red", "dimgrey", "lightgrey"])
axes[1].set_ylabel("15-NN label purity")
axes[1].set_xlabel("")
axes[1].set_title("Do neighbours share the marker-based label?")
axes[1].tick_params(axis="x", rotation=30)
fig.tight_layout()
out_png = OUT_DIR / "v2_cell_embeddings.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
