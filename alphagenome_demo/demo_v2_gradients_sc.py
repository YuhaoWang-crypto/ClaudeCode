"""AlphaGenome v2 demo: input gradients of SC_GENE_EXPRESSION_DENOISED.

Combines two v2 features: expression conditioning on real cells
(`cell_conditioning=`) and input gradients (`predict_interval_gradients`).
The question: for the same gene and the same DNA, does the model attribute
expression to different bases depending on which cells it is conditioned on?

Cells: 10x PBMC 3k, grouped by marker genes into monocytes, B cells, T cells.
Genes: LYZ (monocyte marker), MS4A1 (B cell marker), CD3E (T cell marker).

Requires ALPHAGENOME_API_KEY.

    python alphagenome_demo/demo_v2_gradients_sc.py [data_dir]
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
from alphagenome.models.v2 import dna_client

OUT_DIR = pathlib.Path(__file__).parent
DATA_DIR = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else OUT_DIR / "data"
PBMC_URL = "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
MTX_DIR = DATA_DIR / "filtered_gene_bc_matrices" / "hg19"
pd.set_option("display.width", 200)

# ---------------------------------------------------------------- cells
if not MTX_DIR.exists():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(PBMC_URL) as resp:
        tarfile.open(fileobj=io.BytesIO(resp.read())).extractall(DATA_DIR, filter="data")
counts = scipy.sparse.csr_matrix(scipy.io.mmread(MTX_DIR / "matrix.mtx").T)
genes = pd.read_csv(MTX_DIR / "genes.tsv", sep="\t", header=None, names=["gene_id", "symbol"])
adata = anndata.AnnData(counts, var=genes.set_index("gene_id"))

MARKERS = {
    "T cell": ["CD3D", "CD3E", "IL7R"],
    "B cell": ["MS4A1", "CD79A", "CD79B"],
    "Monocyte": ["LYZ", "CD14", "S100A8", "FCGR3A"],
    "NK cell": ["GNLY", "NKG7", "KLRD1"],
    "Platelet": ["PPBP", "PF4"],
}
lib = np.asarray(counts.sum(axis=1)).ravel()
col = {s: i for i, s in enumerate(genes.symbol)}
scores = {}
for name, syms in MARKERS.items():
    x = np.log1p(counts[:, [col[s] for s in syms]].toarray() / lib[:, None] * 1e4).mean(axis=1)
    scores[name] = (x - x.mean()) / x.std()
scores = pd.DataFrame(scores)
labels = scores.idxmax(axis=1).where(scores.max(axis=1) > 0.5, "Unassigned").to_numpy()

rng = np.random.default_rng(0)
GROUPS = ["Monocyte", "B cell", "T cell"]
cells = {}
for g in GROUPS:
    idx = np.flatnonzero(labels == g)
    idx = np.sort(rng.choice(idx, size=min(len(idx), dna_client.MAX_CONDITIONING_CELLS), replace=False))
    cells[g] = adata[idx]
print("Conditioning cells:", {g: a.n_obs for g, a in cells.items()})

TARGETS = {"LYZ": "ENSG00000090382", "MS4A1": "ENSG00000156738", "CD3E": "ENSG00000198851"}
# Observed mean expression in the conditioning cells (CPM), for reference.
obs = pd.DataFrame(
    {g: {s: float(a.X[:, a.var_names.get_loc(TARGETS[s])].sum() / a.X.sum() * 1e6) for s in TARGETS} for g, a in cells.items()}
)
print("\nMeasured expression in the conditioning cells (CPM, pooled counts):")
print(obs.round(1).to_string())

# ------------------------------------------------------------- gradients
model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"], timeout=60)
results = {}
for symbol, gene_id in TARGETS.items():
    for g in GROUPS:
        out = model.predict_interval_gradients(
            gene_id,
            requested_output=dna_client.OutputType.SC_GENE_EXPRESSION_DENOISED,
            sequence_length=dna_client.SEQUENCE_LENGTH_100KB,
            cell_conditioning=cells[g],
            cell_condition_name=g,
        )
        results[symbol, g] = out
    first = results[symbol, GROUPS[0]]
    print(f"done: {symbol} {first.target_gene.interval} window {first.interval}")
first = results["LYZ", "Monocyte"]
print("Gradients shape:", first.gradients.shape, "| objective tracks:", list(first.track_metadata.name))

score = pd.DataFrame({s: {g: results[s, g].scalar_score for g in GROUPS} for s in TARGETS}).T
print("\nscalar_score = log predicted expression of the gene under each conditioning:")
print(score.round(2).to_string())
print("Highest-scoring conditioning per gene:", score.idxmax(axis=1).to_dict())

# Per-position attribution, and how it compares between conditionings.
attr = {k: v.gradients.sum(axis=1) for k, v in results.items()}
BIN = 100
rows = []
for symbol in TARGETS:
    gene = results[symbol, GROUPS[0]].target_gene.interval
    start = results[symbol, GROUPS[0]].interval.start
    pos = np.arange(start, start + len(attr[symbol, GROUPS[0]]))
    tss = gene.end if gene.strand == "-" else gene.start
    for g in GROUPS:
        a = np.abs(attr[symbol, g])
        rows.append({
            "gene": symbol, "cells": g, "total_abs_attr": a.sum(),
            "share_within_1kb_of_tss": a[np.abs(pos - tss) <= 1000].sum() / a.sum(),
            "share_in_gene_body": a[(pos >= gene.start) & (pos < gene.end)].sum() / a.sum(),
            "top_100bp_bin_offset_from_tss": int((np.argmax(a[: len(a) // BIN * BIN].reshape(-1, BIN).sum(axis=1)) + 0.5) * BIN + start - tss),
        })
summary = pd.DataFrame(rows)
print("\nWhere the attribution sits:")
print(summary.round(3).to_string(index=False))

print("\nCorrelation of per-base attribution between conditionings (same gene, same DNA):")
for symbol in TARGETS:
    c = pd.DataFrame({g: attr[symbol, g] for g in GROUPS}).corr()
    print(f"  {symbol}: " + ", ".join(f"{a} vs {b} r={c.loc[a, b]:.3f}" for a, b in [(GROUPS[0], GROUPS[1]), (GROUPS[0], GROUPS[2]), (GROUPS[1], GROUPS[2])]))

# Cell-type-specific hotspots: bins where |attribution| differs most.
print("\nBins (100 bp) with the largest difference in |attribution| between the expressing cell type and the others:")
OWNER = {"LYZ": "Monocyte", "MS4A1": "B cell", "CD3E": "T cell"}
diff_rows = []
for symbol, owner in OWNER.items():
    gene = results[symbol, owner].target_gene.interval
    start = results[symbol, owner].interval.start
    tss = gene.end if gene.strand == "-" else gene.start
    n = len(attr[symbol, owner]) // BIN * BIN
    b = {g: np.abs(attr[symbol, g][:n]).reshape(-1, BIN).sum(axis=1) for g in GROUPS}
    others = np.mean([b[g] for g in GROUPS if g != owner], axis=0)
    d = b[owner] - others
    for i in np.argsort(-d)[:3]:
        diff_rows.append({"gene": symbol, "owner": owner, "bin_start": start + i * BIN, "offset_from_tss": int(start + (i + 0.5) * BIN - tss),
                          "owner_attr": b[owner][i], "others_attr": others[i]})
print(pd.DataFrame(diff_rows).round(3).to_string(index=False))
summary.to_csv(OUT_DIR / "v2_gradients_sc_summary.csv", index=False)

# ------------------------------------------------------------------ figure
colors = {"Monocyte": "tab:red", "B cell": "tab:orange", "T cell": "tab:blue"}
fig, axes = plt.subplots(len(TARGETS), 1, figsize=(14, 3.0 * len(TARGETS)))
SM = 200
for ax, symbol in zip(axes, TARGETS):
    gene = results[symbol, GROUPS[0]].target_gene.interval
    start = results[symbol, GROUPS[0]].interval.start
    n = len(attr[symbol, GROUPS[0]]) // SM * SM
    x = (start + (np.arange(n // SM) + 0.5) * SM) / 1e6
    for g in GROUPS:
        y = np.abs(attr[symbol, g][:n]).reshape(-1, SM).sum(axis=1)
        ax.plot(x, y, color=colors[g], lw=0.9, label=f"{g} (score {results[symbol, g].scalar_score:.2f})")
    ax.axvspan(gene.start / 1e6, gene.end / 1e6, color="grey", alpha=0.15)
    ax.set_ylabel(f"{symbol}\n|attribution| / {SM} bp")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlabel(f"{gene.chromosome} position (Mb); grey = gene body")
axes[0].set_title("Input gradients of denoised single-cell expression, conditioned on PBMC cell groups (AlphaGenome v2)")
fig.tight_layout()
out_png = OUT_DIR / "v2_gradients_sc.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"\nSaved {out_png}")
