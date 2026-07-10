# Computing a signature from raw single-cell perturbation data

Source: the single-cell best-practices book, *Differential gene expression*
chapter. This is the ✅-rigorous path from counts to a `Signature`.

## The one rule: pseudobulk, don't test cells

Cells from the same sample/donor are **not independent replicates**
(pseudoreplication). A per-cell Wilcoxon or t-test treats thousands of
correlated cells as thousands of samples and **inflates false positives**. The
sample-level (pseudobulk) approach has been shown to outperform cell-level tests
for scRNA-seq. So:

> aggregate counts to **(sample × cell-type/condition)** pseudobulk, then run a
> **bulk DE method — DESeq2 / PyDESeq2 / edgeR / limma-voom**.

`MAST` with a per-donor random effect is the only cell-level method the book
keeps on the table; plain Wilcoxon/t on cells is **not recommended**.

## Workflow (maps 1:1 to `perturbomics.pseudobulk`)

```python
import scanpy as sc, decoupler as dc
import pertpy as pt   # wraps PyDESeq2 / edgeR

# 1. QC on raw counts; keep counts in a layer
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()

# 2. Pseudobulk: SUM counts within sample × group  (the unit of replication)
pb = dc.pp.pseudobulk(adata, sample_col="sample", groups_col="cell_type",
                      layer="counts", mode="sum")
dc.pp.filter_samples(pb, min_cells=10, min_counts=1000)

# 3. Look before modelling: PCA on pseudobulk reveals covariates for the design
sc.pp.normalize_total(pb, target_sum=1e6); sc.pp.log1p(pb); sc.tl.pca(pb)

# 4. Gene filter PER cell type (different populations express different genes)
dc.pp.filter_by_expr(pb, group="condition", min_count=10, min_total_count=15)

# 5. DE with PyDESeq2 (recommended engine)
ds = pt.tl.PyDESeq2(adata=pb, design="~ condition")
ds.fit()
res = ds.test_contrasts(ds.contrast("condition", baseline="ctrl",
                                    group_to_compare="perturbed"))
```

`res` has `log2FoldChange`, `stat` (Wald), `pvalue`, `padj`, `baseMean`.

## Turn it into a Signature

```python
from perturbomics import Signature
sig = Signature.from_deseq2(res, stat_col="stat",
                            name="perturbed_vs_ctrl", modality="scrna_dge")
```

Use the **Wald `stat`**, not raw `log2FoldChange`: it already folds effect size
*and* its uncertainty into one signed, ranked number — exactly the ranking a
connectivity score wants. (`Signature.from_deseq2` defaults to `stat`.)

The package's `signature_from_pseudobulk(pb, design, "perturbed", "ctrl")` does
steps 2→5→Signature in one call, using **PyDESeq2 when installed** and a
transparent **CPM-log + Welch-t fallback** otherwise (fine for the demo /
teaching; install `pydeseq2` for real data).

## Complex designs (perturbation × cell type)

To ask whether a perturbation's effect **differs by cell type** — the
interaction — fit `~ cell_type * condition` and test the interaction contrast
(the book's example subtracts the (stim−ctrl) effect of one cell type from
another). Do this when your combination hypothesis is cell-type-specific.

## Pitfalls (all ✅ deterministic to avoid)

- Filter lowly-expressed genes **separately per cell type** — populations
  express distinct gene sets.
- Put **every major covariate** (batch, donor sex, etc.) surfaced by the PCA
  into the design matrix, or confounding leaks into the signature.
- Keep **raw counts** for DE (DESeq2 models counts); never feed log/normalised
  values to PyDESeq2.
