---
name: single-cell-analysis
description: >-
  End-to-end single-cell (scRNA-seq / multiomic / spatial) analysis with
  benchmark-grounded method choices: QC → normalization → feature selection →
  embedding → clustering → cell-type annotation → integration → differential
  expression / compositional analysis → trajectories, plus a task→method map of
  the deep-learning single-cell literature. Use when analyzing or planning a
  single-cell dataset, annotating cell types (manual markers, SingleR,
  CellTypist, AUCell, scArches/scANVI label transfer, LLM annotation),
  correcting batch effects or mapping to a reference atlas, testing
  differential expression or cell-type abundance between conditions, choosing
  between a deep-learning method and a simple baseline, or reviewing someone
  else's single-cell pipeline for the standard failure modes. Distilled from
  sc-best-practices (Theis lab), OSCA (Bioconductor), and the OmicsML
  deep-learning-single-cell literature. Labels every default as
  ✅-benchmark-backed vs ⚠️-convention.
---

# Single-cell analysis: the benchmark-grounded spine

Three sources fused into one operating manual:

| Source | What it contributes | Ecosystem |
|---|---|---|
| [sc-best-practices](https://www.sc-best-practices.org/preamble.html) (Heumos/Theis et al., *Nat Rev Genet* 2023) | the workflow + which method to pick, chosen from external benchmarks out of >1700 tools | Python / scanpy / scverse |
| [OSCA](http://bioconductor.org/books/3.15/OSCA.basic/) (Amezquita/Lun/Marioni et al.) | rigorous annotation + marker statistics; the `SingleCellExperiment` world | R / Bioconductor |
| [awesome-deep-learning-single-cell-papers](https://github.com/OmicsML/awesome-deep-learning-single-cell-papers) (OmicsML; survey: Molho et al., *ACM TIST* 2024) | task → deep-learning method map, and where DL does *not* beat a linear baseline | mixed |

## The spine — run in this order

| # | Step | Default | Label | Detail |
|---|---|---|---|---|
| 1 | Raw processing | CellRanger / alevin-fry / STARsolo → count matrix | ⚠️ convention | `reference/preprocessing.md` |
| 2 | **QC** | **MAD-based** outliers (5 MAD counts/genes/top-20%, 3 MAD mito **+** hard mito > 8%) — *never* fixed cutoffs | ✅ | `reference/preprocessing.md` |
| 3 | Ambient RNA | SoupX with preliminary clusters | ✅ | ″ |
| 4 | Doublets | scDblFinder, **per batch, never on pooled batches** | ✅ | ″ |
| 5 | **Normalization** | shifted log for DR + DE; scran pooling when count depth varies hard; analytic Pearson residuals for gene selection / rare types | ✅ | ″ |
| 6 | Feature selection | deviance on **raw counts** (`scry`); 2000–4000 genes | ✅ | ″ |
| 7 | Embedding | PCA (≈50 PCs) → kNN graph → UMAP **for display only** | ✅ | ″ |
| 8 | **Clustering** | Leiden (not Louvain), **several resolutions** (0.25/0.5/1.0), `n_iterations=2` | ✅ | ″ |
| 9 | **Annotation** | manual markers **and** an automated method, then reconcile | ✅ | `reference/annotation.md` |
| 10 | Integration | only if batch effect is real; pick by scIB metrics, not by UMAP | ✅ | `reference/integration.md` |
| 11 | **DE between conditions** | **pseudobulk** + DESeq2/edgeR — not cell-level Wilcoxon | ✅ | `reference/conditions.md` |
| 12 | Abundance shifts | scCODA / Milo — not a proportion t-test | ✅ | ″ |
| 13 | Trajectories | DPT/Palantir; velocity only if the timescale is right | ⚠️ | `reference/trajectories.md` |

`reference/deep-learning-map.md` is the cross-cutting one: for any step above,
which deep-learning method exists, and whether it's worth it.

## The six decisions that actually change the answer

Everything else is bookkeeping. These are where analyses go wrong.

1. **QC by MAD, not by round numbers.** A fixed `n_genes > 500` filter silently
   deletes small resting lymphocytes and every low-RNA cell type. Use lenient
   MAD cutoffs, then check *what you deleted* by clustering the discarded cells.
2. **Cell-level DE across conditions is pseudoreplication.** Cells from one
   donor are not independent samples. Bulk methods on pseudobulk beat
   single-cell-specific DE tools in benchmarks; cell-level tests
   systematically mislabel highly-expressed genes as DE. Power comes from
   **more donors**, not more cells per donor.
3. **Cell-type proportions are compositional.** Sum-to-one forces negative
   correlations: one population genuinely expanding makes every other
   population "shrink". Use scCODA (with a declared reference cell type) or
   Milo on kNN neighbourhoods.
4. **Automated annotation is a starting point, never an endpoint.** Every
   annotation ships with marker-gene evidence and a confidence/uncertainty
   column; high-uncertainty cells become `Unknown`, not a guess.
5. **Integration is optional and destructive.** It can erase the biology you
   came for. Look at the unintegrated data first; if you integrate, judge it
   with scIB batch-vs-bio metrics, because *any* method makes a pretty UMAP.
6. **UMAP distances mean nothing.** Cluster adjacency, inter-cluster distance
   and velocity arrows on a UMAP are all projection artifacts. Never draw a
   biological conclusion from a 2-D layout alone.

## When deep learning actually wins

From the DL literature map + its own benchmarks — the honest version:

| Task | Use DL when | Otherwise use |
|---|---|---|
| Batch integration / atlas building | many datasets, heterogeneous protocols, labels available → **scANVI/scVI** | Harmony (fast, strong on simple tasks) |
| Reference mapping onto a big atlas | always → **scArches** (+ uncertainty) | kNN on shared PCA |
| Denoising / imputation | rarely; imputation inflates gene-gene correlation | just don't impute before DE |
| Cell-type annotation | cross-tissue, atlas-scale → **CellTypist**, scBERT-class models | markers + SingleR |
| Perturbation prediction | ⚠️ **benchmarks show DL does not yet beat simple linear baselines** on unseen-perturbation prediction | linear/additive baseline first, always |
| Spatial deconvolution / domains | yes → cell2location, DestVI, SpaGCN | NMF-based (SPOTlight) |
| Multiomic translation / integration | yes → GLUE, BABEL, Multigrate | matched-feature CCA |

**Always run the dumb baseline.** The single most cited negative result in this
field is that deep perturbation-response models fail to outperform linear
methods; the same pattern recurs in imputation and clustering benchmarks.
A DL method with no baseline comparison is an unvalidated claim.

## Ecosystem choice

- **Python/scverse** (scanpy, anndata, scvi-tools, squidpy, pertpy, decoupler):
  scale, deep-learning models, atlas mapping. Default for anything > 100k cells.
- **R/Bioconductor** (SingleCellExperiment, scran, scater, SingleR, AUCell,
  DESeq2/edgeR): the rigorous statistics — marker detection effect sizes,
  reference-based annotation diagnostics, bulk DE engines.
- They interoperate (`anndata2ri`, `zellkonverter`, `.h5ad`/`.rds`). Mixing is
  normal and correct: cluster in scanpy, annotate with SingleR, test in DESeq2.
- GPU: `rapids-singlecell` mirrors the scanpy API for 10–100× on large atlases.

## Run the template

```bash
pip install scanpy leidenalg igraph          # + scvi-tools celltypist decoupler pydeseq2 as needed
python3 assets/sc_workflow.py --demo         # synthetic data, no download, end-to-end
python3 assets/sc_workflow.py --h5ad my.h5ad --sample-key sample --condition-key condition
```

`assets/sc_workflow.py` implements steps 2–12 with the defaults above and
prints a per-step audit line (how many cells each filter removed, how many
clusters at each resolution, how many DE genes). It is the reference
implementation of this skill — read it before writing a new pipeline.

## Honesty labeling (same discipline as `network-biomarker`)

Every recommendation carries:

- **✅ benchmark-backed** — an external benchmark or the book's explicit
  recommendation supports it (scIB/Luecken 2022 for integration; Squair 2021 /
  the DE chapter for pseudobulk; the QC chapter for MAD).
- **⚠️ convention** — community practice with no decisive benchmark
  (n_PCs = 50, resolution = 1.0, 2000 HVGs, most trajectory choices).

Never present a ⚠️ default as if it were ✅. When a parameter is a convention,
say so and show the sensitivity (e.g. cluster at three resolutions and report
whether the conclusion survives).

## Hard-won gotchas

- **Filter cells, keep genes.** Gene-level filtering before QC gives no
  measured benefit and breaks ambient-RNA estimation, which needs the raw
  droplet matrix.
- **Doublet detection on pooled batches invents doublets.** Run per sample.
- **HVG selection must be batch-aware** before integration, or you select the
  batch effect as your biological signal.
- **Feature selection on raw counts**, normalization for everything else — the
  log transform of exact zeros is the reason deviance-based selection exists.
- **`rank_genes_groups` on Leiden clusters is circular** — the same data
  defined the clusters and tests them, so p-values are anticonservative. Use
  the ranking to *find* markers, and effect sizes (Cohen's d, AUC / log2FC +
  detection fraction) to *report* them. Never quote those p-values as evidence.
- **A strong condition effect makes clustering split a cell type by
  condition** — and then that "cell type" is collinear with your contrast, so
  DE inside it is undefined (PyDESeq2 fails with a singular matrix; a naive
  script reports nonsense). Reproduced by
  `sc_workflow.py --h5ad` on unannotated clusters. Fix: annotate to *cell
  type*, not to raw cluster ids, and integrate over the condition covariate
  before comparing.
- **Reference and query gene spaces must match exactly** (same IDs, same
  order, missing genes zero-filled) before any label transfer, or scArches/
  CellTypist output is silently garbage.
- **Protein-level markers often fail at RNA level** (CD4 is the canonical
  example: barely detected in CD4⁺ T cells). Validate markers *in your data*
  before trusting a textbook panel.
