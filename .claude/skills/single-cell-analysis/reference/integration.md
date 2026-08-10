# Integration, batch correction and reference mapping

Source: sc-best-practices *Data integration*; benchmark of record —
**Luecken et al., *Nat Methods* 2022, "Benchmarking atlas-level data
integration in single-cell genomics"** (scIB, 16 methods × 5 RNA tasks).

## First: do you need it? ✅

"Batch effect correction is not always required and it might mask the
biological variation of interest." Integration is a *destructive* operation —
it removes variance you cannot get back.

Look at the unintegrated data first:

```python
sc.pp.neighbors(adata); sc.tl.umap(adata)
sc.pl.umap(adata, color=[batch_key, "celltype"])
pd.crosstab(adata.obs["leiden"], adata.obs[batch_key], normalize="index")
```

If clusters contain cells from every batch in roughly the sample proportions,
you have no batch effect to remove. If a cluster is 100% one batch, decide
whether that is technical or a real donor-specific population **before**
correcting it away.

**Choosing the batch covariate** is itself a decision: sample < donor < 10x
chemistry < lab < protocol. Finer batch resolution removes more technical
effect but "fine batch variation is also more likely to be confounded with
biologically meaningful signals." Correcting at the sample level in a
case/control design can delete the disease effect outright.

## Method taxonomy

| Class | Methods | Output |
|---|---|---|
| Global / linear models | ComBat | corrected matrix |
| Linear embedding | **MNN**, fastMNN, **Seurat CCA/RPCA anchors**, **Scanorama**, **Harmony** | embedding (+ matrix for some) |
| Graph-based | **BBKNN** | corrected kNN graph only |
| Deep generative | **scVI**, **scANVI**, scGen, scPoli, trVAE | latent embedding (+ decoded counts) |

## Which one ✅ (from the scIB benchmark)

| Situation | Pick |
|---|---|
| Simple task: few batches, one protocol, one tissue | **Harmony** or **Seurat** — consistently strong, fast, few cells needed |
| Complex task: many datasets, different protocols/technologies | **scVI**, **scGen**, **Scanorama** |
| Cell labels available (even partial) | **scANVI** — "keeps the differences between cell labels while removing batch effects" |
| Only a graph is needed downstream | BBKNN (cheapest) |
| ATAC / multiome | see the benchmark for the modality — RNA rankings do not transfer |

"An optimal method for all scenarios does not exist." Run 2–3 and measure.

```python
# scVI — raw counts, batch as covariate
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
model = scvi.model.SCVI(adata); model.train()
adata.obsm["X_scVI"] = model.get_latent_representation()

# scANVI — initialise FROM the trained scVI model, add labels
sca_model = scvi.model.SCANVI.from_scvi_model(
    model, labels_key=label_key, unlabeled_category="unlabelled")
sca_model.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = sca_model.get_latent_representation()

# Harmony / BBKNN — on PCA
sc.external.pp.harmony_integrate(adata, key=batch_key)   # -> X_pca_harmony
bbknn.bbknn(adata, batch_key=batch_key, neighbors_within_batch=3)

sc.pp.neighbors(adata, use_rep="X_scANVI"); sc.tl.umap(adata)
```

**Label harmonization first.** scANVI (and any label-aware method) breaks when
the same cell means `"T cell"` in one dataset and `"CD8+ T cell"` in another.
Map to a common ontology (Cell Ontology / CL terms) before training.

## How to judge it ✅

**Not by UMAP.** "It is tempting to select an integration based on the UMAPs,
but this does not fully represent the quality of an integration." Every method
produces a well-mixed-looking UMAP; that is what they optimize.

Use `scib-metrics`, which scores the two competing objectives separately:

| Batch correction | Bio conservation |
|---|---|
| kBET | NMI / ARI vs. labels |
| graph iLISI | cell-type ASW |
| batch ASW | isolated-label F1 (rare types!) |
| PCR comparison | graph cLISI, trajectory conservation |

```python
from scib_metrics.benchmark import Benchmarker
bm = Benchmarker(adata, batch_key=batch_key, label_key=label_key,
                 embedding_obsm_keys=["X_pca", "X_scVI", "X_scANVI", "X_pca_harmony"])
bm.benchmark(); bm.plot_results_table(min_max_scale=False)
```

The overall score in the benchmark weights **bio conservation 0.6 / batch
removal 0.4** — deliberately, because over-correction is the worse failure.

## Over-correction: the tell-tales

- A rare population present in only one batch vanishes after integration
  (watch the isolated-label F1 metric specifically).
- Cell types that should be distinct merge (the CD4/CD8 boundary is a good probe).
- The condition effect you are studying disappears — check by running DE
  before and after and comparing.
- Integrating across a covariate that is **confounded with your biology**
  (all controls run in batch 1, all cases in batch 2) removes the effect by
  construction. No method can fix a confounded design; only experimental
  design can. Say so in the report rather than integrating anyway.

## Reference mapping (query → existing atlas)

Different problem from de-novo integration: the reference is fixed and you map
onto it. Use **scArches** (architecture surgery on a pre-trained scVI/scANVI/
totalVI model), or Seurat's `FindTransferAnchors`/`MapQuery`, or Symphony.

Non-negotiable preconditions:
1. Query gene space = reference gene space, **same identifiers, same order**,
   missing genes zero-filled.
2. Raw counts as input for scVI-family models.
3. Query preprocessing must mirror what the reference model was trained on.

Then transfer labels with uncertainty and threshold it → `annotation.md`
Strategy E.

Available references: CellxGene Census, Human Lung Cell Atlas, Human Cell
Atlas projects, Azimuth references, scvi-hub pretrained models.
