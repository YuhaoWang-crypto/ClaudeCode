# Cell-type annotation

Sources: OSCA.basic ch. 7 *Cell type annotation* (R/Bioconductor) +
sc-best-practices *Annotation* chapter (Python/scverse). These two disagree in
tooling and agree completely in discipline, which is what makes the
intersection trustworthy.

> "Cell type" has no rigorous definition. Practitioners recognise types
> intuitively, borders are partly subjective and move over time. Manual
> interpretation is the acknowledged bottleneck of the whole workflow.

## The discipline (non-negotiable)

1. **Run two independent strategies** — one reference/automated, one
   marker-based — and reconcile them. OSCA: "apply both strategies to examine
   the agreement." Disagreement is *information*, not failure: the classic case
   is reference-based methods cleanly separating CD4⁺/CD8⁺ T cells while
   unsupervised clustering groups them together, because activation state
   dominates the transcriptome more than the CD4/CD8 distinction does.
2. **Automated annotation is a starting point, never an endpoint.**
   Every published label needs marker-gene support.
3. **Ship a confidence column.** Low-confidence / high-uncertainty cells become
   `Unknown`. An honest `Unknown` beats a confident wrong label.
4. **Reference quality is the ceiling.** A reference-based method can only
   return labels its reference contains — it cannot discover a novel state, and
   it will confidently mislabel one.

---

## Strategy A — manual markers (cluster level)

The baseline everyone still trusts. "Expression of known marker genes is still
the most accepted support" for an annotation.

```python
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1")
sc.pl.umap(adata, color=marker_genes, vmax="p99")          # p99: kill outlier scaling
sc.pl.dotplot(adata, var_names=marker_dict, groupby="leiden_1", standard_scale="var")
adata.obs["celltype"] = adata.obs.leiden_1.map(annotation_dict)
```

Find candidates first, then check them:

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_1", method="wilcoxon")
sc.tl.filter_rank_genes_groups(adata, min_in_group_fraction=0.2,
                               max_out_group_fraction=0.2)
```

Or, per-cell set scoring: `sc.tl.score_genes(adata, gene_list, score_name=...)`.

**Caveats**
- Labor-intensive and subjective; results depend on cluster resolution — test
  several.
- "Markers that work in one dataset often underperform in others."
- **Protein markers frequently fail at RNA level** — CD4 is barely detected in
  CD4⁺ T cells. Validate any textbook panel in your own data before trusting it.
- Dropout means a marker-negative cell is not evidence of absence; this is why
  annotation happens at cluster level, not cell level.

Marker sources: CellMarker, PanglaoDB, Azimuth references, the tissue atlas
papers themselves, MSigDB/GO for programs rather than identities.

## Strategy B — reference-based: SingleR (R) ✅

The OSCA workhorse. Spearman correlations between each test cell and reference
samples, restricted to markers that distinguish each pair of labels, then a
fine-tuning round using only the top-scoring labels' markers.

```r
library(SingleR)
ref  <- celldex::BlueprintEncodeData()          # or HumanPrimaryCellAtlasData(), etc.
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
table(pred$labels)

plotScoreHeatmap(pred)                           # per-cell scores, scaled to [0,1]
plotDeltaDistribution(pred)                      # gap between best and median label
pruned <- pred$pruned.labels                     # NA where the call is unreliable
table(Assigned = pruned, Cluster = colLabels(sce))   # reconcile with clustering
```

**Diagnostics that matter**
- `plotScoreHeatmap`: a cell scoring high on *many* labels is ambiguous, not
  confidently multi-labeled.
- **delta** = best score − median score. A small delta means the call is
  low-information; `pruneScores`/`pruned.labels` encodes this as `NA`.
  Keep the NAs.
- Cross-tabulate against your unsupervised clusters. A label that scatters
  across every cluster is noise; a cluster that splits cleanly into two labels
  is a real sub-structure your resolution missed.

**Custom single-cell reference** (usually better than a bulk reference):
```r
pred <- SingleR(test = sce, ref = ref_sce, labels = ref_sce$celltype,
                de.method = "wilcox")   # pairwise markers for a single-cell ref
```

**Limitation, stated:** "restricted by the diversity and resolution of the
available labels." Novel states are invisible to it.

## Strategy C — marker-set enrichment per cell: AUCell (R) ✅

Reference-free: needs only marker *identities*, not reference expression.

```r
library(AUCell)
rankings <- AUCell_buildRankings(counts(sce), plotStats = FALSE)
aucs     <- AUCell_calcAUC(marker_sets, rankings)
results  <- t(assay(aucs))
AUCell_exploreThresholds(aucs, plotHist = TRUE)   # expects BIMODAL distributions
```

Ranks genes within each cell, computes the AUC of the recovery curve for each
marker set, assigns the max-AUC label.

**Diagnostic:** the AUC histogram must be **bimodal** (a high-scoring
population + everything else). Unimodal ⇒ the gene set is uninformative in this
dataset, and the "threshold" is meaningless.

**Limitations:** relative expression is discarded (only ranks survive);
scoring is *competitive*, so activity in one set depresses the others; not
valid for overlapping sets or sets containing opposing-direction genes.

## Strategy D — classifiers: CellTypist (Python) ✅

Logistic-regression models trained on thousands of genes across large atlases —
much broader signal than a hand-written marker panel.

```python
sc.pp.normalize_total(adata, target_sum=1e4)   # CellTypist REQUIRES CP10K + log1p
sc.pp.log1p(adata)
model = celltypist.models.Model.load(model="Immune_All_Low.pkl")
pred  = celltypist.annotate(adata, model=model, majority_voting=True)
adata = pred.to_adata()      # .obs: predicted_labels, majority_voting, conf_score
```

- `majority_voting=True` over-clusters and takes the cluster majority — much
  more stable than per-cell calls.
- `conf_score` is the confidence column: threshold it, don't ignore it.
- Note the input contract: CellTypist wants CP10K-normalized log1p data
  specifically. Feeding it raw or differently-scaled data fails silently.
- Models: `Immune_All_High/Low`, Lung/Gut/Brain atlases, `models.models_description()`.

## Strategy E — reference mapping + label transfer: scArches / scANVI ✅

The atlas-scale approach: map your query into a pre-trained reference latent
space, transfer labels by weighted kNN, and get **uncertainty per cell**.

```python
# 1. Gene space MUST match the reference exactly: same IDs, same order,
#    missing genes zero-filled.  This is the #1 silent failure.
adata_q = adata_q[:, ref_genes].copy()          # after zero-filling missing ones
# 2. Load the reference model with the query attached
model = sca.models.SCVI.load_query_data(adata_q, ref_model_path)
model.train(max_epochs=500, plan_kwargs={"weight_decay": 0.0})
adata_q.obsm["X_scVI"] = model.get_latent_representation()
# 3. Weighted kNN label transfer with uncertainty
knn = sca.utils.knn.weighted_knn_trainer(adata_ref, "X_scVI", n_neighbors=15)
labels, uncert = sca.utils.knn.weighted_knn_transfer(
    adata_q, "X_scVI", adata_ref.obs, label_keys="cell_type", knn_model=knn)
adata_q.obs["celltype"] = labels["cell_type"].where(uncert["cell_type"] < 0.2, "Unknown")
```

- Uses **raw counts**; needs a GPU to be pleasant.
- scANVI (label-aware scVI) is the recommended model when reference labels
  exist — it keeps label differences while removing batch effects.
- **Uncertainty is the point**: high-uncertainty regions are either novel
  biology or bad mapping. But the book's own caveat — "uncertainty scores are
  often imperfect and sometimes fail to highlight new cell types or states."
  Do not treat low uncertainty as proof.

## Strategy F — cluster-level gene-set testing (what is this cluster *doing*)

Identity vs. activity. For programs rather than types:

```r
go.out <- limma::goana(unique(entrez[is.de]), species = "Hs", universe = entrez_all)
# then inspect the genes behind an interesting term
aggregated <- scuttle::sumCountsAcrossFeatures(sce, by.go,
                                exprs_values = "logcounts", average = TRUE)
```
Python equivalents: `decoupler` (with PROGENy/DoRothEA/MSigDB), `gseapy`.

**Caveats:** all conclusions are *relative to the other clusters* — an outgroup
result, not an absolute identity. Filter GO to BP terms with ≤200 genes or the
top hits are uselessly general. Per-set averaging is for **visualization**, not
discovery: non-DE members add noise and opposing members cancel.

## Strategy G — LLM-assisted annotation ⚠️

Recent work (e.g. GPTCelltype; mLLMCelltype 2025, which uses multi-LLM
consensus) shows LLMs annotate competitively **from marker-gene lists**.
Treat as a fast hypothesis generator:
- feed the top-N markers per cluster, not the expression matrix;
- require multi-model or multi-run consensus, and treat disagreement as
  low confidence;
- verify every accepted label against marker expression in your data.

Not validated enough to stand alone. Same rule as every other automated method.

---

## Reconciliation checklist

```python
pd.crosstab(adata.obs["leiden_1"], adata.obs["automated_label"])  # where do they disagree
sc.tl.dendrogram(adata, groupby="celltype")                       # is the hierarchy sane
sc.pl.dotplot(adata, marker_dict, groupby="celltype", standard_scale="var")
adata.obs.loc[uncertainty > 0.2, "celltype"] = "Unknown"
```

- [ ] Two independent methods run; disagreements listed and explained.
- [ ] Every label has ≥2 supporting markers visible in a dotplot **of this dataset**.
- [ ] Coarse and fine annotations are mutually consistent (fine nests in coarse).
- [ ] Dendrogram groups related types together — if myeloid sits inside T cells,
      something upstream is wrong.
- [ ] Uncertain cells labelled `Unknown`, counted, and reported.
- [ ] Annotation repeated at ≥2 cluster resolutions; conclusions survive both.
- [ ] Doublet-derived clusters (co-expressing two exclusive programs) removed,
      not annotated as a "transitional" population.
