# Preprocessing: raw counts → clusters

Source: sc-best-practices, *Preprocessing and visualization* + *Identifying
cellular structure* parts; OSCA *Basic* chapters 1–6.

Labels: ✅ = benchmark-backed / explicitly recommended · ⚠️ = convention.

---

## 1. Raw data processing

| Assay | Tool | Note |
|---|---|---|
| 10x 3'/5' | CellRanger, STARsolo, **alevin-fry** (fast, low-memory) | alevin-fry/kb-python give near-identical counts at a fraction of the cost |
| plate-based (Smart-seq) | STAR + featureCounts / salmon | full-length ⇒ no UMIs ⇒ no deduplication, and TPM-style length correction matters |
| velocity | **velocyto** or alevin-fry `--use-velocity` | you must decide *before* quantification: spliced/unspliced layers cannot be recovered later |

Empty-droplet calling: `DropletUtils::emptyDrops` (OSCA) rather than the
knee-point cutoff — the knee removes real low-RNA cells. ✅

## 2. Quality control ✅

Three covariates, considered **jointly**, never alone:

1. count depth (`total_counts`)
2. number of detected genes (`n_genes_by_counts`)
3. mitochondrial fraction (`pct_counts_mt`)

Optional: ribosomal fraction, hemoglobin fraction (RBC contamination),
`pct_counts_in_top_20_genes` (library complexity).

**Interpretation is joint.** High mito + low counts = dying cell. High mito +
high counts = a metabolically active cell type (cardiomyocytes, hepatocytes,
proximal tubule) — filtering on mito alone deletes it.

### MAD thresholds (the recommendation)

```
MAD = median(|x_i - median(x)|)
outlier if x < median - k*MAD or x > median + k*MAD
```

| Metric | k |
|---|---|
| `log1p_total_counts` | 5 |
| `log1p_n_genes_by_counts` | 5 |
| `pct_counts_in_top_20_genes` | 5 |
| `pct_counts_mt` | 3 **plus** a hard `> 8%` filter |

Why not fixed cutoffs: "filtering should be based on median absolute
deviations with lenient cutoffs to avoid bias against smaller subpopulations."
A `n_genes > 500` rule deletes resting lymphocytes, neutrophils and every
low-RNA type in the tissue.

**The check nobody does:** cluster the *discarded* cells. If they form a
coherent cluster with real markers, you deleted a cell type, not debris.

Feature/gene filtering before QC shows **no measured downstream benefit** ✅ —
and it breaks ambient-RNA estimation, which needs the unfiltered droplet matrix.

## 3. Ambient RNA ✅

Free-floating mRNA from lysed cells contaminates every droplet, and it is
tissue-specific (a pancreas dataset where every cell "expresses" insulin).

```r
sc <- SoupChannel(raw_matrix, filtered_matrix, calcSoupProfile = TRUE)
sc <- setClusters(sc, cluster_vector)     # improves the estimate a lot
sc <- autoEstCont(sc, doPlot = FALSE)
out <- adjustCounts(sc, roundToInt = TRUE)
```

Alternatives: `CellBender` (deep generative, also removes empty droplets),
`DecontX`. Run **before** doublet detection and QC-filtering, on the raw matrix.

## 4. Doublets ✅

`scDblFinder` (best-in-benchmark, R) or `scDblFinder`-equivalents
`DoubletFinder`, `scrublet`, `solo` (Python).

```r
sce <- scDblFinder(SingleCellExperiment(list(counts = mat)))
# scDblFinder.score / .ratio / .class
```

**Hard rule:** never run doublet detection on data pooled across batches — the
simulated doublets become cross-batch chimeras and the calls are meaningless.
One sample at a time.

Genetic multiplexing (`souporcell`, `demuxlet`, `vireo`) identifies
cross-donor doublets directly and is strictly better when donors are pooled.

## 5. Normalization ✅

| Method | Use it for | Call |
|---|---|---|
| **Shifted logarithm** `log(y/s + 1)` | dimensionality reduction, DE — "outperforms other methods for uncovering latent dataset structure" | `sc.pp.normalize_total(adata, target_sum=None)` then `sc.pp.log1p` |
| **scran pooling** deconvolution size factors | strongly varying count depth, batch correction | `scran::computeSumFactors(sce, clusters=quickCluster(sce), min.mean=0.1)` |
| **Analytic Pearson residuals** | gene selection, rare cell types; no pseudo-count, no log | `sc.experimental.pp.normalize_pearson_residuals(adata)` |

`target_sum=None` (median depth) rather than `1e4`: CP10K is an arbitrary
scale that exaggerates the pseudo-count's effect on low-depth cells. ⚠️ but
preferred by the book.

**Always keep the raw counts** in `adata.layers["counts"]`. Every count-model
method downstream (scVI, DESeq2, deviance selection, velocity) needs them, and
they cannot be reconstructed after normalization.

The book's own caveat: no method wins for every downstream task; the
22-algorithm benchmark did not settle it. Choose by task, and say which you chose.

## 6. Feature selection ✅

**Binomial deviance on raw counts** (`scry::devianceFeatureSelection`,
Townes et al. 2019). Closed-form, no normalization or pseudo-count needed —
which matters because you cannot log-transform exact zeros without distorting
the variance estimate that variance-based HVG methods depend on.

```r
sce <- devianceFeatureSelection(sce, assay = "X")
```

Number of genes: 2000 is the convention ⚠️; the book's example takes 4000 of
~20000; 500–2000 is typical in standard pipelines. Known bias: high-mean genes
tend to be selected as highly deviant.

Alternatives: `seurat_v3` (variance-stabilizing, on raw counts),
`cell_ranger`/`seurat` (dispersion on log data).

**Batch-aware selection before integration:**
`sc.pp.highly_variable_genes(adata, batch_key="batch")` selects genes variable
*within* batches. Skipping this means you select the batch effect and then ask
an integration method to remove what you just enriched for.

## 7. Dimensionality reduction

- **PCA**, ~50 components ⚠️, on scaled HVGs. Everything downstream (kNN,
  clustering, most integration) runs on the PCA space, so this is the real
  representation; UMAP is decoration.
- **UMAP/t-SNE for visualization only.** ✅ Do not measure distances, do not
  read cluster adjacency as relatedness, do not draw trajectories on it.
  The book states this explicitly; it is the most-violated rule in the field.
- Diffusion maps for continuous/trajectory structure.
- Deep alternatives (scVI latent space, hyperbolic embeddings for hierarchies)
  → `deep-learning-map.md`.

## 8. Clustering ✅

**Leiden**, not Louvain (Louvain can produce internally-disconnected
communities and is unmaintained); community detection on the kNN graph beats
k-means/hierarchical on scRNA-seq.

```python
sc.tl.leiden(adata, resolution=0.25, key_added="leiden_res0_25",
             flavor="igraph", n_iterations=2, directed=False)
```

- Resolution 1.0 is a starting point ⚠️ — **run several** (0.25/0.5/1.0) and
  keep them all in `.obs` under explicit keys.
- `n_iterations=2` standard; `-1` forces convergence (slow on big data).
- **Sub-cluster** to resolve states within a type — but sub-clustering finds
  structure in noise if you push it; require marker evidence for every split.
- Clustering is a *prerequisite* for annotation, not annotation itself.

### The circularity trap

`sc.tl.rank_genes_groups` on Leiden clusters tests the same data that defined
the clusters. The p-values are anticonservative by construction — inflated,
sometimes wildly. Use the ranking to **find** candidate markers; report
**effect sizes** (Cohen's d, AUC, log2FC + fraction detected) or validate on
held-out data. OSCA makes the same point and provides effect-size-first
marker detection (`scoreMarkers`).

## 9. GPU

`rapids-singlecell` mirrors the scanpy API (`rsc.pp.*`, `rsc.tl.*`) with
10–100× speedups on PCA/neighbors/UMAP/Leiden. Worth it above ~500k cells.
