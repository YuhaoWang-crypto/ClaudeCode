# Comparing conditions: DE, composition, pathways, perturbation

Source: sc-best-practices *Dealing with conditions* part. This is where most
published single-cell analyses are statistically wrong, so the rules here are
strict.

---

## 1. Differential expression — use pseudobulk ✅

**The rule:** aggregate counts per (sample × cell type), then run a bulk method.

**Why:** cells from one donor are not independent replicates. Testing them as
if they were is *pseudoreplication*: "inferential statistics is applied to
biological replicates which are not statistically independent", and the FDR
blows up. Benchmarks found bulk methods on pseudobulk beat single-cell-specific
DE tools, and that single-cell methods were "especially prone to wrongly
labeling highly expressed genes as differentially expressed" — because with
10 000 cells any effect is p < 10⁻³⁰⁰, including the effect of one outlier donor.

### Workflow

```python
import decoupler as dc, pertpy as pt

# 1. aggregate raw counts by SUM per sample x cell type
pdata = dc.pp.pseudobulk(adata, sample_col="sample", groups_col="cell_type",
                         layer="counts", mode="sum")
#    (decoupler < 2.0: dc.get_pseudobulk with the same arguments)

# 2. drop underpowered pseudobulks
dc.pp.filter_samples(pdata, min_cells=10, min_counts=1000)

# 3. PER CELL TYPE: drop lowly expressed genes -- different cell types
#    express different gene sets, so this must not be done globally
sub = pdata[pdata.obs["cell_type"] == "Monocyte"].copy()
dc.pp.filter_by_expr(sub, group="condition", min_count=10, min_total_count=15)

# 4. inspect the design BEFORE testing
sc.pp.normalize_total(sub); sc.pp.log1p(sub); sc.pp.pca(sub)
sc.pl.pca(sub, color=["condition", "donor", "batch", "sex"])   # what actually drives PC1?

# 5. test
dds = pt.tl.PyDESeq2(adata=sub, design="~batch + condition")
dds.fit(); res = dds.test_contrasts(dds.contrast("condition", "stim", "ctrl"))
```

`assets/sc_workflow.py` implements steps 1–5 in ~40 lines of numpy + PyDESeq2
if you would rather not add decoupler/pertpy.

### Requirements and caveats

- **Replicates, not cells.** "The best way to increase statistical power is to
  increase the number of independent experimental samples." More cells per
  donor improves precision within a donor and does almost nothing for
  between-donor power. ~8 samples per condition is a reasonable floor.
- **n < 3 per group: do not report p-values.** Describe effect sizes and call
  it exploratory.
- **Model your covariates.** Run PCA on the pseudobulks first; if PC1 is
  `batch` or `sex`, put it in the design matrix. Failing to "inflates the FDR."
- Sum vs. mean aggregation: sum is what the count models expect; the book
  notes the comparison "requires further investigation."
- BH correction always; report `padj`, not `pvalue`.
- **Cell-type-level DE within one sample** (marker detection) is a different
  question and cell-level tests are fine there — the pseudoreplication problem
  is specific to generalising across *individuals*.

### The confounding trap you will actually hit

Test the design **before** testing genes: every group in the contrast needs
≥2 pseudobulk replicates.

A strong condition effect makes Leiden split one cell type into
condition-specific clusters. If you then run DE "per cluster", the cluster is
collinear with the contrast — the design matrix is singular and the question is
undefined. (`assets/sc_workflow.py --h5ad` reproduces this on purpose: with raw
cluster ids, cluster 4 is 100% control and cluster 5 is 100% stimulated, both
of them the same underlying monocyte population.)

Fixes, in order of preference: annotate to **cell type** rather than cluster
id; integrate over the condition covariate before clustering; or merge the
split clusters after confirming with markers that they are one type. Silently
reporting DE for a population present in one arm only is the failure mode this
guard exists to prevent.

### When cell-level methods are defensible

MAST with a random effect for donor, or a proper mixed model (glmmTMB, muscat's
mixed mode), when you genuinely need cell-level resolution (e.g. a continuous
covariate measured per cell). They are slower and still need the donor random
effect. Never a plain Wilcoxon across conditions.

## 2. Compositional analysis ✅

**The rule:** cell-type proportions sum to 1, so they are compositional data
and cannot be tested independently.

If enterocytes double in absolute number and nothing else changes, every other
population's *proportion* falls. A per-type Wilcoxon "will falsely perceive
cell-type population shifts as statistically sound effects, although they were
induced by inherent negative correlations."

| Method | Level | Notes |
|---|---|---|
| **scCODA** | labelled cell types | Bayesian Dirichlet-multinomial; designed for few replicates |
| **tascCODA** | cell-type *tree* | detects effects on aggregated branches; for fine-grained atlases |
| **Milo / miloR** | kNN neighbourhoods | no clustering needed — catches transitional states clusters hide |
| **DA-seq** | per cell | scores each cell by neighbour condition prevalence |

```python
import pertpy as pt
sccoda = pt.tl.Sccoda()
mdata  = sccoda.load(adata, type="cell_level", generate_sample_level=True,
                     cell_type_identifier="cell_type", sample_identifier="sample",
                     covariate_obs=["condition"])
mdata  = sccoda.prepare(mdata, formula="condition",
                        reference_cell_type="automatic")
sccoda.run_nuts(mdata); sccoda.summary(mdata)
```

**Parameters that decide the answer:**
- **Reference cell type** — scCODA needs one assumed-unchanged population for
  identifiability. `"automatic"` picks a low-variance one; a manual choice
  "substantively influences interpretations", so state which you used.
- **FDR** — default 0.05 is conservative for this model; up to 0.2 is
  acceptable practice and reveals more effects. Report the value you used.
- **MCMC acceptance rate** should sit in 0.4–0.9; outside that the sampler is
  broken and the credible intervals are meaningless.

Limitations: log-linear covariate relationship assumed; detects mean shifts
only, not variance changes; cannot recover correlation structure.

Also confounded by **FACS sort gates** and dissociation protocol — if samples
were sorted under different gates, whole-tissue composition is a property of
the gate, not the biology. Restrict to a gate present in every sample.

## 3. Gene set / pathway analysis

Per-cell or per-pseudobulk scoring rather than a hypergeometric test on a gene
list, when you have continuous data:

```python
import decoupler as dc
dc.mt.ulm(adata, net=dc.op.progeny(organism="human"))     # pathway activities
dc.mt.ulm(adata, net=dc.op.collectri(organism="human"))   # TF activities
```

- **PROGENy** (14 pathways, footprint-based — uses downstream response genes,
  not pathway members: much better signal than KEGG membership).
- **CollecTRI/DoRothEA** for TF activity.
- **MSigDB Hallmark** for a broad first pass; GO for detail.
- Over-representation on a DE list is fine but ranks small sets artificially
  high — filter to sets with ≤200 genes; for a gene-list ORA against Reactome/
  GO/InterPro with BH FDR, the `omics-ppi-pathway` skill does this offline.

Interpret enrichment on **pseudobulk** results between conditions; per-cell
scores are for visualization and are competitive (one program's activity
depresses others).

## 4. Perturbation modelling ⚠️

CRISPR screens (Perturb-seq), drug panels, cytokine stimulation. Tooling:
`pertpy` (mixscape for guide-assignment QC, augur for perturbation-responsive
cell-type ranking, scGen/CPA/GEARS wrappers), `scPerturb` for harmonized data.

**The honest state of the art:** a 2024 benchmark found that *deep-learning
predictions of gene-perturbation effects do not yet outperform simple linear
methods*. Before reporting that a deep model predicts unseen perturbations:

1. Run the trivial baselines — predict the training mean; predict the control
   profile unchanged; linear additive model of single perturbations.
2. Evaluate on genuinely unseen perturbations, not random cell splits.
3. Report the metric on the DE genes, not on all genes (all-gene correlation is
   dominated by the unchanged majority and looks great for a null model).

Mixscape first: cells that received a guide but show no transcriptional effect
("escapers") must be removed, or they dilute every downstream estimate.
