# Deep learning in single-cell: task → method map

Distilled from [awesome-deep-learning-single-cell-papers](https://github.com/OmicsML/awesome-deep-learning-single-cell-papers)
(OmicsML) and its survey — Molho et al., *Deep learning in single-cell
analysis*, **ACM TIST** 15(3), 2024. Companion list for foundation models:
[awesome-foundation-model-single-cell-papers](https://github.com/OmicsML/awesome-foundation-model-single-cell-papers).

The list is a bibliography; this file is the part you need at the bench —
**which task has a DL method worth using, and what it has to beat.**

---

## The gate: does this task need deep learning?

Run this before adopting any model.

1. **What is the non-DL baseline?** PCA + Harmony; a linear/additive model;
   marker scoring; NMF. If the paper has no baseline comparison, the claim is
   unvalidated.
2. **Does an external benchmark cover this task?** scIB for integration,
   BEELINE for GRN inference, the deconvolution and CCC benchmarks below. Prefer
   the benchmark's ranking to the method paper's self-report — a meta-analysis
   of single-cell benchmarks (*Genome Biol* 2023) exists precisely because
   self-reported rankings do not replicate.
3. **Is the evaluation split biologically meaningful?** Random cell splits leak;
   held-out *donors*, *datasets*, or *perturbations* do not.
4. **Is the metric dominated by the trivial majority?** All-gene correlation in
   perturbation prediction looks excellent for a model that predicts "no
   change". Evaluate on the DE genes.
5. **Do you need what DL adds?** Its genuine advantages are: a count noise
   model, conditioning on covariates, an amortized encoder that generalizes to
   new cells, and uncertainty. If none apply, the linear method is better
   because it is auditable.

**The field's own cautionary result:** *"Deep learning-based predictions of
gene perturbation effects do not yet outperform simple linear methods"*
(bioRxiv 2024). Similar patterns recur in imputation and clustering
benchmarks. This is the default expectation, not an outlier.

---

## Task → method

### Representation learning / dimensionality reduction
- **scVI** — *Deep generative modeling for single-cell transcriptomics*, Nat Methods 2018. The workhorse: ZINB/NB VAE, batch covariates, latent space for everything downstream.
- **scvis** — *Interpretable dimensionality reduction ... deep generative models*, Nat Commun 2018.
- **SIMBA** — Nat Methods 2023: cells *and* features (genes, peaks, TFs) in one embedding space.
- Hyperbolic embeddings for hierarchical/developmental structure — Nat Commun 2021; *Complex hierarchical structures ... deep hyperbolic manifold learning*, Genome Res 2023.
- **SATURN** — cross-species embedding via protein-language-model gene representations.

### Batch integration / atlas building → also `integration.md`
- **scVI / scANVI** (scANVI = *Probabilistic harmonization and annotation ... deep generative models*, Mol Syst Biol 2021) — top tier in the scIB benchmark for complex tasks.
- **scArches** — *Mapping single-cell data to reference atlases by transfer learning*, Nat Biotech 2021. Architecture surgery: extend a pretrained reference model with query-specific weights.
- **trVAE** (Bioinformatics 2020), **scDREAMER** (Nat Commun 2023), **CLAIRE** (Bioinformatics 2023) — contrastive, balances mixing vs. heterogeneity.
- Non-DL competitors that win on simple tasks: **Harmony** (Nat Methods 2019), **Seurat CCA** (Nat Biotech 2018), **MNN** (Nat Biotech 2018).
- Benchmarks: **scIB** (Nat Methods 2022); batch-correction benchmark (Genome Biol 2020); **kBET** metric (Nat Methods 2018).

### Cell-type annotation → also `annotation.md`
- **CellTypist** — *Cross-tissue immune cell analysis*, Science 2022 + *Scaling cross-tissue single-cell annotation models*, 2023. Logistic regression at atlas scale; the pragmatic default.
- **scBERT** — Nat Mach Intell 2022, transformer for annotation. But read the *Reusability report: learning the transcriptional grammar ... using transformers* (Nat Mach Intell 2023) before believing the headline numbers.
- **TOSICA** — *Transformer for one stop interpretable cell type annotation*, Nat Commun 2023 (attention → interpretable pathway/TF tokens).
- **scDeepSort** (NAR 2021, weighted GNN), **scGCN** (Nat Commun 2022), **scIAE** (Brief Bioinform 2022), **ACTINN**, **SciBet** (Nat Commun 2020), **CHETAH**, **SingleCellNet**.
- **STELLAR** — Nat Methods 2022: annotated reference cell *graph* + query graph, for spatial data with neighbourhood context.
- **mLLMCelltype** — 2025: multi-LLM consensus annotation from marker lists.

### Foundation models ⚠️
scGPT, **Geneformer** (*Transfer learning enables predictions in network
biology*, Nature 2023), scFoundation, CellPLM, scPRINT, xTrimoGene,
Cell2Sentence, tGPT, GET.

Real value today: a pretrained embedding for **few-shot / low-data** settings
and for tasks with no good supervised baseline. Persistent caveats: zero-shot
performance often below a tuned task-specific model; heavy compute; evaluation
leakage between pretraining corpora and benchmark datasets is common. Always
compare against scVI + a linear probe. *Evaluating the Utilities of Large
Language Models in Single-cell Data Analysis* (bioRxiv 2023) is the sober read.

### Clustering
**scDeepCluster** (Nat Mach Intell 2019), **DESC** (Nat Commun 2020, clustering
+ batch removal), **scGAC** (graph attention), **scDSC**, ZINB graph-embedding
autoencoder (AAAI 2022), **SC3s** (scalable consensus).
→ Leiden on a good representation (PCA or scVI latent) remains the baseline to beat.

### Imputation / denoising ⚠️
**MAGIC** (Cell 2018), **scImpute** (Nat Commun 2018), **DeepImpute** (Genome
Biol 2019), **scGNN** (Nat Commun 2021), **VIPER**, **G2S3**, SAVER-X
(*Data denoising with transfer learning*, Nat Methods 2019).

**Do not impute before DE or gene-gene correlation analysis** — imputation
manufactures correlation structure and inflates false positives. Use a count
model (scVI, NB GLM) that handles sparsity natively instead.

### Perturbation & drug response ⚠️ (weakest DL claims)
**scGen** (Nat Methods 2019), **CPA** (Mol Syst Biol 2023), **chemCPA**
(NeurIPS 2022), **GEARS** (multi-gene perturbations), **CellOT** (*neural
optimal transport*, Nat Methods 2023), **PerturbNet**, **CellOracle**
(*Dissecting cell identity via network inference and in silico gene
perturbation*, Nature 2023). Data: **scPerturb** (Nat Methods 2024),
benchmark: **CausalBench**.
See the gate above — linear baselines are competitive.

### RNA velocity & dynamics → also `trajectories.md`
**veloVI** (*Deep generative modeling of transcriptional dynamics*, Nat Methods
2023, adds uncertainty), **scTour** (Genome Biol 2023), **PRESCIENT**,
**scNODE**, **scDiffEq** (neural SDE), **moscot** (optimal transport across
time and space).

### Spatial: domains, deconvolution, segmentation
- Domains: **SpaGCN** (Nat Methods 2021), **STAGATE** (Nat Commun 2022, graph attention AE), **CCST** (Nat Comput Sci 2022), **BayesSpace** (Nat Biotech 2021, subspot), **SpatialDE** (SVGs), **CellCharter**, **SpiceMix** (Nat Genet 2023), **stLearn**, **Giotto**.
- Deconvolution: **cell2location** (Nat Biotech 2022), **RCTD** (Nat Biotech 2021), **DestVI** (Nat Biotech 2022), **CARD** (Nat Biotech 2022), **SPOTlight**, **spatialDWLS**, **BayesPrism** (Nat Cancer 2022), **CytoSPACE**, **DSTG**. Benchmark: Nat Commun 2023 + **Spotless**.
- Segmentation: **Cellpose** (Nat Methods 2021), **Baysor** (Nat Biotech 2021), **Mesmer** (Nat Biotech 2021), **BIDCell** (Nat Commun 2024), **JSTA** (joint segmentation + annotation).
- Reviews: *Museum of spatial transcriptomics* (Nat Methods 2022); *Cell segmentation in imaging-based spatial transcriptomics* (Nat Biotech 2021).

### Cell-cell communication
**CellPhoneDB** (*Single-cell reconstruction of the early maternal–fetal
interface*, Nature 2018), **NicheNet** (Nat Methods 2020, ligand → target
genes), **NATMI**, **NCEM** (*Modeling intercellular communication using
spatial graphs of cells*, Nat Biotech 2022), Scriabin (*Inferring cell–cell
communication at single-cell resolution*, Nat Biotech 2023), **scTensor**.
Benchmark/consensus: **LIANA** (*Comparison of methods and resources for
cell–cell communication inference*, Nat Commun 2022) — different methods
disagree substantially, so use the consensus framework, and treat every
prediction as a hypothesis for co-localization or perturbation follow-up.

### Gene regulatory networks
**STGRNS** (transformer, Bioinformatics 2023), **DynGFN** (GFlowNets),
**DeepTFni** (scATAC → TF network, Nat Mach Intell 2022), **SCENIC+/scMEGA**
(multiomic enhancer-based), **Inferelator 3.0**, **SIGNET**.
Benchmark: **BEELINE** (Nat Methods 2020) — the sobering finding is that
inference algorithms have low reproducibility and modest accuracy. Downstream
network-dynamics analysis of an *inferred* GRN inherits that uncertainty
(the `network-biomarker` skill in this repo labels exactly this boundary).

### Multiomic integration & cross-modality translation
**GLUE** (*graph-linked embedding*, Nat Biotech 2022), **Cobolt**, **MultiVI**
(*mixture-of-experts deep generative model*, Cell Rep Methods 2021),
**Multigrate**, **scJoint** (Nat Biotech 2022), **totalVI** (RNA + surface
protein), **BABEL** (PNAS 2020), **sciPENN** (Nat Mach Intell 2022),
**scButterfly** (Nat Commun 2024), **SCIM**, **scDART**, **SMILE**.
Benchmarks: RNA+ATAC integration (bioRxiv 2023), joint paired/unpaired
integration (Genome Biol 2023), multi-omics fusion for cancer (Genome Biol 2022).

### Simulation (for benchmarking your own method)
**scDesign3** (Nat Biotech 2023, multimodal + spatial), **scDesign2**,
**scMultiSim** (Nat Methods 2025, GRN- and CCI-guided), **scReadSim**
(read-level), **GRouNdGAN**. Benchmark of simulators: Nat Commun 2021.

Use these rather than a hand-rolled negative binomial when you need ground
truth — they preserve gene-gene correlation, which naive simulators do not,
and correlation is exactly what most methods exploit.

---

## Reporting standard for a DL result

- [ ] Named non-DL baseline, run on the same split, reported.
- [ ] Split held out at the biologically meaningful level (donor / dataset / perturbation).
- [ ] Metric computed on the signal, not the constant majority.
- [ ] Seeds varied; variance across runs reported.
- [ ] External benchmark consulted where one exists.
- [ ] Compute cost stated — a 3-day GPU run that ties Harmony is a negative result.
