# scRNA-seq of melanoma before/after checkpoint immunotherapy — reanalysis

**Dataset:** Sade-Feldman et al., *Defining T Cell States Associated with Response
to Checkpoint Immunotherapy in Melanoma*, **Cell 2018** (GEO series **GSE120575**,
raw data dbGaP **phs001680**). CD45⁺ immune cells from 48 melanoma tumor biopsies
taken **before or on/after** immune-checkpoint blockade (anti–PD‑1 ± anti–CTLA‑4),
full-length Smart-seq2.

This is a fully reproducible reanalysis run from the raw GEO matrix
(`analysis/download_data.sh` → `analysis/run_scrna_analysis.py`). Data files are
large and referenced by name, not committed.

## What was run

- Loaded the genes × cells TPM matrix (**55,737 genes × 16,291 cells**) and the
  per-cell clinical annotation (patient, timepoint Pre/Post, Responder/Non-responder).
- QC (median **2,151 genes/cell**), log1p, 2,000 HVGs, PCA → kNN → **Leiden
  (24 clusters)** → UMAP.
- Annotated clusters to lineages by canonical markers.
- Scored two CD8 T-cell programs and derived a per-lesion responder signature.

Cells resolve into clean immune populations:

![cell types](figures/umap_celltypes.png)

## 1. Populations that expand / contract

### With response (Responder vs Non-responder, per-lesion fractions)

| population | mean frac R | mean frac NR | log2FC (R/NR) | p (Mann–Whitney) |
|---|---|---|---|---|
| **B cells** | 0.160 | 0.041 | **+1.94** | **0.0029** |
| **CD4 T** | 0.127 | 0.052 | **+1.28** | **0.0008** |
| Plasma | 0.027 | 0.019 | +0.48 | 0.61 |
| Treg | 0.164 | 0.159 | +0.04 | 0.90 |
| CD8 T | 0.373 | 0.447 | −0.26 | 0.15 |
| NK | 0.098 | 0.137 | −0.47 | 0.20 |
| **pDC** | 0.007 | 0.022 | **−1.46** | **0.0033** |
| **Myeloid** | 0.043 | 0.123 | **−1.49** | **0.0088** |

**Expanded in responders:** B cells and CD4 T cells (both significant).
**Expanded in non-responders:** myeloid cells and pDCs (both significant).
The B-cell enrichment reproduces the original study and the later
tertiary-lymphoid-structure findings (Helmink et al., *Nature* 2020). Note the
*quantity* of CD8 T cells is **not** what separates responders — their *state* is
(section 3).

### With treatment (Post vs Pre)

Aggregated across all patients, **no population shifts significantly** with
treatment timepoint (all p > 0.05; CD4 T trends down, p≈0.06). This reproduces
the paper's observation that per-cluster frequencies show no consistent global
Pre→Post change — the informative axis is responder vs non-responder, and (in the
original TCR analysis) per-lesion transitions between memory and exhausted states.
Full tables: `population_enrichment_Post_vs_Pre.csv` and
`..._Post_vs_Pre_responders.csv`.

![composition](figures/composition_by_group.png)

## 2. Marker genes (Wilcoxon, top per lineage)

| lineage | top markers |
|---|---|
| CD8 T | CD8A, CD8B, CCL5, NKG7, GZMA, CST7, TRAC |
| CD4 T / Treg | CD4, IL32, TNFRSF18, ICOS, CD28, TRAC |
| B | MS4A1, CD79A, BANK1, CD22, IRF8, HLA-DRA |
| NK | GNLY, NKG7, PRF1, KLRK1, TRDC |
| Myeloid | LYZ, IFI30, CST3, TYROBP, FCER1G, AIF1 |
| pDC | LILRA4 / IL3RA-defined cluster |

Full table: `cluster_marker_genes.csv`.

## 3. Responder signature from the two CD8 T-cell states

Two CD8 programs (gene sets curated from the paper's reported markers):

- **CD8_G — memory-like** (`TCF7, IL7R, SELL, CCR7, CD28, LEF1, CD27, GZMK, …`)
- **CD8_B — exhausted** (`PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, ENTPD1, CD38, TOX, …`)

![cd8 signatures](figures/umap_cd8_signatures.png)

Per-lesion mean scores over 6,699 CD8 T cells (47 lesions with ≥10 CD8 cells;
17 responder, 30 non-responder):

| | CD8_G (memory) | CD8_B (exhaustion) | memory − exhaustion |
|---|---|---|---|
| Responder | 0.196 | 0.095 | **+0.101** |
| Non-responder | 0.047 | 0.346 | **−0.299** |

**Signature = mean(CD8_G) − mean(CD8_B) per lesion** separates responders with
**AUC = 0.827** (Mann–Whitney p = 2.3×10⁻⁴). The single-marker **TCF7⁺ fraction**
proxy reaches **AUC = 0.731**, confirming TCF7 alone carries much of the signal —
consistent with the paper's TCF7 protein-staining result.

![signature](figures/signature_boxplot.png)

### Proposed stratifier (deployable on bulk RNA too)

```
score = mean(TCF7, IL7R, SELL, CCR7, CD28, LEF1, CD27, GZMK)
      − mean(PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, ENTPD1, CD38, TOX)
      (+ optional B-cell / TLS module: MS4A1, CD79A, BANK1)
```

High score → predicted responder. The gene sets are lineage-interpretable and
deconvolvable from bulk RNA-seq, so the signature can be tested beyond scRNA-seq
cohorts (e.g. Riaz, Van Allen anti–PD‑1 cohorts).

## 4. External validation — Riaz anti–PD-1 cohort (bulk, GSE91061)

Applied to bulk RNA-seq of the independent Riaz nivolumab cohort (109 samples;
Responder = CR/PR, Non-responder = PD):

| visit | score | AUC | p |
|---|---|---|---|
| On | **CD8 abundance** | **0.78** | **0.006** |
| On | exhaustion genes | 0.75 | 0.014 |
| On | memory genes | 0.71 | 0.036 |
| On | **state difference (mem − exh)** | 0.39 | 0.27 |
| Pre | exhaustion genes | 0.69 | 0.10 |
| Pre | memory genes | 0.68 | 0.11 |
| Pre | state difference | 0.54 | 0.71 |

**Key result (honest):** the single-cell *state* signature (memory − exhaustion)
does **not** translate to bulk (AUC ≈ 0.4–0.54). In bulk tissue, exhaustion genes
(PDCD1, LAG3…) track overall **T-cell infiltration**, which is favorable — the
opposite of their within-CD8 meaning — so the difference cancels the predictive
signal. What validates in bulk is CD8/immune **abundance**, strongest
on-treatment (AUC 0.78). The state signature is a single-cell property; validating
it properly requires single-cell (not bulk) external data — motivating section 6.
Tables: `riaz_validation_stats.csv`, `riaz_signature_scores.csv`.

## 5. CD8 sub-states (GSE120575)

Sub-clustering the 6,699 CD8 T cells yields five states. Enrichment mirrors outcome:

| CD8 state | n | CD8_G | CD8_B | log2FC (R/NR) |
|---|---|---|---|---|
| Cytotoxic/effector | 636 | 0.02 | 0.27 | **+0.48** |
| Naive/memory | 1777 | 0.19 | 0.07 | **+0.42** |
| Effector-memory | 2785 | 0.09 | 0.42 | +0.04 |
| Exhausted | 809 | 0.01 | 0.49 | **−0.63** |
| Proliferating | 692 | −0.15 | 0.39 | **−1.43** |

Memory scores peak in naive/memory; exhaustion peaks in exhausted/proliferating.
The **proliferating** (cell-cycle + exhaustion) and **exhausted** states are the
most strongly non-responder-enriched — matching the paper's G11 cluster.
Table: `cd8_substate_summary.csv`; figures: `umap_cd8_states.png`,
`cd8_state_enrichment.png`.

## 6. Cross-cancer gene-set validation + TCR clonal dynamics (Yost, GSE123813)

Basal/squamous cell carcinoma, paired scRNA + TCR, pre/post anti–PD-1.

**(a) Gene sets generalize across tumor type.** Scoring CD8_G/CD8_B on Yost's CD8
states: memory score peaks in `CD8_mem` (0.30) and `Naive` (0.19); exhaustion score
peaks in `CD8_ex_act` (1.70) and `CD8_ex` (1.11). The signature's cell-state
discrimination transfers to a different cancer. Table: `yost_geneset_by_state.csv`.

**(b) Clonal dynamics.** Using TCR CDR3 clonotypes (28,371 T cells with a TCR):
- Exhausted CD8 states are the **most clonally expanded** (CD8_ex 0.76, CD8_ex_act 0.80).
- Activated-exhausted CD8 clonal expansion rises **0.40 → 0.80** from pre to post therapy.
- On average **54% of post-treatment expanded clones are novel** (not detected
  pre-treatment) — direct evidence of **clonal replacement** rather than
  reinvigoration of pre-existing clones, reproducing Yost et al.

Tables: `yost_clonal_expansion_by_state.csv`, `yost_clonal_replacement.csv`;
figures: `yost_geneset_validation.png`, `yost_tcr_dynamics.png`.

## Shareable page

A self-contained visual summary is generated by `analysis/make_artifact.py`
(`results/artifact.html`).

## Outputs

- `SUMMARY.txt` — headline numbers
- `population_enrichment_R_vs_NR.csv`, `population_enrichment_Post_vs_Pre*.csv`
- `cluster_marker_genes.csv`, `cluster_lineage_scores.csv`
- `per_lesion_cd8_signature.csv`, `per_lesion_TCF7_fraction.csv`
- `figures/*.png`

## Reproduce

```bash
python3 -m venv .venv && .venv/bin/pip install -r analysis/requirements.txt
bash analysis/download_data.sh          # fetches GSE120575 into data/ (not committed)
.venv/bin/python analysis/run_scrna_analysis.py
```
