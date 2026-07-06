---
name: scrna-immunotherapy-signatures
description: >-
  Reanalyze single-cell RNA-seq of tumor biopsies before/after immunotherapy (immune
  checkpoint blockade) to find cell populations that expand/contract with response,
  derive marker genes, build a CD8 memory-vs-exhaustion responder signature, validate
  it on external cohorts (bulk or single-cell), and analyze TCR clonal dynamics. Use
  when the user wants to: pull a scRNA-seq immuno-oncology dataset from GEO (e.g.
  GSE120575, GSE123813, GSE91061), cluster/annotate the immune compartment, compare
  responders vs non-responders or pre vs post treatment, score response signatures,
  do CD8 sub-state or TCR/clonal-replacement analysis, or turn such an analysis into
  a shareable report/Artifact. Triggers: "scRNA-seq immunotherapy", "checkpoint
  blockade single cell", "responder signature", "TCR clonal dynamics", "tumor
  biopsies before and after immunotherapy", "GEO reanalysis".
---

# scRNA-seq immunotherapy response signatures

A reproducible pipeline for single-cell tumor immunotherapy datasets. It was built
and verified end-to-end on Sade-Feldman et al. *Cell* 2018 (melanoma, GEO
**GSE120575**), with external validation on Riaz (bulk, **GSE91061**) and Yost
(BCC/SCC scRNA + TCR, **GSE123813**).

## What it does

1. **Load** a genes × cells expression matrix + per-cell clinical metadata (response,
   timepoint) from GEO.
2. **Cluster & annotate** the immune compartment (Leiden + canonical lineage markers).
3. **Population dynamics** — per-lesion fraction contrasts: Responder vs Non-responder
   and Post vs Pre, with Mann-Whitney p-values and log2 fold-changes.
4. **Marker genes** per population (Wilcoxon).
5. **Responder signature** — score a CD8 memory set (`CD8_G`) vs an exhaustion set
   (`CD8_B`); the per-lesion `mean(CD8_G) − mean(CD8_B)` stratifies responders (ROC AUC).
6. **CD8 sub-states** — sub-cluster CD8 cells into naive/memory, effector-memory,
   cytotoxic, exhausted, proliferating; relate to response.
7. **External validation** — apply the gene sets to an independent cohort.
8. **TCR / clonal dynamics** — clone sizes from CDR3, clonal expansion by state,
   pre→post clonal replacement.
9. **Artifact** — a self-contained HTML report with embedded figures.

## Files (`scripts/`)

| script | purpose |
|---|---|
| `download_data.sh` | fetch GSE120575 (+ GSE91061, GSE123813) into `data/`; extracts Yost signature-gene rows + per-cell library sizes |
| `requirements.txt` | scanpy, leidenalg, igraph, scikit-learn, scipy, pandas, matplotlib |
| `run_scrna_analysis.py` | main pipeline: load → cluster → annotate → dynamics → signature → figures (caches a processed `.h5ad`) |
| `validate_riaz.py` | apply signature components to Riaz bulk RNA-seq; ROC AUC pre/on-treatment |
| `cd8_states_tcr.py` | CD8 sub-states (Part A), cross-dataset gene-set validation (Part B), TCR clonal dynamics (Part C) |
| `make_artifact.py` | build `results/artifact.html` (figures embedded as base64 data URIs) |

## How to run

```bash
python3 -m venv .venv && .venv/bin/pip install -r scripts/requirements.txt
bash scripts/download_data.sh
.venv/bin/python scripts/run_scrna_analysis.py    # ~4 min (matrix parse dominates)
.venv/bin/python scripts/validate_riaz.py
.venv/bin/python scripts/cd8_states_tcr.py
.venv/bin/python scripts/make_artifact.py
```
Outputs land in `results/` (CSVs + `figures/*.png` + `REPORT.md` + `artifact.html`).
Raw matrices and the processed `.h5ad` are large — keep them gitignored; reference
datasets by GEO accession name, not as committed files or links.

> The bundled scripts hard-code `DATA`/`OUT` to `/home/user/ClaudeCode/{data,results}`.
> When reusing in another project, edit those two constants at the top of each script
> (or copy the scripts to `analysis/` and point them at the new project root).

## Signature gene sets (curated from Sade-Feldman markers)

- **CD8_G (memory-like, "good"):** `TCF7, IL7R, SELL, CCR7, CD28, LEF1, CD27, GZMK`
- **CD8_B (exhaustion, "bad"):** `PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, ENTPD1, CD38, TOX`
- Portable stratifier: `mean(CD8_G) − mean(CD8_B)` per lesion; `TCF7` alone is a
  strong single-marker proxy.

## Adapting to a NEW dataset

- **Metadata parsing** is dataset-specific — locate the per-cell columns for
  `response` (map to Responder/Non-responder) and `timepoint` (Pre/Post). For GEO,
  the labels are often in the series matrix (`!Sample_characteristics_ch1`) or a
  supplementary sample table, not the counts file.
- **Gene identifiers vary:** symbols vs Entrez vs Ensembl. Riaz uses Entrez
  (`hg19KnownGene`) — `validate_riaz.py` holds a symbol→Entrez map; extend it if you
  add genes. Always report which signature genes were found vs missing.
- **Normalization:** GSE120575 ships TPM (just `log1p`). 10x UMI counts need
  CP10k + `log1p`. Bulk: `log2(FPKM+1)` then z-score per gene across samples.
- **Response binarization:** RECIST CR/PR = Responder; PD (± SD) = Non-responder.
  Report the exact grouping; SD is ambiguous — try both including/excluding it.

## Hard-won gotchas (read `references/gotchas.md` before debugging)

- **OOM on load:** never build AnnData from a dense DataFrame with several copies.
  Convert to `scipy.sparse.csr_matrix`, `del` intermediates, filter genes early, and
  do **not** keep a `.raw`/dense `tpm` layer. A 55k×16k dense matrix ×4 copies blows
  past 15 GB.
- **Wide-TSV parse quirks:** GSE120575 gene rows have a **trailing tab** (one extra
  column) and the metadata file is **latin-1** (a `µ` byte). GEO 10x counts often
  have an R-style header **off by one** vs data columns — verify field counts before
  aligning library sizes.
- **`dtype='float32'` in `read_csv`** must be applied per-column (exclude the
  gene-name column) or it tries to cast gene symbols to float.
- **Categorical groupby NaN:** after an `.h5ad` round-trip, string obs columns become
  categorical; `groupby(...).unstack()` then injects empty-combination NaN rows that
  make `scipy.mannwhitneyu` return NaN. Use `observed=True` (and `nan_policy='omit'`).
- **matplotlib:** `Axes.boxplot(labels=...)` → use `tick_labels` or
  `set_xticklabels`.

## The single most important scientific caveat

The CD8 **state** signature (memory − exhaustion) is a *within-CD8-cell* property and
does **not** transfer to **bulk** RNA-seq: in bulk, exhaustion genes track overall
T-cell *infiltration* (favorable), the opposite of their single-cell meaning, so the
difference cancels. In bulk, CD8/immune **abundance** predicts response (Riaz
on-treatment AUC ≈ 0.78). Validate a state signature on **single-cell** external data,
or deconvolve CD8 state from CD8 quantity before applying to bulk. Always report this
distinction rather than fishing for a positive bulk AUC.
