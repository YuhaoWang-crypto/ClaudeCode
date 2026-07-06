---
name: biomarker-discovery
description: Discover and prioritize novel protein/gene biomarkers for a disease indication across the FDA-NIH BEST categories (risk, diagnostic, predictive/companion-diagnostic, response/PD, monitoring, prognostic). Combines a single-cell foundation model (Geneformer) for de novo cell-state discovery, protein/regulatory variant scoring (ESM-2 VEP plus PROTO/Evo2/AlphaGenome), mechanistic backup (Boltz/EDEN/ChEMBL) and clinical evidence grounding (PubMed/bioRxiv/ClinicalTrials). Use when the user wants to find or rank biomarkers for an indication, run Geneformer de novo cell-state discovery (locally or on a GPU via Modal), score or track a candidate panel, or port the whole pipeline to a new disease. Triggers include biomarker discovery, find biomarkers, predict drug response, companion diagnostic, Geneformer de novo, cell-state discovery, candidate panel, and port to new indication.
---

# Biomarker Discovery

End-to-end, evidence-grounded workflow to mine **new protein/gene biomarkers** for one indication and map them to clinical use. Only the **verified-working** components are packaged here (Geneformer de novo discovery, GPU run on Modal, candidate scoring/tracking) plus the methodology.

## When to use
- User asks to discover / prioritize biomarkers for a disease (diagnostic, predictive/CDx, response, monitoring, prognostic).
- User wants to run **Geneformer de novo cell-state discovery** on scRNA (CPU demo or real data on GPU).
- User wants to **score/track a candidate panel** or **port the pipeline to a new indication**.

## The method (BEST framework)
Every clinical use maps to an FDA-NIH BEST category. Do not ask "is this a biomarker?" — ask "which BEST category, versus what comparison group, measurable in what sample?". Full playbook + the 6-step process + "swap these 6 things" porting table: `references/workflow.md`. Per-model use/input/output: `references/models-io-reference.md`.

**Core loop:** define question + comparison group → ground in evidence (PubMed/bioRxiv/ClinicalTrials/ChEMBL via MCP tools if available) → generate candidates on 3 omic legs (Geneformer scRNA / ESM-2 coding variants / PROTO regulatory variants) → mechanistic filter (Boltz/EDEN/ChEMBL) → clinical triage (novelty/accessibility/alignment) → combine into per-mechanism scorecards. Multi-modal cross-corroboration is the key false-positive control.

## Scripts (all verified to run)

### 1. Geneformer de novo cell-state discovery — `scripts/geneformer_denovo_discovery.py`
Discovers cell states, markers, an inflammation/disease module, and condition enrichment; has an in-silico perturbation hook.
- **CPU demo (runs anywhere, no data/GPU needed):**
  `python3 scripts/geneformer_denovo_discovery.py --synthetic --backend pca`
  Needs: `numpy scipy pandas scikit-learn`. Verified: recovers an IL13RA2/CXCL13/TNFRSF11B inflammatory-fibroblast state, ~1.9x inflamed enrichment.
- **Real Geneformer embeddings (GPU):** `--h5ad DATA.h5ad --backend geneformer` (needs `geneformer`, `torch`, weights). Adapt the `INFLAMMATION_MODULE` list to the target disease.
- Input: scRNA counts (cells×genes, Ensembl IDs + n_counts). Output: `denovo_discovery_report.json` + `denovo_markers.csv`.

### 2. GPU Geneformer on Modal — `scripts/modal_geneformer.py`
Runs the **real pretrained Geneformer-V1-10M** (+ official gc30M vocab, rank-value encoding) on a Modal T4, returns cell embeddings + de novo clustering.
- `python3 -m modal run scripts/modal_geneformer.py` (needs `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET`; `pip install 'modal[api-proxy-support]'`).
- Verified working end-to-end (device=cuda, BertModel, 256-d embeddings). IMPORTANT lesson (see `references/models-io-reference.md`): synthetic random-gene data yields near-chance separation because Geneformer encodes *real* learned co-expression — genuine discovery needs **real scRNA** (GSE134809 Martin CD; SCP259 Smillie UC — see `references/scRNA_datasets.md`). Use the smallest general **pretrained** model, not a fine-tuned classifier, and match the vocab to the model (V1↔gc30M, V2↔gc104M).

### 3. Candidate scoring / tracking — `scripts/score_candidates.py`
Reads a candidate panel + PROTO variant list, computes priority (evidence grade, novelty, accessibility, contradiction penalty, gap-opportunity bonus), routes each variant to ESM-2/PROTO/Boltz engine skeletons, writes a tracking CSV.
- `python3 scripts/score_candidates.py --plan-variants` (needs `pyyaml`). Defaults read `references/candidates.yaml` + `references/proto_variants_il23r_tyk2.tsv`; override with `--candidates/--variants`.

## Reference files
- `references/workflow.md` — BEST framework, 6-step process, "swap 6 things" porting guide, worked IBD example.
- `references/models-io-reference.md` — every model's use/input/output + data-flow.
- `references/candidates.yaml` — machine-readable candidate schema (example: IBD use case A, 10 candidates with cited evidence cards).
- `references/proto_variants_il23r_tyk2.tsv` — PROTO variant-scoring input format.
- `references/scRNA_datasets.md` — open scRNA datasets + how to run on real data.

## Porting to a new indication
Keep the 6 steps; swap only: (1) indication + comparison groups, (2) scRNA/omics dataset (GEO / Single Cell Portal / CELLxGENE), (3) mechanism/drug classes, (4) disease-specific gene module (`INFLAMMATION_MODULE`), (5) candidate variant list (GWAS Catalog / Open Targets), (6) validation cohorts (ClinicalTrials NCTs). Everything else (schema, scoring, pipelines) is indication-agnostic. Details in `references/workflow.md`.

## Guardrails
- Attribute PubMed and include DOI links for every cited article.
- Discovery ≠ validation: candidates need analytical + independent prospective clinical validation. Flag evidence-grade-1 items as genetics-supported hypotheses.
- Record contradictions explicitly (e.g. discovery-cohort positive vs RCT-negative) and stratify by subtype/age/sample; watch batch/center effects.
- Do not fabricate rsIDs or citations; mark unconfirmed ones `to_validate`.
