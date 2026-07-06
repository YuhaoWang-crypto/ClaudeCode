---
name: target-discovery
description: >-
  Discover novel drug targets for a disease indication by pairing Geneformer
  single-cell in-silico perturbation (on Modal GPU) with an MCP evidence layer
  (ClinicalTrials, PubMed, ChEMBL, Boltz). Use when the user wants to find,
  mine, or nominate new/first-in-class drug targets for an indication; run
  in-silico gene knockout / perturbation on single-cell data to rank genes;
  assess target druggability or novelty; or design a small-molecule/biologic
  starting point. Runs the full funnel: landscape -> hypothesis engine ->
  novelty gate -> druggability -> ranked target dossier.
---

# Target discovery (foundation model + evidence layer)

A counterfactual pipeline: the engine (Geneformer) ranks genes by how much their
in-silico knockout shifts diseased cells toward the healthy state; the evidence
layer (MCP databases) keeps only what is novel, druggable, and tractable.
Every stage below was verified end-to-end on IPF. See `perturb_app.py` (bundled)
for the Modal engine.

## Prerequisites
- **Modal** for GPU: `pip install modal python-socks` (the `python-socks` extra
  is required behind a proxy) and `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` set.
- **MCP servers** (search with ToolSearch): ClinicalTrials, PubMed/bioRxiv,
  ChEMBL, Boltz_API. Load schemas via ToolSearch before calling.

## Two layers

### A. Evidence layer — MCP, no GPU (do this first; it's fast and free)
1. **Landscape / novelty gate** — `Clinical_Trials.search_trials(condition=<disease>,
   phase=[PHASE2,PHASE3])`. Extract the mechanisms already in the clinic; these
   become the crowded set to DOWN-weight (novelty gate).
2. **Disease map** — `PubMed.search_articles` + `bioRxiv` for pathology, cell
   states, and candidate genes; note the single-cell atlases mentioned.
3. **Known targets** — `ChEMBL.drug_search(indication=<disease>)`.
Output a `crowded_targets.csv` (gene, mechanism, phase) and a candidate gene list
(whitespace + known effectors as positive controls + housekeeping as negatives).

### B. Hypothesis engine — Modal GPU (`perturb_app.py`)
Run these in order (a gene LIST is a *combined* deletion in geneformer, so the
targeted stage loops one gene at a time; genome-wide uses `genes_to_perturb=all`):

```bash
modal run perturb_app.py::inspect_census      # 1. find the exact Census disease
                                              #    label, dataset_id, cell lineages
# 2. EDIT perturb_app.py: PF_DATASET, LINEAGE, STATES (start=disease label,
#    goal='normal'), CANDIDATES  -- these are the only per-indication changes
modal run perturb_app.py::fetch_ipf_data --max-per-state 8000   # 3. raw counts -> Volume
modal run perturb_app.py::tokenize                              # 4. rank-value tokens + HVGs
modal run perturb_app.py::finetune                             # 5. disease CellClassifier
modal run perturb_app.py::perturb_targeted \
    --model-type CellClassifier --model-path /data/ft_model \
    --out-name perturb_ft --max-ncells 200                     # 6. in-silico KO
modal run perturb_app.py::read_shifts --out-name perturb_ft    # 7. per-gene goal-shift ranking
```

Interpret: positive shift = KO pushes diseased cells toward normal. **Validate**
by checking known effectors rank high and housekeeping genes rank low/negative —
if not, the run is noisy (use the fine-tuned classifier, more cells).

### C. Druggability + dossier — MCP
For top hits: `ChEMBL.target_search`/`get_bioactivity` (ligandability, IC50);
`Boltz` `structure_and_binding` (pocket), `small_molecule_screen` (rank a
chemotype panel — set `boltz_smarts_catalog_filter_level: disabled` for
halogen-rich series), or `protein_design` (`boltz_nanobody`/`boltz_antibody`
for secreted targets); `PubMed` for genetic evidence. Rank into a dossier.

## Adapting to a NEW indication — only 4 things change
1. ClinicalTrials / PubMed queries (the disease name).
2. In `perturb_app.py`: `PF_DATASET`, `LINEAGE` (cell types), from `inspect_census`.
3. `STATES` start_state = the Census disease label (goal usually `'normal'`).
4. `CANDIDATES` gene list.
Everything else — tokenizer, model, Boltz calls, scoring — is reused as-is.

## Gotchas (already fixed in perturb_app.py)
- Proxy: Modal client needs `python-socks`.
- Geneformer install: full `git clone` then `pip install ./dir` (HF git promisor
  breaks pip's partial clone); pull weights with `huggingface_hub.snapshot_download`,
  not git-lfs.
- Pin `transformers==4.46.3` (newer drops top-level `SpecialTokensMixin`).
- V2 checkpoint needs `emb_mode='cls'`.
- CellClassifier path: `nproc=1` and move `state_embs` to CPU (CUDA-cannot-
  reinit-in-fork during geneformer's `.map`).
- MLM decoder OOMs at high batch → A100 + `forward_batch_size` 16.
- Census `disease` labels differ per disease (e.g. IPF = `'pulmonary fibrosis'`);
  always confirm with `inspect_census` first.

## Cost / compute
Modal T4 (smoke) / A100 (perturbation, fine-tune); checkpoint + data cached in
Volumes so reruns are fast. Boltz jobs are ~$0.02–0.05 each.
