# IPF Target Discovery — Full Report

A foundation-model + MCP pipeline that nominates novel drug targets for a chosen
indication. This report documents (1) every model/tool used — what it is, what it
does, its inputs and outputs; (2) the real results produced for Idiopathic
Pulmonary Fibrosis (IPF); and (3) how to run the same pipeline on a new
indication.

---

## 1. The funnel

```
 Step 1  Disease map & landscape   PubMed/bioRxiv + ClinicalTrials + ChEMBL   (MCP)
 Step 2  Hypothesis engine         Geneformer in-silico perturbation          (Modal GPU)
 Step 3  Novelty gate              crowded_targets.csv (from ClinicalTrials)
 Step 4  Druggability / structure  ChEMBL bioactivity + Boltz-2.1             (MCP)
 Step 5  Target dossier            Claude synthesis  ->  pilot_dossier.md
```

---

## 2. Models & tools — role / input / output

### 2.1 Geneformer-V2-104M  (from NVIDIA-BioNeMo/bionemo-recipes)  — the engine
- **What it is**: a transformer foundation model pretrained on ~104M single cells
  (rank-value "sentences" of genes ordered by expression). 104.4M params, vocab
  20,275 genes, context 4096.
- **Role**: *counterfactual* target discovery. It answers "what does deleting
  gene X *do* to a diseased cell's state" — not mere correlation with disease.
- **Input**: single-cell raw counts (Ensembl IDs) → tokenized to rank-value
  encoding; a `cell_states_to_model` spec (start = disease, goal = healthy).
- **Output**: per-gene **goal-shift score** — how far an in-silico knockout moves
  diseased cells toward the healthy-cell embedding. Ranked list of candidate
  genes.
- **Where it runs**: Modal GPU (`modal/geneformer_app.py`, `modal/perturb_app.py`).
  Verified loading + forward pass on a T4; perturbation on an A100.

### 2.2 ESM-2  (from bionemo-recipes)  — protein readout  *(scaffolded, not yet run)*
- **What it is**: protein language model (embeddings from amino-acid sequence).
- **Role**: characterize a nominated target protein — family, likely fold,
  variant-effect scoring, druggable-domain signal.
- **Input**: protein sequence (one-letter). **Output**: per-residue / per-protein
  embeddings; variant effect scores. *(Planned layer; the config references it.)*

### 2.3 Boltz-2.1  (Boltz MCP)  — structure, binding, screening  — RUN, real data
- **What it is**: a structure + binding-affinity prediction model (AlphaFold-class
  co-folding with a binding readout).
- **Role**: turn a "target" into an "attackable target" — predict the fold, find
  the pocket, dock ligands, rank chemotypes, design binders.
- **Input**: protein sequence(s) (+ optional ligand SMILES / a SMILES library).
- **Output**: 3D structure (CIF), confidence metrics (pTM, iPTM, pLDDT,
  ligand_iPTM), a `binding_confidence` / `optimization_score`, and per-molecule
  ADME for screens.

### 2.4 ChEMBL MCP  — druggability & bioactivity  — RUN, real data
- **Role**: is the target druggable, and with what chemical matter?
- **Input**: gene symbol / target / compound IDs.
- **Output**: target class + PDB count, known ligands with IC50/Ki/Kd + pChEMBL,
  ligand efficiency, ADMET.

### 2.5 PubMed / bioRxiv MCP  — literature  — RUN
- **Role**: disease mechanism map; functional/genetic evidence per candidate.
- **Input**: query strings. **Output**: ranked articles (PMIDs / DOIs), abstracts.

### 2.6 ClinicalTrials.gov MCP  — competitive landscape / novelty  — RUN
- **Role**: is this target already crowded in the clinic? (the novelty gate)
- **Input**: condition / intervention / phase. **Output**: trials with
  sponsor, phase, mechanism → `crowded_targets.csv`.

### 2.7 Claude  — orchestration & synthesis
- **Role**: pick the indication, define cell states, write/debug the Modal jobs,
  weigh evidence across all sources, and assemble the ranked dossier. Not a
  replacement for the models — the conductor.

---

## 3. What we actually produced for IPF

**Indication**: IPF — orphan disease, only 2 approved drugs (both slow, don't
reverse), rich single-cell data.

- **Landscape (ClinicalTrials, 55 Ph2/3 trials)**: crowded mechanisms = PDE4B,
  LPA1, integrin αvβ6, ROCK2, CSF1R, TG2 → `data/crowded_targets.csv` (novelty
  gate). Whitespace = aberrant-basaloid drivers, fibroblast-fate TFs, senescence.
- **Evidence-triaged shortlist** (`data/seed_shortlist.csv`): PTGES, MDK, SFRP2,
  CDKN2A, PRRX1, TWIST1 with ChEMBL druggability.
- **PTGES / mPGES-1** (ChEMBL): best IC50 **1 nM**, LE 0.40–0.46 → strong
  small-molecule target. **Boltz**: monomer + inhibitor conf 0.89; **homotrimer**
  run confirmed the inhibitor sits in a well-defined inter-subunit pocket
  (ligand_iPTM 0.78→0.85, interface error 3.1→1.35 Å).
- **PTGES inhibitor screen** (Boltz, `data/ptges_screen_results.csv`):
  `optimization_score` cleanly separates 5 actives (0.28–0.42) from decoys
  (0.00–0.01); indole-acetic-acid chemotype best on binding + ADME.
- **MDK / midkine**: ChEMBL nM binders only as sulfated-glycan mimetics; Boltz apo
  flexible/no deep pocket → **biologic** route.
- **Geneformer goal-shift** (`data/geneformer_goalshift.csv`): first real
  in-silico deletion ranking on 200 fibrosis cells. **SFRP2 (a whitespace
  candidate) ranks #1** for shifting fibrotic cells toward normal.

### Upgrade: fine-tuned classifier (verifiable-grade)
Re-run on **16k cells** with a **fine-tuned CellClassifier** (disease vs normal,
**94.9% eval accuracy**) — `data/geneformer_goalshift_finetuned.csv`:
- **Known fibrosis effectors now top the ranking** — FAP, POSTN, CTHRC1, TGFB1
  (deleting canonical drivers pushes cells toward normal — biologically correct).
- **Housekeeping controls sink to the bottom** (ACTB, B2M) — controls behave
  correctly, which the pretrained pass failed. Signal ~10× stronger.
- All 6 whitespace candidates (TWIST1, PTGES, CDKN2A, PRRX1, MDK, SFRP2) score
  positive; SFRP2 remains a positive hit though the classifier ranks the dominant
  ECM/TGF-β effectors higher.
This is the "noisy → verifiable" upgrade: fine-tuning gives the model a learned
disease boundary, so goal-shift is measured in classifier space.

### Honest read on the original (pretrained) Geneformer result
This was a **working but noisy first pass**, not a validated screen:
- Magnitudes are tiny (~1e-5); single-gene KO on a **pretrained** (not fine-tuned)
  model over only **200 cells** gives weak signal.
- Controls are imperfect: housekeeping (ACTB/B2M) and collagens land at the
  *bottom*, not cleanly neutral — so ranking is **directional, not definitive**.
- SFRP2 #1 is biologically sensible (pathologic-fibroblast Wnt modulator, Tsukui
  2020) and encouraging, but should be treated as a hypothesis to confirm.
- **To sharpen**: fine-tune a Geneformer CellClassifier on the disease/normal
  labels, use thousands of cells, and average more perturbation samples;
  optionally run the full genome-wide sweep (needs a longer GPU budget than the
  targeted run used here).

---

## 4. Playbook — applying this to a NEW indication

The pipeline is indication-agnostic. To retarget it:

**A. Landscape & whitespace (MCP, no GPU)**
1. ClinicalTrials `search_trials(condition=<disease>)` → build `crowded_targets.csv`.
2. PubMed/bioRxiv → disease mechanism map, cell states, candidate genes.
3. ChEMBL `drug_search(indication=<disease>)` → approved/known targets.

**B. Data (Modal `perturb_app.py`)**
4. `inspect_census` → find the exact Census `disease` label + `dataset_id` and the
   relevant `cell_type` labels for your disease (labels differ per disease!).
5. Edit `PF_DATASET`, `LINEAGE`, and `STATES` (start = disease label, goal =
   `normal`) in `perturb_app.py`. Run `fetch_ipf_data`.

**C. Engine (Modal GPU)**
6. `tokenize` → rank-value dataset (unchanged).
7. `perturb_targeted` with a `CANDIDATES` list for your disease (whitespace +
   known effectors as positive controls + housekeeping as negatives), or
   `perturb` (genome-wide) if you have the GPU budget. Then `read_shifts`.

**D. Druggability & dossier (MCP)**
8. Top hits → ChEMBL `target_search`/`get_bioactivity`, Boltz
   `structure_and_binding` / `small_molecule_screen`, PubMed evidence.
9. Claude assembles the ranked `pilot_dossier.md`.

**Only these change per indication**: the ClinicalTrials/PubMed queries, the
Census `disease` label + `dataset_id` + `LINEAGE`, the `STATES` start label, and
the `CANDIDATES` gene list. Everything else (tokenizer, model, Boltz calls,
scoring) is reused as-is.

---

## 5. Compute, cost, artifacts

- **Modal**: T4 (smoke), A100 (perturbation). Image + 418 MB checkpoint cached in
  a Volume; data cached in `gf-ipf-data` Volume.
- **Boltz**: ~$0.35 total this session (structure preds + inhibitor screen).
- **Key files**: `modal/geneformer_app.py`, `modal/perturb_app.py`,
  `data/*.csv`, `pilot_dossier.md`, `src/` (the MCP-orchestrated funnel scaffold).
- **Setup gotchas solved** (see `modal/README.md`): proxy needs python-socks;
  install via full clone (HF git promisor breaks pip partial clone);
  `huggingface_hub.snapshot_download` not git-lfs; pin `transformers==4.46.3`;
  V2 needs `emb_mode='cls'`; gene LIST = combined deletion (loop for per-gene);
  MLM decoder OOMs at high batch (use A100 + small batch).
