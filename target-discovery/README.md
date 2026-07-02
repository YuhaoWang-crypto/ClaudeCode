# IPF Target Discovery Pipeline

An end-to-end, foundation-model-driven pipeline to nominate **novel drug targets**
for **Idiopathic Pulmonary Fibrosis (IPF)**, combining:

- **bionemo-recipes** foundation models (Geneformer, ESM-2) as the *hypothesis engine*
- **Claude + MCP tools** (PubMed, ChEMBL, ClinicalTrials.gov, Boltz) as the *evidence & triage layer*

The core idea: don't look for *correlation* with disease — do **counterfactual
in-silico perturbation** with Geneformer to find genes whose virtual knockout
pushes pathological cell states back toward health, then gate them against the
competitive landscape (novelty) and druggability before committing.

---

## Why IPF

- Orphan / rare disease, prevalence ~13–20 / 100k.
- Only 2 approved drugs (nintedanib, pirfenidone) — both only *slow* decline.
  Median survival 3–5 years. Huge unmet need → high value for first-in-class.
- Rich single-cell atlases (Habermann 2020, Adams 2020) → Geneformer perturbation
  has data to work with.

## The funnel

```
 Step 1  Disease map            PubMed / bioRxiv  (Claude-orchestrated)
 Step 2  Hypothesis engine      Geneformer in-silico perturbation   <-- this repo
         + protein readout      ESM-2 embeddings / variant effect
 Step 3  Novelty gate           ClinicalTrials.gov + ChEMBL  (crowded_targets.csv)
 Step 4  Druggability / triage  ChEMBL + Boltz structure/binding
 Step 5  Target dossier         Claude synthesis  ->  results/
```

`crowded_targets.csv` and `seed_shortlist.csv` under `data/` are the concrete
outputs of Steps 1 & 3 already run via MCP tools (see repo history / chat).
They seed the novelty gate and validate the triage template end-to-end.

---

## Pipeline scripts (`src/`)

| Script | Stage | Needs GPU | Needs data download |
|---|---|---|---|
| `01_download_data.py`        | fetch IPF scRNA atlas (CELLxGENE / GEO)     | no  | yes |
| `02_tokenize.py`             | tokenize to Geneformer rank-value encoding  | no  | — |
| `03_insilico_perturbation.py`| virtual KO, shift-toward-healthy scoring     | **yes** | — |
| `04_rank_and_gate.py`        | rank by shift, apply novelty gate           | no  | — |
| `05_evidence_triage.py`      | ChEMBL druggability + PubMed evidence (MCP) | no  | — |

Run in order; each writes to `results/` and reads the previous stage's output.

## Quickstart

```bash
pip install -r requirements.txt
python src/01_download_data.py   --config config.yaml
python src/02_tokenize.py        --config config.yaml
python src/03_insilico_perturbation.py --config config.yaml   # GPU node
python src/04_rank_and_gate.py   --config config.yaml
python src/05_evidence_triage.py --config config.yaml         # emits the dossier
```

## Model provenance

- **Geneformer** — from `NVIDIA-BioNeMo/bionemo-recipes` (`recipes/geneformer`).
  Fine-tune / load per the recipe; this pipeline uses the checkpoint for
  in-silico perturbation via the `geneformer` package's `InSilicoPerturber`.
- **ESM-2** — from `bionemo-recipes` (`models/esm2`), HF checkpoint
  `facebook/esm2_t33_650M_UR50D` or the 15B TE-accelerated variant, for protein
  representation and variant-effect scoring of nominated targets.
