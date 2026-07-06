---
name: lipid-library-screen
description: >-
  Design and rank ionizable-lipid libraries for LNP mRNA delivery (Track A). Use
  when the user wants to enumerate a combinatorial lipid library (amine heads ×
  ester/amide linkers × tails, incl. degradable disulfide/ester tails), score it
  with the LiON D-MPNN model, and get a ranked/confidence-aware shortlist —
  optionally for a specific target organ (lung / liver / spleen). Triggers:
  "design lipids", "enumerate lipid library", "screen ionizable lipids", "rank
  lipid candidates", "spleen/liver/lung-targeted lipid".
---

# Lipid library design & screening (Track A)

Pipeline: **enumerate → score on Modal (LiON) → rank (ensemble mean ± std)**.
Scoring uses the reproduced LiON model (see the `lion-modal` skill); this skill
covers enumeration and ranking and calls the Modal `screen` function.

## Prerequisites
- `pip install -r requirements.txt && pip install -e .` (RDKit, pandas, sklearn).
- A trained LiON model on the Modal volume `lion-models` (train it via the
  `lion-modal` skill: `train_lite` for a fast model, `train_resilient` for the
  5-fold ensemble). `MODAL_TOKEN_ID/SECRET` must be set.

## Steps

1. **Enumerate** a library (aza-Michael chemistry; tail count = head reactive N-H):
   ```bash
   python analysis/enumerate_library.py      # v1: 14 heads × 2 × 16 tails → ~323 lipids
   python analysis/enumerate_library_v2.py   # v2: + degradable S–S / ester tails, +heads → ~538
   ```
   To customise blocks, edit the `HEADS` / `TAILS` dicts in those scripts, or call
   `lipidlib.lion_library.enumerate_michael_lipids(heads, tails, linkers, mw_range)`
   directly. Output: `data/libraries/<lib>/<lib>.csv` (+ `_enumerated.csv` with
   head/linker/tail/n_tails/mw annotations).

2. **Score on Modal** with the trained model. Pass `folds=5` for the ensemble (must
   match the trained folds), and a `context` for the organ:
   ```bash
   SM=$(python -c "import pandas,json;print(json.dumps(pandas.read_csv('data/libraries/combo_v1/combo_v1.csv').smiles.tolist()))")
   modal run modal_app/lion.py::screen --smiles-json "$SM" --folds 5 \
       --context spleen_IV --name combo_v1_spleen
   ```
   Contexts: `KK_HeLa` (in-vitro baseline), `lung_epithelium_IT`, `liver_IV`,
   `spleen_IV`. Results land on the volume as `screen__<name>.csv`; fetch with
   `modal volume get lion-models screen__<name>.csv <local>`.

3. **Rank (confidence-aware)**: merge predictions with `_enumerated.csv`, compute
   ensemble mean and cross-fold std, and prioritise **high mean AND low std**
   (folds agree). See `analysis/*ranked_ensemble* / *shortlist_top50*` for the
   pattern, or reuse the ranking cells in the docs.

## Interpreting results
- Recovered SAR sanity checks: unsaturated tails and 3-tail lipids score well
  in vitro; **the SAR shifts by organ** (spleen → 2-tail + small heads + degradable
  S-tails; lung/liver → 3-tail + H7). See `docs/COMBO_LIBRARY.md`,
  `docs/ORGAN_SCREEN.md`, `docs/COMBO_V2_SPLEEN.md`.
- Prefer **low-ensemble-std** picks; the lite (1-fold) model is a rough pre-filter,
  the 5-fold ensemble is the trustworthy ranking.

## Caveats
Predictions are ranked hypotheses; in-vivo/organ training data is sparser, and the
model conditions on context. Validate top picks experimentally.
