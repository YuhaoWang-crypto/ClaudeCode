---
name: lion-modal
description: >-
  Train and run the LiON ionizable-lipid potency model (Witten 2024 D-MPNN,
  Chemprop) on Modal GPUs, and score lipid libraries with it. Use when the user
  wants to (re)train the LNP-delivery model, pick a lightweight vs full run, or
  screen SMILES against a trained model. Handles the Python-3.8/Chemprop-1.7 vs
  Modal-3.10 conflict and preemption-resilient training. Triggers: "train LiON",
  "train the lipid model", "screen lipids on Modal", "reproduce Witten/LiON",
  "lite vs full training".
---

# LiON on Modal (train / screen infrastructure)

`modal_app/lion.py` runs the upstream LNP_ML (LiON) Chemprop pipeline on Modal.
Chemprop 1.7 needs Python 3.8 but Modal's builder needs ≥3.10, so the Modal
function runs in 3.10 while Chemprop lives in a `uv`-built 3.8 venv (invoked via
subprocess). LiON's hardcoded 5-fold loop is made configurable via `LION_CV_NUM`.

## Prerequisites
- `pip install modal 'modal[api-proxy-support]'`; `MODAL_TOKEN_ID/SECRET` set (or
  `modal token new`). Behind a proxy, `python-socks` is required.
- The upstream data is cloned inside the image at build (no local clone needed);
  for local inspection use `scripts/fetch_lnp_ml.sh`.

## Two training tiers
```bash
# lite: 1 fold, 15 epochs, cheap T4 — ~6 min, RMSE ~0.90 (rough pre-filter model)
modal run modal_app/lion.py::train_lite

# full ensemble: 5 folds, 30 epochs, A10G — RMSE 0.78–0.86, preemption-RESILIENT
# (trains fold-by-fold, persists each, resumes after worker preemption)
modal run modal_app/lion.py::train_resilient
```
`train` (non-resilient) also exists; prefer `train_resilient` for the full run
because a single 40-min job gets preempted and restarts from scratch otherwise.

## Smoke test / analyze / screen
```bash
modal run modal_app/lion.py::demo_screen          # toy split train+predict, end-to-end check
modal run modal_app/lion.py::analyze              # held-out Pearson/Spearman/Kendall
modal run modal_app/lion.py::screen --smiles-json '["<smi>", ...]' \
    --folds <1|5> --context <KK_HeLa|lung_epithelium_IT|liver_IV|spleen_IV> --name <run>
```
`--folds` MUST match the trained model (1 for `train_lite`, 5 for the ensemble).
Predictions persist to the `lion-models` volume as `screen__<name>.csv` and the
per-fold `pred_file.csv` (columns: `cv_i_pred_delivery`, `avg_pred_delivery`).

## Data schema (for building screen inputs)
LiON conditions on the ionizable-lipid SMILES + 36 context features (formulation
molar ratios + one-hot delivery-target / helper-lipid / route / cargo / cell). Build
valid screen inputs with `lipidlib.lion_library.build_library(smiles, out, name,
formulation)` — schema is locked by `tests/test_lion_library.py`.

## Build gotchas already handled in the image
Python-3.8 venv via `uv`; `setuptools`/`wheel` added (hyperopt needs `pkg_resources`);
apt `libxrender1 libxext6 libsm6 libglib2.0-0` for RDKit Draw; screen output path is
`<split>_preds/<library>/pred_file.csv`. See `docs/PLAN.md`.
