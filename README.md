# lipidlib — AI-guided lipid / targeting-ligand library screening

A working repository for using machine learning to design and screen molecular
libraries for **LNP (lipid nanoparticle) mRNA/RNA delivery**.

The project is deliberately split into **two distinct design problems** that are
often conflated (see [`docs/PLAN.md`](docs/PLAN.md) for the full rationale):

| | **Problem A — ionizable-lipid design** | **Problem B — active-targeting ligand design** |
|---|---|---|
| What you optimize | the ionizable lipid *inside* the LNP | a separate ligand *conjugated to the LNP surface* |
| Governs | encapsulation, endosomal escape, potency, **passive** organ tropism (via protein corona) | **active** receptor targeting / receptor-mediated endocytosis |
| Right model class | structure→property QSAR (graph NN / descriptors), MD-derived shape | protein–ligand binding / docking / contrastive retrieval |
| Reference method | **LiON** (Witten 2024) & spatial-conformation ML (Su 2026) | **DrugCLIP** (Gao/Lan 2023–2025), Boltz, docking |
| Endocytosis role | shapes corona → *which organ* takes up the particle | *which cell-surface receptor* triggers uptake |

Both papers you supplied address **Problem A**. The "small molecule binds a CD
ectodomain → conjugate to the liposome → targeted lipid" idea is **Problem B**.
The ionizable lipid is *not* the recognition element and is *not* extracted from
the particle — that concern is unfounded (details in the plan).

## Status

**Phase 1 in progress — track A (LiON reproduction) + Modal, with a GLP1R Problem-B scaffold.**

Delivered:
- **Featurization** (`lipidlib/featurize.py`): SMILES → Morgan/MACCS/RDKit-2D + CLI.
- **LiON reproduction (Problem A)**: LNP_ML cloned; data schema mapped (13,331 pts,
  smiles + 36 X-features → `quantified_delivery`); **Modal training/screen app**
  (`modal_app/lion.py`) faithful to the upstream `split→train→analyze→predict`
  pipeline (Chemprop 1.7.0, GPU); screen-library builder (`lipidlib/lion_library.py`)
  matching the exact `predict` input schema (schema-locked by `tests/`).
- **GLP1R pilot (Problem B)**: 1,422 real GLP1R ligands pulled from ChEMBL
  (`scripts/fetch_glp1r_ligands.py` → `data/targets/GLP1R/`); strategy +
  Boltz-validation plan in [`docs/PROBLEM_B_GLP1R.md`](docs/PROBLEM_B_GLP1R.md).
- Reference lipids, resource inventory, full plan under `docs/`.

## Quickstart

```bash
pip install -r requirements.txt

# Problem A — featurize + reproduce LiON
python scripts/fetch_reference_lipids.py
python -m lipidlib.cli data/reference_lipids.csv --fp morgan -o out.npy
bash scripts/fetch_lnp_ml.sh                         # clone the LiON repo
# then on a Modal-enabled machine:
#   modal run modal_app/lion.py::train --split all_random_split_for_paper --epochs 30

# Problem B — GLP1R target data
python scripts/fetch_glp1r_ligands.py --min-pchembl 6

python -m pytest tests/ -q
```

See [`docs/PLAN.md`](docs/PLAN.md), [`docs/RESOURCES.md`](docs/RESOURCES.md),
[`docs/PROBLEM_B_GLP1R.md`](docs/PROBLEM_B_GLP1R.md).
