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

Phase 0 (scaffolding) — in progress. This turn delivers:
- Featurization module (`lipidlib/featurize.py`): SMILES → ECFP/Morgan, MACCS, RDKit 2D descriptors.
- Reference ionizable-lipid dataset builder (`scripts/fetch_reference_lipids.py`, pulls canonical SMILES from PubChem).
- Integration notes + fetch helpers for LNP_ML, DrugCLIP, Boltz, ChEMBL, Modal.
- Full technical plan and resource inventory under `docs/`.

## Quickstart

```bash
pip install -r requirements.txt
python scripts/fetch_reference_lipids.py            # build data/reference_lipids.csv from PubChem
python -m lipidlib.cli data/reference_lipids.csv --smiles-col smiles --fp morgan -o out.npy
```

See [`docs/PLAN.md`](docs/PLAN.md) and [`docs/RESOURCES.md`](docs/RESOURCES.md).
