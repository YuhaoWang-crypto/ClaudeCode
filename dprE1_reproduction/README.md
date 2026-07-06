# Open-Source Reproduction — Generative-AI + Structure-Based DprE1 Inhibitor Design

An open-source reproduction of the computational workflow in:

> Chikhale, Cabayé, Goupil-Lamy, Popovic. *Generative AI and Structure-Based
> Workflow for the De Novo Design and Optimization of DprE1 Inhibitor
> Candidates.* chemRxiv 2026, DOI: 10.26434/chemrxiv.15004861/v2

## Why "reproduction" and not "re-run"

The original study runs entirely on **commercial, closed-source software** —
BIOVIA Generative Therapeutic Design (GTD), GOLD docking, Discovery Studio and
CHARMm. Those cannot be re-executed without expensive licences, and the GTD
generative algorithm itself is proprietary. This project rebuilds the
**scientific core** with an open-source stack and checks the reproduced numbers
against the paper's reported values. Everything needed is public:

* training data + drawn structures — the authors' Zenodo release
  (DOI 10.5281/zenodo.20340949)
* candidate SMILES (TCA1, GTD_9.1–9.10) — the preprint Appendix-2
* protein structures — PDB 4KW5 (DprE1 WT), 5OEL (Y314C mutant), 5W0C (CYP2C9)

| Paper component | Original tool | Open-source substitute here |
|---|---|---|
| Activity model (DprE1 v2) | GTD Random Forest | RDKit descriptors + ECFP4 + scikit-learn RF |
| De novo generation (GFSP) | BIOVIA GTD | *not reproduced — proprietary algorithm* |
| Docking | GOLD (ChemPLP) | AutoDock Vina / smina (Part 2) |
| MD 500 ns + MM-GBSA | CHARMm / Discovery Studio | OpenMM + AMBER on Modal, via *making-it-rain* (Part 3) |

## Part 1 — Activity model + candidate analysis  ✅ reproduced

```
pip install -r requirements.txt
python src/01_prepare_dataset.py      # Zenodo xlsx + .cdx  -> dataset
python src/02_build_activity_model.py # RF model, 3-iteration CV ROC AUC
python src/03_candidate_analysis.py   # Table-2 property check + RF scoring
```

### Results vs paper

| Quantity | Paper | Reproduced |
|---|---|---|
| IC50 molecules used | 406 | **396 parsed / 366 after structure resolution + dedup** |
| Actives (pIC50 ≥ 5.75) | 192 | **161** |
| DprE1 v2 ROC AUC | 0.92 | **0.902 ± 0.005** (5-fold CV, 3 iterations) |
| GTD_9.x MolWt | Table 2 | **mean abs error 0.4 Da** |
| GTD_9.x MolPSA (≈TPSA) | Table 2 | near-exact (e.g. 94.3/94.3, 116.0/116) |
| GTD_9.x rotatable bonds | Table 2 | 82% exact match |

The ROC AUC is reported here under 5-fold cross-validation (more conservative
than a single held-out split), landing within 0.02 of the paper. The candidate
property recomputation confirms the transcribed SMILES are correct.

**Observation:** the reproduced RF gives the GTD_9.x candidates only moderate
P(active) ≈ 0.34–0.47. This is consistent with the paper's own DprE1 v2 sub-score
(0.54–0.66, not near 1.0): the reported "overall desirability ≈ 1.0" is driven by
the intestinal-absorption and hepatotoxicity terms, not by predicted potency.

## Part 2 — Docking (see `src/` docking scripts)
## Part 3 — MD on Modal (see `modal_md/`)

## Layout
```
data/      Zenodo downloads + converted structures (gitignored where large)
src/       pipeline scripts
results/   dataset, metrics, model, ROC plot, candidate table
modal_md/  Modal MD pipeline (Part 3)
```
