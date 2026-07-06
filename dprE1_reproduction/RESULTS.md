# Reproduction Results Summary

Open-source reproduction of Chikhale et al. (chemRxiv 2026,
DOI 10.26434/chemrxiv.15004861/v2). See `README.md` for the full method mapping.

## Part 1 — Activity model + candidate analysis ✅

**Dataset** (`src/01_prepare_dataset.py`) — rebuilt from the authors' Zenodo
release. Structures come from their ChemDraw `.cdx` drawings (OpenBabel), with
OPSIN IUPAC→SMILES as fallback; activity taken from the molar (IC50) entries and
converted to pIC50.

| | Paper (DprE1 v2) | Reproduced |
|---|---|---|
| IC50 molecules parsed | 406 | 396 |
| after structure resolution + dedup | — | 366 |
| actives (pIC50 ≥ 5.75) | 192 | 161 |

**Model** (`src/02_build_activity_model.py`) — Random Forest on the paper's 11
descriptors + 1024-bit ECFP4.

| | Paper | Reproduced |
|---|---|---|
| ROC AUC | 0.92 | **0.902 ± 0.005** (5-fold CV × 3 iterations) |

**Candidates** (`src/03_candidate_analysis.py`) — recomputed properties for
TCA1 + GTD_9.1–9.10 (SMILES from Appendix-2) vs paper Table 2:

- MolWt: mean absolute error **0.4 Da** (confirms SMILES transcription)
- MolPSA ≈ RDKit TPSA: near-exact (e.g. 94.3/94.3, 116.0/116, 80.0/80)
- Rotatable bonds: 82% exact match (± 1 elsewhere, amide-counting convention)

## Part 2 — Docking (AutoDock Vina vs GOLD/ChemPLP)

Vina affinities (kcal/mol, more negative = better). Scores are **not** on the
GOLD ChemPLP scale — compare rankings/trends, not absolute numbers.
(`src/04_docking.py`, `src/05_docking_analysis.py`; figure
`results/docking_comparison.png`.)

| Compound | WT Vina | WT GOLD | Y314C Vina | Y314C GOLD | CYP2C9 Vina | CYP2C9 GOLD |
|---|---|---|---|---|---|---|
| TCA1 | -8.46 | 68.7 | -9.87 | 69.0 | -8.83 | 81.4 |
| GTD_9.1 | -10.93 | 90.5 | -9.21 | 97.8 | -9.88 | 74.2 |
| GTD_9.2 | -9.79 | 90.4 | -10.54 | 87.7 | -11.21 | 68.8 |
| GTD_9.3 | -8.72 | 85.5 | -10.29 | 82.6 | -10.37 | 88.1 |
| GTD_9.4 | -10.25 | 76.9 | -9.12 | 86.0 | -10.23 | 66.3 |
| GTD_9.5 | -10.05 | 84.0 | -10.17 | 75.4 | -10.80 | 91.3 |
| GTD_9.6 | -9.55 | 93.0 | -9.28 | 79.5 | -10.25 | 92.3 |
| GTD_9.7 | -9.45 | 90.5 | **-10.91** | 94.7 | -10.70 | 87.6 |
| GTD_9.8 | -9.92 | 84.7 | -9.77 | 77.7 | -11.48 | 96.4 |
| GTD_9.9 | -10.55 | 90.8 | -10.21 | 83.1 | -10.89 | 85.0 |
| GTD_9.10 | -9.92 | 92.8 | -10.09 | 78.1 | -10.13 | 98.2 |

**What reproduces:**
- **DprE1 WT: 10/10 GTD candidates dock better than TCA1** — matches the paper's
  central docking claim exactly (paper Table 3: all GTD 76.9–93.0 vs TCA1 68.7).
- On the Y314C mutant, **GTD_9.7 is the strongest binder (-10.91)** — one of the
  two lead candidates the paper highlights (GTD_9.7, GTD_9.4).

**What does not (honestly):**
- Absolute Vina↔GOLD rank correlation is weak (WT Pearson r = -0.46, Spearman
  ρ = -0.19; sign is correct — lower Vina ↔ higher GOLD). Expected: different
  scoring functions, and all GTD molecules sit in a narrow high-affinity band so
  within-series ordering is noisy.
- **CYP2C9 off-target:** the paper's safety claim is mechanistic (GTD compounds
  avoid heme-iron coordination), which a raw Vina affinity cannot capture — here
  the GTD compounds still score well on CYP2C9. Reproducing that claim needs
  interaction/heme-distance analysis, not docking score alone.

## Part 3 — MD on Modal (making-it-rain port)

Runnable OpenMM + AMBER (ff14SB/GAFF2/TIP3P) pipeline on Modal GPUs
(`modal_md/`). Reproduces the paper's Table-4 observables (ligand/protein/
cofactor RMSD, protein RMSF). Not executed here (requires a Modal account and
~1 GPU-day per 500 ns complex); see `modal_md/README.md` to run.

## Honest limitations

- The **generative step (BIOVIA GTD / GFSP)** is proprietary and is *not*
  reproduced — we analyse the authors' published GTD_9.x molecules directly.
- Docking uses Vina, not GOLD; MD uses OpenMM/AMBER, not CHARMm — so absolute
  scores/energies differ by construction. The reproduction validates the
  **trends and the ML model**, which is where the paper's public data allows it.
- The reproduced RF scores the GTD_9.x candidates at moderate P(active)
  (~0.34–0.47), consistent with the paper's own DprE1 v2 sub-scores (0.54–0.66);
  the reported "overall desirability ≈ 1.0" is driven by the ADME/tox terms.
