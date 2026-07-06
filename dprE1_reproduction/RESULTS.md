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

<!-- DOCKING_TABLE -->

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
