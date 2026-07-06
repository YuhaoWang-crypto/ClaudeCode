# Script templates

Battle-tested templates from a full reproduction of Chikhale et al. (DprE1).
They are **starting points to adapt**, not a turnkey CLI — edit paths, PDB codes,
HET codes, SMILES, descriptor lists and paper reference numbers for your target.

| Script | Part | Purpose |
|---|---|---|
| `build_activity_model.py` | 1a+1b | Zenodo dataset → structures (CDX/OPSIN) → RF activity model + ROC. Split into two files in a real project. |
| `convert_cdx.py` | 1a | Crash-isolated ChemDraw `.cdx` → SMILES via the OpenBabel CLI (one subprocess per file). |
| `docking_vina.py` | 2 | Receptor prep (meeko + rigid cofactor) and AutoDock Vina docking into WT / mutant / off-target. Other scripts import this. |
| `offtarget_metal_analysis.py` | 2b | Ligand→metal / →cofactor-core geometry to test a mechanistic safety claim. Imports `docking_vina`. |
| `prepare_md_inputs.py` | 3 | Build docked starting complexes (full-ligand reconstruction via meeko) + extract cofactor SDF. Imports `docking_vina`. |
| `modal_md_app.py` | 3 | **Working** OpenMM+AMBER MD on Modal (ff14SB/GAFF2/TIP3P, FAD cofactor, MM-GBSA). Ported from making-it-rain. |

## Dependencies between scripts
`offtarget_metal_analysis.py` and `prepare_md_inputs.py` do
`import_module("docking_vina")` (originally `04_docking`) — keep them in the same
dir, or update the import name. `modal_md_app.py` is self-contained (runs on Modal).

## Cross-references to a full worked example
The complete, runnable project these came from lives at
`dprE1_reproduction/` in this repo (src/, modal_md/, results/, RESULTS.md) — use
it as the reference implementation.

## Read `../references/gotchas.md` before running `modal_md_app.py`.
