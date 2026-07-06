---
name: reproduce-genai-drug-design
description: Reproduce a computational drug-discovery paper (generative-AI + structure-based hit-to-lead, e.g. BIOVIA GTD / GOLD / Discovery Studio / CHARMm workflows) with an open-source stack. Use when the user hands you a chemistry/drug-design preprint or PDF and asks to reproduce it, or to rebuild the ML activity model, docking, off-target/mechanistic analysis, or MD of such a study. Substitutes RDKit+scikit-learn (activity model), AutoDock Vina (docking), and OpenMM+AMBER on Modal (MD, ported from making-it-rain) for the proprietary tools, and checks reproduced numbers against the paper.
---

# Reproduce a generative-AI + structure-based drug-design paper (open source)

Papers like Chikhale et al. *"Generative AI and Structure-Based Workflow for the
De Novo Design and Optimization of DprE1 Inhibitor Candidates"* run entirely on
**proprietary software** — BIOVIA Generative Therapeutic Design (GTD), GOLD,
Discovery Studio, CHARMm. You cannot re-run those, and the generative algorithm
is closed. But the **scientific core is reproducible** with open tools, and these
papers usually publish everything you need (dataset on Zenodo, candidate SMILES
in an appendix, PDB codes). This skill rebuilds and *validates* that core.

## First: set expectations honestly

Tell the user up front what can and cannot be reproduced:

| Paper component | Proprietary tool | Open-source substitute | Reproduces |
|---|---|---|---|
| Activity model | GTD Random Forest | RDKit descriptors + ECFP4 + scikit-learn | ✅ the metric (ROC) |
| De novo generation | BIOVIA GTD (GFSP) | — (algorithm is closed) | ❌ analyse published molecules instead |
| Docking | GOLD (ChemPLP) | AutoDock Vina / smina | ✅ trends/ranking, not absolute scores |
| Off-target safety | Discovery Studio | direct geometry analysis of docked poses | ✅ the mechanism |
| MD + MM-GBSA | CHARMm / Discovery Studio | OpenMM + AMBER on Modal (making-it-rain) | ✅ Table-4 observables |

Different scoring functions / force fields mean **absolute** scores and energies
differ by construction — you reproduce **trends, the ML model, and mechanisms**.

## Workflow (each part is independent; scope to what the user asks)

Read the PDF first (`references/read_pdf.md`), extract: the target + PDB codes,
the hit compound, the TPP/objectives, the model type + descriptors, the candidate
SMILES (usually an appendix), and the data-availability link (Zenodo/figshare).

### Part 1 — ML activity model + candidate analysis (always start here; fully runnable anywhere)
1. Download the dataset (Zenodo API → find the real file `links.self`).
2. Resolve structures: prefer the authors' ChemDraw **`.cdx`** via OpenBabel CLI
   (run **one file per subprocess** — the C++ CDX reader can SIGABRT and kill the
   whole run); fall back to **OPSIN** on IUPAC names. Match `.cdx` to activity
   rows by `(DOI-folder, label)` — filenames are inconsistent.
3. Parse activity: keep molar (µM/nM) = IC50 entries → pIC50; ug/mL are usually
   MIC (separate model). Label active at the paper's pIC50 threshold.
4. Train RF on the paper's descriptors + 1024-bit ECFP4; report ROC AUC over 3
   iterations of stratified CV (the paper averages 3). Target: paper's ROC.
5. Recompute candidate properties (MolWt, TPSA≈MolPSA, rot-bonds) from the
   appendix SMILES to *validate the transcription* (MolWt should match ±1 Da).
   → `scripts/build_activity_model.py`

### Part 2 — Docking (AutoDock Vina)
1. Fetch PDBs (`files.rcsb.org/download/XXXX.pdb`). Identify native ligand +
   cofactor HET codes.
2. Clean receptor: chain A, strip altlocs (keep A, blank column), keep cofactor
   (FAD/HEM) in the **rigid** receptor so π-stacking is represented; extract the
   native ligand to define the box centre.
3. Receptor PDBQT via meeko `mk_prepare_receptor.py` (`--allow_bad_res
   --default_altloc A --box_enveloping <native> --padding 6`); append the
   cofactor PDBQT (obabel `-xr`).
4. Ligand: RDKit embed + `mk_prepare_ligand.py`. Dock with the `vina` Python API.
   Validate by re-docking the native ligand (RMSD vs crystal ≈ paper's).
5. Compare rankings to the paper's Table-3 (Spearman/Pearson; sign flips because
   Vina lower-is-better vs GOLD higher-is-better). The robust signal is usually
   categorical ("N/N candidates beat the reference").
   → `scripts/docking_vina.py`

### Part 2b — Off-target / mechanistic analysis (reproduces claims scores can't)
If the paper makes a *mechanistic* safety/selectivity claim (e.g. "does not
coordinate the heme iron"), measure the geometry directly on the docked poses:
min ligand→metal distance, min ligand→cofactor-core distance, coordination
(<2.5 Å) yes/no. Anchor on the native co-crystal ligand (its crystal pose gives
a ground-truth distance to check against the paper). → `scripts/docking_vina.py`
(reuse the receptor prep) and `references/gotchas.md`.

### Part 3 — MD + MM-GBSA on Modal (making-it-rain port)
`scripts/modal_md_app.py` is a **working** OpenMM+AMBER pipeline on Modal GPUs,
ported from making-it-rain's Protein_ligand notebook. It parametrises the ligand
**and the FAD cofactor** (GAFF2), solvates (tleap/TIP3P), runs NPT MD, and
computes RMSD/RMSF (Table 4) + MM-GBSA (MMPBSA.py, igb=5).
1. `scripts/prepare_md_inputs.py` builds docked starting complexes (reconstruct
   the full ligand from the Vina pose with **meeko `RDKitMolCreate`** — obabel
   drops nonpolar H and antechamber then fails).
2. `modal run scripts/modal_md_app.py --system <name> --ns 1` for a fast test;
   `--all --ns 500` for the full protocol (~1 GPU-day/complex).
3. Needs `MODAL_TOKEN_ID`/`MODAL_TOKEN_SECRET` and `pip install
   'modal[api-proxy-support]'` (the proxy needs python-socks).

**The tleap/parametrisation pipeline has many sharp edges — read
`references/gotchas.md` before touching `modal_md_app.py`.** Every one of these
cost a debugging round: openff can't write MOL2 (feed SDF to antechamber);
PDBFixer H-names clash with ff14SB (let tleap protonate via pdb4amber);
`addIons`/`addIonsRand` reject the two-ion `0 0` form (use two single-ion calls);
MMPBSA namelists need name+`/` on their own lines; ante-MMPBSA skips the complex
topology when there's no solvent (use `complex_dry.prmtop` as `-cp`).

## Environment setup

```bash
pip install rdkit scikit-learn pandas numpy openpyxl joblib matplotlib scipy \
            py2opsin openbabel-wheel vina meeko gemmi 'modal[api-proxy-support]'
```
`py2opsin` needs Java. If a broken `cryptography`/`_cffi_backend` blocks pypdf,
`pip install --force-reinstall cffi`. PDB rendering (`pdftoppm`) is often missing;
extract PDF text with pypdf instead (`references/read_pdf.md`).

## Deliverables

Reproduce into a `<paper>_reproduction/` project: `src/` scripts, `results/`
(dataset JSON, model metrics, ROC + comparison figures, docking/heme tables,
MD metrics), `modal_md/`, and a `RESULTS.md` that tabulates **reproduced vs paper**
for every claim, ending with an explicit **honest-limitations** section.
