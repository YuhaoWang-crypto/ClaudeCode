# mdscreen — MD-based small-molecule screening on Modal GPUs

Run [making-it-rain](https://github.com/pablo-arantes/making-it-rain)-style
molecular dynamics on **your own [Modal](https://modal.com) GPUs**, then use the
trajectories to judge whether a small molecule is **active** (binds and stays
bound), whether it acts as an **allosteric modulator**, and to rank a series of
candidates for **enzyme inhibition/activation**.

It reproduces the making-it-rain physics (OpenMM engine, ff14SB protein force
field, **GAFF2** small-molecule parameters, TIP3P water, Langevin/PME) but
replaces the Colab + AmberTools/tleap plumbing with the modern
**OpenMM + OpenFF** stack and a serverless, parallel Modal backend.

---

## What it does

| Workflow | Command | Purpose |
|---|---|---|
| Protein-only MD | `modal run modal_app.py::protein` | stability / apo dynamics |
| Small-molecule MD | `modal run modal_app.py::ligand` | conformational stability of a candidate |
| Protein–ligand MD | `modal run modal_app.py::complex` | pose stability + **MM/GBSA ΔG** |
| Activity/allostery screen | `modal run modal_app.py::screen` | rank a ligand series, flag allostery |

Each run produces trajectories (`.dcd`), a topology (`.pdb`), state logs, plots
(`.png`), per-metric CSVs, and a JSON summary.

### How activity is judged
- **Binding free energy** — single-trajectory end-state **MM/GBSA** (GBn2
  implicit solvent). More negative ΔG ⇒ stronger predicted binder. Entropy is
  omitted, so treat it as a **ranking** score for a congeneric series, not an
  absolute affinity.
- **Pose stability** — ligand RMSD vs. the docked/initial pose. A favourable ΔG
  with a drifting pose is flagged as *borderline*.
- **Engagement** — protein–ligand contact count, H-bonds, and the persistent
  binding-site residues.
- **Verdict** — combined into `likely active` / `borderline` / `likely inactive`
  with a confidence level (thresholds in `ScreenConfig`).

### How allosteric effects are detected
- An **apo** (ligand-free) reference trajectory is run alongside the **holo**
  (bound) one.
- Per-residue **RMSF** and **dynamic cross-correlation (DCCM)** changes from
  apo→holo are computed. Substantial dynamic changes at residues **distal** to
  the binding site are the fingerprint of an **allosteric** (vs. purely
  orthosteric) modulator.
- A Cα interaction **network** (betweenness centrality) highlights communication
  relay residues.
- Output: `allosteric_verdict` + ranked distal perturbed residues.

> Enzyme **inhibition vs. activation** is inferred from *where* and *how* the
> ligand perturbs dynamics: strong, stable binding in the orthosteric/active
> site ⇒ likely competitive inhibitor; a distal allosteric signal that rigidifies
> or loosens catalytic-loop dynamics ⇒ candidate allosteric modulator
> (activator or inhibitor). The reported metrics give you the evidence; final
> mechanistic calls should be confirmed with assays.

---

## Quick start

```bash
pip install modal pyyaml          # local orchestrator only
modal token new                    # authenticate to YOUR Modal account (uses your GPUs/credits)

# one protein–ligand run (2 ns) with MM/GBSA
modal run modal_app.py::complex --pdb 3poz --ligand "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1" --ns 2

# a full screen
modal run modal_app.py::screen --config examples/screen_config.yaml

# pull all artifacts locally
modal volume get mdscreen-outputs <run_id> ./results
```

Pick the GPU with an env var (default `A10G`):

```bash
MDSCREEN_GPU=A100 modal run modal_app.py::screen --config examples/screen_config.yaml
```

The first run builds the conda image (OpenMM/OpenFF/AmberTools) — a few minutes,
cached thereafter.

### Inputs
- **Receptor**: a local `.pdb` file **or** a 4-letter RCSB id (auto-downloaded,
  cleaned and protonated with PDBFixer).
- **Ligand**: a **SMILES** string **or** an `.sdf`/`.mol2` file. Charges are
  assigned with AM1-BCC (AmberTools); a 3D conformer is generated if absent.

---

## Architecture

```
modal_app.py          Modal app: GPU image + remote functions + entrypoints
mdscreen/
  config.py           SimConfig / ScreenConfig (fully serialisable)
  prepare.py          PDBFixer + OpenFF + SystemGenerator (GAFF2) + solvation
  simulate.py         OpenMM: minimise -> NVT -> NPT -> production
  analyze.py          MDAnalysis: RMSD / RMSF / Rg / contacts / H-bonds
  binding.py          end-state MM/GBSA (GBn2) binding free energy
  allostery.py        DCCM + network + apo/holo perturbation
  pipeline.py         chains the above; screening campaign + classification
  cli.py              local CLI (same pipeline, local GPU/CPU)
examples/
  screen_config.yaml  worked example
```

The `mdscreen` package is **engine code with no Modal dependency** — the exact
same functions run locally (`python -m mdscreen.cli ...`) or on Modal. Modal
only adds the GPU image, parallel fan-out (`complex_md.starmap`), and a
persistent output Volume.

### Scaling a screen
`modal run ...::screen` runs the apo reference once, then fans every ligand out
to its own GPU container in parallel via `starmap`, then runs allostery
comparisons. Ten ligands finish in roughly the wall-clock time of one.

---

## Configuration

Everything is in `mdscreen/config.py`. Common knobs:

| Field | Default | Meaning |
|---|---|---|
| `production_ns` | 5.0 (2.0 for screens) | production trajectory length |
| `timestep_fs` / `hydrogen_mass_amu` | 4.0 / 1.5 | HMR-enabled 4 fs steps |
| `temperature_k` | 300 | thermostat set point |
| `small_molecule_forcefield` | `gaff-2.11` | or `openff-2.1.0` |
| `solvent_padding_nm` | 1.0 | box padding |
| `ionic_strength_molar` | 0.15 | NaCl |
| `platform` | `CUDA` | falls back to CPU if unavailable |

---

## Local run (no Modal)

If you have the MD stack installed on a workstation with a GPU:

```bash
mamba install -c conda-forge openmm openmmforcefields openff-toolkit \
    ambertools pdbfixer rdkit mdanalysis numpy scipy matplotlib networkx pyyaml

python -m mdscreen.cli complex --pdb 3poz --ligand "CCO" --ns 1
python -m mdscreen.cli screen --config examples/screen_config.yaml
```

---

## Notes, scope & honest limitations

- **Not yet executed on a GPU here.** This repository is the full, reviewed
  implementation; it has been syntax-checked but the actual MD runs on your
  Modal GPUs. Expect to tune force-field/residue edge cases for your specific
  proteins (non-standard residues, metals, cofactors, covalent ligands are not
  auto-handled).
- **MM/GBSA is approximate** (single-trajectory, no entropy) — good for ranking
  a congeneric series, not for absolute ΔG. For rigor, use alchemical free
  energy (a natural extension point).
- **Short screening trajectories** (≤ a few ns) sample local relaxation, not
  large conformational change. Lengthen `production_ns`, add replicas, or seed
  from docked poses for production decisions.
- Membrane systems, CHARMM-GUI inputs, GLYCAM carbohydrates, and the AlphaFold2
  front-ends from making-it-rain are out of scope for v0.1 but fit the same
  `prepare -> simulate -> analyse` structure.

## Attribution
Physics and defaults follow the making-it-rain notebooks
(Arantes et al.). This is an independent reimplementation for Modal.
