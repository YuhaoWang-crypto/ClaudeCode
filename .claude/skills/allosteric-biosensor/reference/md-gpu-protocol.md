# GPU apo/holo MD + ligand parameterization protocol

The rigorous open-source route to the switch mechanism: run production MD of the
chimera **apo** and **holo** (ligand bound), then compare **conformational
entropy** ΔS(holo−apo) and **active-site RMSF**. The paper's mechanism predicts
ligand binding **lowers** conformational entropy (ΔS < 0) and pre-organizes the
reporter active site — the structural basis of the ON state / dynamic range.

Implemented in `biosensor_pipeline/md_entropy.py` (production functions) and
`biosensor_pipeline/modal_app.py` (cloud-GPU runner).

## What runs

1. **Ligand parameterization** — `make_system_generator()` uses
   `openmmforcefields.SystemGenerator` with an **OpenFF SMIRNOFF** small-molecule
   force field (`openff-2.2.0`, pure-python, no ambertools) matched to the
   analyte SMILES. GAFF-2 is a drop-in alternative if AmberTools is present.
2. **System build** — `build_complex()` loads the protein (Boltz model, OXT-
   capped), adds the parameterized ligand (from Boltz's `sample_0_ligands.sdf`),
   solvates (TIP3P + 0.15 M ions, PME) or uses GBn2 implicit for cheap scans.
3. **Production MD** — `run_leg()`: minimize → 200 ps equilibration (Langevin +
   Monte-Carlo barostat, 4 fs with HMR) → `prod_ns` production, sampling Cα
   frames.
4. **Analysis** — Kabsch superposition (`_superpose`) → per-residue + active-site
   RMSF; **Schlitter absolute entropy** (`schlitter_entropy`, mass-weighted
   covariance, J/mol/K). `compare_apo_holo()` returns ΔS and Δ(active-site RMSF).

## Environment (conda-forge; not pip-installable)

```bash
mamba create -n biomd -c conda-forge python=3.11 openmm=8.1 openmmforcefields \
      openff-toolkit mdtraj numpy biotite
mamba activate biomd
```

## Run on a cloud GPU with Modal (recommended)

`modal_app.py` defines a micromamba GPU image with the stack above and runs the
two legs in parallel (`gpu="A10G"`).

```bash
pip install modal && modal token new

# 1) prepare inputs from Boltz models (local; strips the ligand for the protein PDB)
python3 -c "from biosensor_pipeline.modal_app import cif_to_pdb as c; \
  c('biosensor_out/vitd_apo.cif','biosensor_out/vitd_apo.pdb'); \
  c('biosensor_out/vitd_holo.cif','biosensor_out/vitd_holo_protein.pdb')"
#    the holo ligand SDF is Boltz's sample_0_ligands.sdf (download from the prediction)

# 2) run apo + holo on GPU
modal run biosensor_pipeline/modal_app.py \
    --apo-pdb  biosensor_out/vitd_apo.pdb \
    --holo-pdb biosensor_out/vitd_holo_protein.pdb \
    --ligand-sdf biosensor_out/vitd_ligand.sdf \
    --smiles "CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C" \
    --prod-ns 50
```

Output: `S_apo`, `S_holo`, `ΔS(holo−apo)`, and the active-site RMSF change.

## Sampling & rigor

- **Convergence:** a single 20–50 ns run gives a *trend*, not a converged ΔS.
  For a defensible number use ≥3 replicas × 100–500 ns, block-average, and report
  the spread. Schlitter is an **upper bound**; the quasi-harmonic (Andricioaei–
  Karplus) estimator is a tighter alternative if you swap it in.
- **Ligand match:** the SDF ligand must match the SMILES graph (openff matches by
  graph); Boltz's ligand SDF carries bonds, so it parameterizes cleanly.
- **Honesty labels:**
  - ✅ ΔS, RMSF are *measurements on the trajectories*.
  - ⚠️ ΔS < 0 *supports* the entropic-switch mechanism but is **not** a dynamic
    range. Mapping ΔS (or Δactive-site-RMSF) to a DR number needs calibration
    against measured `kobs(+L)/kobs(−L)` on a training set of constructs.
  - The **only** ground-truth DR is the wet-lab titration.

## Where this sits in the ladder

Boltz (static) → `coupling.py` (apo/holo active-site pLDDT, free) →
**this MD (physical ΔS/RMSF)** → QM/MM (kcat barrier, once ON/OFF ensembles are
credible) → bench kobs. Each rung is more physical and more expensive; MD is the
first rung that actually sees the entropic driver of the switch.
