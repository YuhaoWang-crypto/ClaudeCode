---
name: 2d-desalination-dft-mlp
description: >-
  End-to-end multiscale workflow for studying 2D materials (graphene, MoS2, …)
  as reverse-osmosis desalination membranes: DFT → machine-learning potential →
  non-equilibrium MD. Use when the task is to build a nanoporous 2D membrane,
  compute water/ion pore-crossing barriers and selectivity, generate QE
  (Quantum ESPRESSO) SCF/CI-NEB inputs, train a DeePMD/MACE/NequIP potential
  (with DP-GEN active learning), or run pressure-driven NEMD to get water flux
  and salt rejection. Also covers the MLP-vs-ReaxFF decision and running the
  heavy steps on Modal cloud (CPU for DFT, GPU for MLP/NEMD). Every physical
  number is honesty-labelled: surrogate/order-of-magnitude vs. converged DFT.
---

# 2D-material desalination: DFT → MLP → NEMD

A reusable, runnable pipeline (bundled under `pipeline/`) that takes a 2D
membrane from geometry to water-flux / salt-rejection numbers along the standard
multiscale route. It runs **fully open-source** (ASE, Quantum ESPRESSO / CP2K,
DeePMD-kit / MACE / NequIP, LAMMPS, DP-GEN) and everything that does not need a
DFT/MD engine is validated locally with `numpy + ase`.

```
 Stage 1 (DFT/AIMD) ──▶ Stage 2 (ML potential) ──▶ Stage 3 (NEMD)
 QE/CP2K gold data      DeePMD/MACE/NequIP + DP-GEN   LAMMPS: flux & rejection
```

## When to use this skill
- Build a nanoporous **graphene / MoS2** (or other 2D) membrane + salt-water cell.
- Compute **pore-crossing barriers** and **selectivity** (water vs ions) via DFT.
- Generate ready-to-run **QE SCF / CI-NEB** inputs (incl. **hydrated ions**).
- Train an **MLP** (DeePMD / MACE / NequIP) with **DP-GEN** active learning.
- Run **pressure-driven NEMD** → **water flux + salt rejection**.
- Decide **MLP vs ReaxFF**, or run heavy steps on **Modal** cloud.

## Non-negotiable honesty rule
Label **every** physical number as one of:
- ✅ **converged DFT** (the user's QE/CP2K with a vdW functional, charged
  supercells for ions, explicit water, proper convergence) — publication-grade.
- ⚠️ **surrogate / order-of-magnitude** (universal MLIP like CHGNet, coarse
  scans, bare-ion approximations, unconverged cutoffs) — trends only.

Never present a surrogate barrier or an unconverged energy as a real prediction.
The bundled study (`pipeline/06_dft_surrogate_study/REPORT.md`) is a worked
example of this labelling.

## The critical physics to get right
**Salt rejection comes from the ion DEHYDRATION penalty, not bare-ion sterics.**
A bare Na⁺/Cl⁻ is *smaller* than water and passes pores that stop water — the
opposite of selectivity. Always model **explicitly-hydrated ions** (ion + first
solvation shell) when computing rejection. The bundled `run_hydrated_pmf.py`
demonstrates this: bare Na⁺ shows ~0 eV barrier while Na⁺·(H₂O)₆ shows a real
one. `02_dft/build_hydrated_neb.py` builds the corresponding real-DFT NEB.

## MLP vs ReaxFF (default recommendation)
For *physical* (non-reactive) RO — water and hydrated ions crossing without
breaking bonds — **use an MLP** (DeePMD/MACE/NequIP). Reserve **ReaxFF** only for
*chemistry*: pore-edge protonation/deprotonation, membrane degradation, fouling.
Do not develop both in parallel; ship the MLP path first.

## Workflow (map task → script)

All paths are under `pipeline/`. Run `pip install -r pipeline/requirements.txt`
first (numpy, ase, matplotlib for local steps).

| Goal | Command |
|---|---|
| Build graphene/MoS2 pore + salt-water cell | `01_build/build_graphene_pore.py`, `build_mos2_pore.py`, `build_saltwater_box.py` |
| Ready QE SCF input from a geometry | `02_dft/make_qe_input.py --xyz <geom.xyz> --out scf.in` |
| Cutoff/k-point convergence grid | `02_dft/make_convergence_tests.py --xyz <geom.xyz>` |
| **Bare** ion / water NEB endpoints + input | `02_dft/build_neb_endpoints.py --ion Na` |
| **Hydrated** ion NEB (correct for rejection) | `02_dft/build_hydrated_neb.py --ion Na --n-water 6` |
| Parse a NEB barrier | `02_dft/parse_neb_barrier.py --dat <prefix>.dat` |
| DFT→DeePMD data / npy→extxyz for MACE-NequIP | `03_mlp/prepare_deepmd_data.py` |
| Generate + label training frames (local demo) | `03_mlp/generate_training_frames.py` |
| Train MLP | configs `deepmd_input.json` / `mace_config.yaml` / `nequip_config.yaml`; DP-GEN `dpgen_param.json` |
| Calibrate NEMD piston force for a pressure | `04_nemd/piston_force.py --data <cell.lammps-data> --pressure-mpa 100 --n-piston N` |
| Run NEMD | `04_nemd/in.desalination.lammps` (LAMMPS) |
| Analyze flux + rejection from a dump | `04_nemd/analyze_flux_rejection.py` |
| Validate the analysis pipeline (no engine) | `05_toy_validation/toy_nemd.py --plot out.png` |
| Local surrogate study (both materials) | `06_dft_surrogate_study/run_study.py`, `run_hydrated_pmf.py` |
| One-command local check of everything | `pipeline/run_local_pipeline.sh` |

## Running the heavy steps (Modal cloud / HPC)
- **DFT / NEB are CPU-bound** → `02_dft/modal_run_dft.py`, `modal_run_neb.py`
  (or `slurm/qe_scf.slurm`). Do **not** use a GPU for plane-wave DFT.
- **MLP training + NEMD are GPU-bound** → `03_mlp/modal_train_deepmd.py`,
  `modal_train_mace.py` (or `slurm/deepmd_train.slurm`, `slurm/lammps_nemd.slurm`).
Modal images (QE, DeePMD, MACE) are defined in code; `pip install modal &&
modal token new` then `modal run …`.

## Recommended milestone order
1. Real QE + vdW-DF2 SCF on the r≈2.5 Å pore (both materials); converge cutoffs.
2. **Hydrated-ion** CI-NEB / PMF barrier — the first real selectivity number.
3. A few ps AIMD of confined salt water → MLP training frames.
4. Train + bake-off DeePMD vs MACE vs NequIP; validate against the NEB barrier.
5. 100 MPa NEMD → first water-flux / salt-rejection point; scan P and pore size.
6. Repeat for the next 2D material.

## Deeper reference
See `reference/methodology.md` for the physics, the surrogate-vs-DFT boundary,
units, and common pitfalls. See `pipeline/06_dft_surrogate_study/REPORT.md` for a
complete worked study with honesty labels and figures.
