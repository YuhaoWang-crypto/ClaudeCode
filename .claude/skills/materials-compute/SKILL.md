---
name: materials-compute
description: >-
  Scope, prototype, and scaffold theoretical / computational materials-science
  and molecular-modeling jobs (DFT, AIMD, ML potentials, docking, phonons,
  spectroscopy) for client-style requests. Use when someone asks whether a
  materials/chemistry calculation can be done, wants a quick runnable demo or
  proof-of-concept, needs the workflow and HPC input files for a production run,
  or asks to screen metal centers, predict metal-binding sites, build molten-salt
  / ML-potential training data, model adsorption / SERS / vibrational spectra,
  reproduce a quantum-dot / exciton model, or plan phonon / electron-phonon /
  excited-state (STE) calculations. Splits every request into what runs locally
  (analytical models, small-molecule DFT via PySCF, structure builders,
  geometric predictors) vs. what needs HPC (periodic VASP/QE/CP2K, EPW, GW-BSE),
  and labels demo-proxy vs. production honestly.
---

# Materials & molecular computation — demo + scaffolding toolkit

A reusable capability for turning a client's "can you compute this?" into
(a) a **runnable local demo/proxy** that produces real numbers and a figure in a
lightweight Python environment, and (b) a **production scaffold** (workflow +
HPC input files) for the full job. Everything lives in `materials_demos/`.

## The core discipline: local-proxy vs. HPC-production, labeled

Every materials request splits into two tiers. State which tier each piece is in
and never let a proxy masquerade as the real material.

| Tier | What it is | Where it runs |
|---|---|---|
| ✅ **local demo** | analytical/tight-binding models, small-molecule DFT (PySCF), classical MD (ASE), structure builders, geometric predictors, data packaging | this environment (`pip install` only) |
| ⚠️ **HPC production** | periodic DFT (VASP/QE/CP2K), AIMD at scale, EPW e-ph, GW-BSE, phonons, NEB barriers, ML-potential training | client's cluster / supercomputer |

The local tier proves the workflow, validates the physics on a proxy, and
de-risks the quote. The HPC tier is where real material numbers come from — scaffold
it with `hpc_templates/`, don't pretend to run it here.

## What runs locally — the six worked demos

| Demo | File | Answers | Tech |
|---|---|---|---|
| #1 Quantum-dot FSS | `qd_fss_model.py` | exciton fine-structure splitting vs. uniaxial stress; polarization angle; which dots tune to zero FSS | numpy (analytical, reproduces PRL 106,227401) |
| #2 Metal-center screen | `metal_center_screen.py` | electronic-structure pre-screen of catalytic metals ([M–OH] model): gap, charge, spin | PySCF UKS DFT |
| #3 Method checklist | `demo3_method_checklist.md` | maps reviewer demands (phonon, e-ph, excited-state/STE) to method+software | reference doc |
| #4 Molten-salt MLP | `moltensalt_mlp_pipeline.py` + `moltensalt_builder.py` | liquid MD → g(r)/diffusion; **DeePMD/npy** training-set packaging; real-salt initial structures | ASE, dpdata |
| #5 Cd²⁺ site predictor | `cd_binding_site_predictor.py` | ranks metal-binding sites in a protein + mutation targets (HSAB + clique clustering); validated on 1CA2 | Biopython/numpy |
| #6 SERS chirality | `sers_chirality_dft.py` | proves free-molecule L≡D vibrations ⇒ SERS difference lives on the surface | PySCF DFT + RDKit |

Run any demo directly (`python3 materials_demos/<file>.py`); each prints results
and saves a `*.png`. Install set:
`pip install numpy scipy matplotlib pyscf ase dpdata biopython rdkit pyberny`.

## Scaffolding a production job — `materials_demos/hpc_templates/`

Ready-to-adapt input files keyed to each direction (replace `<<...>>`):
`moltensalt_cp2k/` & `moltensalt_vasp/` (AIMD), `deepmd/` (train.json + workflow),
`cd_docking/` (metal-aware docking + QM/MM refinement), `sers_au_slab/` (Au(111)
adsorption + Raman), `phonon_qe/` (DFPT phonons + Huang-Rhys / STE). See
`hpc_templates/README.md` for the folder→client map and the molten-salt
end-to-end (builder → AIMD → dpdata → dp train → LAMMPS properties).

## How to handle a new client request

1. **Classify** the ask into the tiers above; name the real deliverable
   (e.g. "density & viscosity vs. T", "which residues to mutate", "adsorption
   energy + Raman").
2. **Pick or build a local proxy** that produces a real number + figure. Reuse a
   demo if one fits; otherwise write the smallest honest model (analytical,
   small-molecule DFT, classical MD, or a geometric predictor).
3. **Validate** the proxy against something known (demo #5 recovers the 1CA2
   metal site; demo #1 matches the paper's analytic limits; demo #6 hits machine
   precision on L≡D). A demo that can't be sanity-checked is a liability.
4. **Scaffold the production run** from `hpc_templates/`, marking which settings
   are converged vs. placeholder.
5. **Label everything** demo-proxy vs. production, and say plainly what needs the
   client's HPC.

## Honesty rules (non-negotiable)

- A molecular cluster is **not** the periodic material; a classical/LJ potential
  is **not** DFT; a geometric predictor is **not** docking. Say so each time.
- Template defaults (cutoffs, k-meshes, timesteps, functionals) are **reasonable
  starting points, not converged production settings** — flag them.
- Report what validated and what didn't. If a proxy only demonstrates plumbing
  (e.g. DeePMD packaging with LJ forces), say the physics still needs AIMD.
- This environment has **no** DFT engine beyond PySCF (molecular) and **no** HPC.
  Don't claim a VASP/QE/CP2K/EPW/GW-BSE result was computed here.

## Reference

- `reference/request-playbook.md` — worked mapping from the real client chats
  (quantum dots, Cr-catalyst, scintillator/STE, molten salt, Cd-protein, SERS)
  to the demo + HPC scaffold chosen, and why.
- `materials_demos/README.md` — user-facing index of all demos + capability table.
