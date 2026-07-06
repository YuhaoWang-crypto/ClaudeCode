---
name: fairchem-dft-repro
description: >
  Reproduce quantum-chemistry / DFT calculations with open-source models when the
  user has a study, paper, or service report (often Gaussian/VASP) and wants to
  redo the achievable parts cheaply. Computes adsorption energies of molecules on
  surfaces, binding/interaction energies and ΔG/K_a for molecule–cluster
  complexes, geometry optimizations, frontier molecular orbitals (HOMO/LUMO),
  differential charge density (Δρ) and Bader charges. Use when the user asks to
  reproduce/validate a DFT result, compute an adsorption or binding energy, get
  HOMO/LUMO or charge transfer with fairchem UMA / PySCF, or replace VASP/Gaussian
  with open tools. Knows exactly what a machine-learning potential can and cannot
  do and routes the electronic-structure parts to PySCF/QE.
---

# Open-source DFT/QC reproduction (fairchem UMA + PySCF + Bader)

Reproduce the computable parts of a DFT study using **fairchem UMA** (energies &
forces), **PySCF** (orbitals, electron density, redox), and the **Henkelman
`bader`** code — with an honest map of what an MLIP can and cannot do.

## What an MLIP (UMA) can vs cannot do
- ✅ **Can**: geometry optimization, adsorption energy, binding/interaction energy,
  reaction/relative energies, approximate ΔG via harmonic vibrations, fast
  screening at scale.
- ❌ **Cannot** (no wavefunction/electron density): molecular orbitals (HOMO/LUMO),
  differential charge density, Bader/atomic charges, DOS, work function. Route
  these to **PySCF** (finite molecules/clusters) or **Quantum ESPRESSO / GPAW**
  (periodic surfaces).

## Setup (once)
```bash
bash scripts/setup.sh     # torch(CPU) + fairchem-core + ase + pyscf + rdkit + matplotlib
export HF_TOKEN=...        # access to facebook/UMA (gated)
# for Bader: curl -O http://theory.cm.utexas.edu/henkelman/code/bader/download/bader_lnx_64.tar.gz && tar xzf
```

## Building blocks (scripts/)
| Script | Does | Engine/task |
|---|---|---|
| `00_build_structures.py` | build molecules (SMILES→3D), surfaces (slabs), H-saturated clusters | RDKit/ASE |
| `01b_fe110_adsorption.py` | adsorption energy `E(slab+ads)−E(slab)−E(mol)` via a docking search | UMA `oc20` |
| `02_rsv_tiox_binding.py` | binding energy ΔE, K_a for molecule–cluster complexes | UMA `omol` |
| `02b_vib_dG.py` | approximate ΔG via harmonic vibrations (ASE HarmonicThermo) | UMA `omol` |
| `03_pyscf_fmo.py` | HOMO/LUMO, gap, χ/η/ω, dipole; HOMO/LUMO cubes | PySCF B3LYP/6-31G* |
| `04_charge_density_bader.py` | Δρ + Bader charge, validated on a small analog | PySCF PBE + `bader` |
| `05_full_systems_inputs.py` | full-size drivers + Quantum ESPRESSO inputs for periodic systems | PySCF / QE |

## How to run for a NEW system
1. **Edit `00_build_structures.py`** — put your molecule SMILES, surface
   (element/facet/size), and/or cluster cut. Run it to generate `structures/*.xyz`.
2. **Adsorption energy on a surface** → `01b_*`: point it at your slab + adsorbate;
   it relaxes the bare slab to its true minimum, docks the molecule in several
   orientations, and reports the lowest `Eads`.
3. **Binding energy / ΔG / K_a of a complex** → `02_*` then `02b_*` (ΔG).
4. **HOMO/LUMO** → `03_*` on the relaxed molecule.
5. **Δρ + Bader** → `04_*` (finite clusters, on this machine) or `05_*` to emit QE
   inputs for periodic/metal systems to run on GPU/HPC.

## Critical correctness rules (learned the hard way)
- **Reference-state discipline.** Every reference (bare slab, isolated molecule,
  fragment) must be relaxed to the SAME standard as the complex. An under-relaxed
  slab once produced a spurious −10.8 eV adsorption energy; relaxing it properly
  gave −1.8 eV. Compare only energies from the **same UMA task** (shared reference).
- **Functional offsets.** UMA task heads ≈ different DFT (oc20≈RPBE, omol≈ωB97M-V,
  omat/omc/odac per name). Absolute numbers differ from a paper's PBE/B3LYP; trust
  **signs, magnitudes, and relative trends**.
- **Fragments for Δρ** must be frozen at their in-complex geometry (else ∫Δρ ≠ 0).
- **TiOₓ / small-gap clusters**: add Fermi smearing + `density_fit()` in PySCF or
  SCF won't converge.
- **Task choice**: `oc20` = molecule on catalyst/metal surface; `omol` = finite
  molecules/clusters (needs `charge`,`spin`); `omat`/`omc`/`odac` = bulk materials
  /molecular crystals/MOF-DAC.

## Outputs
`results/*.json` (all numbers), `structures/relaxed_*.xyz`, `results/*.cube`
(orbitals, Δρ), `figures/*.png`. See the bundled `REPORT_TEMPLATE.md` for how to
write up the comparison to the original study.
