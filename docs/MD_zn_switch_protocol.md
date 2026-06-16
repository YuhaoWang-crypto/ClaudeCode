# Distal His-pair / Zn-switch → substrate-binding ΔΔG: MD protocol

Goal: test whether a **distal engineered Zn site** (His-pair grafted into a
large-motion region) perturbs the **active-site ↔ substrate binding** of the
enzyme — i.e. a computational proxy for "does the switch disturb activity".

We do this as a **relative** MM-GBSA calculation:

```
ΔΔG = ΔG_bind(substrate | +Zn)  −  ΔG_bind(substrate | no Zn)
```

A large |ΔΔG| (or a clear shift in the substrate Kd) ⇒ the distal metal site is
allosterically coupled to the active site = predicted activity perturbation.
Sign tells you the *likely* direction (tighter vs weaker substrate binding), but
— as the Soumillion 2010 cpTEM-1 paper shows — activate-vs-inhibit is not
reliably predictable and must ultimately be measured (kcat/Km ± metal).

## Models built here (`data/models/`, demo = GCK)

| file | contents | role in MM-GBSA |
|------|----------|-----------------|
| `GCK_WT_glc.pdb`        | WT GCK (1V4S:A) + glucose            | control |
| `GCK_mutHis_glc.pdb`    | G193H/F195H mutant + glucose         | **state 1: no Zn** |
| `GCK_mutHis_Zn_glc.pdb` | G193H/F195H mutant + glucose + Zn²⁺  | **state 2: +Zn** |

- His-pair `G193H / F195H` = top hit from `results/GCK.md` (distal region
  186–199, ~24 Å from the active centre, Cα–Cα 7.06 Å).
- Glucose (GLC) is the **substrate**, already crystallised in the active site of
  1V4S → no docking needed, the cleanest of the targets.

> ⚠️ As-built, the two His are in pdbfixer default rotamers (ND1–ND1 ≈ 7.6 Å),
> too far to coordinate Zn directly. The Zn is placed at the midpoint as a
> *starting* position. You MUST refine the metal site with distance restraints
> (below) during minimisation/equilibration so the His rotate in to coordinate.

## Running on making-it-rain `Protein_ligand.ipynb`

That notebook's standard path = **1 protein (ff14SB) + 1 small-molecule ligand
(GAFF2/AM1-BCC via antechamber) + TIP3P, OpenMM MD, MM-PBSA/GBSA**. Two gaps for
our case — neither is handled out of the box:

### Gap 1 — the ligand is the *substrate* (glucose), not a drug
Feed `glucose` as the notebook's ligand (extract `GLC` to its own PDB/MOL2).
Glucose is neutral and small → antechamber GAFF2 + AM1-BCC parametrises it fine.
The protein input is the **mutant** PDB with glucose removed.

### Gap 2 — Zn²⁺ (the notebook does NOT parametrise metals)
Pick one of:

1. **Nonbonded (12-6) ion + restraints (quickest).** Add Zn²⁺ as an ion in
   `tleap` (`addIons`/`loadAmberParams frcmod.ions234lm_126_tip3p`), then during
   minimisation + early equilibration impose harmonic distance restraints
   Zn–N(His) ≈ 0.21 nm so the His coordinate. Simple, but nonbonded Zn drifts;
   keep weak restraints through production or accept it as approximate.
2. **Bonded model — ZAFF (recommended for a real number).** Use AmberTools
   `MCPB.py` / ZAFF to derive bonded Zn–His parameters, giving a stable
   tetrahedral site. More setup, far more physical. Do this once the geometry
   from (1) looks reasonable.
3. **Cationic dummy-atom model.** Good middle ground for tetrahedral Zn.

In the notebook's `tleap` cell add, e.g.:
```
loadAmberParams frcmod.ions234lm_126_tip3p   # 12-6 divalent ion params
# (for ZAFF: loadAmberPrep / loadAmberParams of the MCPB-generated files)
ZN = loadpdb zn.pdb
complex = combine { protein glucose ZN }
addIonsRand complex Na+ 0 Cl- 0
solvatebox complex TIP3PBOX 12.0
```

### His coordination restraint (the "约束 His 配位" step)
In OpenMM, before production, add:
```python
from openmm import CustomBondForce
r = CustomBondForce("0.5*k*(d-d0)^2")
r.addPerBondParameter("k"); r.addPerBondParameter("d0")
# d0=0.21 nm, k=200000 kJ/mol/nm^2 for each Zn–N(His193/His195) pair
for n_idx in (his193_NE2, his195_NE2):
    r.addBond(zn_idx, n_idx, [200000.0, 0.21])
system.addForce(r)
```
Minimise + 0.5–1 ns restrained NPT equilibration, then relax/keep restraints for
production (50–100 ns recommended; ≥3 replicas).

## MM-GBSA / ΔΔG analysis
Use AmberTools `MMPBSA.py` (single-trajectory) on each state with glucose as the
ligand:
```
&general startframe=1, endframe=5000, interval=10 /
&gb igb=8, saltcon=0.15 /
```
Then `ΔΔG = ΔG_bind(+Zn) − ΔG_bind(no Zn)`. Also report, per state:
- active-site pocket RMSD / volume vs WT (geometric distortion),
- catalytic-residue distances to glucose,
- domain/hinge angle distribution (the metal's actual mechanical effect).

## Per-target "intermediate" (the substrate/ligand differs)
| target | conformer w/ ligand | what to put in the active site for MM-GBSA |
|--------|--------------------|--------------------------------------------|
| **GCK**   | 1V4S (glucose bound) | glucose (present) ± ATP/AMP-PNP |
| **PTP1B** | build from 1T49      | **phospho-Cys215 covalent intermediate** (model the thiophosphate; or pTyr-peptide Michaelis complex) — needs covalent/QM-MM setup |
| **AdK**   | 1AKE (Ap5A bound)    | **Ap5A** transition-state mimic (present; large, charged → careful GAFF charges) |
| **β-lac (cpTEM1)** | 1BTL + model | **acyl-enzyme intermediate** on Ser70 (covalent) — needs covalent parametrisation |

PTP1B and β-lactamase involve **covalent** intermediates → not a plain GAFF
ligand; use a covalent/QM-MM or modelled phosphointermediate. GCK and AdK have
non-covalent ligands already in the PDB → start with those.

## Honest scope
- MM-GBSA ΔΔG is a **ranking/triage** signal, not an experimental verdict; ±3–5
  kJ/mol is within noise. Use ≥3 replicas and report spread.
- Direction (activate/inhibit) is geometry-dependent and emergent — confirm with
  the wet-lab kcat/Km ± Zn²⁺/Ni²⁺/Co²⁺ assay (the paper's actual validation).
- Zn parametrisation choice changes the number; state which model you used.
