---
name: allosteric-pocket-finder
description: >-
  Find distal (allosteric-candidate) conformational hotspots from an active/inactive
  protein structure pair and propose His-pair/His-triplet metal-switch engineering sites.
  Use when the user wants to download PDB structure pairs from RCSB, run structural
  alignment to locate remote conformational responses to active-site binding, rank
  5-20 distal candidate regions, or design Zn/Ni/Co metal-switch sites on hinge/interface
  loops. Triggers: allosteric pocket, hinge, conformational change, apo vs bound / active
  vs inactive alignment, Ca displacement, distal hotspot, His-pair, His-triplet, metal switch.
---

# Allosteric-pocket / metal-switch finder

Reproducible pipeline: **download an active/inactive structure pair → robust Cα
superposition → rank 5–20 distal candidate regions → flag His-pair/His-triplet metal
sites**. Built for the allosteric-pocket / Zn-switch demo but works on any protein pair.

## When to use
- The user gives (or you can find) two PDB structures of the same protein in different
  states (apo vs ligand-bound, open vs closed, inactive vs active) and wants to know
  **where the protein moves away from the active site** (allosteric hotspots).
- The user wants to design a **His₂/His₃ metal-binding site** on a mobile hinge/loop to
  make catalysis metal-switchable, then screen Zn²⁺/Ni²⁺/Co²⁺ on kcat/Km.

## Requirements
```bash
pip install biopython numpy
```

## Workflow

### 1. Download structures + sequences
```bash
# explicit PDB IDs:
python3 scripts/download_structures.py --pdb 1V4S 1V4T --out data
# or the bundled 7-system manifest (GCK, PTP1B, AdK, PFK, TEM-1, ATCase, GP):
python3 scripts/download_structures.py --manifest assets/targets.json --out data
```
Writes `data/pdb/*.pdb` and `data/fasta/all_targets.fasta`. Idempotent + retrying.

### 2. Find distal candidate regions
Define the active site with **either** a ligand (HETATM resname in the *bound* structure)
**or** catalytic residue numbers.
```bash
# single pair, active site = a bound ligand:
python3 scripts/find_candidates.py --bound 1V4S --unbound 1V4T --active-ligand GLC --name GCK
# single pair, active site = catalytic residue(s):
python3 scripts/find_candidates.py --bound 1T49 --unbound 1SUG --active-residues 215 --name PTP1B
# batch mode from a JSON config (see assets/pairs.json for the schema):
python3 scripts/find_candidates.py --pairs assets/pairs.json
```
Options: `--chain A` `--pdb-dir data/pdb` `--out results` `--top 12`.

Outputs per system in `results/`:
- `<name>_candidates.csv` — per-residue Cα displacement, distance-to-active-site, core flag.
- `<name>_regions.json` — ranked distal regions + `his_pair_candidates` (Cα–Cα 4.5–12 Å).
- `summary.json` — RMSD stats across all systems.

### 3. Prioritise and design (human/wet-lab step)
Pick 3–6 of the ranked regions that are **distal + high-displacement + at a hinge/interface
+ surface loop + contain a plausible His-pair**. Turn the flagged Cα–Cα pairs into His
substitutions (validate side-chain rotamer geometry before ordering — spacing is a coarse
filter). Build a small His-pair/His-triplet library, purify, strip adventitious metal
(EDTA + Chelex buffer), then measure **kcat, Km, kcat/Km** across a Zn²⁺/Ni²⁺/Co²⁺ series.

## Method notes / caveats
- **Robust core superposition**: iteratively drops the worst ~30% of residues each round so
  the reference frame locks onto the rigid domain — essential for hinge motions, otherwise
  the fit averages across the moving parts and washes out the signal.
- The retained core fraction is **imposed** by `CORE_REJECT_FRAC` (top of the script), not
  discovered; a low core-RMSD confirms the retained set is genuinely rigid.
- **Two static crystal snapshots** capture endpoint displacement, not dynamics. For
  subtle-motion systems (e.g. PTP1B, overall RMSD < 1 Å) raw Cα displacement under-reports
  coupling — complement with loop RMSF, contact-map differencing, or normal-mode analysis.
- **Single chain only.** Oligomeric/interface allostery (PFK, ATCase, GP) needs multi-chain
  superposition + inter-subunit contact-map differencing — not yet implemented here.
- Cα–Cα spacing is a **coarse** His-site proxy: the Cβ vectors must point toward a shared
  coordination point (His–metal ≈ 2.0–2.2 Å, His–metal–His ≈ 90–120°).

## Tunable thresholds (top of `scripts/find_candidates.py`)
`DISPLACEMENT_MIN=2.0` (moved), `DISTAL_MIN=15.0` (distal from active site),
`CORE_REJECT_FRAC=0.30`, `CORE_ITERS=5`, `HIS_PAIR_MIN/MAX=4.5/12.0` Å.

## Bundled assets
- `assets/targets.json` — 7-system manifest (PDB pairs, UniProt IDs, roles, references).
- `assets/pairs.json`   — ready-to-run batch config for the three monomeric positive
  controls (GCK, PTP1B, AdK). GCK recovers its known activator pocket in the top regions.
