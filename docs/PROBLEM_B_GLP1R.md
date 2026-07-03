# Problem B pilot — GLP1R-targeted LNP surface ligand

Goal: find/validate a ligand that binds the **GLP1R ectodomain** so it can be
conjugated to the LNP surface (via PEG-lipid) and drive receptor-mediated uptake
into GLP1R⁺ cells (pancreatic β-cells; also GI, some CNS). This is the
active-targeting track — **separate from the ionizable-lipid model** (Problem A).

## Target facts (verified via ChEMBL)

| | |
|---|---|
| ChEMBL target | `CHEMBL1784` (single protein, confidence 9) |
| UniProt | `P43220` (AlphaFold model available) |
| Class | Class B (secretin-like) GPCR — large extracellular domain (ECD) |
| Structures | ~60 PDB entries incl. ECD-only (e.g. `3IOL`, `3C59`) and full 7TM cryo-EM |
| Internalization | Agonist binding → β-arrestin-mediated endocytosis ✅ (good for delivery) |

The large **ECD is the accessible, surface-facing binding site** — the right
target for a conjugatable ligand, versus the 7TM pocket buried in the membrane.

## Data we already have (`data/targets/GLP1R/GLP1R_ligands.csv`)

1,422 unique ligands at pChEMBL ≥ 6, pulled + deduped from ChEMBL:

- **729 peptide/large** — GLP-1 analogs (exenatide, semaglutide-like). The most
  potent binders (pChEMBL up to 11, i.e. sub-pM). **A GLP-1-derived peptide is
  the natural, highest-affinity targeting ligand** and is directly conjugatable.
- **693 small molecules** — modern oral agonists (danuglipron/orforglipron class,
  MW ~600, pChEMBL ~10.9). Smaller, cheaper to conjugate, but fewer and mostly
  bind the 7TM/interface rather than the ECD.

Rebuild anytime: `python scripts/fetch_glp1r_ligands.py --min-pchembl 6`.

## Key design decision

Because the top binders are **peptides** and GLP1R is a peptide-hormone receptor,
this pilot should probably lead with a **peptide targeting ligand** (GLP-1(7-37)
or a stabilized analog) rather than a small molecule. This also cleanly matches
your "small molecule OR peptide library" question — for GLP1R, peptide wins.

## Pipeline (Problem B)

```
known binders (ChEMBL ✅)  ─┐
                            ├─►  actives + decoys  ─►  retrieval model (DrugCLIP)
molecule/peptide library  ─┘                            rank candidates
                                                              │
                                                              ▼
                                        structural validation of top-k
                                        Boltz co-fold (ECD + ligand) → pose + affinity
                                                              │
                                                              ▼
                                        shortlist → PEG-lipid conjugation design
```

### Step 1 — actives/decoys (ready now)
`GLP1R_ligands.csv` = actives. Generate property-matched decoys (DUD-E style) for
a retrieval/classifier baseline.

### Step 2 — retrieval (DrugCLIP, GPU/Modal)
Embed the GLP1R ECD pocket + a candidate library; rank by contrastive similarity.
Code: `bowen-gao/DrugCLIP`. Uni-Mol based, needs GPU → Modal (Phase 5).

### Step 3 — structural validation (Boltz, hosted MCP — no local GPU)
For each top candidate, co-fold with the GLP1R ECD sequence to get a predicted
pose + binding-affinity readout. The Boltz MCP tools are live in this session:
`boltz_start_structure_and_binding` → `boltz_get_job_status/results`. Validate the
peptide GLP-1 analog first as a positive control, then novel candidates.

## Immediate next actions
- [ ] Pull the GLP1R **ECD sequence** (UniProt P43220, residues ~24–145) for Boltz.
- [ ] Boltz positive-control: co-fold GLP-1(7-37) ↔ ECD; confirm sensible pose/affinity.
- [ ] Generate decoys + a DrugCLIP-ready actives file.
- [ ] Decide peptide-first vs. small-molecule-first for the screening library.
