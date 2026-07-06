---
name: target-evidence
description: >-
  Mine and fuse binder/interaction evidence for a protein target (Track B —
  active-targeting ligand design for LNP surfaces). Use when the user names a
  receptor / CD marker / UniProt / gene and wants: measured actives (ChEMBL),
  predicted small-molecule binders (DrugCLIP genome-wide), predicted protein
  partners (humanPPI), a cross-comparison, and a cell-surface enrichment analysis —
  in one command. Triggers: "find binders/ligands for <protein>", "targeting ligand
  for <CD/receptor>", "mine DrugCLIP/humanPPI for X", "is <target>'s interactome
  surface-enriched".
---

# Target evidence mining & fusion (Track B)

Reverse-engineered public APIs (all free, network + CPU) fused into one report per
target. Engine: `lipidlib/targetpipe.py`; one-command CLI: `analysis/run_target.py`.

## One command
```bash
python analysis/run_target.py --uniprot P43220 --name GLP1R \
    --drugclip-decoy P00533 --ppi-decoys P04406,P00533
```
Resolve the UniProt from a gene via UniProt/ChEMBL if needed. Degrades gracefully
when a source has no coverage (e.g. ASGR1 has no DrugCLIP). Outputs:
- `data/targets/<NAME>/` — `*_chembl_ligands.csv`, `*_drugclip_predicted.csv`,
  `*_humanppi_partners.csv`, `*_surface_partners.csv`
- `results/figures/<NAME>_molecular.png`, `<NAME>_humanppi.png`
- `results/reports/<NAME>.md` — one-page summary

## What each source gives (also callable individually)
- `scripts/fetch_glp1r_ligands.py --uniprot <id>` → ChEMBL measured actives (pChEMBL).
- `scripts/fetch_drugclip.py --uniprot <id>` → DrugCLIP predicted hits
  (drugclip_score, docking_score, pocket residues, SMILES). API:
  `dtwgapi.yanyanlan.com/complexes/{uniprot}` + `POST /get_smiles`.
- `scripts/fetch_humanppi.py --uniprot <id>` → predicted PPI partners (AF/RF/DCA
  scores, subcellular locality). API: `prodata.swmed.edu/humanPPI/data/{uniprot}`.

## Analyses (in the engine / run_target)
- **molecular_enrichment**: nearest-neighbour Tanimoto of predicted hits to measured
  actives vs a decoy target. NB: DrugCLIP is a *fragment* library (MW ~250), so it
  usually won't match drug-sized actives — expect low overlap; that's real, not a
  bug (validate the machinery with `analysis/validate_pipeline.py`, which shows a
  1674× positive control on held-out actives).
- **fragment_enrichment**: substructure containment of fragments in actives.
- **ppi_surface_enrichment**: is the target's interactome cell-surface enriched vs
  decoys + genome baseline (Fisher / binomial)? A positive result flags surface
  co-targets/handles (e.g. GLP1R → DLK1, PAM).

## Then (optional) structural validation
Validate a shortlisted ligand/peptide ↔ receptor-ectodomain with the hosted **Boltz**
MCP tools (structure + binding affinity; no local GPU). For GPCRs/peptide-hormone
receptors, a peptide ligand to the extracellular domain is often the right handle.

See `docs/PIPELINE.md`, `docs/PROBLEM_B_GLP1R.md`, `docs/GLP1R_crosscompare.md`,
`docs/RESOURCES.md`.
