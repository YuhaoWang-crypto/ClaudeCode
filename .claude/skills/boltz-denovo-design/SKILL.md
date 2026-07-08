---
name: boltz-denovo-design
description: Run a de-novo, multi-modality drug-design campaign against a protein target using the Boltz-2 API (MCP tools). Use when you need to design small molecules (Enamine REAL space), antibodies, nanobodies, or peptide binders against a target, then triage them - independent affinity/interface confirmation, ADME prediction, and ranking. Covers the full pilot→confirm→ADME→scale loop with cost guardrails. Triggers - "design a binder/inhibitor/antibody/nanobody for <target>", "de novo drug design", "small molecule design", "screen molecules", "affinity confirmation", "ADME", "developability triage".
---

# Boltz-2 de-novo drug design — pipeline & recipes

Uses the Boltz API MCP tools (load via ToolSearch `select:mcp__Boltz_API__…`). Every `*_start_*` tool costs credit; every `*_estimate_*` tool is free. Always estimate before starting, and gate on cost.

## Core loop
1. **Design** (generate candidates against the target)
2. **Confirm** (independently re-predict top hits → keep real binders, drop geometry-only)
3. **ADME / developability** (small molecules) or **immunogenicity** (biologics)
4. **Rank & scale** the survivors

## Tools by step
| Step | Tool | Notes |
|---|---|---|
| Small-molecule design | `boltz_estimate/start_small_molecule_design` | `target.type="no_template"` + protein entity; omit `pocket_residues` to auto-detect, or pass `reference_ligands` (SMILES) to bias the pocket/chemotype; `chemical_space="enamine_real"`; `molecule_filters` (Lipinski / rdkit descriptors for solubility). |
| Antibody / nanobody / peptide design | `boltz_estimate/start_protein_design` | target `no_template` + mature-domain entity; `binder_specification=boltz_curated` (`boltz_nanobody`/`boltz_antibody`) is the proven default; optional `epitope_residues` to steer to a functional interface. |
| Structure + affinity confirm | `boltz_estimate/start_structure_and_binding` | re-predict binder+target complex (`binding=protein_protein_binding`, `num_samples=3`); the decisive filter. |
| ADME | `boltz_start_small_molecule_adme` | adme-v1: solubility class, permeability, lipophilicity. |
| Screen a library | `boltz_start_small_molecule_screen` / `protein_screen` | score a provided set against the target. |
| Status / results | `boltz_get_job_status`, `boltz_get_job_results` | results are URL-heavy → they overflow to a file; use `jq` on the saved file, never paste full payloads. |

## Reading the metrics
- Protein design: rank by **ipTM** (>0.9 excellent interface) → **structure_confidence** → **min_interaction_pae** (lower, <2 Å good). `binding_confidence` is a conservative affinity proxy.
- Small-molecule design: rank by **binding_confidence** → **optimization_score** → **structure_confidence**; also `iptm`, `smiles`, embedded ADME.
- **Key lesson:** high ipTM can be geometry-only — the independent `structure_and_binding` confirm is what separates real binders (affinity corroborated) from artifacts. In a real run, only 2 of 6 high-ipTM biologics survived affinity confirmation.

## Cost guardrails (observed rates)
- Small molecule design ≈ $0.025/molecule; nanobody ≈ $0.025/design; antibody ≈ $0.05–0.10/design; ADME ≈ $0.005–0.01/mol; structure+binding ≈ $0.05–0.10/complex.
- ALWAYS estimate first; if a run exceeds budget or account credit, scale `num_*` down. Boltz has **no job-stop/refund** — a job that outruns credit just **stops partial** (results still retrievable). Check credit before large runs; production scale (1k+) needs a top-up.
- Pilot at ≤100–200 to prove the flow; scale winners only.

## Recipe: solubility-optimized analogs
Re-run `small_molecule_design` with `reference_ligands=[<best hit SMILES>]` + `rdkit_descriptor_filter` (e.g. `mol_logp≤3, tpsa≥60, mol_wt 300–460`). Trades potency for solubility — loosen logP→3.5 if binding drops too far.

See `references/workflow.md` for the full worked obesity-program example (5 targets, pilot→confirm→ADME→scale).
