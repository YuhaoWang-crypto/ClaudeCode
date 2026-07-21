# Designing metabolite-discriminating receptors (specificity engineering)

How to build receptors that tell **closely-related metabolites apart** — e.g.
the vitamin-D series D3 → 25(OH)D3 → 1,25(OH)₂D3 → 24,25(OH)₂D3, which differ by
a single hydroxyl. Implemented in `biosensor_pipeline/specificity.py`.

## Why this is a *negative-design* problem

Distinguishing 25(OH)D3 from 1,25(OH)₂D3 means recognizing **one extra A-ring
1α-OH**. Maximizing on-target affinity is not enough — a tight pocket for
25(OH)D3 usually also fits 1,25(OH)₂D3. Specificity is a **margin**:
`specificity = affinity(on-target) − max affinity(off-targets)`.
You have to *design against* the off-targets (counter-selection).

## Route A — mine a naturally-selective receptor from the PDB

`discover.py` finds proteins already co-crystallized with a specific metabolite,
and nature has *already solved* several of these discriminations:

| target metabolite | selective template (PDB) | recognizes |
|---|---|---|
| **1,25(OH)₂D3** (active hormone) | **VDR** nuclear-receptor LBD (1DB1, 7QPI, 259–287 aa) | the A-ring **1α-OH** |
| **25(OH)D3** (clinical status marker) | **DBP / GC** vitamin-D binding protein (1J78, ~458 aa); anti-25(OH)D antibody Fv | 25-OH, discriminates vs 1,25 |
| **Vitamin D3** | designed **CDL2.x** binders (5IEN/O/P, 137 aa); β-lactoglobulin (2GJ5) | the parent |

Extract the ligand-binding domain (VDR-LBD, DBP domain, an scFv, or a designed
CDL), then run it through the biosensor pipeline (chimera → reporter). VDR-LBD is
large (~260 aa) → prefer the CP-reporter/terminal-fusion or split topology, or a
smaller antibody-fragment / designed binder.

## Route B — generate a selective binder de novo (no good template)

For a metabolite with no selective PDB binder (e.g. 24,25(OH)₂D3):

1. **Scaffold** — RFdiffusion-AllAtom (RFdiffusionAA) or Boltz protein-design
   against the target metabolite, building a pocket complementary to the
   **discriminating group**.
2. **Sequence** — LigandMPNN designs pocket residues that H-bond the
   distinguishing –OH and leave no room for the off-target substituents.
3. **Positive filter** — co-fold vs the target (Boltz/AF3); keep high interface
   confidence with correct contact to the discriminating –OH.
4. **Counter-select (the key step)** — co-fold survivors vs **every** off-target;
   keep only those with a large specificity margin. This is what single-target
   design skips and what actually buys discrimination.
5. Feed the selective binder into the biosensor pipeline.

Available here via the Boltz protein-design MCP tools and the sibling
`boltz-denovo-design` / `rfantibody-epitope-campaign` skills.

## The shared scoring core — cross-panel co-folding matrix

Whichever route, the deciding number is a **receptor × metabolite matrix**:
co-fold each candidate against **every** metabolite (Boltz holo,
ligand_protein_binding), read `ligand_iptm` / binding probability, and rank by
`specificity_score` (on-target − best off-target).

```bash
python3 -m biosensor_pipeline.specificity     # panel, templates, matrix framework, plans
```

`boltz_specificity_plan()` emits the jobs (e.g. 3 receptors × 5 metabolites = 15
holo predictions, ~$0.75) → fill `specificity_score()` → `specificity_matrix()`
ranks and labels each receptor selective / weakly-selective / non-selective.

## What the real matrix showed (3×3 run, ✅ Boltz-2.1)

A 9-job run (VDR-LBD 1DB1 · CDL2.2 · DBP 1J78 × D3 / 25(OH)D3 / 1,25(OH)₂D3)
delivered the honest verdict this section predicts: **co-folding does not resolve a
single hydroxyl here.** Every specificity margin was within method noise (±0.05–0.08),
and in two of three receptors a *wrong* metabolite scored highest. Two useful signals
survived: (1) VDR ranked its on-target hormone 1,25(OH)₂D3 highest (ligand_iptm 0.981,
correct affinity order D3<25(OH)<1,25) — a tiny-margin positive; (2) the *designed D3
binder CDL2.2 was not D3-selective* — it scored the hydroxylated metabolites far higher
(0.92/0.84 vs 0.46), a concrete counter-selection red flag. Full numbers:
`biosensor_out/specificity_vitd_matrix.json` + `specificity_vitd_RESULTS.md`.

Takeaway: no PDB template cleanly separated 25(OH) from 1,25(OH)₂ by co-folding alone,
which is exactly what pushes a real design to the **generative + counter-selection**
route (`generative_design_plan`).

## Directional-contact test (borrowed from Baker de-novo design, ✅ Boltz)

Two Baker papers (de-novo luciferases, Nature 2023 `s41586-023-06735-9`; modular
peptide binders, `s41586-023-05909-9`) teach that specificity comes from a **specific
directional polar contact** to the discriminating group — not global affinity. Retested
on VDR by pinning the ligand to VDR's native **1α-OH anchors Ser237+Arg274** (our
construct residues 68+105) with a forced Boltz `pocket` constraint, then comparing to
the unconstrained co-folds:

- Focusing the score on the 1α-subsite moved every number the right way: **1,25 binding-
  confidence +0.009** (satisfied by its 1α-OH), **25(OH)D3 −0.004** (penalized — no 1α-OH);
  the hard 1,25-vs-25(OH) margin corrected from **−0.014 (wrong) → −0.000 (tied)**.
- ⚠️ but still **not a usable margin** — a single –OH stays at/below co-folding resolution
  even with a mechanism-focused constraint. Scoring *focuses* the metric; only **designing
  the contact into the pocket** (RifDock/LigandMPNN + counter-selection) manufactures the
  margin. Full write-up: `biosensor_out/specificity_vitd_pocket1a_RESULTS.md`.

## Honesty

- ✅ the matrix numbers are real Boltz metrics; the SMILES and PDB templates are real.
- ⚠️ resolving a **single –OH** by co-folding interface confidence is at the edge
  of the method's resolution (confirmed empirically above). Treat the ranking as a
  **prioritization**, and confirm selectivity with a **competition / cross-reactivity
  assay** (the wet-lab ground truth) — the same discipline as the rest of the skill.
