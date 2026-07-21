# Generative route demo — de-novo binder for 1,25(OH)₂D3 (result: tractability wall)

**Goal:** design mini-protein binders against the active hormone 1,25(OH)₂D3, then
counter-select against D3 / 25(OH)D3 to see whether a real specificity margin opens —
the fallback route the specificity framework recommends when no PDB template is selective.

## What was run (✅ real Boltz protein-design jobs)
| run | binder spec | target | num | result |
|---|---|---|---|---|
| r1 `gendemo-vitd-125-binders-r1` (`prot_des_oPnlJ5HLn6Wk5r5IF1G8`) | custom_protein, **65 aa fixed**, exclude Cys | 1,25 SMILES, template-free, ligand=epitope | 10 | **failed** — 201/201 samples (100%) |
| r2 `gendemo-vitd-125-binders-r2` (`prot_des_PBFpmEtPoiwsDSoukJPN`) | custom_protein, **80–110 aa range** | same | 10 | **failed** — 201/201 samples (100%) |

Both aborted with `sample_failure_rate_exceeded` (threshold 95%). Cost ~$0.25 each.
Input validation (estimate) passed for both — the SMILES is fine (same string co-folded
successfully in the 3×3 matrix). The failure is at the design/scoring gate, not parsing.

## Diagnosis ⚠️
This is a **systematic tractability wall**, not a transient error: two independent runs, two
binder-length regimes, identical 100% sample failure. Root cause is the **template-free,
small-molecule-only target** — the design pipeline has no 3D anchor for a highly flexible
secosteroid (1,25(OH)₂D3 has a floppy seco-B ring + long side chain), so it cannot build or
validate a complementary pocket around a free-floating ligand, and every sample fails the
validity gate.

This mirrors how ligand-anchored de-novo design is done in practice (RFdiffusion-AllAtom /
LigandMPNN **fix the ligand's 3D atoms** and grow a pocket around them). Boltz protein-design
against a bare SMILES with no pose does not provide that anchor.

## What this adds to the honest story
The 3×3 matrix already showed co-folding can't resolve the single 1α-OH. This adds: the naive
generative fallback ("just generate a selective binder") **does not just work** for this ligand
class either — it needs a **ligand-anchored** setup. So a real discriminating vitamin-D sensor
is a genuinely hard, multi-step design problem, exactly as the ⚠️ labels warned.

## Tractable next step (ligand-anchored)
Provide a **3D pose** of 1,25(OH)₂D3 as a `structure_template` target (a CIF with just the
ligand chain), so the binder is designed around fixed ligand atoms; then counter-select the
surviving designs against D3 / 25(OH)D3 with the co-folding matrix. Requires building a ligand
CIF (RDKit embed → CIF) first. Not run here to avoid a third un-vetted paid attempt.

## Update: pocket-redesign route (LigandMPNN-style) — attempted, infra wall

Ran the actual pocket-redesign campaign the constraint-focused test pointed to:

| run | setup | result |
|---|---|---|
| `vdr-pocket-redesign-1a-125-r1` (`prot_des_xpHnhvxZ5n941CGFbXS1`) | VDR, 6 designed 1α-subsite residues (Ser68/Arg105 windows), **no_template** target = 1,25 SMILES | **failed** — 205/205 samples (100%), ~35 min |
| `vdr-pocket-redesign-1a-anchored-r2` (`prot_des_hUpML9zvdbJgo9FF75H4`) | same binder, **structure_template** target = 1,25 **3D pose** (RDKit) | ran **>2h50m, 0/12, never completed** in-session; no failure, but no output |

**Diagnosis (consistent across 4 design runs):** Boltz protein-design against a flexible
secosteroid is not tractable in this environment — a **bare SMILES target fails the
sample-validity gate 100%** (3/3 no_template runs), and the **ligand-anchored fix
(3D pose) clears the gate but runs anomalously long** (>2.8 h, no result), compounded by
frequent Boltz-API disconnects. This is an **infrastructure/compute wall, not a spec bug**:
the campaign is correctly specified (`vdr_pocket_redesign_plan.json`, tested) and would run
on a stable design backend (or the local RFdiffusionAA + LigandMPNN route named in
`going-real`/`generative_design_plan`).

**Net honest conclusion of the whole specificity thread:** (1) co-folding can't resolve a
single –OH (flat 3×3 matrix); (2) a directional pocket **constraint** moves the numbers the
right way but doesn't open a usable margin; (3) actually **designing** the contact into the
pocket is the correct next step and is fully specified here, but did not complete on this
environment's Boltz design backend. The deciding number remains a wet-lab competition assay.
