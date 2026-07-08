# Worked example — obesity multi-target campaign (pilot → confirm → ADME → scale)

A real end-to-end run across 5 targets, showing the pattern and the numbers.

## Targets & lead modality
| Target | Layer | Lead modality |
|---|---|---|
| GPR75 (O95800) | orphan GPCR | small molecule / siRNA |
| INHBE (P58166) | ligand/transcript | siRNA (+ Ab trap) |
| ACVR1C/ALK7 (Q8NER5) | receptor kinase | ECD antibody / SMKI |
| ActRIIB / myostatin-GDF8 (O14793) | muscle | anti-MSTN nanobody |
| GDF15 (Q99988) | metabolic | agonist antibody |

## Step 1 — pilot design (≤100–200 each)
- SM design vs GPR75 (auto-pocket + 20-HETE `reference_ligands`) and ALK7 ATP-site (SB-431542 reference), 200 mol, ~$5 each.
- Protein design (boltz_curated): anti-MSTN nanobody (100), ALK7-ECD antibody (100), GDF15 antibody (100).
- Targets are the **mature domain** (slice off signal peptide / pro-domain): GDF8 res 267–375, ALK7 ECD 22–113, GDF15 197–308.

## Step 2 — full re-rank
Paginate `boltz_get_job_results` over ALL results (jq on the saved overflow files). Top hits found deeper than page 1: MSTN nanobody ipTM 0.955 / binding 0.649; ALK7 SM binding 0.513; GPR75 SM binding 0.491.

## Step 3 — affinity confirmation (the decisive filter)
`structure_and_binding` (num_samples=3) on the top 2 binders per program. Result: interface geometry reproduced for all, but affinity corroborated for only 2 (MSTN nanobody ipTM 0.90/affinity 0.74; ALK7 antibody ipTM 0.93/affinity 0.21). GDF15 antibody + several ALK7 designs = geometry-only → deprioritised.

## Step 4 — ADME
`small_molecule_adme` on top 20 molecules ($0.20). ALK7_1 (`Cc1cc(-c2nc(-c3cncc(Br)c3)no2)cnc1O`, binding 0.513) = clean ADME → top SM lead; GPR75 best binder flagged solubility-high-risk → analog redesign.

## Step 5 — scale a winner
`protein_design` num_proteins=500–2000 on the confirmed nanobody. Watch credit: a 2000-run ($50) 402'd on a $15.68 balance; a 500-run stopped partial at 441/500 when credit ran out — partial results still retrievable.

## Lessons
- Auto-pocket + reference_ligands works even for orphan/no-structure GPCRs (GPR75).
- ADME embeds in design results — the separate ADME job mostly confirms.
- Solubility-vs-potency: aggressive descriptor filters fix solubility but cost binding; rebalance logP.
- Always check credit before scale-up; jobs stop partial, no refund.
