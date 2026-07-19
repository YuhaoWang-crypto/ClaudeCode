# Worked campaign: cardiac Troponin I (cTnI) diagnostic aptamer — specificity-first SELEX

**Use case:** point-of-care myocardial-infarction (MI) diagnostic. The whole clinical value of a
cTnI assay is *not* cross-reacting with the skeletal-muscle troponin I isoforms (TNNI1 slow,
TNNI2 fast) that spill into blood after any muscle injury. This is the textbook analytical-
specificity problem, so it is the ideal real target for the specificity-first workflow.

- **On-target:** cTnI / TNNI3 (UniProt P19429), cropped to residues 1–90.
- **Counter-targets (paralogs):** TNNI2 fast-skeletal (P48788) 1–70; TNNI1 slow-skeletal (P19237) 1–70.
- **Unrelated / electrostatic control:** hen-egg lysozyme (pI ~11).

## Stage 1 — MSA: locate the cardiac-specific (divergent) epitope
Pairwise-aligned cTnI vs both skeletal isoforms (`tni/msa.py`, BLOSUM62):
- cTnI core is **58 % identical to TNNI2, 64 % to TNNI1** — the conserved core is exactly what a
  naive binder would hit and cross-react on, so it must be **avoided**.
- The alignment head shows cTnI residues 2–31 align to a **gap in both skeletal isoforms** — the
  cardiac-unique N-terminal extension.
- **Cardiac-divergent block (differs from BOTH skeletal): residues 4–32 `GSSDAAREPRPAPAPIRRRSSNYRAYATE`.**
  This carries the PKA phospho-site (RRSS) — a bona fide cardiac signature. → **epitope = cTnI N-term (1–32).**

## Stages 2–3 — pool aimed at the divergent epitope
Six folded 33-nt RNA candidates (`tni/design.py`, ViennaRNA fold-checked, GC-clamped 5-bp stem +
recognition loops biased toward the Arg/Pro epitope) plus a scramble decoy.

## Stages 4–5 — Boltz 4-target screen + calibrated gate  (25 co-folds)
`examples/ctni_diagnostic.json`, gated with `scripts/specificity_gap.py` (decoy baseline 0.749):

| cand | cTnI | TNNI2 fast | TNNI1 slow | Lyso | gate |
|---|---|---|---|---|---|
| scramble (decoy on cTnI) | **0.749** | — | — | — | (baseline) |
| **tni1** | **0.899** | 0.798 | 0.701 | 0.768 | ✅ **PASS** |
| tni2 | 0.777 | 0.754 | 0.775 | 0.474 | ❌ no decoy gap + TNNI1 |
| tni3 | 0.615 | 0.680 | 0.659 | 0.698 | ❌ below decoy |
| tni4 | 0.812 | 0.873 | 0.786 | 0.603 | ❌ TNNI2 wins |
| tni5 | 0.629 | 0.779 | 0.661 | 0.573 | ❌ below decoy + TNNI2 |
| tni6 | 0.767 | 0.256 | 0.822 | 0.740 | ❌ TNNI1 wins |

**Note the scramble scores 0.749 on cTnI** — the Arg-rich cardiac N-term is an electrostatic sink,
so absolute cTnI ipTM is nearly meaningless and the decoy gate is what does the work.

## Result
**`tni1 = GGCAGAUGUCGAUCAAGCCUGUUCCUGUCUGCC` is the single survivor** — it beats its matched
scramble decoy by +0.15 AND beats fast-skeletal (+0.10), slow-skeletal (+0.20), and the unrelated
protein (+0.13). This is a **cleaner** outcome than the CTLA-4 demo (where the real WT aptamer
failed outright), because here Stage 1 pointed the pool at a genuinely cardiac-unique epitope
instead of a conserved motif — the MSA step earns its keep.

## Honest scope (unchanged)
`tni1` is a **provisional, well-prioritised starting candidate**, not a validated specific binder:
single Boltz sample; the labile N-terminal epitope is proteolytically clipped in real serum (a wet-
lab assay would capture on a stable epitope or full-length antigen); and specificity is only *proven*
by counter-SELEX (pre-clear on TNNI1+TNNI2+serum → positive-select on cTnI, toggle 8–15 rounds) plus
SPR Kd on cTnI ≪ Kd on each skeletal isoform. The in-silico loop **prioritised**; the bench **confers**.
