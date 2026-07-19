# Specificity-first SELEX: MSA-divergent epitope + counter-SELEX + calibrated in-silico

For diagnostic/therapeutic aptamers where the target has close paralogs (CTLA-4/CD28,
EGFR/HER2/HER3, PD-L1/PD-L2, GFRα1/α2/α3). **Specificity is EVOLVED by counter-selection
and TARGETED at the divergent surface — it is not something an in-silico score can confirm.**

## The lesson that motivates this (CTLA-4 case, verified in this repo)
A real RNA aptamer (`aptamerd6`) against CTLA-4 was cross-validated against paralogs and an
unrelated protein with BOTH scorers:

| scorer | CTLA-4 (target) | CD28 | PD-L1 | Lysozyme (unrelated) |
|---|---|---|---|---|
| Boltz ipTM | 0.514 | 0.851 | 0.848 | 0.814 |
| HDOCK score | −313.8 | — | — | **−382.1** (stronger!) |

It scored **better on unrelated proteins than on its own target** on both methods →
`specificity_gap.py` verdict: **PROMISCUOUS**. Mechanism: a polyanionic RNA sticks to any
cationic small protein (lysozyme pI ~11) by electrostatics. **Takeaways:**
- Absolute ipTM/HDOCK scores are meaningless here; a high on-target score + positive
  decoy gap (aptamerd6 had +0.25) is NOT specificity.
- Both AI (Boltz) and physics (HDOCK) can be non-discriminative in the SAME direction —
  "consensus" does not rescue a promiscuous binder. Only a **paralog + unrelated panel** does.

## Stage 1 — MSA: locate the specificity-determining epitope (target the DIFFERENCE)
Align the target with its paralogs; the **divergent, surface-exposed** columns are the
specificity handle. AVOID conserved functional motifs (they give pan-family cross-reactivity).
- Tooling: `Bio.Align.PairwiseAligner` (BLOSUM62) for a pair; MAFFT/Clustal for ≥3 members.
- Worked example (CTLA-4 vs CD28 IgV, this repo): **27% identity**; the **MYPPPY B7-binding
  motif (CTLA-4 pos 97–102) is CONSERVED → do NOT target it** (targeting it is likely why
  aptamerd6 cross-reacts). CTLA-4-specific epitope blocks to target instead include
  `49–57 AATYMMGNE` (CD28 `VVYGNYSQQ`), `24–30 ASPGKAT`, `103–107 YLGIG` (CD28 `LDNKS`).
- Output: a target-epitope residue list (feed to Boltz `pocket`/`contact` constraints and to
  the wet-lab immobilization/epitope-masking strategy).

## Stage 2 — Database seed + biased pool
- Search **Apta-Index / UTexas / PubMed** for aptamers vs the target or its family (CTLA-4:
  27 PubMed hits) as scaffold seeds — but treat family binders as pan-reactive starting
  points to be specialized, not final answers.
- Or generate a motif-biased pool with `model/generate_lm.py`. Fold-filter (ViennaRNA).

## Stage 3 — Doped mutation (not full random)
Partially randomize the recognition loops of each seed (15–30 %/position); keep the folded
stem/scaffold. This is affinity/specificity maturation space, not a naive N40 library.

### Worked maturation loop (this repo, verified) — 8 doped mutants of aptamerd6
Doped `aptamerd6` at ~20 %/position on the recognition loops (keeping the two stems), aimed
the diversified positions at the CTLA-4-divergent epitope, and co-folded all 8 against the
3-target panel (CTLA-4 / CD28 / Lysozyme) with Boltz, then gated with `specificity_gap.py`
(decoy baseline 0.267). Full 8×3 ipTM matrix (`examples/ctla4_maturation.json`):

| mut | CTLA-4 | CD28 | Lyso | gate |
|---|---|---|---|---|
| WT (aptamerd6) | 0.514 | 0.851 | 0.814 | ❌ CD28+Lyso win |
| m1 | 0.717 | 0.789 | 0.785 | ❌ CD28 wins |
| m2 | 0.380 | 0.289 | 0.531 | ❌ Lyso wins |
| m3 | 0.714 | 0.768 | 0.540 | ❌ CD28 wins |
| m4 | 0.265 | 0.786 | 0.610 | ❌ no decoy gap + CD28 |
| m5 | 0.717 | 0.727 | 0.484 | ❌ CD28 wins (−0.01, marginal) |
| **m6** | **0.800** | **0.551** | **0.735** | ✅ **PASS** (+0.25 vs CD28, +0.07 vs Lyso) |
| m7 | 0.751 | 0.848 | 0.800 | ❌ CD28+Lyso win |
| m8 | 0.148 | 0.687 | 0.620 | ❌ no decoy gap + CD28 |

**Result: 7/8 screened OUT, 1 survivor (m6 = `UGUACAGAGGGCUGGUAACGAUCCGGAUAAGAA`).** This is
exactly what the gate is for — it kills the promiscuous majority. But note m6's margin over the
*unrelated* protein is thin (+0.07, one Boltz sample), so it is a **provisional** in-silico
survivor, NOT a proven specific binder: it earns a seat in the next counter-SELEX round, nothing
more. The loop **prioritises**; the bench (Stage 4) **confers** specificity.

### Round 2 — re-dope the survivor (m6 → 8 grandchildren), verified
Took m6 as the new seed, folded it (ViennaRNA: one stem-loop, everything else loop), kept the
stem fixed and doped the loops ~25 %/position → 8 grandchildren, same 3-target Boltz gate:

| mut | CTLA-4 | CD28 | Lyso | gate | worst-off margin |
|---|---|---|---|---|---|
| m6 (parent) | 0.800 | 0.551 | 0.735 | ✅ | **+0.065** (vs Lyso) |
| r2_1 | 0.467 | 0.465 | 0.259 | ❌ | +0.002 (vs CD28) |
| r2_2 | 0.633 | 0.849 | 0.701 | ❌ | −0.22 |
| r2_3 | 0.464 | 0.787 | 0.674 | ❌ | −0.32 |
| **r2_4** | 0.644 | 0.533 | 0.433 | ✅ | **+0.111** (vs CD28) |
| r2_5 | 0.731 | 0.828 | 0.741 | ❌ | −0.10 |
| r2_6 | 0.363 | 0.777 | 0.569 | ❌ | −0.41 |
| r2_7 | 0.248 | 0.375 | 0.743 | ❌ | −0.49 |
| r2_8 | 0.635 | 0.738 | 0.712 | ❌ | −0.10 |

**Again 7/8 screened out, 1 survivor — but the survivor's worst-off margin nearly doubled
(+0.065 → +0.111): `r2_4 = UGCUCACCCAGCUGGCAACUAUCCGGAAAUGAC`.** Two honest caveats: (1) the
gain came at the cost of *on-target* confidence (0.800 → 0.644) — the loop traded absolute binding
for selectivity, exactly the affinity/specificity frontier real SELEX walks; (2) +0.11 is still
one Boltz sample, so r2_4 is a *better-prioritised* provisional survivor, not a validated specific
binder. Takeaway: **iterated in-silico maturation demonstrably improves the selectivity margin but
plateaus fast (7/8 still die each round) — which is exactly why the decisive selectivity gain must
come from wet-lab counter-SELEX, with the gate used only to rank.**

## Stage 4 — ★ Positive + counter (toggle) SELEX — where specificity is EVOLVED ★
Each round:
1. **Negative / pre-clear:** flow the pool over immobilized **paralog(s) (e.g. CD28) + serum
   proteins + bare matrix**; DISCARD binders.
2. **Positive:** collect binders to the **target (CTLA-4)**, ideally presented so the
   divergent epitope (Stage 1) is accessible and shared motifs are masked.
3. **Toggle:** alternate targets across rounds and add **soluble paralog competitor**;
   increase stringency late (lower target conc., more washes).
8–15 rounds. This is the only reliable source of specificity — a naive library + positive-only
selection gives family cross-reactors (as the in-silico panel predicts).

## Stage 5 — NGS-guided tracking + CALIBRATED in-silico triage
- Deep-sequence each round; keep families **enriched on target AND depleted in the CD28
  counter-selection**.
- In-silico triage with `scripts/specificity_gap.py` — the ONLY trustworthy readout:
  a candidate PASSES only if, on **every** scorer, on-target beats (a) its matched scramble
  decoy on the same target, and (b) **every** off-target (paralog + unrelated). Any off-target
  win = PROMISCUOUS → send back to counter-SELEX. Never rank by absolute ipTM/HDOCK.

## Stage 6 — Validation
- Synthesize survivors (+ 2'-F/2'-OMe for RNA, 3'-invdT). SPR/BLI Kd on target AND on each
  paralog; require target Kd ≪ paralog Kd. Cell assays on target⁺ vs paralog⁺ lines.

## Files
- `scripts/specificity_gap.py` — the calibrated specificity gate (example: `examples/ctla4_crossval.json`).
- MSA: `Bio.Align.PairwiseAligner`; see the CTLA-4/CD28 worked numbers above.

## Honest scope
In-silico here **screens out** promiscuous candidates and points the library at the right
epitope; it does **not** confer or prove specificity. Specificity is won at the bench by
counter-SELEX. Every score is a relative triage signal, never an affinity.
