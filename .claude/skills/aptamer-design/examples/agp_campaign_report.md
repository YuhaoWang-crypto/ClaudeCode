# Worked campaign: α1-acid glycoprotein (AGP / ORM1, PDB 3KQ0) aptamer library

**Target:** human α1-acid glycoprotein (orosomucoid, ORM1; UniProt P02763; PDB 3KQ0 =
human AGP crystal structure). Lipocalin β-barrel, major acute-phase inflammation biomarker
and plasma drug-binder. **Two hard features:** (1) a near-identical paralog ORM2 (89% id),
and (2) ~45% carbohydrate — **5 N-glycans** shield much of the barrel surface.

## Stage 1 — MSA + glycan mapping (`agp/msa.py`)
- ORM1 vs ORM2 = **89% identical** → very hard paralog-discrimination problem.
- N-glyc sequons (Asn, mature #): **15, 38, 54, 75, 85** → glycan-shielded, must avoid.
- One clean ORM1-divergent, non-glycosylated block: **residues ~112–117 `LAFDVN`**
  (ORM2 has `FGSYLD`). Mildly acidic/hydrophobic loop (not an Arg cluster) → chosen epitope.

## Stage 2–3 — library design (`agp/design.py`)
10 fold-checked RNA candidates (33 nt, diverse GC-clamped stem-loop scaffolds, ViennaRNA
MFE −7 to −16), epitope-biased loops, + scramble decoy. Library in `examples/agp_lib.json`.

## Stage 4 — Boltz co-fold screen (on-target ORM1 + scramble decoy)
**Decoy (scramble × ORM1) ipTM = 0.862 — extremely high** (AGP's surface is a strong
electrostatic/sticky sink, even more than cTnI's 0.749). On-target ipTM, ranked:

| rank | cand | sequence | ORM1 ipTM | vs decoy 0.862 | fold |
|---|---|---|---|---|---|
| 1 | agp1 | GGCAGACGGAAUGGCUGCAGCGGAAUGUCUGCC | 0.845 | −0.02 | −14.2, 11bp |
| 2 | agp9 | GGCAGGAGGUCAGGUAGUGUAGAACAAUCUGCC | 0.840 | −0.02 | −7.0, 9bp |
| 3 | agp5 | GGCAGGGCAUGGAACCUGGAGUACGGUGCUGCC | 0.823 | −0.04 | −10.3, 9bp |
| 4 | agp3 | GCGGCGGACUAUCGGCACCCCUAUGUACGCCGC | 0.796 | −0.07 | −12.5, 12bp |
| 5 | agp6 | GCGGCUGGACGACUGGACUGUAUUGGGAGCCGC | 0.771 | −0.09 | −10.7, 8bp |
| 6 | agp10| CGGAGGGCAUGAUGGGAGGUAGUGUACACUCCG | (pending) | — | −8.3, 12bp |
| 7 | agp2 | GCGGCAUGGCAUGGCCGUCACGAAUUGUGCCGC | 0.660 | −0.20 | −11.9, 10bp |
| 8 | agp4 | GGGACAUGGCCUGGACUUGGCCCUUUACGUCCC | 0.656 | −0.21 | −12.6, 9bp |
| 9 | agp7 | GGGACCUGGAUGACGGGCAUUUAAUAUUGUCCC | 0.627 | −0.24 | −7.7, 11bp |
| 10| agp8 | GCGGCCAGGUGAGGUCGGUGUUUUCACGGCCGC | 0.567 | −0.30 | −15.8, 12bp |

## Verdict — honest, and consistent with the whole skill
**NO candidate beats its scramble decoy (max on-target 0.845 < decoy 0.862).** By the
calibrated gate, the entire library **FAILS** on-target — before even bringing in the ORM2
paralog and lysozyme counter-screens (those jobs were run; they can only lower a candidate
further, never rescue a sub-decoy one). The ORM1-vs-ORM2 and unrelated margins refine the
*order* but not the *verdict*.

**Why:** AGP presents a sticky, glycan-decorated surface; a structured RNA adheres about as
well scrambled as designed, so Boltz ipTM carries no specificity signal here. This is the same
electrostatic-sink failure mode seen with cTnI (decoy 0.749) and CTLA-4 (aptamerd6), only more
extreme. It is exactly the case the skill is built to *catch*, not paper over.

**Deliverable status:** `agp1` / `agp9` are the least-promiscuous **triage-ranked SELEX
starting sequences** (top of a design-diverse pool aimed at the ORM1-divergent epitope), NOT
predicted binders. Real AGP aptamers must come from wet-lab SELEX with (a) counter-selection vs
ORM2 + serum albumin, (b) attention to the glycan shield (select against glycoform or target
the deglycosylated core), and (c) the divergent 112–117 loop presented accessibly. The in-silico
pass here **narrowed the pool and flagged the target as electrostatically hard** — that is the
honest, useful output, not a ranked list of "hits."
