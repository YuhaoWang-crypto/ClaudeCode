---
name: mhc-assay-feasibility
description: >-
  Decide what an MHC/HLA wet-lab study actually needs to buy, before it is
  commissioned — and say honestly which parts computation cannot touch. Use when
  reviewing an RFP, CRO quote, or study design involving peptide-MHC assays; when
  asked whether an immunology validation can be done in silico; to pre-screen
  which peptide-MHC complexes are worth custom tetramer/multimer synthesis; to
  work out how many donors must be HLA-typed to find a required allele
  combination; to check a T-cell hybridoma or restriction-element system against
  its documented epitope; or to score peptides on specific class I (HLA-A/B/C) or
  mouse class II (H2-IA/IE) alleles. Sorts every requested item into
  replace / de-risk / cannot, and refuses to score molecules outside what the
  predictors represent.
---

# Feasibility triage for peptide-MHC assay programmes

Turns "can we do this in silico?" into an item-by-item answer with evidence
attached, and — where the answer is yes — into a calibrated pre-screen that
changes what gets commissioned.

For **immunogenicity risk assessment of a protein** (impurities, ligands, HCP,
biologics) use the sibling skill `hla-dr-immunogenicity` instead. This one is
about **de-risking an assay programme**, not scoring a molecule's risk.

## The only honest framing: three buckets, and the first is usually empty

Sort every line item of the study into one of these. Say which bucket each lands
in *before* showing any number, because the bucket is the answer.

| Bucket | Meaning | Typical members |
|---|---|---|
| **REPLACE** | Computation answers it outright; the assay is not needed for that answer | Almost always empty. Say so out loud — it is the most useful sentence you will write |
| **DE-RISK** | Cannot answer it, but changes what you commission, how much, or whether it is feasible at all | Binding pre-screens, reagent build lists, donor cohort feasibility, panel design |
| **CANNOT** | A physical measurement with no sequence-level surrogate, or a molecule outside what the models represent | Cytokine output, proliferation, precursor frequency, refolding yield, anything with non-natural residues |

A proposal that claims REPLACE for a functional readout is wrong, and a
feasibility review that does not name the CANNOT items is not a review.

## Hard boundaries — check these first, they are disqualifying

Run these checks before any prediction, because each one voids results rather
than degrading them.

**Non-natural residues void the prediction entirely.** D-amino acids,
N-methylation, β-amino acids, cyclisation, non-standard capping: every MHC
predictor in use is trained on L-peptides of the 20 canonical residues and has
no representation of stereochemistry. Submitting a D-substituted sequence
returns the score of its all-L twin, which is a **different molecule**. This is
not a caveat to note in a limitations section — it is a refusal to score. When a
study's core comparison is "all-L control vs D-substituted test", computation
cannot perform that comparison at all; it can only characterise the L arm.

**Class II has a length floor.** The IEDB class II endpoint refuses input below
**11 residues**, measured. An 8-mer or a bare 9-mer core cannot be submitted.
Padding it to clear the floor is not a workaround — see `reference/boundaries.md`
for the measured size of that artefact.

**Class II binding needs flanking residues.** The bound core is 9 residues but
the groove is open-ended and affinity depends on the peptide flanking residues
either side. A bare core is a different molecule from the same core in context,
and usually a much weaker binder. This is normally *why* a study includes a
truncated core as a control.

**Class I %Rank is not stability.** Predictors report binding, not off-rate or
complex half-life. Tetramer survival through staining is governed by stability.
When an RFP asks for "affinity/stability", computation addresses the first word
only.

**Precursor frequency is not predictable.** It is a property of one donor's
repertoire, set by their thymic selection and exposure history. No sequence model
reports it. Say this plainly rather than gesturing at "T-cell epitope prediction".

## Calibrate on the alleles in front of you, never on the average

Class I prediction is good *on average*, and the average is carried by
well-studied alleles. A*02:01 performance says nothing about A*32:01. Before a
pre-screen is allowed to say "do not build that complex":

1. Pull every labelled IEDB record for **those specific alleles** (`mhc_search`).
2. Score them and report ROC AUC per allele with a bootstrap CI.
3. Choose the threshold from the measured curve, not from convention.

**The class-imbalance trap.** IEDB class I data is overwhelmingly eluted ligands,
which are positive by construction — one allele in the worked example comes back
**96 % positive**. PPV and NPV computed at that prevalence describe IEDB's
collection policy, not your decision. Recompute them at a *stated* prior for
"a designed peptide proposed for this allele actually binds", and sweep it.

**The decision is asymmetric, so pick the threshold accordingly.** A complex
wrongly built costs money; a complex wrongly skipped is an epitope nobody goes
back for. So the governing quantity is **NPV at the skip threshold**, and the
rule is: take the *widest* cut that still keeps ~95 % of real binders on the
build side. Do not optimise MCC or balanced accuracy here — they price the two
errors equally and this decision does not.

## Donor feasibility is arithmetic, and it is often the finding

When a protocol names required HLA alleles, work out the cohort burden before
anything else — it is cheap, and it can invalidate the design.

```
carrier frequency  P = 1 - (1 - f)^2        # single-locus Hardy-Weinberg
donors to type     n / P                    # for n carriers, in expectation
```

Allele frequencies come from the IEDB population-coverage tables (see
`reference/boundaries.md` for extraction). **Two loci are not independent** —
HLA-A and -B are ~1.3 Mb apart and in strong linkage disequilibrium, so a
product of two carrier frequencies is a planning figure with a stated direction
of error, never a quotable number. Check the specific pair against haplotype
data before costing.

The question that usually needs asking back: **must each donor carry *both*
named alleles, or does either suffice?** In the worked example that distinction
moves the screening burden from ~89 donors to ~2,900 — the difference between a
feasible study and an impossible one, and the RFP did not say which it meant.

## Check the system's own control first

If a study rests on a restriction element and a reference epitope — a hybridoma,
a transgenic line, a reference tetramer — reproduce that pairing before anything
else. The DO11.10 / I-A(d) / cOVA323-339 system in the worked example returns the
documented `VHAAHAEIN` register at %Rank 0.17. If the control does not reproduce,
nothing downstream is worth reading, and that is a finding worth an email on its
own.

## Run it

```bash
pip install pyyaml numpy
python scripts/c1_benchmark.py            # labelled IEDB data for the named alleles
python scripts/c3_calibrate.py            # per-allele AUC + skip threshold; ~30-60 min
python scripts/c4_prescreen.py            # the build/skip list
python scripts/c5_donor_feasibility.py    # cohort burden; seconds, run it first
python scripts/c6_nanopeptide_classII.py  # system control + boundary demonstration
python scripts/c7_feasibility_matrix.py   # the item-by-item verdict
```

`c5` costs nothing and can change the study — run it before the long jobs.
`c3` is the long pole and is resumable. Sponsor sequences drop into
`data/candidate_peptides.tsv`; absent that, `c4` draws a stand-in panel from
benchmark peptides **held out of calibration**, so the demo makes real calls
against known answers rather than scoring peptides it was tuned on.

## Reference

- `reference/boundaries.md` — the measured limits: class II length floor, the
  padding artefact, non-natural residues, IEDB class I API mechanics, allele
  frequency extraction, and the class-imbalance handling.

## Writing the answer

Lead with the bucket counts and the sentence that nothing is replaceable. Then
the two or three items that *are* de-riskable, each with its number. Then the
questions the document left ambiguous — in the worked example, "both alleles or
either?" was worth more than any prediction. Close with what the wet-lab study
should still be, which is normally the same study, smaller and better targeted.

Never let a pre-screen's ranking be read as a grade. It says which complexes to
build first, not how immunogenic anything is.
