---
name: epitope-gated-binder-triage
description: >-
  Run and triage a de novo binder campaign whose goal is FUNCTIONAL BLOCKADE,
  not just binding — nanobody/VHH, minibinder or antibody designs against a
  structurally defined epitope, using Boltz-2 (or any co-folding design API).
  The central finding this skill encodes: co-folding confidence scores
  (binding_confidence, ipTM, PAE, ipSAE) tell you whether a complex is real,
  NOT whether it sits on your blocking epitope — so selection must gate on
  computed interface geometry and then rank on diversity, never on score
  alone. Use this whenever someone wants to design a blocker/antagonist/
  neutralizing binder, asks how many designs to order, asks whether a
  co-folding score can be trusted for ranking, wants epitope occupancy or
  paratope footprints computed from predicted complexes, is planning
  predict-and-redesign iteration rounds, or is about to send designed protein
  sequences for synthesis. Also use it when a campaign "looks great by score"
  and you need to check whether that means anything. Enforces per-file
  structure verification, honest stratified statistics, and a final
  developability QC that catches fixes which introduce new defects.
---

# Epitope-gated binder triage

A campaign discipline for designing binders that must **block**, not merely
bind. It exists because the obvious pipeline — generate designs, rank by the
design API's confidence score, order the top N — quietly fails at the step
that matters, and fails in a way that produces confident, well-formatted,
wrong answers.

## The one result that drives everything else

Measured on 1,600 de novo VHH designs against a single target, across two
generation modes (undirected and parent-seeded CDR3 redesign), with every
design's epitope occupancy computed deterministically from its predicted
complex:

| score | Spearman ρ vs epitope recall, informative range |
|---|---|
| `binding_confidence` | +0.067 (p=0.27) round 1; −0.152 round 2 |
| `ipTM` | null |
| `min_interaction_pae` | null |
| `ipSAE` (all variants) | +0.078 (p=0.18) round 1; +0.071 (p=0.28) round 2 |

Every one of these shows a healthy-looking positive correlation **over the full
score range** (+0.32 to +0.48) that vanishes once you restrict to designs that
are actually docked. The full-range correlation is an artifact: the bottom of
the distribution is designs with essentially no interface, so "low score → low
recall" is true but vacuous.

The mechanism is worth understanding rather than memorising. Confidence scores
answer *"is this complex real?"*. They correlate strongly with how much surface
is buried (ρ ≈ +0.28 to +0.44 with contacting atom pairs) and not at all with
*where* it is buried. ipSAE specifically corrects ipTM's length confounding —
but in a VHH library every binder is 112–132 aa (sd ≈ 5), so there is no length
artifact left to correct and ipSAE degenerates to ipTM (ρ(ipSAE, ipTM) = +0.75).

**Consequence for the pipeline:** epitope occupancy is something you *compute*,
not something you *infer from a score*. Geometry gates; score does not.

A second, non-obvious consequence: in this campaign 69% of designs meeting the
blocker criterion had `binding_confidence` ≤ 0.1. Using score to decide *which
structures are worth analysing* discards most of the hits. Score every design.

## Pipeline

Seven phases. Phases 0 and 1 are gates — they can end the project cheaply, which
is their whole value. Do not skip them because the compute is fun.

### Phase 0 — Does this need to exist?

Before any compute, establish what an existing validated reagent already does.
A campaign that duplicates an off-the-shelf tool is a year lost, and the
reviewer will ask. Search for: species-matched blocking antibodies, ligand
traps, conditional-knockout lines, and prior art (patents, not just papers —
single-domain binders are heavily claimed).

Write one sentence: *"Existing reagents X, Y, Z demonstrably cannot do ___."*
If that sentence will not write, stop.

Also grade every reagent claim by evidence source — peer-reviewed vs patent vs
vendor datasheet. Vendor cross-reactivity claims for receptor ectodomains are
frequently unvalidated, and a vendor datasheet that itself reports <5%
cross-species reactivity is a finding, not a footnote.

### Phase 1 — Define and verify the epitope

Derive the epitope from an experimental complex where one exists
(`references/epitope-occupancy.md` has the contact computation). Then verify
three things before designing anything:

1. **Numbering.** Map the PDB's numbering onto your target's canonical numbering
   by aligning sequences, not by assuming an offset. Record the offset you
   derived and the residue identities at the epitope positions.
2. **The model you will design against.** Superpose the predicted/AlphaFold
   model onto the experimental structure. Sub-Ångström agreement at the epitope
   is what licenses designing against the model; report the number.
3. **Confidence and crop.** Exclude low-pLDDT termini — designing against
   unreliable coordinates produces designs that bind noise.

Then check the epitope's *neighbourhood*: glycosylation sequons within ~10 Å of
an epitope residue are invisible in an aglycosylated model but present on the
real receptor. This determines how the screening antigen must be expressed
(mammalian, glycosylated) and is a classic silent failure — a binder selected
against bare protein that cannot reach the glycosylated target.

### Phase 2 — Cost the campaign before running it

Ground the order size in experimentally-labelled data rather than intuition.
Published labelled binder sets show hit rates that vary far more by target
(0–80%) than by generator or scoring threshold, and pooled rates flatter the
planning estimate. Plan at the pessimistic end.

Separate two numbers that are easy to conflate: how many designs to **generate**
(cheap — generate widely) and how many to **order** (expensive — this is where
saturation applies). Generating 100 and ordering 96 is a selection ratio of
1.04×, which throws away the only computational lever you have.

### Phase 3 — Generate

See `references/design-spec-recipes.md` for the Boltz-2 spec patterns:
epitope-constrained targets, framework locking, developability rules, and
anchored CDR3 maturation for iteration rounds.

Two spec-level choices pay for themselves:
- **Lock the framework** by making it fixed sequence rather than designed. In
  one campaign this took framework anomalies from 11–14% of gate-passers to
  0/600. A spec-level fix beats post-hoc screening because it cannot be forgotten.
- **Exclude developability liabilities at design time** (glycosylation sequons,
  free cysteines, oxidation- and isomerisation-prone motifs) — cheap now,
  expensive after you have a lead.

### Phase 4 — Score every design's epitope occupancy

Compute, per design, from the predicted complex: contacted target residues,
epitope recall, purity, atom-weighted purity, hotspot coverage, and whether
off-epitope contacts are rim (adjacent to the epitope) or a genuinely separate
patch. `references/epitope-occupancy.md` has the definitions and the per-file
verification that keeps this honest.

Define the blocker criterion explicitly and justify each term from the
distribution you actually observe, not from an imported threshold. Absolute
cutoffs borrowed from another target are not portable; within-pool percentiles are.

### Phase 5 — Iterate, and test whether iteration worked

Iteration is the one lever that published data supports, but the headline
"blocker rate went up" is usually not the number that matters. Check:

- the **joint bar that fills a plate** (blocker AND whatever quality threshold
  you'd actually order at) per design, round over round;
- whether the **ceiling** moved — did any child beat the best parent;
- **per-lineage** rates. In one campaign the entire gain came from 1 of 11
  lineages (56.5% plate-grade); the other ten returned 1.1%, *worse than the
  undirected baseline*. A pooled improvement can hide ten lineages being made
  worse.

State a falsifiable prediction before the round and check it after. Ours was
"bc-driven optimisation will push designs off the epitope"; it was falsified
(2.5% showed the signature, sign test p=0.20) and saying so plainly was more
useful than the prediction would have been if it had held.

### Phase 6 — Select the plate, then QC it

When no available score discriminates within the gate, ranking by score is
selection on noise. The defensible objective is the **number of independent
shots on goal**: cluster the gate survivors on paratope-relevant features and
take one representative per cluster, with a quality floor so that "diverse"
does not become "worse". `references/plate-selection.md` has the method and the
comparison that justifies it.

Then run a **full QC of the final sequences as a batch** — see
`references/failure-modes.md`. The reason this is a separate step and not a
per-edit check: fixing one liability can create another. In this campaign,
reverting a free cysteine C→S turned `N-T-C` into `N-T-S`, a new N-glycosylation
sequon next to the paratope. It was caught only by re-running the whole QC on
the final sequences.

## Non-negotiable verification gates

These are the places where being wrong is silent — the pipeline runs, the
numbers format nicely, and they mean nothing.

**Verify structure numbering per file, not once.** Assert that the target chain
in every predicted complex has the expected residue range *and* the expected
sequence. A numbering shift makes every occupancy number confidently wrong. Two
campaigns, 1,600 files, zero failures — the assertion cost nothing and the
alternative was unfalsifiable output.

**Relocate reported positions in the actual sequence.** A report saying "free
cysteine at position 32" may be off by one. Find the anomaly by pattern in the
sequence itself before editing; the edit is then correct regardless of whether
the report was.

**Stratify every correlation.** A full-range correlation on a score whose low
end is undocked designs is an artifact. Report the stratified number as the
result and the full-range number as the artifact it is. If you report a positive
correlation, show it survives restriction to the informative range.

**Run a control before trusting a filter.** Before using a score to decide what
to analyse, score a random sample of what the filter would discard. Our
"bc > 0.1 is a safe prefilter" assumption was falsified this way (21.7% blocker
rate in the discarded stratum, not ~0).

**Check edited positions against the paratope.** Framework reversions are safe
because they are outside the interface — verify that per design rather than
assuming it.

## Honest reporting

Two distinctions carry most of the weight:

- **Computational confidence is not a hit rate.** A statement like "10.5% of
  designs scored above 0.5" and a published "26.8% experimental hit rate" are
  different quantities that happen to share the interval [0,1]. Keep planning
  estimates anchored to experimental data.
- **Within-campaign null is not universal null.** Our scores were null for
  epitope location on one target, one modality, at this n. Say that, rather
  than that the scores are worthless.

Label assumptions as assumptions. The expected-independent-modes calculation
that justifies a diversity plate depends on an independence assumption that was
not measured; it is still the right decision procedure, and it should say so.

## Reference files

- `references/epitope-occupancy.md` — contact and occupancy metric definitions,
  the blocker criterion, per-file verification, CIF parsing conventions.
- `references/scoring-is-not-location.md` — the full empirical result, the
  ipSAE formula and its validation, and how to re-run this test on a new campaign.
- `references/design-spec-recipes.md` — Boltz-2 spec patterns: epitope
  constraints, framework locking, developability rules, anchored CDR3 maturation.
- `references/plate-selection.md` — diversity-first selection, the comparison
  against score ranking, and the expected-modes argument.
- `references/failure-modes.md` — every trap encountered, with the symptom, the
  cause, and the check that catches it. Read this before any long run.
