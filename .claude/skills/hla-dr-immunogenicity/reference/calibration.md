# Calibrating the decision rule against measured T-cell outcomes

Everything upstream of this layer is a *claim* about accuracy. A rule like
"require the affinity head to agree, because eluted-ligand scoring over-calls"
is a reasonable argument and no evidence at all. Without ground truth there is
no way to know whether it raises specificity, costs more sensitivity than it
gains, or does nothing.

## Building the benchmark (M10)

Pull every HLA-DR-restricted T-cell assay result IEDB holds for the panel's
molecules, in human hosts, and label each **(peptide, allele)** pair:

- **1** — at least one Positive T-cell assay for that peptide on that molecule
- **0** — only Negative assays for that pair
- **excluded** — pairs with both; genuinely contested, counted and reported

This yields ~9,600 labelled pairs across ~5,800 peptides for a 25-molecule
panel.

## The three biases that manufacture accuracy that does not exist

**Redundancy.** Overlapping peptides from one protein and one study are not
independent observations. Cluster peptides sharing any 9-mer (union-find over a
9-mer index) and report per-cluster as well as per-pair. Every bootstrap
resamples **clusters**, not rows — a paired cluster bootstrap on the AUC
*difference* between two rules, since both rules see the same peptides.

**Source-organism confounding.** Positives skew toward pathogens and allergens,
negatives toward self. A predictor that merely recognised "foreign" would score
well for the wrong reason. Record the self/non-self split of each label and
report metrics stratified. (Here the best rule holds up in both strata, so the
discrimination is not an artefact of foreignness — but that has to be shown,
not assumed.)

**Assay-context loss.** IEDB negatives include peptides that were never going
to respond in that donor for reasons unrelated to binding. This inflates
apparent specificity. Never present the negative set as a clean non-binder set.

## Ascertainment bias — the one that cannot be corrected away

**Peptides are tested because someone already thought they were interesting**,
and self peptides enter IEDB largely through autoimmunity studies. Two
consequences, both measured here:

- Test count is itself a strong classifier: **"how many DR molecules IEDB
  tested this peptide on" reaches AUC 0.575**, higher than every
  sequence-derived predictor in this benchmark, with the positive rate climbing
  45.8 % → 80.0 % → 95.5 % across 1, 2 and 3 molecules tested.
- Any metric that counts molecules (breadth, population-weighted presentation)
  inherits this. **The only honest comparison is inside a single test-count
  stratum**, which is why M13's headline is the 1,069-peptide single-allele
  layer rather than the pooled number.

When a bias cannot be corrected, say what it prevents you from concluding. M12
could not bound the tolerance down-weight at all for this reason, so it reports
a sensitivity sweep and an explicit statement that the test lacks power —
rather than an optimum it cannot support.

## Leakage: report absolute numbers as an upper bound

NetMHCIIpan is trained on IEDB binding and mass-spec data. The labels here are
T-cell assay outcomes — a different endpoint — but many of these peptides also
carry binding measurements in IEDB, so **partial training-set overlap is
certain and cannot be excluded from outside**.

The comparison *between rules* is far more robust: every rule uses the same
possibly-leaky predictors on the same peptides, so leakage inflates both arms
and largely cancels in the difference. Frame conclusions as rule-vs-rule
wherever possible.

## Report the operating point, not just the AUC

An AUC does not tell a wet-lab scientist anything actionable. Report, at the
threshold actually in use: sensitivity, specificity, MCC, and **the PPV a flag
carries at a realistic scan prevalence** (5 % here). At this pipeline's
operating point that is sensitivity 0.16, specificity 0.96, and a flagged
peptide is real ~16 % of the time.

Note also that MCC is the wrong objective for this application: the
MCC-maximising threshold flags ~3× as many peptides and drops flag PPV. When
the deliverable is a shortlist to send into an assay, **precision is what is
being bought**. Choose the objective deliberately and say which one you chose.

## The four rejected improvements — worked examples

Each of these was proposed on a sensible argument and killed by measurement.
They are the reason the discipline exists.

### 1. The EL×BA consensus gate (M11)

Argument: EL over-calls, so require BA agreement. Measured against 5,795
labelled outcomes it removed **23 true positives and 6 false positives** —
about four real epitopes per false one — and MCC went 0.188 → **0.184**. The
AUC difference between the consensus score and EL alone was +0.013, 95 % CI
[−0.001, +0.025]: not distinguishable from zero.

Gate switched off. **This moved the headline from 2.92× the Protein A benchmark
to 1.05×** — the gate had been suppressing the benchmark's epitopes harder than
the test article's, so an unvalidated rule was setting the published answer.

### 2. "3× more promiscuous than the universal epitopes" (M6)

A bug, not a finding. Control epitopes were matched to scanned 15-mers by
**string equality**, so the 13-residue HA306-318 matched nothing and scored
breadth 0. The fix is position-overlap matching with a minimum 9-residue
overlap:

```python
MIN_OVERLAP = 9
if min(e, w1) - max(s, w0) + 1 < MIN_OVERLAP:
    continue
```

With overlap matching the universal epitopes reach 2–6/25 and the test article
6/25 — **comparable to, not above**. Any control shorter than the scan window
will silently score zero under equality matching; check this first whenever a
control comes back at exactly zero.

### 3. Bounding the tolerance down-weight (M12)

The sweep is flat and the test has no power, because of the ascertainment bias
above. Reported as a sensitivity sweep with a `why` field explaining the
flatness, plus an explicit three-branch recommendation that includes "the data
cannot decide this". Do not manufacture an optimum from a flat sweep.

### 4. Panel-wide breadth vs best single-allele rank (M13)

Argument: a peptide presented by many molecules is the more credible epitope,
so requiring breadth should raise specificity. This is the last obvious
specificity lever, and it is not there.

| vs best single-allele %Rank | ΔAUC | 95 % CI |
|---|---|---|
| breadth (molecules at %Rank < 1) | +0.0005 | [−0.0156, +0.0169] |
| population-weighted presentation | +0.0013 | [−0.0197, +0.0221] |
| breadth, single-allele stratum (n=1069) | +0.0054 | [−0.0113, +0.0221] |
| breadth, multi-allele stratum (n=131) | **−0.0877** | [−0.1599, −0.0145] |

A breadth ≥ 2 gate lifts specificity 0.61 → 0.79 but drops sensitivity
0.47 → 0.28 with MCC essentially flat (0.074 → 0.086): it buys specificity with
sensitivity at par and adds no information.

**Consequence for the pipeline: nothing changes, deliberately.** Flagging stays
per molecule on %Rank; population weighting stays an aggregation for reporting,
never a criterion for calling. What the result closes off is the tempting
addition of a "presented by ≥ N molecules" or "≥ 20 % of the population" gate.

## Validating the tolerance filter itself (M4)

The filter must be validated in the same run that uses it, against a
**shuffled-sequence null**. Real cores hit the human proteome at 24.8 %; the
null hits 0.9 % at the same cut — 27× enrichment, so the cut carries
information. If the null rate is not far below the real rate, the run should
report the filter as uninformative rather than apply it.

A **5-of-9 "TCR-face" screen was tested and rejected**: five specified positions
match the human proteome by chance ~3.5 times per query, so it flags nearly
every core and carries no information. Whole-9-mer best identity (9/9 exact,
8/9 near) is used instead. Do not report statistics the index cannot actually
compute — a one-mismatch masked index cannot give you 7/9 or 6/9 counts.
