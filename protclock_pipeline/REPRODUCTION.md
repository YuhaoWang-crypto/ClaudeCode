# Reproducing doi:10.1038/s41587-026-03286-y

Zhavoronkov, Galkin, Chen, Ren, Aliper, Durymanov, Sidorenko, Cui, Han, Xu,
Liu, Xu, Kuppe, Argentieri, Ying, Goeminne, Moqri, Tyshkovskiy, Gladyshev.
*Integration of proteomic aging clocks in a phase 2a clinical trial supports
simultaneous geroprotective assessment.* Nature Biotechnology, 7 Sep 2026.

Six proteomic aging clocks applied to serum from a 12-week phase 2a trial of
rentosertib, a TNIK inhibitor, in idiopathic pulmonary fibrosis. 42 patients,
four arms, four visits, Olink Explore 3072. All six clocks read the treated
arms as biologically younger. 326 circulating proteins altered, and those
proteins 1.74 times more likely than background to be aging-linked.

## Verdict

**The paper's numbers cannot be reproduced from public material.** The
obstacle is data access, not code. What this package does instead is
reimplement the analysis end to end and validate it against a synthetic
cohort with a known answer, so that the moment the controlled data arrives
the pipeline produces the real numbers with no further work.

| Requirement | Status |
|---|---|
| Trial serum proteome (OMIX008341) | Controlled access, request required |
| UK Biobank reference (n=55,319) | Separate application, months |
| ProtAge trained model | Not released, authors must be emailed |
| ipfP3GPT (IPF fine-tune) | Not released; only the base model is public |
| PAOPAC | Repository exists, no README, no weights confirmed |
| OrganAge chrono + mortality | Coefficients published in supplementary tables |
| PAC | Scoring script public |
| The paper's own pipeline library | **Not findable** |

That last row is worth stating plainly. The press release says the analysis
pipeline "has been released as an open-source Python library on GitHub". The
Insilico Medicine GitHub organisation carries 22 repositories and none of
them concerns aging clocks, proteomics, or biological age. As of September
2026 the library could not be located.

So three of the six clocks have obtainable weights, and neither of the two
datasets does.

## What this package is

A working reimplementation, validated against planted ground truth.

```
m1_cohort.py      synthetic Olink-shaped cohorts, known answer written in
m2_preprocess.py  LOD censoring, QC, imputation, plate centring, panel harmonisation
m3_clocks.py      the six clocks; published weights if supplied, surrogates otherwise
m4_ageaccel.py    age acceleration, all six clocks on one scale
m5_trialstats.py  within-patient change, mixed model, concordance, negative control, power
m6_enrichment.py  differential abundance, aging set, Fisher enrichment
m7_pathways.py    over-representation, mean shift, aging-aligned shift
real_data.py      adapter for OMIX008341 and a real reference cohort
```

The simulation writes a specific answer into the data: 326 drug-responsive
proteins, 72 of them aging-associated, giving a true enrichment odds ratio of
1.730; the 30 mg BID arm carries the largest effect. The pipeline then has to
find those numbers back without being told. That tests the code. It says
nothing whatsoever about rentosertib.

## What the validation showed

Running it surfaced five things, four of them problems in the first draft.

**The pipeline is calibrated.** With the drug effect set to zero, the
false-positive rate across replicates is 0.056 in the power sweep and 0.067
in the standalone negative control, against a nominal 0.05. Mean arm shift is
-0.042 SD and clock agreement sits at chance. A pipeline that reported a
result there would be manufacturing it.

**Replicates have to be independent.** The first power sweep showed the null
arm shifted -0.585 SD with 75% clock concordance, which is impossible under
no effect. The cause was a fixed seed inside the cohort simulators: every
"replicate" redrew the same patients and the same noise, so one patient-level
fluke repeated in all of them. Threading the seed through fixed it. Any
resampling analysis is exposed to this.

**Converting a mortality clock to years inflates its noise.** A
mortality-trained clock outputs a log hazard whose slope against age is
around 0.01 logit per year. Dividing the residual by that slope to express it
"in years" multiplies the noise by a hundred. On the years scale
OrganAge_mortality has a within-patient standard deviation near 20 years
against 2 to 5 for the chronological clocks. Standardising by the reference
residual SD instead removes the problem, so all inference runs on that scale
and years are kept only for interpretation.

**A plain mean shift is the wrong pathway test.** Proteins in an aging set
move in both directions with age, so a drug that rejuvenates all of them
produces shifts that cancel. Senescence came out at Cohen's d = -0.07,
p = 0.91 despite being built entirely from aging proteins. Projecting each
effect onto its own aging direction recovers it at q = 0.0001, while the
immune-activation negative control stays null.

**Thresholded enrichment is biased upward.** The odds ratio depends on which
contrast is used and inherits a selection effect, because proteins with
larger effects are easier to detect and in this simulation those are
disproportionately the aging-associated ones. From identical data:

| Contrast | Proteins altered | Odds ratio | 95% CI |
|---|---|---|---|
| Treated vs placebo | 37 | 1.22 | 0.53 to 2.79 |
| Paired within treated | 178 | 2.22 | 1.58 to 3.13 |
| Planted truth | 326 | 1.730 | |

Both intervals contain the truth. Neither is precise. The treated-vs-placebo
contrast recovers 11% of the responsive proteins because it spends its power
on an 11-patient control group; the paired contrast recovers 53% but cannot
separate a drug effect from disease progression or assay drift. An order of
magnitude separates the two protein counts, from one dataset, on an analyst's
choice.

## The power problem

The clearest result is about sample size. Six replicates per effect size:

| Planted effect (yr) | Tests at p<0.05 | Replicates with 5+/6 clocks agreeing | Strongest arm identified |
|---|---|---|---|
| 0 (null) | 0.056 | 0.33 | 0.17 |
| 15 | 0.241 | 0.67 | 0.17 |
| 30 | 0.519 | 0.83 | 0.50 |
| 60 | 0.759 | 1.00 | 0.50 |
| 120 | 0.833 | 1.00 | 0.50 |

Clock concordance saturates at 1.00 while the dose ranking never rises above
0.50, even at four times the planted effect. The pipeline detects that
something happened long before it can say which dose did most, and more data
per arm, not a larger effect, is what would close that gap.

This bears directly on how the paper's dose ranking should be read. A
42-patient trial split four ways leaves 10 to 11 patients per arm. Reading a
dose-response ordering off that is not well supported, whatever the analysis.
The six-clock concordance argument has a related weakness: the clocks are
trained on overlapping proteins from overlapping cohorts, so they are not
independent, and a binomial p-value over six agreeing clocks is
anticonservative by an unknown factor. `m5_trialstats.concordance` reports it
and labels the column `binom_p_ANTICONSERVATIVE`.

None of this says the paper is wrong. It says these specific quantities are
imprecise at this sample size, and the pipeline quantifies by how much.

## Switching to real data

1. Obtain `OMIX008341` and a reference cohort, per `data/README.md`.
2. Drop clock coefficient CSVs into `data/weights/`.
3. Build the input through `real_data.build_preprocessed` instead of
   `m2_preprocess.run`, then run M3 through M7 unchanged.

`real_data.py` has never been executed against the real spreadsheet, because
it cannot be. Its column names are inferred from the accession metadata and
the Olink export format. Expect to adjust it on first contact, and check the
coefficient match count against each clock's expected protein count before
trusting any output.

Once real data is in, everything that compares against planted truth becomes
meaningless and must be skipped: M1's truth table, M6's recovery numbers,
M5's power sweep. The measurements themselves remain valid.

## Honesty labels

- **Rigorous:** the preprocessing, calibration, statistics and enrichment
  code, and every claim above about the pipeline's own behaviour, all
  computed and reproducible with `python3 -m protclock_pipeline.run_all`.
- **Rigorous:** the weight-availability table, checked against each
  repository in September 2026.
- **Hypothesis:** every biological number the demo prints. The drug effect,
  the dose ordering and the enrichment were injected by `m1_cohort` and
  measured back out.
- **Not established:** anything about rentosertib. This package contains no
  evidence for or against the paper's conclusions.
