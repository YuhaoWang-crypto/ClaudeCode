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

**Partly reproducible, and more than it first appears.** The trial-level
numbers cannot be recomputed, because the proteome is controlled-access. But
the paper releases its own clock library with real weights, and its
supplementary workbook is openly downloadable and contains the per-sample
clock predictions. So a real, non-simulated piece of the paper is reproduced
here, and it changes how the headline claim reads.

| Requirement | Status |
|---|---|
| Supplementary tables S1-S15 | **Open.** Downloaded and analysed in M8 |
| `proteoclock`, the paper's library | **Open.** Real weights for 3 of 6 clocks |
| Trial serum proteome (OMIX008341) | Controlled access, request required |
| UK Biobank reference (n=55,319) | Separate application, months |
| PAC, OrganAge chrono + mortality | **Real weights**, shipped in proteoclock |
| PAOPAC | Weights downloadable, but a Windows-only .pyd for Python 3.9 |
| ipfP3GPT | Weights usable only inside the UK Biobank RAP |
| ProtAge trained model | Not released, authors must be emailed |

An earlier draft of this document said the released library "could not be
located". That was wrong. It is at
https://github.com/Insilico-org/proteoclock, under the `Insilico-org`
organisation, not the `insilicomedicine` one that carries the company's other
repositories. Searching the wrong organisation produced a confident negative.

Three of the six clocks therefore run here on their real published weights,
and M8 analyses real trial-derived data rather than simulation.

The adapter that bridges this pipeline to the library is verified, not
assumed. `proteoclock_backend.validate_against_package` runs GSE169148, the
31-sample real Olink dataset the package bundles, through both a direct
proteoclock call and this pipeline's wrapper, and requires the two to agree
bitwise. All three clocks return a maximum absolute difference of exactly
zero. The wrapper reshapes wide NPX into long form, subsets each clock's own
proteins and reindexes the output back to caller order, and any of those
steps could silently reorder or drop samples while still returning plausible
ages, so the check runs every time the backend module is executed.

One caveat about that dataset. The package also ships
`new_clock_res_GSE169148.tsv`, a table of expected results. It is NOT used as
the validation target, because it does not match what the current code
produces: its values correlate about 0.96 with the present API output but sit
on a different scale, roughly 2 against 63, which reads as a raw log hazard
against an age-converted one. It appears to predate the current code. A
correlation of 0.96 rather than 1.0 means it is not merely a rescaling of the
same numbers, so the shipped table cannot serve as a golden test of the
library as released. The live package is the correct comparison for an
adapter regardless.

## What the real data shows

Supplementary Table S2 carries 168 rows, 42 patients times 4 visits, with
predictions from all 43 clock variants. Each row is keyed by a SHA256 of its
own NPX values, so predictions cannot be tied to an arm or visit without the
controlled data. Per-arm statistics are therefore still out of reach.

What S2 does give is the joint distribution of the six clocks over the real
trial samples, and that settles a question simulation cannot.

**The six clocks are not six independent votes.** Mean pairwise correlation
is +0.622. PAC and OrganAge-mortality correlate at 0.940, ProtAge and PAOPAC
at 0.917. The eigenvalues of the correlation matrix are 4.16, 1.37, 0.25,
0.09, 0.09, 0.05: two components carry almost everything.

| Estimator | Effective independent clocks |
|---|---|
| Li and Ji 2005 | 3.00 of 6 |
| Nyholt 2004 | 3.79 of 6 |

That matters directly for the concordance argument. Six of six clocks
agreeing has a two-sided binomial p of 0.031 if they are independent. At
three effective clocks it is 0.25. The agreement of the six clocks is a
weaker piece of evidence than its face value.

**The paper's own tally is more measured than its coverage.** Supplementary
Table S5 reports 21 of 54 tests passed, 39%. Chronological clocks passed 13
of 36, mortality clocks 8 of 18, and the mortality clocks passed 0 of 18 in
the 60 mg QD arm.

**The enrichment is heterogeneous and includes depletion.** Supplementary
Table S8 splits clock-feature enrichment by response trajectory rather than
reporting one number:

| Arm | Trajectory | Odds ratio | p |
|---|---|---|---|
| 30 BID | Sustained | 4.65 | <0.001 |
| 30 BID | Delayed | 0.38 | <0.01 |
| 30 BID | Transient | 0.46 | 0.29 |
| 60 QD | Sustained | 2.19 | <0.05 |
| 60 QD | Delayed | 0.60 | 0.16 |

Sustained responders are strongly enriched; delayed ones are significantly
DEPLETED. The single 1.74-fold figure that appears in press coverage is not
what this table says, and no single odds ratio summarises it.

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
m8_supplementary.py  REAL data: the paper's own tables S2, S5, S8
proteoclock_backend.py  the paper's released library, real weights for 3
                     clocks, plus a bitwise adapter check on real GEO data
real_data.py      adapter for OMIX008341 and a real reference cohort
```

M1 through M7 run on simulation. M8 does not; it is real published data.

The simulation writes a specific answer into the data: 326 drug-responsive
proteins, 72 of them aging-associated, giving a true enrichment odds ratio of
1.730; the 30 mg BID arm carries the largest effect. The pipeline then has to
find those numbers back without being told. That tests the code. It says
nothing whatsoever about rentosertib.

## What the validation showed

Running it surfaced six things, five of them problems in the first draft.

**A namespace mismatch sent every clock silently to a surrogate.** The panel
symbols were read from proteoclock's `feature_order.txt` a line at a time.
The file is two tab-separated columns with the symbol repeated, so every
identifier came out as `A1BG\tA1BG`, matched nothing, and all three real
clocks fell back to surrogates while reporting success. This is exactly the
failure `real_data.py` warns about, and it happened here. The guard that
caught it was checking matched-protein counts against each clock's expected
total, which is why that count is now printed on every load.

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

**Real clocks cannot be validated on arbitrary synthetic biology.** Once the
published clocks were wired in they read the synthetic cohort as noise,
because its age slopes had been assigned at random over real gene symbols and
a published clock is keyed on which proteins actually move with age in
humans. M1 now derives its aging axis from a published clock's own
coefficients. That makes the recovery circular for the clocks sharing that
axis and is not evidence any clock works; it exists so the downstream
statistics operate on realistic covariance rather than noise.

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

The six-clock concordance argument has a related weakness, and M8 measures it
rather than asserting it: the clocks correlate at +0.622 on average over the
real trial samples, giving three effective independent clocks, so the
binomial p for six agreeing rises from 0.031 to 0.25.

None of this says the paper is wrong. It says these specific quantities are
imprecise at this sample size, and the pipeline quantifies by how much.

## Switching to real data

1. Obtain `OMIX008341` and a reference cohort, per `data/README.md`.
2. Install `proteoclock` for the three clocks whose real weights ship with
   it, and drop coefficient CSVs into `data/weights/` for any others.
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
- **Rigorous:** everything in M8. The cross-clock correlations, the effective
  clock count, and the S5 and S8 figures are computed from the paper's own
  published tables, not simulated.
- **Rigorous:** the weight-availability table, checked against each
  repository and against the paper's code-availability statement, September
  2026.
- **Hypothesis:** every biological number the demo prints. The drug effect,
  the dose ordering and the enrichment were injected by `m1_cohort` and
  measured back out.
- **Not established:** whether rentosertib has a geroprotective effect. The
  M8 findings bear on how strongly the paper's statistical argument supports
  its claim, not on whether the claim is true. Settling that needs the
  controlled data.
