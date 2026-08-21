# Benchmark: is the ranking better than NetMHCpan alone?

A neoantigen pipeline that never gets scored against ground truth is a
confident-looking sort. This is the only ground truth obtainable from open data
without a data-access agreement, and this document is explicit about what it is
and is not.

## Positives — mined, not hand-picked

`benchmark.build_positive_set()`:

1. Query the IEDB `query-api` for T-cell assay records with
   `qualitative_measure != Negative`, `host = Homo sapiens`, `mhc_class = I`,
   `source_organism = Homo sapiens`, linear peptides.
2. Keep those with a four-digit HLA-A/B/C restriction and a length the pipeline
   predicts (9, 10 by default).
3. Keep only peptides that are **exactly one substitution** away from some
   k-mer in the UniProt reviewed human proteome, and that do **not** themselves
   occur in the proteome.

Step 3 is the whole trick. A peptide that is one residue off a human protein
and is not itself human is mutation-shaped — the signature of a somatic
missense neoepitope — and the search hands back the wild-type counterpart, so
`agretopicity` is computed *exactly* as the main pipeline computes it, not
approximated.

The search is a two-sided prefix/suffix bucket index over all self k-mers,
so a substitution in either half of the peptide is still found. Ambiguous cases
(more than one distinct one-mismatch self k-mer) keep the first match.

**What this is not:** a curated neoantigen database. Some mined positives will
be post-translationally modified self peptides, minor histocompatibility
antigens, or curation artifacts rather than tumor neoantigens. The rule is
mechanical and stated, which is better than a hand-assembled list that is
selected on the same intuitions the score encodes.

## Negatives — real, unlabelled, and run two ways

`benchmark.build_decoy_set()` takes real missense mutations from *other*
TCGA-SKCM patients (the demo patient is excluded) and tiles them through the
same `peptides.py` code path the pipeline uses.

These decoys are **unlabelled, not verified non-immunogenic.** Some fraction
are presumably immunogenic and were simply never tested. That biases every AUC
*downward*. Reported numbers are therefore lower bounds on separability and are
labelled `AUC_lower_bound` in the output table — not decoration, the column is
literally named that so the number cannot be quoted without the caveat.

### Benchmark A — matched on allele and length only

Decoys sampled to match each positive's (allele, length) at a 5:1 ratio.

**This benchmark is a trap, and the demo report says so.** IEDB epitopes were
largely *discovered* because they bind well — many were found by running a
binding predictor in the first place. Random mutant peptides mostly do not
bind. So the two classes differ chiefly in predicted binding, and a binding
predictor separates them almost perfectly without saying anything about
immunogenicity. On the demo run, NetMHCpan %rank alone reaches AUC ≈ 0.97
here. Anyone reporting that number as evidence that their pipeline works has
measured NetMHCpan, not their pipeline.

### Benchmark B — presentation-controlled

The question the selection layer actually exists to answer is: *among peptides
that are all presented, which does a T cell see?* So benchmark B removes
binding as an explanatory variable:

1. Restrict to alleles with at least `--bench-min-allele-n` validated positives
   (default 10) — the rest have too few positives for any number to mean
   anything, and dropping them is reported, not silent.
2. Draw decoys from a much larger pool (`--bench-pool-ratio`, default 250 per
   positive), because strong binders are rare among random peptides: by the
   definition of a percentile rank, only ~0.5% of random peptides fall below a
   0.5% rank, so a small pool simply contains no decoy that binds as well as
   the positives do.
3. `balance_by_stratum()` bins everything by predicted %rank and subsamples
   decoys so each stratum carries the same class mixture. After this the two
   classes have approximately the same binding distribution, so an overall AUC
   above 0.5 cannot be explained by presentation alone.
4. `stratified_evaluate()` additionally reports the AUC of each score *within*
   each stratum, with the per-class counts next to it, and NaN wherever a
   stratum has fewer than 8 of either class — a number computed from four
   positives is not a result.

Benchmark B is harder and the numbers are lower. That is the point.

## Leakage control

`tcr_prior` is computed by aligning to IEDB positives, and the benchmark's
positives come from IEDB. That is circular. Two mitigations:

1. The query peptide and exact duplicates are removed from the reference set.
2. More importantly, the composite is reported **three ways**:
   `score_netmhcpan_only`, `score_composite_no_tcr`, `score_composite`.

Read `score_composite_no_tcr` vs `score_netmhcpan_only`. That comparison is
leak-free and answers the actual question: does the selection layer add
anything over the binding predictor?

## What is deliberately excluded

`expression` and `clonality` are set to zero for every benchmark row. An IEDB
epitope has no tumor RNA value and no CCF. Assigning positives a plausible TPM
and decoys a random one would manufacture a large AUC out of nothing. So the
benchmark evaluates only the **peptide-intrinsic half** of the score; the
expression/clonality half is defensible on first principles (an unexpressed
gene makes no protein) but is not evaluated here, and the report says so.

## Metrics

| metric | meaning |
|---|---|
| `AUC_lower_bound` | Mann-Whitney AUC with tie-corrected ranks, biased down by unlabelled decoys |
| `precision@34` | fraction of true positives in a 34-slot payload drawn from the pooled candidates — the operationally relevant number, since 34 is the real budget |
| `fold_enrichment` | `precision@34` / base positive rate |
| per-stratum AUC | the same AUC computed inside one predicted-%rank band, with `n_pos` / `n_dec` shown so a small-sample number can be discounted on sight |

## What the mined positives actually are

Worth looking at before quoting any number. The mining rule is mechanical
("experimentally T-cell-positive, human source, exactly one substitution from
self"), and what it returns is a mixture:

* genuine tumor neoantigens (e.g. a melanoma-associated CHSA/`Q9Y5P2` variant),
* **post-translationally modified** self peptides — the classic tyrosinase
  369-377 `YMDGTMSQV` is a deamidated `YMNGTMSQV`, one residue from self and
  genuinely T-cell-reactive, but not a somatic mutation,
* **minor histocompatibility antigens** — allelic variation between donor and
  recipient, e.g. HA-1-like peptides, T-cell-reactive in transplant settings,
* HLA and other polymorphic-locus peptides, which are population variation
  rather than tumor-specific change.

All of them are "one substitution from self and demonstrably immunogenic",
which is the property the score is meant to detect, so they are not noise. But
they are not all neoantigens, and the demo report says which is which where the
antigen name makes it obvious.

## Running it

```bash
python -m neoantigen_pipeline.run_demo --benchmark --out demo_out
```

Writes `benchmark_scored.csv` (every row with its features) and
`benchmark_metrics.csv` (the table above), and appends a labelled section to
`REPORT.md`.

## Refitting the weights

`benchmark_scored.csv` has `label` and all `feat_*` columns. Fit whatever you
like on it — logistic regression is the honest default given the sample size —
and write the coefficients back into `config.Weights`. If you do, say in the
report that the weights are fitted and on what, because the literature defaults
and a fitted set are different claims.

## What the demo run actually found

Run on 2026-08-21: 119 mined positives, TCGA-SKCM decoys, NetMHCpan-4.1 EL
through the IEDB cloud API.

| benchmark | NetMHCpan %rank alone | composite (no TCR prior) | composite |
|---|---|---|---|
| A — allele/length-matched decoys | **0.966** | 0.935 | 0.939 |
| B — presentation-controlled (86 positives, 348 decoys, 3 alleles) | **0.599** | 0.578 | 0.576 |

Within each binding stratum of benchmark B:

| %rank stratum | n pos | n dec | NetMHCpan alone | composite (no TCR) |
|---|---|---|---|---|
| 0.0–0.1% | 33 | 83 | 0.570 | 0.476 |
| 0.1–0.5% | 18 | 90 | 0.563 | 0.523 |
| 0.5–2.0% | 22 | 110 | 0.584 | 0.554 |
| 2.0–100% | 13 | 65 | 0.860 | 0.697 |

**Read this honestly.** Three things follow, and none of them is "the composite
score is validated":

1. Benchmark A's 0.966 is an artifact of how positives were discovered. The
   drop to 0.599 when binding is controlled is the size of that artifact.
2. Once presentation is equalized, **nothing here separates validated
   immunogenic neoepitopes from unlabelled decoys much better than chance** on
   this data — not agretopicity, not hydrophobicity, not dissimilarity, and not
   the weighted composite. The 0.860 in the weakest-binding stratum is residual
   binding signal inside a very wide bin, not a fourth feature working.
3. The composite scores *slightly below* NetMHCpan alone in benchmark B. On
   this evidence the extra peptide-intrinsic features are not earning their
   weight, and the literature defaults should not be presented as tuned.

Caveats that cut in both directions: decoys are unlabelled (biases AUC down),
n = 86 positives is small, and the mined positives are a mixture that includes
PTM variants and minor histocompatibility antigens as well as true neoantigens.

What this does **not** measure, and what is still doing real work in the
pipeline: the gates (expression, tumor-specificity vs the self proteome,
clonality) and the payload constraints (allele spread, gene caps, junction
control). None of those can be evaluated with a peptide-intrinsic benchmark,
and all of them change which 34 mutations get made.

If you want the composite to beat the binding predictor, the honest routes are
labelled outcome data to fit the weights on, or better features — not a
benchmark whose decoys make the problem look easier than it is.

---

# The TESLA mirror: real labels, real negatives (`tesla.py`)

Everything above works around one weakness — the decoys are *unlabelled*, so
every AUC is a lower bound and the true negative rate is unknown. The public
TESLA mirror removes that weakness where it can be removed.

`data/tesla_deepimmuno_public.csv` holds **522 peptide-HLA pairs, 35
experimentally immunogenic, across 6 patients**, with the per-model scores
published alongside the DeepImmuno evaluation. A label of 0 here means
**assayed against patient T cells and found negative** — not "never tested". It
is a processed public mirror, not the full 608-pMHC supplemental table, and
every number inherits that scope.

Two properties make it the right dataset to argue over:

* **6.7% base rate.** Accuracy is meaningless; average precision (PR-AUC) and
  top-N recovery are what move. `tesla.evaluate()` reports AP, AP relative to
  baseline, AUC, and top-20 / top-34 recovery.
* **It is already presentation-filtered.** 100% of positives and 91% of
  negatives sit at %rank ≤ 2; at ≤ 0.5 it is 97% vs 68%. The negatives are
  mostly predicted binders that were tested and came back negative — which is
  exactly the presentation-controlled comparison benchmark B tries to construct
  synthetically, here for free.

Wild-type counterparts are recovered the same way `build_positive_set` mines
its positives — the self k-mer one substitution away — so agretopicity is
computed identically to the main pipeline. On this mirror 487 of 522 peptides
have one; the rest keep the neutral value and are counted, not dropped.

## Results (2026-08-21, NetMHCpan-4.1 EL via the IEDB cloud API)

Random baseline AP = **0.067**.

| score | AP | ×baseline | AUC | positives in top-34 / patient |
|---|---|---|---|---|
| **NetMHCpan-4.1 EL %rank alone** | **0.207** | 3.08 | 0.791 | **31 / 35** |
| this package's composite | 0.149 | 2.22 | 0.729 | 24 / 35 |
| `cnn_regress` (best published column) | 0.132 | 1.96 | 0.654 | 19 / 35 |
| composite without the TCR prior | 0.127 | 1.90 | 0.738 | 24 / 35 |
| `rf_regress` | 0.108 | 1.62 | 0.619 | 19 / 35 |
| `cnn_classify` | 0.108 | 1.61 | 0.550 | 14 / 35 |
| DeepImmuno `immunogenic score` | 0.083 | 1.23 | 0.477 | 13 / 35 |
| `IEDB` immunogenicity score | 0.070 | 1.04 | 0.523 | 17 / 35 |
| agretopicity alone | 0.066 | 0.99 | 0.449 | 10 / 35 |
| dissimilarity alone | 0.065 | 0.97 | 0.494 | 13 / 35 |
| hydrophobicity alone | 0.057 | 0.84 | 0.405 | 15 / 35 |

Three things follow.

**1. A current presentation predictor beats every published model column on
this mirror.** NetMHCpan-4.1 EL %rank alone: AP 0.207 against a best published
column of 0.132, and 31 of 35 experimentally confirmed positives inside a
34-slot budget against 19. Note what is being compared: NetMHCpan re-run today
versus model scores *as published* in the mirror. It is a fair comparison of
what a designer would actually use, not a re-training of those models.

**2. The composite scored below plain presentation, again.** Same direction as
the presentation-controlled IEDB benchmark, now on real labels with real
negatives. Two independent datasets agreeing is enough to act on, so the
defaults changed: presentation 0.30 → 0.45, and agretopicity, dissimilarity,
TCR prior and hydrophobicity all cut. Weight settings compared per patient:

| weights | pooled AP | mean AP per patient | positives in top-20 |
|---|---|---|---|
| presentation only | **0.207** | **0.266** | 17 / 35 |
| presentation .85 + TCR prior .15 | 0.178 | 0.250 | **18 / 35** |
| presentation-dominant .85 | 0.176 | 0.226 | 15 / 35 |
| presentation-dominant .70 | 0.165 | 0.211 | 16 / 35 |
| old literature-balanced (.30) | 0.149 | 0.195 | 15 / 35 |

Monotone in the weight on presentation. `config.PRESENTATION_ONLY` ships as a
preset for anyone whose only objective is hit rate; the shipped default keeps a
little weight on the others plus the full weight on expression and clonality,
which this dataset cannot evaluate at all.

**3. Per-patient variance dwarfs the model differences.** Ranked by
presentation, in a 20-slot budget:

| patient | pMHC | positives | recovered in top-20 |
|---|---|---|---|
| 1 | 82 | 9 | 4 |
| 2 | 97 | 4 | 3 |
| 3 | 84 | 12 | 5 |
| 10 | 59 | 3 | 3 |
| 12 | 79 | 4 | 2 |
| 16 | 121 | 3 | **0** |

Patient 16 gets nothing. Any claim about a pipeline's hit rate that is not
reported per patient is hiding this.

**Caveats**: 35 positives across 6 patients is small; the mirror is a subset of
the published supplement; and TESLA's tested peptides were nominated by
prediction pipelines in the first place, so the *candidate* set is not an
unbiased sample of all mutant peptides — though the labels, which are assay
outcomes, are unbiased with respect to that.

```bash
python -m neoantigen_pipeline.run_demo --out demo_out --tesla
```

---

## The number that matters more than the AUC

Published prospective work — the TESLA consortium's blinded comparison of
neoantigen prediction pipelines (Wells 2020, *Cell* 183:818) — found that the
fraction of top-ranked predicted neoantigens that were actually recognized by
patient T cells was in the low tens of percent, across every participating
pipeline. Any pipeline claiming near-certainty about a 34-peptide payload is
overclaiming, this one included. The correct posture is: the ranking makes a
34-slot budget much better spent than random, and the readout is still the
patient's T cells.
