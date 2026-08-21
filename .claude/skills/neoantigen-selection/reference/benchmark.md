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

## The number that matters more than the AUC

Published prospective work — the TESLA consortium's blinded comparison of
neoantigen prediction pipelines (Wells 2020, *Cell* 183:818) — found that the
fraction of top-ranked predicted neoantigens that were actually recognized by
patient T cells was in the low tens of percent, across every participating
pipeline. Any pipeline claiming near-certainty about a 34-peptide payload is
overclaiming, this one included. The correct posture is: the ranking makes a
34-slot budget much better spent than random, and the readout is still the
patient's T cells.
