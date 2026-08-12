# Results

*Numbers in this file are rendered from `results/*.json` by
`python -m virtualcell.report`; the tables live in `results/tables.md`. Nothing
here is hand-transcribed.*

**Status: benchmark running.** This document holds the design, the reading
protocol, and the findings already established. The four-fold tables land when
the run completes.

---

## The question

Virtual Cell Challenge 2 asks for predictions of *"how multiple cell lines
respond to specified gene knockdowns, given only non-targeting control
profiles"*, scored against new Arc data *"without your model having been trained
on any of that data"*.

That is not something a model can be checked against before submission, because
the scoring data does not exist publicly. So the question this benchmark answers
is the one that can be answered:

> Given only a cell line's non-targeting controls, can a model predict its
> knockdown responses better than the challenge's own baseline — and better than
> simply copying the effect measured in other cell lines?

The second half matters more than the first. Beating a trivial baseline shows
nothing; naive cross-line transfer is the honest comparator.

## The benchmark

Four CRISPRi Perturb-seq contexts, one library design, one processing pipeline,
so a knockdown means the same thing in all of them:

| Cell line | Origin | Source | Perturbations | Control replicates | Median cells/perturbation |
|---|---|---|---|---|---|
| K562 | chronic myeloid leukaemia | Replogle 2022 | 2,057 | 109 | 110 |
| RPE1 | retinal pigment epithelium | Replogle 2022 | 2,393 | 130 | 69 |
| HepG2 | hepatocellular carcinoma | Nadig 2025 | 2,393 | 56 | 44 |
| Jurkat | T-cell leukaemia | Nadig 2025 | 2,393 | 55 | 81 |

**6,642 genes and 2,053 knockdowns shared by all four.** Each fold hides one
line completely: the model gets that line's non-targeting controls and nothing
else from it. Hyperparameters come from a nested leave-one-out over the *source*
lines only.

Two regimes:

- **context** — the knockdown was measured in the source lines, never in the
  target. The challenge's own setting.
- **double-blind** — the knockdown is also deleted from every source line.
  Neither context nor perturbation has been seen.

Models compared: return the control profile unchanged (`Δ=0`); predict the mean
perturbation response for everything (`cell-eval`'s reference model, which the
leaderboard normalises against); unweighted cross-line transfer; and
`ContextTransfer`.

## How to read the scores

The three headline metrics **conflict by construction**, so a single number
cannot summarise them and any report that offers one is hiding something.

- **Perturbation discrimination** is an L1 *retrieval* metric: it asks whether a
  prediction is closer to its own truth than to 2,052 others. It rewards large,
  well-spread effects. 0.5 is chance.
- **MAE** rewards small effects. The `Δ=0` prediction is very hard to beat, and
  in the 2025 challenge most submissions were worse than baseline on it.
- **DE overlap** sits in between.

Two scorings are therefore always reported:

- **VCC score** — `cell-eval`'s aggregation, which clips a metric worse than
  baseline to zero. This is the leaderboard number.
- **balanced** — the same without the clip, so a metric worse than baseline
  counts against the model. This is what hyperparameters were selected on,
  because under the clipped form a metric already below baseline is free to get
  arbitrarily worse, and a search will sell all of MAE for a little
  discrimination. (Observed directly: the clipped objective returned a model
  with discrimination 0.823 and MAE 0.1023 against a baseline MAE of 0.0459.)

## Established findings

These were measured during development and do not depend on the final tables.

### Cross-context transfer works, and has a ceiling

| Quantity | Value |
|---|---|
| Across-line mean effect, as a fraction of total effect variance | **43%** |
| Per-perturbation cross-line effect correlation | median **0.11–0.26** |
| Relative effect *magnitude* between lines (Spearman) | **0.62–0.69** |
| On-target residual expression after CRISPRi | **9–18%** median |
| On-target efficiency conservation between lines (Spearman) | **0.38–0.54** |

Effect *magnitude* and *knockdown efficiency* transfer well; effect *shape*
transfers weakly. A model that transfers the shared component and gets its
magnitude right captures roughly the 43%. The remaining 57% is context-specific
or noise, and nothing here predicts it.

### ❌ Target-gene expression does not predict effect magnitude

The intuitive mechanism for context-specificity — *a gene expressed less in the
new cell line should respond less to being silenced* — is **not supported**:

| Test | Result |
|---|---|
| Relative target expression vs relative effect magnitude, across contexts | Spearman **−0.014** (n = 7,080) |
| Target expression vs effect magnitude, within a line | Spearman 0.05–0.15 |

Effect magnitude is a property of the perturbation, not of the context. This
idea was dropped rather than shipped.

### ❌ Soft-thresholding does not resolve the MAE/discrimination tension

Zeroing small predicted components to cut MAE costs far more than it recovers:
discrimination 0.714 → 0.653 for an MAE gain of 0.0035. Every threshold tested
scored worse overall. L1 retrieval uses precisely the components being zeroed.

### ✅ Magnitude renormalisation is what makes denoising affordable

Smoothing, shrinkage and low-rank projection all improve the *shape* of a
predicted response and all pull effects toward each other, flattening the
across-perturbation magnitude spread that discrimination depends on. Restoring
each row's pre-denoising norm keeps the improved shape at the original spread.

### ✅ Gene-response similarity reaches knockdowns seen nowhere

In the double-blind regime a knockdown has no consensus row to transfer, and
both baselines sit at chance. Routing through the silenced gene's behaviour
across the rest of the knockdown panel recovers real signal — DE overlap 0.060
(chance) → 0.132 on RPE1 held out. This uses no measurement of the gene *as
perturbed*, only of the gene *as observed*, so it is leakage-free.

## Limitations

- **Depth.** 44–110 cells and ~11–14k UMI per perturbation, against the
  challenge's ~1,000 cells and >50k UMI on 10x Flex. Absolute scores here are
  not comparable to leaderboard numbers.
- **Pseudobulk, not single cells.** The challenge computes DE from submitted
  cells. Significance for the measured data here comes from a control-replicate
  empirical null that separates batch from per-cell sampling variance. A real
  submission needs a single-cell sampling layer.
- **Essential genes only.** All four screens target essential genes; behaviour
  on regulatory or non-essential targets is untested.
- **Four contexts is few.** Context weights are fitted from three source lines
  per fold — barely more than a heuristic.
- **One lab, one design.** What makes these four lines comparable also means
  shared systematic error is invisible here in a way it would not be across
  independent studies.
- **Gene universe and evaluation set** are intersections taken across all four
  lines, which touches the target's gene and perturbation lists. Neither leaks a
  perturbation *response*, and in the real challenge both are known anyway.

## What a real Challenge 2 entry would add

In descending order of expected value:

1. **More contexts.** The VCC 2025 H1 hESC data as a fifth line (Challenge 2
   explicitly permits its reuse), and Tahoe-100M's ~50 lines for a context
   encoder that is more than a three-point heuristic. The context ablation
   measures what each added context is worth.
2. **A pretrained gene/protein embedding.** ESM-2 was the cheapest consistent
   gain among 2025 entrants, and would replace the data-derived gene-similarity
   routing used for unseen knockdowns.
3. **A single-cell output layer**, so DE scoring is native rather than
   approximated.
4. **The K562 genome-wide screen** (already downloaded) to widen the
   perturbation panel well beyond essential genes.
