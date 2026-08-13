# Results

*Numbers in this file are rendered from `results/*.json` by
`python -m virtualcell.report`; the tables live in `results/tables.md`. Nothing
here is hand-transcribed.*

---

## The answer, in one table

Context zero-shot — the challenge's own setting. 2,053 knockdowns, mean over
four held-out cell lines, each scored with no perturbation data from it in
training.

| model | discrimination | DE overlap@100 | MAE | VCC score | balanced |
|---|---|---|---|---|---|
| control (Δ=0) | 0.5010 | 0.0047 | 0.0570 | 0.009 | −0.020 |
| global mean *(challenge baseline)* | 0.5011 | 0.0716 | 0.0576 | 0.000 | +0.000 |
| naive cross-line transfer | **0.6495** | **0.2489** | 0.0613 ✗ | **0.163** | +0.139 |
| **ContextTransfer** | 0.6235 | 0.2432 | **0.0551** | 0.158 | **+0.158** |

*✗ = worse than the challenge baseline.*

Both transfer models are far above the baseline on discrimination (0.50 is
chance) and DE overlap. Between them the result is **close, and which one wins
depends on which scoring is quoted** — so here is the honest reading:

**ContextTransfer is the only model that clears the baseline on all three
headline metrics, on every fold.** Naive transfer is worse than baseline on MAE
in all four (0.053/0.069/0.061/0.062 against 0.046/0.068/0.059/0.058). Since the
challenge "enforce[s] minimum thresholds on all metrics … discouraging models
that perform well on one metric at the expense of the others", a model that
fails one threshold on every fold is in a materially different position from one
that fails none.

**ContextTransfer predicts a more accurate response; naive transfer predicts a
more distinguishable one.** On every accuracy measure the ordering is
consistent:

| | naive transfer | ContextTransfer |
|---|---|---|
| DE direction agreement | 0.8191 | **0.8410** |
| DE log2FC Spearman | 0.4664 | **0.4921** |
| Pearson (effect) | 0.2781 | **0.3278** |
| MAE (effect) | 0.0613 | **0.0551** |

Naive transfer buys its discrimination lead by over-predicting effect
magnitude — which is also why its error is worse than doing nothing at all.

**On the leaderboard's own aggregation naive transfer edges ahead (0.163 vs
0.158); on the unclipped aggregation ContextTransfer leads (+0.158 vs +0.139).**
The gap between those two readings *is* the finding: the clipped score forgives
naive transfer's MAE failure entirely, because a metric below baseline
contributes zero rather than a penalty.

So: **comparable in aggregate, better on every measure of response accuracy,
and the only model meeting every threshold** — not a sweep, and not nothing.

### Per-fold

| held out | | discrimination | DE overlap@100 | MAE |
|---|---|---|---|---|
| **K562** | naive | **0.726** | 0.261 | 0.053 ✗ |
| | ContextTransfer | 0.676 | **0.281** | **0.043** |
| **RPE1** | naive | 0.568 | **0.255** | 0.069 ✗ |
| | ContextTransfer | 0.567 | 0.207 | **0.065** |
| **HepG2** | naive | **0.665** | 0.270 | 0.061 ✗ |
| | ContextTransfer | 0.616 | **0.279** | **0.056** |
| **Jurkat** | naive | **0.639** | **0.209** | 0.062 ✗ |
| | ContextTransfer | 0.636 | 0.206 | **0.057** |

RPE1 is the hardest context for every model (discrimination 0.57 against K562's
0.68–0.73), and DE overlap is the one metric whose winner flips between folds.

Full tables, including supplementary metrics and per-fold hyperparameters, are
in [`../results/tables.md`](../results/tables.md).

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

### ⚠️ Magnitude renormalisation works, and was still rejected

Smoothing, shrinkage and low-rank projection all improve the *shape* of a
predicted response and all pull effects toward each other, flattening the
across-perturbation magnitude spread that discrimination depends on. That
mechanism is real: restoring each row's pre-denoising norm lifted discrimination
from 0.714 to 0.761 on K562 held out.

**It costs more MAE than it buys.** That configuration's error rose from 0.0464
to 0.0550, and cross-validation set `renorm=0.0` on all four folds. What
survived selection was restraint rather than compensation — blended denoising
(`rank_mix` 0.5–0.75), light smoothing (`smooth` ≤ 0.15), shrinkage off on three
of four folds, and a single global scale `beta` between 0.8 and 1.25 doing the
magnitude work.

Recorded here because an earlier draft of this report called renormalisation
"the fix that matters" on the strength of a single-fold sweep, before the full
cross-validation contradicted it.

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
