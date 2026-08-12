# Benchmark protocol

The challenge scores against data nobody outside Arc has. Everything below
reconstructs the *task shape* on public data so a model can be judged before
submission.

## Leave-one-cell-line-out

Each fold hides one cell line completely. The model receives that line's
non-targeting control profiles and nothing else from it. With four contexts
this gives four folds, each with three source lines.

Run two regimes:

**`context`** — the knockdown was measured in the source lines, never in the
target. This is the challenge's own setting, and the one that matters.

**`double`** — the knockdown is *also* deleted from every source line. Neither
context nor perturbation has been seen. Almost every model collapses to the
baseline here, which is exactly why it is worth running: it separates models
that transfer a lookup from models that generalise.

## Nested tuning, source lines only

Hyperparameters are chosen by a leave-one-out *inside the source lines*. For
target `c*` with sources `S`, hold out each `s ∈ S` in turn, fit on `S \ {s}`,
score on `s`, average. The held-out line never influences model selection.

In the `double` regime the inner folds must drop their evaluation perturbations
from the inner sources too, or hyperparameters get tuned for an easier task than
the one they will face.

Coordinate ascent over a small grid is sufficient; two passes captures almost
all of it (one observed run: 0.1902 after pass 0, 0.1909 after pass 1). **Check
whether any switch saturated at a grid edge** — if it did, the grid is wrong,
not the model. A global scale pinned at the maximum is the usual sign that the
objective is broken (see `metrics.md`).

## Leakage audit

Grep the harness for every use of the target line and confirm `fit()` receives
only `target.mu`. In a correct implementation `target.pert` and `target.delta`
appear *only* inside the evaluation function, as ground truth.

Two things are legitimate and should still be disclosed:

- **The gene universe** is the intersection across all four lines, which uses
  the target's gene list. In the real challenge the assay panel is known
  anyway.
- **The evaluation perturbation set** is the intersection across lines. In the
  real challenge you are told which genes to predict.

Neither leaks a perturbation *response*.

## Making the comparison meaningful

- Evaluate on all shared perturbations, not a favourable subset. If you subsample
  for cost, say so — discrimination depends on how many candidates it ranks
  against, so a score over 100 candidates is not comparable to one over 2,053.
- Select the evaluation set on technical criteria (cell count) or not at all.
  Selecting on measured effect size in the *target* line biases toward
  perturbations that happen to be strong there.
- Stratify the report by measured effect strength. Roughly half an essential-gene
  panel produces almost no transcriptional response at 45–120 cells per
  perturbation, and whole-panel averages are dominated by those. Separate
  "predicts strong responses" from "predicts that most things do nothing".

## The context ablation

The single most informative extra experiment. Fit on each possible *subset* of
source lines — 1, 2, then 3 — with hyperparameters held fixed, and plot score
against the number of contexts.

Because all four lines share nearly the same perturbation panel, adding a source
line adds almost no new *knockdowns*, only new *contexts*. Any improvement is
therefore attributable to contextual diversity rather than to scale — which is
the claim in *Virtual Cells Need Context, Not Just Scale* (2026), stated there
from data-generation arguments and testable here directly.

## Reporting

Report per-fold and averaged, and always include:

- all three headline metrics as raw numbers, not only the aggregate
- both the clipped (leaderboard) and unclipped (balanced) aggregate
- naive cross-line transfer in the table as a first-class comparator
- the number of perturbations evaluated and the DE-gene distribution
- cells per perturbation and sequencing depth relative to the challenge data,
  since public Perturb-seq is far thinner (~45–120 cells, ~11–14k UMI) than the
  challenge's 10x Flex data (~1,000 cells, >50k UMI), and absolute scores are
  not comparable across that gap

A single-number claim of "better" is not reportable here. The metrics conflict
by construction, and which model looks best depends on which scoring is quoted.
