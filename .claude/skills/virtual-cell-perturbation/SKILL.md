---
name: virtual-cell-perturbation
description: >-
  Build and honestly benchmark a virtual cell model that predicts single-cell
  transcriptional responses to gene knockdowns in a cell line it was never
  trained on, given only that line's non-targeting control profiles — the
  Virtual Cell Challenge 2 task. Use when asked to build a virtual cell or
  perturbation-response model, to transfer CRISPRi/Perturb-seq effects across
  cell lines or contexts, to predict knockdown effects zero-shot or few-shot, to
  score predictions with the challenge metrics (perturbation discrimination,
  differential expression, MAE) or Arc's cell-eval, to assemble multi-cell-line
  Perturb-seq data (Replogle K562/RPE1, Nadig Jurkat/HepG2), or to design a
  leave-one-cell-line-out benchmark without leakage. Encodes which modelling
  ideas measurably worked and which measurably did not.
---

# Zero-shot cross-context perturbation prediction

A methodology plus a working `virtualcell` package for the task the Arc
Institute's Virtual Cell Challenge 2 poses: **given only non-targeting control
profiles of a cell line, predict how it responds to specified gene knockdowns**,
with no perturbation data from that line in training.

## Start here

1. **Know what you are being scored on.** The three headline metrics are not
   interchangeable and they actively conflict. Read `reference/metrics.md`
   before writing any model — the discrimination score is an L1 *retrieval*
   metric, and that single fact determines most design decisions.
2. **Get four contexts, not one.** `reference/data-access.md` has exact URLs,
   file sizes and the pseudobulking recipe for four matched CRISPRi cell lines
   sharing ~2,000 knockdowns. One cell line cannot test context generalisation.
3. **Build the model from `reference/model-design.md`.** It documents the seven
   stages, and — more usefully — the ideas that *failed* with the numbers that
   killed them.
4. **Benchmark with `reference/benchmark-protocol.md`.** Leave-one-cell-line-out
   with nested tuning on source lines only. The leakage rules are specific.

## The one thing that decides everything

Perturbation responses split into a part that transfers between cell lines and
a part that does not. On four matched CRISPRi lines, **the across-line mean
effect carries 43% of total effect variance**, and per-perturbation cross-line
correlation has a median of only 0.11–0.26. So:

- A model that transfers the shared component and gets its magnitude right will
  beat the challenge baseline comfortably.
- Anything claiming to predict the context-specific 57% needs evidence, because
  the obvious mechanisms for it do not work (see below).

## What measurably worked

| Idea | Effect |
|---|---|
| Cross-line consensus (average the effect over source lines) | discrimination 0.50 → 0.72 vs the challenge baseline |
| Per-gene on-target knockdown transferred from sources | knockdown efficiency is a guide-pair property, conserved across lines at Spearman 0.38–0.54, and ranges from <5% to >30% residual — one global constant wastes that |
| Program-basis denoising, *blended* not applied outright | selected on all four folds (`rank_mix` 0.5–0.75); improves agreement with the measured effect |
| Keeping smoothing and shrinkage light | cross-validation drove `smooth` to ≤0.15 and `shrink` to 0 on three of four folds — both flatten the magnitude spread discrimination needs |
| Gene-response-similarity routing for unseen knockdowns | the only thing that works at all when the gene appears in no source line; DE overlap 0.060 (chance) → 0.132 |

⚠️ **Magnitude renormalisation** — restoring each row's pre-denoising norm — is
implemented and available as a switch, and it does move discrimination up
(0.714 → 0.761 on one fold). But it costs more MAE than it buys, and
cross-validation **rejected it on all four folds** (`renorm=0.0` everywhere).
Reach for it only if you are optimising discrimination specifically and are
willing to fail an MAE threshold. The global scale `beta` is the better-behaved
control on the same trade-off.

## What measurably did not work — do not re-derive these

| Idea | Why it failed |
|---|---|
| Scaling the effect by the target gene's expression in the new line ("a gene expressed less should respond less") | Spearman **−0.01** between relative target expression and relative effect magnitude across contexts. There is no signal. Effect magnitude is a property of the perturbation (conserved at ρ≈0.65), not of the context |
| Soft-thresholding the prediction to cut MAE | discrimination 0.714 → 0.653 for an MAE gain of 0.0035. Zeroing small components removes exactly what L1 retrieval uses |
| Per-gene expression-ratio context modulation | real but marginal: Pearson 0.3003 → 0.3029 |
| Reliability shrinkage toward the generic response, at full strength | destroys discrimination — shrunk perturbations become mutually indistinguishable |
| Heavy smoothing / low-rank projection alone | same failure: both blur perturbations together |

## Honesty requirements

Every claim in an output must be labelled, in the style this repository already
uses:

- ✅ **measured here** — a number this pipeline computed, with the fold or
  dataset it came from.
- ⚠️ **hypothesis** — a mechanism proposed but not demonstrated.
- ❌ **tested and rejected** — with the number that rejected it.

Two specific disclosures are mandatory whenever results are reported:

1. **Which scoring is quoted.** Arc's aggregation clips a metric worse than
   baseline to zero. That is right for a leaderboard and wrong as a search
   objective — under it, a metric already worse than baseline is free to get
   worse, and tuning will sell all of MAE for a little discrimination. Report
   both the clipped (leaderboard) and unclipped (balanced) score.
2. **Where the DE significance came from.** At pseudobulk there are no single
   cells to test, so significance comes from a control-replicate empirical null
   that separates batch from per-cell sampling variance. Say so.

## Reference files

| File | Contents |
|---|---|
| `reference/data-access.md` | Exact download URLs, sizes, schemas, pseudobulking, disk strategy |
| `reference/metrics.md` | The cell-eval port, the metric conflict, the pseudobulk DE null |
| `reference/model-design.md` | The seven stages, and every rejected idea with its number |
| `reference/benchmark-protocol.md` | LOCO design, nested tuning, leakage audit, the two regimes |

The runnable implementation lives in `virtualcell/` at the repository root:
`data.py`, `metrics.py`, `model.py`, `benchmark.py`, `run.py`, `analysis.py`,
`figures.py`, with `prep_nadig.py` for the one-off GEO pseudobulking.
