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

## Before anything else: does your source cover the target panel?

Two lines, run against every candidate source, before harmonising a single file:

```python
panel = np.loadtxt("pert_counts.csv", dtype=str, skiprows=1)   # the target list
print(sum(p in set(line.names) for p in panel), "/", panel.size)
```

The four matched CRISPRi lines this skill recommends share an *essential-gene*
panel, and the Virtual Cell Challenge 2026 validation panel contains **0 of
300** of it — essential-gene screens target what a cell needs to survive, which
is exactly what a challenge testing regulatory prediction leaves out. Replogle's
genome-wide K562 arm covers **272 of 300** and becomes the source instead, at the
cost of dropping from four source contexts to one. That is an architecture
decision, and finding it late means rebuilding.

Keep the four-line atlas as the **benchmark** — four contexts is what makes
context generalisation measurable — and check panel coverage separately for the
**training corpus**. They are not the same choice. See `reference/data-access.md`.

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
| **A frontier language model as the context-specific term** | loses to the *context-free* baseline in both contexts (DE overlap@100 0.105/0.089 against 0.123/0.138) and to a transfer model by 2.5×. Blending its ranking in is monotonically worse at every weight. See below — the failure is not what you would expect |
| **A DNA sequence model (Evo2 and kin) for cross-context transfer** | structural, not empirical: it has no cell-context input, so its score for a gene is identical in every context and it cannot express a cross-context term at all. The context-free quantity it *could* supply — how strongly a knockdown acts — is already measured in the source and conserved across lines at ρ 0.62–0.69 |

### Why the language model fails, since the obvious diagnosis is wrong

Test it under the most favourable conditions: name the cell lines instead of
anonymising them, and ask explicitly for lineage-specific answers.
`virtualcell/llm_prior.py` does this on 24 knockdowns with strong measured
effects in both Jurkat and RPE1.

The intuitive explanation — *it gives a generic answer that ignores the
context* — is **false**. Its two answers share only a quarter of their genes
(Jaccard 0.248), and its reasoning is specifically lineage-aware: T-cell
identity genes, Jurkat's inactive p53 demoting CDKN1A, epithelial programmes in
RPE1. It differentiates strongly and is still wrong.

The likely reason, and the thing to check before trying this again: a language
model's knowledge is **pathway-level causal narrative** (*knock down NEDD8 →
cullins lose neddylation → IκBα persists → NF-κB output falls*), while measured
pseudobulk DE is dominated by **global state** — growth rate, ribosome content,
stress response. The biology is right and is not what the assay measures. Any
plan to inject literature knowledge has to bridge that gap explicitly rather
than assume the two vocabularies align.

## Cross-checked against two independent efforts

Two independent implementations of this task tested this methodology directly.
Their corrections are folded in above and in `reference/model-design.md`; the
full reconciliation with verification numbers is in `virtualcell/RESULTS.md`.
The load-bearing points:

- **Discrimination pays for across-perturbation magnitude spread, and this is
  the single biggest lever on it.** Measured truth has a magnitude-spread CV of
  ~0.45. A model that compresses toward ~0.15 scores at chance for unseen
  perturbations *by construction*, whatever its biology. Check the CV of your
  predicted effect norms before concluding anything about mechanism.
- **Which operating point cross-validation picks is decided by how the objective
  prices MAE.** Price it at a third and the search compresses spread; clip it to
  zero and the search expands spread. Same code, opposite conclusion. Report
  which you used.
- **The seen-target advantage of any transfer model is largely a scale
  artifact.** Independently found twice: at matched output scale, a plain
  cross-context consensus reproduces almost all of it. Always give every
  baseline its own best scale before claiming an architectural win.
- **`cell-eval` computes DES by Wilcoxon on the cells you submit.** Replicated
  pseudobulk means have zero within-group variance and manufacture significance
  — a predictor of nothing scores 0.083 that way. Submit cells with realistic
  dispersion, and never compare a pseudobulk DE-overlap number to a leaderboard
  DES.
- **External knowledge is what carries the unseen-perturbation regime.** Purely
  data-derived gene similarity reaches ~0.57 discrimination there; STRING +
  Reactome + ESM-2 embeddings reach ~0.67. Those same blocks contribute nothing
  measurable where transfer is possible.

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
