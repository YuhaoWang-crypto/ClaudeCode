# The metrics, and why they fight each other

Definitions here are ported from Arc's `cell-eval` (v0.8.2, PyPI). Read the
source rather than the README — the README lists metric names without formulas.

```bash
pip download --no-deps cell-eval && unzip -q cell_eval-*.whl -d src
# src/cell_eval/metrics/_anndata.py  _de.py  _impl.py ; src/cell_eval/_score.py
```

## The three headline metrics

**Perturbation discrimination (`discrimination_score`).** For each
perturbation, rank *every measured* perturbation effect by L1 distance to the
*predicted* effect for that one, then report `1 - rank / n_perturbations`. 1.0
is perfect, **0.5 is chance**. The knocked-down gene is dropped from the feature
set (`exclude_target_gene=True`), so nailing on-target knockdown cannot win it.

This is a **retrieval** metric, and that governs model design. It does not ask
whether the prediction is close to the truth; it asks whether the prediction is
closer to *its own* truth than to 2,052 others. Consequences:

- The *spread of effect magnitudes across perturbations* carries most of the
  signal. If every prediction comes out the same size, the nearest measured
  effect is whichever one happens to be that size.
- Therefore every operation that pulls effects toward each other — smoothing,
  shrinkage, low-rank projection — costs discrimination even while it improves
  agreement with the truth. Renormalise magnitude afterwards.
- It is not scale-invariant. Scaling a prediction changes the ranking.

**Differential expression (`overlap_at_k`).** Real genes are those passing
FDR, ranked by |log2FC|; predicted genes ranked by predicted |log2FC|; score is
`|intersection| / k_eff` with `k_eff = min(k, n_real_significant)`.

**MAE.** Mean absolute error on *absolute expression*, not on the delta. The
`Δ=0` prediction — just return the control profile — is very hard to beat here.
In the 2025 challenge most submissions were worse than baseline on MAE.

## The aggregation, and the trap in it

`cell_eval/_score.py` renormalises every metric against a reference model:

```
best-is-zero (MAE, MSE):  1 - user / base
best-is-one  (others):    (user - base) / (1 - base)
then: nan -> 0, clip(lower=0), mean
```

The reference model is `build_base_mean_adata` — predict the *mean perturbation
response* for every perturbation.

**The clip at zero is correct for a leaderboard and wrong as a search
objective.** Once a metric is worse than baseline it contributes zero, so
making it arbitrarily worse is free. A tuner given this objective will discover
that it can multiply its predicted effects by 3, destroy MAE completely, and
buy discrimination with the proceeds. Observed directly: a search on the clipped
objective returned a model with discrimination 0.823 and MAE 0.1023 against a
baseline MAE of 0.0459.

Fix: **select on the unclipped form, report the clipped one.** The challenge
itself guards this differently — it "enforce[s] minimum thresholds on all
metrics to promote a balanced performance, discouraging models that perform
well on one metric at the expense of the others" — but the thresholds are not
published.

**Second trap: metric-family imbalance.** Averaging a full metric suite puts
four DE-derived metrics against one error metric, so MAE is already only an
eighth of the objective before any clipping. Score on the three the commentary
actually names, which keeps MAE at a third.

## Differential expression without single cells

The challenge computes significance from submitted single cells. At pseudobulk
there are none, so significance for the *measured* data comes from an empirical
null built out of the non-targeting control replicates.

The variance model matters. A pseudobulk profile from `n` cells has variance
`batch + sampling / n`: the batch term is gem-group drift that does not shrink
with more cells, the sampling term does. Control replicates span a range of cell
counts, so regressing their squared deviations on `1/n` separates the two.

Without that split, a perturbation profiled from 45 cells gets judged against a
null built from control replicates of ~1,000 and nearly every gene looks
significant.

Then moderate each variance component toward its **expression-matched** trend,
limma-style. Moderating toward a global constant instead inflates the variance
of lowly expressed genes by orders of magnitude and silences them entirely —
this produced *zero* significant genes across every perturbation before it was
caught.

Sanity check the result: median 11–38 DE genes per perturbation, 54–68% of
perturbations with ≥10, maxima in the low thousands for the strongest
knockdowns. Zero, or "everything is significant", both mean the variance model
is broken.

## Sanity checks that catch real bugs

Run these before trusting any number:

| Check | Expected |
|---|---|
| discrimination, prediction = truth | 1.00 |
| discrimination, random prediction | 0.50 |
| discrimination, same prediction for every perturbation | 0.50 |
| DE overlap, `Δ=0` prediction | ≈ k/n_genes (chance) |
| DE overlap, global-mean prediction | small but above chance |

A model scoring 0.5 on discrimination is not broken — it is predicting the same
thing for everything, which is what the `Δ=0` and global-mean baselines do by
construction.
