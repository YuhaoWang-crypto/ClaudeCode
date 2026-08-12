# Designing the model, and the ideas that failed

## Start from the 2025 evidence, not from architecture fashion

Arc's own wrap-up of Virtual Cell Challenge 2025: *"purely AI-based approaches
did not consistently outperform statistical baselines."* Three of the four
prize-winning entries were summary-level statistical transfer with context
conditioning — the winner combined deep learning with classical statistics and
protein embeddings, second place used pseudobulk plus a plain fully-connected
network and ESM-2, third place used gene-level summaries with similarity-aware
aggregation and a global linear scale.

So the default should be a structured statistical model whose components are
switchable and cross-validated, not an end-to-end network. A pretrained gene or
protein embedding is the highest-value addition (second place's ESM-2), not a
bigger decoder.

## The decomposition

Predict the **effect**, never the absolute profile:

```
y(context c, knockdown p) = μ_c + δ(c, p)
```

`μ_c`, the control mean, is given at test time — that is the entire content of
"only non-targeting control profiles". So the whole problem is `δ`, and MAE is
mostly determined by how well you resist over-predicting it.

## The seven stages

1. **Context-weighted consensus.** Average `δ(c, p)` over source lines,
   weighting by control-transcriptome similarity to the target. Take that
   correlation on the genes that *vary between* cell lines — on the full
   transcriptome any two human lines correlate above 0.9 and the weights carry
   no information.
2. **Perturbation-neighbour smoothing.** Blend each effect with its nearest
   neighbours in effect space, weighted by their own reliability. Build the
   mixing matrix sparse; gathering `n × k × G` densely is gigabytes here.
3. **Reliability shrinkage.** Estimate what fraction of a measured effect is
   signal from **cross-line reproducibility** — a knockdown that looks the same
   in two unrelated lines is measuring biology, one that does not is measuring
   sampling noise from a few dozen cells. Fall back to a Wiener factor with one
   source line.
4. **Program-basis denoising.** SVD the source response matrix once, project
   onto the leading directions, blend rather than replace.
5. **Magnitude renormalisation.** Restore each row's pre-denoising norm. See
   below — this is the stage that makes stages 2–4 affordable.
6. **On-target knockdown.** Apply in count space, where CRISPRi is
   multiplicative, at the **per-gene** efficiency measured in the sources.
7. **Global calibration.** One scalar. Averaging several noisy sources shrinks
   toward the shared component, so the transferred effect comes out too small.

Knockdowns absent from *every* source line take a separate route: the silenced
gene is still *measured*, so its behaviour across the rest of the knockdown
panel identifies functionally related genes whose own knockdowns stand in for
it. This uses no measurement of the gene being silenced, only of the gene being
watched, so it is leakage-free.

## Stage 5 is the one that matters

Stages 2–4 all improve the *shape* of the predicted response and all pull
effects toward each other, which flattens the across-perturbation magnitude
spread — the thing discrimination is built on. The result is a model that
predicts a better response and retrieves worse.

Restoring each row's original norm keeps the improved shape at the original
spread. On K562 held out this closed a 0.05 discrimination gap without giving
back the DE or MAE gains.

## Rejected ideas, with the numbers

Do not re-derive these.

**Effect magnitude scales with the target gene's expression in the new line.**
The intuition is compelling — a gene expressed less should respond less to being
silenced. It is false. Relative target expression versus relative effect
magnitude across contexts: **Spearman −0.014** (n=7,080). Meanwhile magnitude
itself transfers well between lines (ρ 0.62–0.69). Effect magnitude is a
property of the perturbation, not of the context.

**Soft-thresholding to cut MAE.** Zero the small components, keep the large
ones. Discrimination 0.714 → 0.653 for an MAE gain of 0.0035; every threshold
scored worse overall. L1 retrieval uses precisely the components being zeroed.

**Reliability shrinkage at full strength.** Pulling unreliable perturbations
toward the generic response makes them mutually indistinguishable, and
discrimination collapses. Cross-validation drives this switch to zero.

**Per-gene expression-ratio context modulation.** Real but marginal: Pearson
0.3003 → 0.3029 at the best exponent, MAE_delta 0.06941 → 0.06894. Keep it as a
tunable switch; do not build a story on it.

**Low-rank projection applied outright.** Improves Pearson and MAE, costs
discrimination. Only useful blended, and only with stage 5.

## Calibrate expectations

On four matched CRISPRi lines:

- across-line mean effect = **43%** of total effect variance
- per-perturbation cross-line correlation, median **0.11–0.26**
- relative effect magnitude between lines, Spearman **0.62–0.69**
- on-target residual expression **9–18%** median, conserved at ρ 0.38–0.54

A model transferring the shared component with correct magnitude should reach
discrimination ~0.68–0.73 against a 0.50 chance level. Claims of predicting the
context-specific 57% need evidence, since the obvious mechanism for it does not
exist.

## Baselines you must include

- **`Δ=0`** — return the control profile. Wins MAE outright; sits at chance on
  discrimination. Its purpose is to show how much of MAE is free.
- **Global mean effect** — `cell-eval`'s own reference model, and what the
  leaderboard normalises against.
- **Naive cross-line transfer** — unweighted mean of source effects. This is the
  honest thing to beat: it already uses cross-context information and reaches
  discrimination ~0.73. A model that beats only the trivial baselines has shown
  nothing.
