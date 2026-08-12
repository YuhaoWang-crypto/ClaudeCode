# virtualcell — a zero-shot, multi-context virtual cell model

Built from the technologies surveyed in
[`SURVEY.md`](SURVEY.md), against the task shape announced for
**Virtual Cell Challenge 2**:

> predict how multiple cell lines respond to specified gene knockdowns, given
> only non-targeting control profiles … without your model having been trained
> on any of that data.

Nothing about that task can be tested on the challenge's own held-out data, so
this package rebuilds the task on public data and answers the question it is
really asking: **does the model transfer into a cell line it has never seen a
perturbation in, and does it beat the baselines while doing so?**

## The benchmark

Four CRISPRi Perturb-seq contexts, one perturbation library, one processing
pipeline, so a knockdown means the same thing in all of them:

| Cell line | Origin | Source | Perturbations | Control replicates |
|---|---|---|---|---|
| K562 | chronic myeloid leukaemia | Replogle 2022, Figshare+ 20029387 | 2,057 | 109 |
| RPE1 | retinal pigment epithelium | Replogle 2022, Figshare+ 20029387 | 2,393 | 130 |
| HepG2 | hepatocellular carcinoma | Nadig 2025, GEO GSE264667 | 2,393 | 56 |
| Jurkat | T-cell leukaemia | Nadig 2025, GEO GSE264667 | 2,393 | 55 |

**6,642 genes and 2,053 knockdowns are shared by all four.** Every fold hides
one line completely: the model receives that line's non-targeting controls and
nothing else from it. Two regimes are run —

- **context** — the knockdown was measured in the source lines, never in the
  target. This is the challenge's setting.
- **double-blind** — the knockdown is *also* deleted from every source line.
  Neither the context nor the perturbation has ever been seen.

Hyperparameters are selected by a nested leave-one-out over the *source* lines
only, so the held-out line never touches model selection.

## The model

`ContextTransferModel` predicts the effect of a knockdown as a transferred
consensus, reshaped for the target context. Six stages, each switchable, all
switches set by cross-validation:

1. **Context-weighted consensus** — average the effect over source lines,
   weighting by control-transcriptome similarity to the target, measured on the
   genes that actually vary between cell lines.
2. **Perturbation-neighbour smoothing** — borrow strength from knockdowns with
   correlated effects. This is the only route to a gene never perturbed
   anywhere, and is what makes the double-blind regime possible.
3. **Reliability shrinkage** — a James-Stein factor per perturbation, estimated
   from cross-line reproducibility. Effects that replicate survive; effects that
   do not are pulled toward the generic response.
4. **Program-basis denoising** — project onto the top components of the source
   response matrix, blended rather than applied outright, because denoising
   sharpens agreement with the measured effect but blurs perturbations together
   and discrimination punishes that.
5. **On-target knockdown** — applied in count space, using the *per-gene*
   residual measured in the source lines. Knockdown efficiency is a property of
   the guide pair, not the cell line: it ranges from under 5% to over 30%
   residual and transfers across contexts at Spearman 0.38–0.54.
6. **Global calibration** — one scalar. Averaging several noisy sources shrinks
   toward the shared component, so the transferred effect is systematically too
   small and needs scaling back up.

Baselines: predicting the control profile unchanged (Δ=0); predicting the mean
perturbation response for everything (this is `cell-eval`'s
`build_base_mean_adata`, the model leaderboard scores are normalised against);
and unweighted cross-line transfer.

## Metrics

`metrics.py` is a port of Arc's [`cell-eval`](https://pypi.org/project/cell-eval/)
0.8.2 — the discrimination score (L1, target gene excluded), `mae`/`mae_delta`,
`pearson_delta`, `overlap_at_k`, direction match, LFC Spearman, and the
baseline-normalised aggregation from `_score.py`.

One deviation is stated rather than hidden: the challenge computes significance
from submitted single cells, and this works on pseudobulk. Significance for the
*measured* data therefore comes from an empirical null built from the
non-targeting control replicates, with variance split into a batch component
and a per-cell sampling component (`control_variance_model`) so that a
perturbation profiled from 45 cells is not judged against a null built from
control replicates of 1,000. Predictions are ranked by predicted |log2FC|
without a significance filter of their own. This is applied identically to
every model compared.

## Running it

```bash
pip install numpy scipy pandas scikit-learn anndata h5py matplotlib

# one-off: pseudobulk the GEO single-cell files (5.6 GB / 9.4 GB downloads)
python -m virtualcell.prep_nadig GSE264667_hepg2_raw_singlecell_01.h5ad hepg2.npz
python -m virtualcell.prep_nadig GSE264667_jurkat_raw_singlecell_01.h5ad jurkat.npz

python -m virtualcell.data                      # check what loaded
python -m virtualcell.run --regime context      # the challenge's setting
python -m virtualcell.run --regime double       # perturbation unseen too
python -m virtualcell.run --ablation            # how much each context buys
python -m virtualcell.figures
```

Results land in `results/*.json` (+ a rendered `.txt`), figures in
`figures/virtualcell/`.

## Files

| File | Role |
|---|---|
| `prep_nadig.py` | streams the GEO h5ad files and pseudobulks them |
| `data.py` | harmonises four cell lines onto one gene axis |
| `metrics.py` | the `cell-eval` port |
| `model.py` | `ContextTransferModel` and the three baselines |
| `benchmark.py` | leave-one-cell-line-out harness, nested tuning, ablation |
| `run.py` | CLI |
| `figures.py` | the two charts |
| `SURVEY.md` | every technology cited by the challenge paper, and its status here |
| `RESULTS.md` | what the benchmark found |
