# Model card — ContextTransfer

## What it does

Given a cell line's **non-targeting control profiles** and a list of genes to
knock down, it predicts the post-perturbation transcriptome — without having
seen any perturbation measured in that cell line. This is the task shape
announced for the Arc Institute's Virtual Cell Challenge 2.

- **Input.** Control pseudobulk replicates for the target context (6,642 genes,
  log1p counts-per-10k), plus perturbation data from one or more *other* cell
  lines, plus the gene symbols to knock down.
- **Output.** A predicted effect per knockdown, and the absolute profile
  `control mean + effect`.
- **Not required.** Any perturbation measured in the target line. Any GPU. Any
  pretrained weights.

## How it works

Seven stages, each individually switchable, with every switch chosen by nested
leave-one-out cross-validation over the *source* lines only:

1. context-weighted consensus over source lines, weighted by control-profile
   similarity on cross-line-variable genes
2. perturbation-neighbour smoothing in effect space
3. reliability shrinkage from cross-line reproducibility
4. program-basis (SVD) denoising, blended
5. **magnitude renormalisation** — restores the across-perturbation magnitude
   spread that stages 2–4 flatten
6. per-gene on-target knockdown, applied in count space
7. one global calibration scale

A knockdown absent from every source line is routed through gene-response
similarity: the silenced gene is still *measured*, so its behaviour across the
rest of the knockdown panel identifies functionally related genes whose own
knockdowns stand in for it. No measurement of the gene-as-perturbed is used.

Roughly 2,000 lines of NumPy/SciPy. No neural network. This is a deliberate
choice, not a limitation of effort: Arc's own 2025 wrap-up concluded that
"purely AI-based approaches did not consistently outperform statistical
baselines", and three of four prize-winning entries were summary-level
statistical transfer with context conditioning.

## Training data

Four matched CRISPRi Perturb-seq contexts — K562 and RPE1 (Replogle 2022),
Jurkat and HepG2 (Nadig 2025) — sharing 6,642 genes and 2,053 knockdowns. One
library design, one processing pipeline, so a knockdown means the same thing in
every line. See `MANIFEST.md`.

Per fold the model trains on three lines and is evaluated on the fourth.

## Evaluation

Leave-one-cell-line-out, two regimes: `context` (knockdown seen in source lines,
never in the target — the challenge's setting) and `double` (knockdown deleted
from every line). Metrics are a port of Arc's `cell-eval` 0.8.2. Numbers are in
`../RESULTS.md` and `results/`.

Baselines: return the control profile unchanged (`Δ=0`); predict the mean
perturbation response for everything (`cell-eval`'s own reference model); and
unweighted cross-line transfer, which is the genuinely hard comparator.

## Known limitations

- **It transfers the shared component.** On these four lines the across-line
  mean effect is 43% of total effect variance; the remaining 57% is
  context-specific or noise, and this model does not predict it. The obvious
  mechanism for context-specificity — scaling by the target gene's expression in
  the new line — has **no signal** here (Spearman −0.01).
- **It predicts pseudobulk, not single cells.** The challenge scores DE from
  submitted cells; this produces a mean profile. Submitting to the real
  challenge needs a sampling layer.
- **The metrics conflict.** Discrimination is an L1 retrieval metric and rewards
  large, well-spread effects; MAE rewards small ones. Which model looks best
  depends on which scoring is quoted, so both the clipped (leaderboard) and
  unclipped (balanced) aggregates are always reported.
- **Thin data.** 45–120 cells and ~11–14k UMI per perturbation, versus the
  challenge's ~1,000 cells and >50k UMI. Absolute scores are not comparable to
  leaderboard numbers.
- **Essential genes only.** All four screens target essential genes. Behaviour
  on non-essential or regulatory targets is untested.
- **Four contexts is few.** Context weights are fitted from three source lines
  per fold, which is barely enough to be more than a heuristic.

## What would improve it

In descending order of expected value:

1. **More contexts.** Add the VCC 2025 H1 hESC data (Challenge 2 explicitly
   permits it) as a fifth line, and Tahoe-100M's ~50 lines for a real context
   encoder. The context ablation in `RESULTS.md` measures what each added
   context is worth.
2. **A pretrained gene/protein embedding.** ESM-2 embeddings were the cheapest
   consistent gain among 2025 entrants, and would replace the data-derived
   gene-similarity routing used for unseen knockdowns.
3. **A single-cell output layer**, to make DE scoring native rather than
   approximated.
4. **Deeper source data** — the K562 genome-wide screen is already downloaded
   and would widen the perturbation panel considerably.

## Intended use

Research benchmarking of cross-context perturbation prediction. Not validated
for, and should not be used for, any clinical or diagnostic purpose. Predicted
knockdown effects are statistical extrapolations from three cell lines, not
experimental measurements, and the model states no uncertainty on individual
gene predictions.

## Reproducing

```bash
./fetch_data.sh
cd ../.. && python -m virtualcell.run --regime context
```

Deterministic given the same inputs and seed: the only randomness is the
evaluation subsample (seeded) and the randomized SVD (fixed `random_state`).
