---
name: enzyme-combinatorial-design
description: >-
  Design improved enzyme variants by recombining measured single mutations:
  FFT-PLSR (AAindex protein spectra + partial least squares) over a combinatorial
  mutation space, with the baselines, held-out tests and numerical guardrails that
  decide whether to trust the ranking. Use when turning measured single-mutant
  activity data into a ranked pick-list of multi-mutant variants; running a
  round-based directed-evolution campaign (nominate sites -> assay singles ->
  recombine -> assay -> retrain); reproducing or extending the PylRS /
  pyrrolysyl-tRNA synthetase FP4COM work; or deciding whether a sequence-encoding
  ML model beats a plain additive model on a small mutant dataset. Enforces
  ✅-validated vs ⚠️-unvalidated labeling on every number it reports.
---

# Combinatorial enzyme design with FFT-PLSR

A reusable methodology and a working `fftplsr` package for the specific, common
situation in enzyme engineering: **you have measured a handful of single mutants,
and you need to decide which combinations to build next.**

Method source: Hu *et al.*, *Nat Commun* **16** (2025),
doi:[10.1038/s41467-025-61952-2](https://doi.org/10.1038/s41467-025-61952-2);
reference code <https://github.com/zjuhaoran/FPFORCOM> (MIT).

## What this does and does not do

| | |
|---|---|
| ✅ **Does** | rank the 2^k combinations of k measured single mutations; pick a diversity-filtered order list; run the round-over-round retrain loop |
| ✅ **Does** | tell you, on your own data, whether the ML model beats a plain additive model — and say so when it doesn't |
| ❌ **Does not** | design enzymes de novo, invent new active sites, or predict which *positions* to mutate — that needs structure/PLM tools (see **Composing**) |
| ❌ **Does not** | predict absolute activity, or transfer across backgrounds or assays — fitness is a fold-change relative to one parent in one assay |
| ❌ **Does not** | work without measured training data. It is a regression model, not a zero-shot predictor |

The hard requirement: **measured activities for single mutants on one parent
background, in one assay, with the parent itself measured as the 1.0 reference.**

## The method in five steps

1. Map each residue to a scalar via an **AAindex** entry (566 available) → `x[0..L-1]`
2. **Subtract the mean** of `x` (this zeroes the DC term — see the trap below)
3. **FFT**, take the magnitude spectrum, normalise by its maximum, keep the first
   `floor(L/2)` bins → the "protein spectrum"
4. **Greedily select** AAindex entries by out-of-fold error (1 or 3 entries,
   concatenated along the feature axis)
5. **PLS regression** (`n_components` chosen by cross-validation) → score the whole
   combinatorial space

The non-obvious part is *why this can model epistasis at all*: steps 1–3 are linear
in the residue values, but the magnitude `|·|` is not, so the features are a
nonlinear function of the mutation pattern. That is the entire mechanism — and it
also bounds what the model can do (see **Is it worth it?**).

## ⚠️ The trap: drop FFT bin 0

Step 2 sets the DC term to **exactly zero**, so bin 0 of the spectrum contains only
floating-point round-off (~1e-16). That is harmless until the regressor standardises
its columns: `PLSRegression(scale=True)` divides each feature by its standard
deviation, and for bin 0 that deviation *is* the round-off — so the column is
amplified to O(1) and PLS fits noise through it.

✅ **Measured consequence** (`python3 -m fftplsr.m2_dc_artifact`), on the paper's
round-1 training set: re-drawing bin 0 anywhere inside its round-off range — which
is all that changing FFT backend, precision or summation order does — moves the
cross-validated error at 10 components across **0.29 to 8.05**, and changes which
of the 566 descriptors wins. With bin 0 dropped, the same perturbation moves it by
**7.7e-14**. The reference implementation keeps bin 0 and stores features as
float32; its published round-1 score (cvMSE 0.359 at k=10) sits inside that
lottery.

`fftplsr.encode.spectrum` therefore defaults to `drop_dc=True`. **Keep it on.**
`drop_dc=False` exists only to reproduce the artifact.

## ⚠️ The second trap: components vs. training-set size

The reference `regscore` reported `predict(X)` scored against the same `y` it was
fit on, i.e. **training-set** MSE/R². With 13 samples and 227 features those
numbers run to R² = 0.999 and mean nothing. `fftplsr` reports out-of-fold scores as
`cvR2`/`cvMSE` and keeps in-sample ones as `train_r2`/`train_mse`.

Related: selection pressure pushes `n_components` toward its ceiling because more
components fit more noise. ✅ The paper's round-1 model used **10 PLS components
from 12 training points**. Cap components well below the sample count; prefer
`k ≤ n/4` unless cross-validation argues hard otherwise.

## Is it worth it? Run the baselines first

In a recombination task only a handful of positions ever change, so the 227-bin
spectrum is a deterministic function of at most k bits. **The spectrum cannot encode
anything the mutation indicators do not.** So always compare against them.

✅ Held-out results on the paper's own prospective splits
(`python3 -m fftplsr.m3_baselines`), R² / mean fitness of the top-8 picks:

| split | FFT-PLSR | one-hot PLS | log-additive | oracle top-8 |
|---|---|---|---|---|
| 13 singles → 25 combos | **0.634** / 5.86 | 0.539 / **6.11** | 0.309 / **6.11** | 6.39 |
| 38 → 64 later variants | 0.750 / 8.49 | **0.827** / **8.76** | −1.652 / 7.08 | 8.88 |

Read this carefully, because it is the operational finding:

- FFT-PLSR wins on R² only in the **smallest** split (13 training points), and even
  there a one-hot model picks **better variants**.
- With 38 training points, a plain ridge/PLS on mutation indicators beats it on both
  metrics, deterministically and in milliseconds.
- `log-additive` (just multiply the measured fold-changes) collapses on
  high-order combinations — multiplying 7 fold-changes overshoots badly. It is a
  fine null for doubles and triples, not for order ≥5.

⚠️ Scope: one enzyme, one assay, held-out n = 25/64/21. This is not a general
verdict on FFT-PLSR — it is the reason to run `evaluate_against_baselines` on
*your* data before trusting a ranking. `design_round` does this automatically and
annotates the report when the ML model fails to beat the best baseline.

### Where the FFT encoding does earn its keep: the cold start

✅ `python3 -m fftplsr.m4_design` reconstructs the position just before the paper's
final round — parent Com1-IFRS, 27 measured saturation singles over 6 positions,
11,492 recombinations (the paper's own 11,520-variant space) — and asks where the
eventual winner `N7Y/H63L/K67N/V74W` (2.75× parent) lands:

| round | training data | FFT-PLSR cvR² | best baseline cvR² | winner's rank |
|---|---|---|---|---|
| A | 27 singles only | **0.336** | −0.075 (`mean`) | 1853 / 11519 |
| B | + 92 measured doubles | **0.634** | 0.405 (`onehot-pls`) | **332 / 11519** |

Three things to take from this:

- **Round A is the encoding's real niche.** With singles only, every substitution
  is seen exactly once, so a one-hot model holding one out has no column for it and
  falls back to the intercept — *every* baseline scores below zero. FFT-PLSR is the
  only model with signal, because the spectrum shares information across
  substitutions. This is the case where the Fourier detour is worth it.
- **Doubles are what makes the winner reachable.** Adding 92 measured doubles
  moves the winner from rank 1853 to 332 — a 5.6× improvement, and the point at
  which epistasis first becomes visible to the model.
- ⚠️ **But top-8 would still have missed it.** Rank 332 of 11,519 is the top 2.9% —
  a ~35× enrichment over chance, not an oracle. Treat the output as a shortlist to
  assay, size the order list to the enrichment, and expect to need a second round.

## What reproduces, and what doesn't

✅ `python3 -m fftplsr.m1_reproduce`

| claim | status |
|---|---|
| Site-triage descriptor triple `QIAN880114_OOBM770105_QIAN880125` | ✅ reproduces exactly (1,698 candidate screens, same 3 in the same order) |
| Round-2 prospective accuracy, paper R² = 0.835 | ✅ R² = **0.833** with the paper's descriptor (`RADA880104`, k=4) |
| Round-2 descriptor choice | ⚠️ `RADA880104` ranks 3rd of 553 here; the top 3 are within 3% cvMSE — a near-tie, not a determination |
| Round-1 descriptor choice + cvMSE 0.359 | ❌ does not reproduce; it is the DC-bin artifact above. Honest round-1 held-out R² is **0.63**, not 0.84 |
| Final-round model `AVBF000109_JUNJ780101_JUKT750101` | ❌ **not rebuildable at all** from current data — see below |

The pattern: **the pipeline reproduces where the model is well-conditioned (few
components, enough samples) and is a coin-flip where it is not.**

### ⚠️ Third trap: the `aaindex` package changed under the method

All 566 AAindex1 entries carry complete residue values in `aaindex` **1.0.5** (the
version the paper pinned). In **1.3.2**, 13 of them return `None` for at least one
standard residue: `AVBF000101`–`AVBF000109`, `GUYH850103`, `ROSM880104`,
`ROSM880105`, `YANJ020101`.

`AVBF000109` is one of them, and it is the first descriptor of the paper's
final-round model. That model therefore cannot be rebuilt from current `aaindex`
data — not approximately, but not at all. Pin `aaindex==1.0.5` if you need it.

`fftplsr` makes this visible instead of silent: `available_indices()` returns the
553 usable entries (not 566), `encode` raises `IncompleteIndexError` naming the
missing residue, and `screen_indices` prints how many candidates it skipped and
why. The reference implementation's bare `except ValueError` scored such a
candidate as `-100` and moved on.

## Run it

```bash
pip install numpy scipy pandas scikit-learn aaindex joblib
python3 -m fftplsr.m1_reproduce    # reproduce the paper's rounds
python3 -m fftplsr.m2_dc_artifact  # the DC-bin finding
python3 -m fftplsr.m3_baselines    # is it worth it, on held-out data
python3 -m fftplsr.m4_design       # design a new variant end to end
python3 -m pytest tests/test_fftplsr.py -q
```

## Designing on a new enzyme

```python
from fftplsr import design

report = design.design_round(
    parent=my_sequence,                  # the background every label is relative to
    measured={"WT": 1.0, "D2N": 3.62, "H62Y": 1.97, ...},
    sites=design.improved_sites(measured, threshold=1.05),
    n_rounds=1,       # AAindex entries to select greedily; 3 for larger training sets
    cv=None,          # leave-one-out; pass an int for k-fold once n > ~50
    pick=8,           # size of the order list
    max_jaccard=0.6,  # reject near-duplicate picks
)
print(report.summary())   # headline + baselines + pick-list + caveats
```

Then read `report.summary()` in this order:

1. **`report.baselines`** — if nothing beats `mean`, you have no signal; stop and
   get more data. If a one-hot model wins, use it and skip the FFT.
2. **`report.selection.cv_r2`** vs the best baseline — the headline states both.
3. **`report.notes`** — flags raised automatically (ML lost to a baseline; scored
   variants that already have measurements).
4. **`report.picks`** — the order list. `epistasis` is predicted minus
   log-additive: large values are where the model claims to add something, and are
   the informative ones to assay whether or not they rank top.

Operational defaults that matter:

- **`improved_sites(threshold=1.05)`** — recombine only singles that beat the
  parent. ⚠️ **Prune positions, never substitutions within a position.** A single
  mutation's solo effect is a poor guide to its value in combination — which is
  the entire reason you are fitting an epistasis-aware model. ✅ Measured on this
  dataset: the best variant in the released panel, Com2-IFRS =
  `N7Y/H63L/K67N/V74W` at 2.75× parent, is built from singles ranking **3rd of 3,
  1st of 7, 4th of 5 and 2nd of 2** at their positions — three of four in the
  bottom half. A `per_position=3` filter drops K67N and makes the answer
  unreachable before the model is fitted. Keeping everything above threshold
  reproduces the paper's own 11,520-variant space exactly.
- **Diversity filter** — a raw top-8 from PLS is typically eight spellings of one
  mutation set. `max_jaccard=0.6` spends the assay budget on distinguishable
  hypotheses. Lower it to diversify harder.
- **Round over round** — after assaying, add the new measurements and re-run.
  ✅ Observed in the paper's data: going from 28 singles to 120 singles+doubles is
  what made the final winner reachable, because doubles are where epistasis first
  becomes visible.

## Composing with other skills

FFT-PLSR only recombines mutations someone already nominated and measured. The
upstream steps belong to other skills in this repo:

| Need | Use |
|---|---|
| Which **positions** to mutate, with no data yet | `enzyme-mutation-ranking` (ESM-2/ESMC + ProteinMPNN consensus) — the paper's own round 3 used ESM-1v, MutCompute and ProRefiner for exactly this |
| Per-residue **dynamics / allostery** to target | `anthropic-skills:dyna1-dynamics-vs-static` |
| **Structure** of a designed variant, or ligand/substrate binding | `anthropic-skills:boltz-denovo-design`, `anthropic-skills:protein-ligand-md` |
| Turning winners into **orderable DNA** | `anthropic-skills:codon-optimize-qc`, `anthropic-skills:synbio-cassette-designer` |

The honest pipeline is: PLM/structure tools nominate positions → **assay singles**
→ FFT-PLSR (this skill) recombines → assay → retrain. The assay steps are not
optional; there is no path from this skill to a variant without them.

## Package layout

| Module | Role |
|---|---|
| `fftplsr/variants.py` | mutation-string parsing, space enumeration, saturation scans; raises on labels inconsistent with the parent |
| `fftplsr/encode.py` | AAindex + FFT spectra. `MutationEncoder` is an exact rank-k spectrum update — same features as the reference path, fast enough for 2^20 variants and 566-descriptor screens |
| `fftplsr/model.py` | PLS with out-of-fold scoring, greedy descriptor selection, `nested_cv_r2` for selection bias, `pls_path_predict` (all component counts from one fit) |
| `fftplsr/baselines.py` | mean, log-additive, one-hot ridge/PLS, pairwise ridge |
| `fftplsr/design.py` | `design_round`, `improved_sites`, diversity filter, automatic baseline comparison |
| `fftplsr/datasets.py` | the paper's data, vendored; see `data/PROVENANCE.md` |
| `fftplsr/m1..m4` | reproduction, the DC artifact, baselines, a worked design round |

## Reporting rules

Carry the repo's honesty convention into every summary:

- ✅ for numbers this pipeline **computed** — always say on what data, with what
  held-out or cross-validated split, and with what n.
- ⚠️ for anything extrapolated: transfer to another background, another assay,
  another enzyme, or a ranking whose model did not beat its baselines.
- ❌ for claims that did not reproduce. Say so plainly and give the number that did.
- Never report an in-sample R² as model accuracy. If a source quotes one (this
  method's reference implementation does), say which it is.
- A pick-list is a hypothesis list. Report predicted fitness as a **ranking**, not
  as an expected fold-change, unless a held-out split supports the calibration.
