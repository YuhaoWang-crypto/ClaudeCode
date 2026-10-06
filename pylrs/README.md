# pylrs — running and auditing the PylRS-libY model

Tooling for [wendao/PylRS-libY](https://github.com/wendao/PylRS-libY), the
AutoGluon model behind *Machine Learning Assisted Pyrrolysyl-tRNA Synthetase
(PylRS) Design*, plus a parser for a curated aaRS × ncAA literature sheet.

| File | What it does |
|---|---|
| [`ENVIRONMENT.md`](ENVIRONMENT.md) | the working environment recipe (Python 3.8 + AutoGluon 0.4.0 via `uv`, ~2 min) and what still needs a Rosetta licence |
| `run_prediction.py` | re-runs the published TCOY/TCOC SSM predictions and **asserts** they match the upstream Colab output to 6 decimals |
| `audit.py` | recovers the dataset's real structure and cross-validates it under four increasingly honest splits |
| `literature.py` | turns the curated `dataset_from_paper_v2.0.xlsx` reading sheet into a tidy campaign table with a chemotype column |
| `data/` | ranked SSM outputs and the parsed literature table |

```bash
uv venv --python 3.8 venv38
uv pip install --python venv38/bin/python "autogluon.tabular[all]==0.4.0"
venv38/bin/python -W ignore pylrs/run_prediction.py --repo /path/to/PylRS-libY

# the audit needs only numpy + scikit-learn, any modern Python
python3 pylrs/audit.py --repo /path/to/PylRS-libY
```

## ✅ It reproduces exactly

All 20 published TCOY values match digit for digit (`306A` axial P=0.263559
rank 1; `305P` equatorial P=0.253156 rank 2; …). The 2022 pickle loads fine under
`uv`-provisioned Python 3.8 despite a wall of scikit-learn unpickle warnings.

Top-ranked single mutants, for reference:

| target | top 5 (position + substitution / TCO ring conformer) |
|---|---|
| TCOY (on the B1 template) | Y306A/ax, L305P/eq, Y306A/eq, C348K/eq, L309W/eq |
| TCOC (on the 106 template) | L309P/eq, N346I/eq, L309I/eq, L309L/eq, Y306A/eq |

⚠️ **No variant scores above P = 0.5** — max is 0.264 (TCOY) and 0.270 (TCOC)
against a 3% training base rate. This model is a *ranker*; its probabilities are
not calibrated recognition likelihoods and should never be quoted as "26% chance
of working".

## ❌ The published accuracy figure is in-sample

`4.prediction/reference_data/leaderboard.info.txt` reports `score_test` = 1.000
for several models and 0.990 for `WeightedEnsemble_L3`. That column is in-sample:
`3.training/analyzing-autogluon.ipynb` calls `predictor.leaderboard(dataset)` with
`dataset` being the training set itself.

✅ Verified: scoring the shipped predictor on its own 936 training rows gives
AUC **0.9905**, matching the published 0.990481 to four decimals.

An independent cross-validation with comparable models gives:

| split | what it holds out | pooled AUC |
|---|---|---|
| random 5-fold | nothing structural | 0.708 ± 0.024 |
| grouped by variant | all 8 rows of a variant | **0.636 ± 0.013** |

So the honest expectation is **AUC ≈ 0.64**, not 0.99. Note this is a weaker
dataset than the headline suggests: **28 positives in 936 rows**.

## The dataset's real structure (not documented upstream)

936 rows = **117 PylRS variants × 8 ncAAs**, in 8 contiguous blocks of 117.
Column layout, read off `analyzing-autogluon.ipynb`, is *not* the order the
README's "256 + 5 + 1 + 18" implies:

```
cols   0..17    18 Rosetta cartesian_ddg terms   <- the ONLY substrate-aware features
cols  18..22     5 ESM-1v scores                 }
col      23      1 ESM-msa-1b score              } variant-only: identical across
cols  24..279  256 eUniRep (PylRS-fine-tuned)    } all 8 rows of a variant
```

The 8 substrates (`3.training/UAA-info.csv`) are **all O-substituted tyrosines** —
A OMe-Tyr, B OEt-Tyr, C O-propargyl-Tyr, D O-cyclohexyl-Tyr, E O-phenyl-Tyr,
F O-benzyl-Tyr, G O-cyclooctenyl-Tyr (the TCO one), H O-(2-nitrobenzyl)-Tyr.
Positives per substrate: 12, 2, **0**, 1, 5, 2, 3, 3.

The screened libraries span **both** *M. mazei* (82 variants) and *M. barkeri*
(43) PylRS, mutated mostly at Y384 (104×), C348 (100×), V401 (77×), N346 (77×),
A302 (71×), W417 (59×), L309 (31×), Y306 (30×), L305 (6×).

## ⚠️ What it does and does not transfer to

Mean per-fold AUC — "can it rank variants *within* one held-out ncAA":

| features | new ncAA, variants already screened | new ncAA **and** new variants |
|---|---|---|
| 18 Rosetta ddG (substrate-aware) | **0.487** (chance) | 0.369 |
| ESM + eUniRep (variant-only) | **0.704** | 0.360 |
| all 280 | 0.697 | 0.398 |

Read that carefully, because it inverts the paper's "sequence + structure beats
sequence-only" claim *for cross-substrate transfer*: the only substrate-aware
block transfers at chance, while the variant-only block carries the signal. What
generalises is **which PylRS variant is promiscuous at all** — a prior on the
scaffold, not a model of fit to a particular ncAA. And for a genuinely new
variant against a new ncAA, nothing here beats chance.

Practical consequence: the model is reusable as a **promiscuity prior** over the
9 hotspot positions. Treating it as a substrate-specific predictor for an ncAA
outside the O-alkyl-tyrosine chemotype is not supported by its own data.
