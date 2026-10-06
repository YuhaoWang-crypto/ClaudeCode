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
| `tyrrs_dataset.py` | builds a labelled MjTyrRS variant × ncAA dataset from the clone tables embedded as **images** in sheet 2 |
| `baseline.py` | the sequence-only baseline on that dataset, with the confounds controlled |
| `data/` | ranked SSM outputs, the parsed literature table, the TyrRS variant table |

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

---

# The MjTyrRS dataset and its sequence-only baseline

```bash
python3 pylrs/literature.py --xlsx dataset_from_paper_v2.0.xlsx --out pylrs/data
python3 pylrs/tyrrs_dataset.py        # -> data/tyrrs_variants.csv
python3 pylrs/baseline.py             # needs fair-esm + torch
```

## ⚠️ A two-character bug that would have invalidated everything

The sequence cells are prefixed with a label and a colon — `WT:MDEFEM…`,
`M.barkeri:MDKK…`. Stripping non-letters without removing the label glues it onto
the N-terminus, shifting every residue number by two (or by eight for
`M.barkeri`). Every mutation the literature cites then lands on the wrong residue,
silently.

`literature.verify_numbering` checks all 17 cited MjTyrRS positions against the
sequence and refuses to continue if any disagree. ✅ With labels stripped, all 17
match (Y32, L65, A67, H70, Y102, V103, E107, F108, Q109, L110, Y114, Q155, D158,
I159, H160, Y161, L162) and the scaffolds come out at their correct lengths —
MjTyrRS 306 aa, MbPylRS 419, McTyrRS 319, EcTyrRS 424.

## The dataset: the clone tables are images

The campaign rows hold prose; the actual clone tables are **21 embedded images**.
`tyrrs_dataset.py` carries a hand transcription with the source image recorded per
campaign, giving **62 positive clones (58 unique mutation sets) across 12 ncAAs**
on one MjTyrRS scaffold:

| ncAA | clones | | ncAA | clones |
|---|---|---|---|---|
| pBpa | 6 | | BipAla | 7 |
| **pAzF** | **7** | | BpyAla | 2 |
| pPRF (alkyne) | 8 | | pBoroPhe | 7 (3 unique) |
| HQ-Ala | 8 | | pCMF | 5 |
| 2-NPA | 6 | | OAY | 4 |
| pIF / pAF | 1 + 1 | | | |

Mutated positions, counted over all positives: 32 (61×), 158 (62×), 162 (50×),
65 (35×), 159 (35×), 107 (27×), 108 (24×), 109 (22×), 70 (20×), then 67, 114,
155, 102, 103, 110, 160, 161.

Two campaigns are deliberately not transcribed: 2.9 (OCF₃Phe) is a dense 14-clone
alignment image where a misread residue would be silent, and 2.15 (Nal) reports
only IC₅₀ values with no sequences.

⚠️ **There are no measured negatives.** These papers report survivors only, so the
data is positive-unlabelled. `decoys()` samples random members of each campaign's
own randomised library as presumed negatives — defensible, but assumed.

## ✅ What a sequence-only model can do

**Task A — ncAA attribution, decoy-free.** Given a selected clone, which of the 12
ncAAs was it selected for? Leave-one-clone-out over the 58 unique variants, using
only measured data.

| scope | features | top-1 | MRR |
|---|---|---|---|
| free | chance (stratified) | 0.105 | — |
| free | PSSM log-odds | 0.517 | 0.679 |
| free | ESM-2 35M Δ-embedding | 0.569 | 0.703 |
| free | one-hot of hotspot residues | 0.621 | 0.720 |
| free | **ESM-2 650M Δ-embedding** | **0.672** | **0.758** |
| controlled | one-hot of hotspot residues | 0.466 | 0.576 |
| controlled | PSSM log-odds | 0.155 | 0.394 |

**Sequence alone names the right ncAA 67% of the time against a 10.5% chance
rate — a 6.4× lift.** So active-site sequence does carry ncAA-specific
information; this is a genuinely better starting point than the PylRS model's
cross-substrate behaviour.

Two caveats that matter:

- ⚠️ **Model size decides whether the PLM is worth it.** ESM-2 **35M scores below
  a plain one-hot** of the hotspot residues (0.569 vs 0.621); only the 650M model
  pulls ahead (0.672). A small PLM here is worse than counting residues.
- ⚠️ **Part of the lift is library design, not biology.** Campaigns randomised
  different position sets, so "which positions are mutated" leaks the ncAA. The
  *controlled* scope uses only the three positions every campaign randomised
  (32, 158, 162); one-hot still reaches 0.466, i.e. **4.4× chance from three
  residues**. Real signal, but about a third of the free-scope lift is design
  artifact.

## ❌ What the impressive-looking number actually measures

**Task B — selected clones vs library decoys**, leave-one-ncAA-out:

| scope | features | mean AUC | mean AP |
|---|---|---|---|
| free | **PSSM log-odds** | **0.980** | 0.865 |
| free | ESM-2 35M | 0.939 | 0.628 |
| free | one-hot | 0.882 | 0.521 |
| controlled | PSSM log-odds | 0.962 | 0.711 |

AUC 0.98 looks like a working classifier. It is not. The simplest possible
model — counting which residues are common among selected clones — **wins**, and
the negatives were never assayed. The task measures "does this look like
something selection would keep", which is a property of the decoy distribution,
not of substrate recognition. Quote Task A, not Task B.

## Where this leaves the plan

The sequence-only baseline is **worth building on**: 6.4× chance on a decoy-free
task, from 58 variants, with no structure and no Rosetta. The next honest step is
more measured data rather than more model — specifically **measured negatives**,
which no published campaign in this sheet provides and which currently cap what
any classifier here can claim.
