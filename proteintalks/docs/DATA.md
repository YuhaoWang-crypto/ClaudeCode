# Data: what the paper used, and what this reproduction can actually run on

## What the paper used

| Asset | Scale |
|---|---|
| DIA-MS runs | 16,311 (after discarding 689 with <1,000 identifications) |
| Protein measurements | >38 million |
| Protein groups quantified | 5,530 groups / 5,143 unique proteins (model input: 5,585) |
| Cell lines | 18 breast cancer lines (16 TNBC, 2 non-TNBC) |
| Drugs | 63 FDA-approved, plus 98 additional anti-cancer compounds for held-out testing |
| Timepoints | 0, 6, 24, 48 h, in triplicate |
| Drug combinations | 914 combination–cell-line tuples |
| Cytotoxicity datapoints | 23,482 |
| Instrument time | ~6,668 h on a TripleTOF 5600+ |
| Missingness | 51.7% overall, imputed at 0.8 × minimum detected intensity |

Modelled conditions after averaging replicates: 1,529 (1,070 train / 305 val /
154 test in Setting 1).

## The reproducibility problem

That corpus is the entire contribution. It took roughly nine months of
continuous mass-spectrometer time to generate, and it cannot be regenerated,
approximated, or substituted. A reproduction that lacks it is not reproducing
the result — it is reproducing the method.

This matters more here than for a transcriptomics foundation model, where public
corpora (CELLxGENE, Human Cell Atlas, Tahoe-100M) let an independent group
retrain from scratch. There is no public perturbation-proteomics corpus at this
scale to fall back on.

## Where the real data lives

- **Raw MS data**: ProteomeXchange via iProX, accession `IPX0007409000`.
- **Processed datasets (PTDS)**: `db.prottalks.com`.
- **Code**: see the availability statement of the Nature version.

See `docs/REPRODUCTION_STATUS.md` for what was and was not retrievable during
this reproduction attempt, with observed HTTP status codes.

## Using the real data with this code

`proteintalks.data.load_real` accepts:

```
protein_matrix.csv     cell_line, drug, time, protein, value
labels.csv             cell_line, drug, label          (1 = effective)
drug_features.csv      drug, f0 ... f934               (881 fingerprint bits +
                                                        54 physicochemical)
```

It applies the paper's preprocessing: imputation at 0.8 × the minimum detected
intensity, then per-sample min–max normalisation. Then:

```python
from proteintalks import load_real, make_splits, train_proteintalks
ds = load_real("protein_matrix.csv", "labels.csv", "drug_features.csv")
tr, va, te = make_splits(ds, setting=1, seed=0)
model, hist = train_proteintalks(ds, tr, va)
```

Nothing else in the pipeline changes. The simulator and the real loader emit the
same `PerturbationDataset` container.

## The simulator, and what it is for

`SyntheticPerturbationProteome` generates perturbation proteomics from a known
ODE:

```
dP/dt = -alpha (P - P0_l) + beta W tanh(P - P0_l) - gamma u_d(P)
```

with `W` a sparse signed protein-interaction matrix with a scale-free degree
distribution, `P0_l` a low-rank cell-line-specific baseline, and `u_d` drug
inhibition on that drug's targets. Efficacy is thresholded on the activity of an
essential module at 48 h, so the label depends on the trajectory rather than on
the baseline alone. Drug fingerprints are generated from a latent
mechanism-of-action class, so chemically similar drugs share targets — the
property that makes leave-one-drug-out generalisation possible at all.

**What it can establish**: whether the architecture trains, whether the neural
ODE inductive bias helps when the ground truth genuinely is an ODE, how the
three evaluation settings behave, and how the reported metrics compare against
trivial controls under a *known* data-generating process.

**What it cannot establish**: any biological claim in the paper. It says nothing
about whether ProteinTalks predicts real drug efficacy at AUROC 0.960, whether
the SHAP-prioritised proteins are real resistance markers, or whether the
predicted synergies validate in cells. Those claims rest on the real corpus and
on the wet-lab validation the authors report, and no simulation substitutes for
either.

The simulator is deliberately generous to the model: the ground truth is exactly
the model class the architecture assumes. Results on it are an upper bound on
the architecture's advantage, not an estimate of it.
