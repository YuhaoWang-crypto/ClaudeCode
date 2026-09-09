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

That corpus is the entire contribution. It took roughly 6,668 hours of
mass-spectrometer time to generate, and it cannot be regenerated, approximated,
or substituted. A reproduction that lacks it is not reproducing the result — it
is reproducing the method.

This matters more here than for a transcriptomics foundation model, where public
corpora (CELLxGENE, Human Cell Atlas, Tahoe-100M) let an independent group
retrain from scratch. There is no public perturbation-proteomics corpus at this
scale to fall back on.

## Where the real data lives

**Open, no login** (all verified retrievable):

- **Code and trained weights**: `github.com/guomics-lab/PTV-1`, MIT licence.
  `ProteinTalks/best_checkpoint.pth` is 10.1 MB and loads cleanly.
- **Nature supplementary tables**, a single 126 MB ZIP, `41586_2026_11001_MOESM3_ESM.zip`
  on `media.springernature.com`. Contains the 1,116-condition efficacy label
  matrix (Table S1 `C_Efficacy`), 63 drug SMILES and targets (Table S1
  `B_Drugs`), a 2,487 × 3,631 multi-timepoint protein matrix across 11
  timepoints and 3 cell lines (Table S2), a 501-patient FFPE matrix with
  survival (Table S13), 3,000 compound SMILES (Table S14), and the complete
  per-cell-line and per-drug benchmark results (Tables S5, S6).

**Gated**:

- **PTDS protein matrix**, 16,311 samples × 5,583 protein groups, 325 MB, at
  `db.prottalks.com`. Requires an application form with an institutional email
  (the site server-side rejects consumer domains), PI name, lab URL and a stated
  intended use, followed by reviewer approval and an emailed link.

**Asserted but not retrievable**:

- **Raw MS data**, ProteomeXchange via iProX `IPX0007409000`. The PROXI endpoint
  returns all-null fields, ProteomeCentral returns zero datasets for the
  accession, and `download.iprox.cn` returns HTTP 403.

**Missing entirely**: the ordering of the 5,585 protein groups that the released
checkpoint's first phenotype layer expects. Without it the trained weights
cannot be applied to any new proteome. The portal also lists 5,583 protein
groups against the checkpoint's 5,585.

See `docs/REPRODUCTION_STATUS.md` §1 for the full table with observed HTTP
status codes.

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

### Calibration against the real label structure

A simulator whose labels are driven by the wrong factor makes every downstream
comparison meaningless. In the real label matrix, **drug identity carries almost
all of the efficacy signal and cell-line identity carries almost none**. The
simulator is tuned to reproduce that asymmetry, and
`scripts/calibrate_simulator.py` measures it:

| Statistic | Real (Table S1, C_Efficacy) | Simulator (defaults) |
|---|---|---|
| Positive rate | 0.334 | 0.318 |
| Drug-mean AUROC | 0.909 | 0.967 |
| Cell-line-mean AUROC | 0.463 | 0.532 |

The asymmetry is reproduced. The drug-mean AUROC is higher than reality, i.e.
the simulator's drug effects are somewhat more deterministic than real drug
responses, which makes it an easier problem than the real one. This is left as
measured rather than tuned further, because tuning a simulator until it matches
a target metric is a good way to build in whatever conclusion you were after.
