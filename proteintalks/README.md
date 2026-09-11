# ProteinTalks reproduction and comparison

An independent reproduction of **ProteinTalks**, the perturbation-proteomics
"virtual cell" model of Sun, Qian, Li et al.
([bioRxiv 2025.02.07.637070](https://www.biorxiv.org/content/10.1101/2025.02.07.637070v1);
published in *Nature* as "An operational perturbation proteomics-based virtual
cell model", [doi:10.1038/s41586-026-11001-9](https://doi.org/10.1038/s41586-026-11001-9)),
plus a comparison against the wider virtual-cell field.

## Headline results

**The architecture reproduces exactly.** `proteintalks/official.py` is a port of
the authors' released `ppODE`. Loaded with the published checkpoint and fed
identical inputs, it matches the reference implementation to **1.8e-7** on the
proteome output and 5.4e-25 on the phenotype output, with all 30 weight tensors
matching in shape.

**The model is 790,306 parameters, and 98.3% of them are in the phenotype
head.** The neural ODE that gives the model its name is 1.7% of it. A single
layer, `drugsens_conv1` (32 × 5585 × 4), is 90.5%.

**The released code differs from the published equations in ten places.** Two
matter: the dynamics module has **no protein-protein coupling whatsoever** (all
its convolutions use `kernel_size=1`, so each protein's trajectory depends only
on its own baseline and its own perturbation entry), and the ODE integrates over
**uniform integer ticks 0,1,2,3** rather than over 0/6/24/48 hours.

**A one-line predictor that uses no proteomics beats every baseline the paper
reports.** On the paper's own openly published label matrix, predicting each
drug's average efficacy rate across training cell lines scores AUROC 0.919 in
the leave-one-cell-line-out setting, ahead of KNN (0.805), DeepSynergy (0.806),
GeneCompass (0.751), Geneformer (0.684) and UCE (0.482). Paired per cell line,
ProteinTalks is the only model significantly *better* than that control, by
0.016 AUROC (p = 0.015), and linear regression is statistically tied with it.

**Removing the dynamics module entirely improves classification on our
simulator**, and the module's proteome predictions are worse than asserting that
nothing changes. That result is from synthetic data at a reduced training budget
and is suggestive rather than conclusive, but the paper reports no ablation of
its own and 98.3% of the model's parameters sit outside the dynamics module.

**Trying to improve the model mostly failed, and the control run found something
bigger.** Five targeted changes (residual decoding, real elapsed time, the
`f(z,t,D)` field the Supplementary Information specifies, protein coupling, a
low-rank head) moved AUROC from 0.918 to at best 0.937 against a no-proteomics
control at 0.908. None of them, and not the released architecture either, learns
the *direction* of the perturbation response: delta correlation stays within
[-0.088, +0.040]. A positive control then showed the response is trivially
learnable on the same inputs and split — **ridge regression reaches delta
correlation 0.949 and recovers 97% of the response variance**. The benchmark is
fine; the architecture does not find signal that a linear map finds almost
perfectly. See [`docs/OPTIMIZATION.md`](docs/OPTIMIZATION.md).

**On unseen drugs the paper's own supplement reports no significant advantage.**
Table S6: ProteinTalks AUROC 0.638 vs random forest 0.648 (p = 0.80) and linear
regression 0.669. That is the setting that matters for prospective drug
discovery.

![Leave-one-cell-line-out AUROC for every model the paper benchmarks, against a
drug-mean control that uses no proteomics](figures/cellline_auroc.png)

Full detail with reproduction commands: [`docs/REPRODUCTION_STATUS.md`](docs/REPRODUCTION_STATUS.md).
Comparison against scGPT, Geneformer, UCE, GEARS, CPA, STATE and others, with
the 2025–26 benchmarking literature: [`docs/COMPARISON.md`](docs/COMPARISON.md).

## What the model is

Two modules. Module 1 takes a baseline proteome plus a perturbation vector and
integrates a learned vector field with RK4 to predict the proteome at 6, 24 and
48 hours. Module 2 reads that predicted trajectory together with 935-dimensional
drug descriptors and classifies drug efficacy or combination synergy. Trained
jointly under a gradient-conflict-aware multi-task objective.

## Layout

```
proteintalks/
  official.py     exact port of the released reference code, checkpoint-verified
  model.py        the paper's Equations 1-11, with ablation backbones
  multitask.py    Equations 12-18: gradient-cosine task weighting
  data.py         dataset container, mechanistic simulator, real-data loader,
                  and the paper's three evaluation splits
  train.py        joint training loop with early stopping
  evaluate.py     AUROC/AUPRC/accuracy, trajectory metrics, bootstrap CIs
  baselines.py    the paper's comparators + the trivial controls it omits
  interpret.py    SHAP prioritisation + a label-permutation null
  improved.py     ProteinTalks-R: five gated changes, each measurement-motivated
scripts/
  verify_checkpoint.py     shape-match the port against the released weights
  equivalence_check.py     numerical equivalence vs the reference, + param census
  real_label_controls.py   no-proteome controls on the real label matrix
  paired_cellline_test.py  paired Wilcoxon vs the paper's per-cell-line results
  calibrate_simulator.py   check the simulator's label structure vs the real one
  run_experiments.py       settings 1/2/3 and architecture ablations
  optimize_model.py        ablation ladder over the five proposed changes
  verify_ladder_baseline.py  proves rung 0 is exactly the released model
  delta_learnability_control.py  positive control: is the response learnable?
docs/
  REPRODUCTION_STATUS.md   what was verified, what was not, with numbers
  COMPARISON.md            ProteinTalks vs the virtual-cell field
  ARCHITECTURE_NOTES.md    paper equations vs what had to be inferred
  DATA.md                  the corpus, and what is retrievable
  OPTIMIZATION.md          the improvement attempt, and why the control matters
  PTDS_ACCESS_REQUEST.md   what the gated-matrix application requires
```

## Install and run

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt

# Reproduce the architecture verification (needs the released repo + weights)
git clone https://github.com/guomics-lab/PTV-1
python scripts/verify_checkpoint.py --checkpoint PTV-1/ProteinTalks/best_checkpoint.pth
python scripts/equivalence_check.py --repo PTV-1 \
       --checkpoint PTV-1/ProteinTalks/best_checkpoint.pth

# Reproduce the control experiments (needs the open Nature supplementary tables)
python scripts/real_label_controls.py  --table-s1 TableS1_drug_cell_0702.xlsx
python scripts/paired_cellline_test.py --table-s1 TableS1_drug_cell_0702.xlsx \
                                       --table-s5 TableS5_cross_cell_final_260823_ESM.xlsx

# Architecture ablations on the simulator
python scripts/run_experiments.py --quick
```

## Data

The full pretraining corpus (16,311 DIA-MS runs, 5,583 protein groups, 325 MB)
is gated behind an application with reviewer approval, and the 5,585-protein
index the released checkpoint expects is not published anywhere we could find,
so the trained weights cannot be applied to a new proteome. Openly downloadable
without login: the model code and weights (MIT), the 1,116-condition efficacy
label matrix, 63 drug SMILES, a 2,487 × 3,631 multi-timepoint protein matrix,
a 501-patient matrix with survival, and the full per-cell-line and per-drug
benchmark results. See [`docs/DATA.md`](docs/DATA.md).

Because the pretraining matrix is unavailable, `data.py` also ships a
mechanistic simulator whose ground truth is a known ODE on a sparse protein
network. It is for testing the architecture, not for any biological claim, and
`scripts/calibrate_simulator.py` checks its label structure against the real
one.

## Licence

MIT. Contains no data from the original study.
