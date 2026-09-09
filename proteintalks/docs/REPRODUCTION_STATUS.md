# Reproduction status

What was actually verified, what was not, and what the checks showed.
Everything below is reproducible from this repository plus openly downloadable
files. Dates and HTTP outcomes are as observed on 2026-09-09, the day the
*Nature* version went live.

---

## 1. Availability: better than expected

| Asset | Status |
|---|---|
| Model code | **Open.** `github.com/guomics-lab/PTV-1`, MIT licence, last commit 2026-09-07 |
| Trained weights | **Open.** `ProteinTalks/best_checkpoint.pth`, 10.1 MB, no login |
| Drug efficacy labels (1,116 conditions) | **Open.** Nature Supplementary Table S1, sheet `C_Efficacy` |
| Drug SMILES + targets (63 drugs) | **Open.** Table S1, sheet `B_Drugs` |
| Multi-timepoint proteome matrix (2,487 × 3,631) | **Open.** Table S2 |
| Per-cell-line and per-drug benchmark results | **Open.** Tables S5, S6 |
| 501-patient FFPE matrix + survival | **Open.** Table S13 |
| **Full PTDS matrix (16,311 × 5,583, 325 MB)** | **Gated.** Application form, institutional email, PI name, reviewer approval, link emailed |
| **Protein index for the 5,585 model inputs** | **Not published anywhere found** |
| Raw MS data, iProX `IPX0007409000` | **Asserted, not retrievable.** PROXI returns all-null fields; ProteomeCentral returns 0 datasets; `download.iprox.cn` returns 403 |

The missing protein index is the sharpest practical obstacle. The released
checkpoint's first phenotype-head layer has shape `(32, 5585, 4)`, so its input
is a specific ordering of 5,585 protein groups. Without that ordering the
trained weights cannot be applied to any new proteome. The portal also lists
**5,583** protein groups against the checkpoint's **5,585**, a discrepancy we
could not resolve.

## 2. Architecture: reproduced exactly

`proteintalks/official.py` is a port of the reference `ppODE` class. Validation:

```
$ python scripts/verify_checkpoint.py --checkpoint best_checkpoint.pth
matched shapes : 30      shape mismatch : 0
in ckpt only   : 0       in port only   : 0
RESULT: port is structurally identical to the released model

$ python scripts/equivalence_check.py --repo PTV-1 --checkpoint best_checkpoint.pth
proteome  max |ref - port| = 1.788e-07
phenotype max |ref - port| = 5.428e-25
RESULT: numerically equivalent to the reference implementation
```

Both models are loaded with the *released trained weights* and fed identical
inputs, so this is an equivalence check against the actual published model, not
against a re-initialised copy.

### The RK4 tableau is load-bearing

Reaching that equivalence required one non-obvious fix. The reference integrates
with `torchdyn`'s `solver='rk4'`, which is the **classical** Runge-Kutta tableau.
`torchdiffeq.odeint(..., method='rk4')` implements the **3/8-rule** tableau
instead. Both are fourth-order and agree as the step size shrinks, but this
model takes three steps of size h = 1, and there the two disagree:

| Quantity | Discrepancy |
|---|---|
| ODE trajectory | 0.19% relative |
| Decoded proteome | 2.6% relative (max 3.6%) |

So "rk4", stated on its own as the paper states it, does not pin down the model.
`proteintalks.official.rk4_classical` implements the correct tableau.

## 3. Where the parameters are

This is the most informative single measurement about the architecture:

| Component | Parameters | Share |
|---|---|---|
| Module 1: neural ODE + encoder/decoder | 13,089 | **1.7%** |
| Module 2: phenotype head | 776,897 | **98.3%** |
| — of which `drugsens_conv1` (32 × 5585 × 4) | 714,880 | 90.5% |
| Unused `convdrug1` | 320 | 0.0% |
| **Total** | **790,306** | |

The model is small — roughly 0.8 M parameters against scGPT's ~50 M and UCE's
~650 M — and the paper's claim of far fewer parameters is correct and verified.

But 90.5% of it is a single linear map from 5,585 proteins into a 32-dimensional
feature, and the differential-equation component that gives the model its name
is 1.7% of it. That ratio should temper how the "dynamical foundation model"
framing is read.

## 4. Ten places where the released code differs from the published equations

Verified line by line against `ProteinTalks/model.py`, `config.py` and
`multi_task_learning.py`. Full detail in the docstring of
`proteintalks/official.py`.

1. **`C1`/`C2` use `kernel_size=1`.** Every operation in module 1 —
   `linear_input`, `conv1`, the ODE field, `conv2`, `layer_final` — acts on one
   protein at a time through a shared map. **Module 1 contains no
   protein-protein coupling at all.** Each protein's trajectory is a function of
   its own baseline value and its own perturbation entry. Cross-protein mixing
   exists only in the phenotype head. Kernel 1 is the *right* choice given that
   protein row order is arbitrary, but it means the protein-network dynamics the
   model is described as learning are not represented in the dynamics module.
2. **The ODE integrates over `linspace(0, 3, 4)`** — ticks 0, 1, 2, 3 — not over
   0, 6, 24, 48 hours. The solver never learns that the 24→48 h gap is four
   times the 0→6 h gap. The continuous-time machinery reduces to a uniform
   three-step recurrence.
3. `C3`/`C4` output **32** channels; Eq. 7–8 say 128.
4. `drugsens_conv1` treats **proteins as channels and time as the spatial axis**,
   which is not what Eq. 7 describes.
5. Default `hidden_size` is **64**; Eq. 1 says C1 raises to 128.
6. Default `dropout_rate` is **0.0**; Eq. 3 says 0.1.
7. **LayerNorm and GroupNorm** are used throughout and appear nowhere in the paper.
8. The second ODE layer has **no activation**; Eq. 3 specifies Softplus.
9. **SWAG** is implemented, and Table S5C reports `ppODE (SWA)` and
   `ppODE (Non-SWA)` separately, but SWAG is absent from the Methods.
10. Task weights are **recomputed from `[1-λ, λ]` every step**, not accumulated
    as Eq. 16–18 imply, and gradients are clipped per tensor *before* the cosine
    similarity is measured.

One ambiguity we flagged before seeing the code was resolved in our favour.
Equation 15 sets `adjustment_factor = 0.01 × clipped_similarity`, which is
negative exactly when the gradients conflict, inverting Eq. 17–18 against the
stated intent. The released code uses `abs(similarity)`, matching the prose. We
had implemented the prose reading; it is correct.

## 5. The control the paper does not report

Under the paper's own splitting protocols, on the paper's own openly published
label matrix, how far can you get with **no proteomics, no chemistry and no
model**?

`scripts/real_label_controls.py`, run on Table S1 `C_Efficacy`
(1,116 conditions, 18 cell lines × 62 perturbations, 373 effective / 743 not):

| Setting | Predictor | AUROC | AUPRC | Accuracy |
|---|---|---|---|---|
| **1** random split | drug-mean | 0.909 ± 0.002 | 0.817 | 0.862 |
| | cell-line-mean | 0.463 ± 0.004 | 0.325 | 0.665 |
| | **ProteinTalks (paper)** | **0.960** | 0.854 | 0.910 |
| **2** leave-one-cell-line-out | drug-mean | 0.919 ± 0.026 | 0.859 | 0.860 |
| | best baseline the paper reports (KNN) | 0.805 | 0.631 | 0.817 |
| | **ProteinTalks (paper, Table S5A)** | **0.953** | 0.891 | 0.902 |
| **3** leave-one-drug-out | cell-line-mean | 0.508 ± 0.032 | 0.509 | 0.562 |
| | best baseline the paper reports (linear regression) | **0.669** | 0.566 | 0.542 |
| | **ProteinTalks (paper, Table S6A)** | 0.638 | 0.573 | 0.775 |

The drug-mean predictor is one line: for a held-out condition, predict that
drug's efficacy rate across the training cell lines.

**Provenance of the ProteinTalks rows.** The Setting 2 and Setting 3 values are
from the *published* Supplementary Tables S5A and S6A, so they are the final
peer-reviewed numbers. The Setting 1 headline (0.960 / 0.854 / 0.910) is from
the **preprint** Results text: the *Nature* body text is paywalled and only its
abstract, figure captions and availability statements are readable without a
subscription, so we could not confirm that this figure is unchanged in the
published version. Treat the Setting 1 row as provisional.

Two things follow. First, that one-line predictor **beats every baseline the
paper benchmarks** in the leave-one-cell-line-out setting — KNN 0.805,
GeneCompass 0.751, Geneformer 0.684, DeepSynergy 0.806, UCE 0.482 — and lands
within 0.034 AUROC of ProteinTalks itself. Second, in the leave-one-drug-out
setting the drug-mean predictor is 0.5 by construction, and ProteinTalks does
show real signal over the cell-line-mean control (0.638 vs 0.508).

## 6. Paired per-cell-line test

Table S5C reports per-cell-line results for every benchmarked model, including
ProteinTalks (`ppODE`), over repeated runs. That permits the paired comparison
the paper does not make. `scripts/paired_cellline_test.py`, over the 15 cell
lines present in both tables, Wilcoxon signed-rank:

| Model | AUROC | drug-mean control | difference | p | wins |
|---|---|---|---|---|---|
| ppODE (Non-SWA) | 0.962 | 0.946 | **+0.016** | **0.015** | 12/15 |
| ppODE (SWA) | 0.951 | 0.946 | +0.004 | 0.454 | 10/15 |
| Linear regression | 0.941 | 0.946 | −0.005 | 0.639 | 7/15 |
| Random forest | 0.917 | 0.946 | −0.029 | 0.004 | 2/15 |
| KNN | 0.822 | 0.946 | −0.125 | 0.0001 | 1/15 |
| DeepSynergy | 0.806 | 0.946 | −0.140 | 6.1e-05 | 0/15 |
| GeneCompass | 0.748 | 0.946 | −0.198 | 6.1e-05 | 0/15 |
| Geneformer | 0.694 | 0.946 | −0.252 | 6.1e-05 | 0/15 |
| UCE | 0.489 | 0.946 | −0.457 | 6.1e-05 | 0/15 |

Reading this fairly, in both directions:

- **ProteinTalks (non-SWA) is the only model that is significantly better than a
  predictor with no inputs**, by 0.016 AUROC (p = 0.015). That is a real result,
  and it is more than most published perturbation models achieve against this
  kind of control.
- Two models are not significantly *worse* than the control: ProteinTalks (both
  variants) and **linear regression** (0.941, p = 0.64). Everything else the
  paper benchmarks is significantly worse, including all three transcriptomic
  foundation models and DeepSynergy.
- The margin over the control is 0.016, not the 0.15–0.28 implied by comparing
  against the reported baselines. Those baselines are weak: three transcriptomic
  foundation models applied zero-shot to a proteomic task, plus a KNN.
- The SWA variant, which is the one the paper's Table S5A headline corresponds
  to most closely, is **not** significantly better than the control (p = 0.45).

**Limitation of this test.** The pairing assumes the paper's held-out set for a
given cell line contains the same conditions as that cell line's rows in Table
S1 `C_Efficacy`. That is the natural reading of leave-one-cell-line-out on a
18 × 62 design, but the per-cell-line test composition is not published, so it
could not be verified. Table S5C covers 17 cell lines and `C_Efficacy` covers
18; 15 are common and only those are used. If the paper's folds exclude
conditions ours include, the control's per-cell-line AUROC would shift. The
direction and size of every gap in the table would have to change a great deal
to alter the conclusion, but the caveat is real.

## 7. The paper's own supplement already reports the negative result for new drugs

Supplementary Table S6B, verbatim significance calls for leave-one-drug-out:

| Comparison | Metric | p | Significant |
|---|---|---|---|
| ProteinTalks vs random forest | AUPRC | 0.626 | n.s. |
| ProteinTalks vs random forest | AUROC | 0.798 | n.s. |
| ProteinTalks vs random forest | Accuracy | 0.015 | * |
| ProteinTalks vs bootstrap | AUPRC | 0.918 | n.s. |
| ProteinTalks vs bootstrap | AUROC | 0.328 | n.s. |

On unseen chemistry — the setting that matters for prospective drug discovery —
ProteinTalks does not outperform a random forest on ranking metrics, and linear
regression has a higher mean AUROC (0.669 vs 0.638). The paper reports this
honestly in its supplement. It deserves to be read alongside the 0.960 headline.

## 8. What was NOT reproduced

- **The headline 0.960 AUROC was not re-derived from data.** That requires the
  gated 325 MB PTDS matrix. We compare against the paper's reported values.
- **No retraining of the published model.** Without the protein index the
  released weights cannot be mapped onto any proteome we can obtain.
- **No biological claim was checked.** The SHAP-nominated resistance proteins
  (AKR1C3, CMPK1), the four validated synergistic combinations, and the patient
  survival analysis all rest on wet-lab work and on data we do not have.
- **The trajectory-prediction quality (Loss1) was not evaluated on real data.**
  Table S5C reports a `Pearson` column for the proteome task; for GeneCompass it
  sits at 0.006–0.025, but we did not establish what the corresponding
  ProteinTalks values mean without the underlying matrix.

## 9. Reproducing these checks

```bash
git clone https://github.com/guomics-lab/PTV-1        # code + weights
# Nature Supplementary ZIP (126 MB, no login):
#   https://media.springernature.com/original/springer-static/esm/
#   art%3A10.1038%2Fs41586-026-11001-9/MediaObjects/41586_2026_11001_MOESM3_ESM.zip

python scripts/verify_checkpoint.py   --checkpoint PTV-1/ProteinTalks/best_checkpoint.pth
python scripts/equivalence_check.py   --repo PTV-1 --checkpoint PTV-1/ProteinTalks/best_checkpoint.pth
python scripts/real_label_controls.py --table-s1 TableS1_drug_cell_0702.xlsx
python scripts/paired_cellline_test.py --table-s1 TableS1_drug_cell_0702.xlsx \
                                       --table-s5 TableS5_cross_cell_final_260823_ESM.xlsx
```
