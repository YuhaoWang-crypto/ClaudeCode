# Trying to improve the model, and what the attempt actually found

> **Correction, added after a further control.** An earlier version of this
> document concluded that "nothing learns the direction of the response" and that
> the five changes largely failed. Both statements were **true only at the
> published loss weighting, lambda = 0.8**, and I did not initially test whether
> that weighting was the cause. It is. Under trajectory-only training the
> architecture does learn the response, and the five changes give a large
> improvement rather than none. The corrected result is in
> "The loss weighting was doing the damage" below. The original lambda = 0.8
> table is kept unchanged, because it is what the published configuration does.

Short version: at the published objective the proteome task is barely trained, so
the ladder looks flat and the model looks broken. Remove that confound and the
changes work, though a linear baseline is still well ahead.

All numbers below are on the mechanistic simulator, not on the real corpus. See
the limits section at the end before quoting any of them.

---

## What was tried, and why

Five changes, each motivated by something this reproduction measured or
something the published Supplementary Information says the model should do but
the released code does not. Implemented in `proteintalks/improved.py`, each
behind its own flag.

| # | Change | Motivation |
|---|---|---|
| A | **Residual decoding** `P(t) = P(0) + g(z(t))` | The released model's trajectory MSE was far worse than asserting the proteome does not move. Zero-initialised so an untrained model *is* the no-change predictor. |
| B | **Real elapsed time** `[0,6,24,48]/48` | The released code integrates over integer ticks 0,1,2,3. The SI's own Taylor justification is over the real hour grid. |
| C | **`f(z,t,D)` conditioning** | The SI writes the field as non-autonomous and perturbation-conditioned. The released field takes neither. |
| D | **Low-rank protein coupling** | Every convolution in the released dynamics uses `kernel_size=1`, so module 1 has no protein-protein coupling at all. Couples through `r` learned factors rather than widening a kernel over an arbitrary protein ordering. |
| E | **Low-rank phenotype head** | `drugsens_conv1` is 90.5% of the model's parameters. |

The ladder starts from the released architecture and adds one change at a time.
That baseline is verified, not assumed: `scripts/verify_ladder_baseline.py`
loads the released weights into the improved model with every flag off and
reproduces the reference outputs to **exactly zero** difference on both heads.

## Results

1,116 conditions, 200 proteins, 18 cell lines, 62 drugs. Hidden width 32,
60 epochs, AdamW 5e-4, fixed λ = 0.8 (the configuration the SI says the reported
results used). Two seeds, mean reported.

`skill` is `1 − MSE/MSE_nochange` on the 6/24/48 h proteome: positive beats "the
proteome does not move". `delta_r` is the correlation between predicted and true
*change* from baseline.

| configuration | skill | delta_r | AUROC | params |
|---|---|---|---|---|
| released architecture | −2.693 | +0.006 | 0.918 | 92,226 |
| + residual decoding | **−0.071** | −0.088 | 0.919 | 92,226 |
| + real elapsed time | −0.259 | −0.026 | 0.929 | 92,226 |
| + `f(z,t,D)` conditioning | −3.410 | −0.034 | 0.919 | 92,754 |
| + protein coupling | −0.381 | **+0.040** | 0.910 | 97,011 |
| + low-rank head (all five) | −3.407 | −0.024 | **0.937** | 92,403 |
| *drug-mean control (no proteome)* | — | — | *0.908* | *0* |

What this says, read honestly:

**A worked, on the metric it targeted.** Residual decoding cuts excess trajectory
error by roughly 38× (−2.693 → −0.071) at zero parameter cost. But it never gets
*above* the no-change baseline, so the model still does not beat "nothing
changed".

**Implementing what the SI describes made things worse.** Change C, the
non-autonomous perturbation-conditioned field the SI's own equations specify,
scored the worst trajectory skill in the ladder (−3.410) and no AUROC gain. On
this benchmark at this budget, the SI's formulation is not an improvement.

**The classification gains are inside the noise floor of the question.** The
whole ladder spans AUROC 0.910 to 0.937 against a no-proteomics control at 0.908.
The best variant is +0.029 over a predictor that uses no protein data at all, on
a 112-condition test set, over two seeds.

**At this weighting, nothing learns the direction of the response.** `delta_r`
never leaves [−0.088, +0.040] anywhere in the ladder, released or improved. The
next two sections establish that this is a property of the objective, not of the
architecture.

## The control that changed the conclusion

`delta_r ≈ 0` everywhere admits two readings: the models fail to learn a response
that is there, or the benchmark has no learnable response and the column measures
noise. Reporting the first when the truth is the second is how benchmark work
goes wrong, so `scripts/delta_learnability_control.py` asks whether *anything*
can predict the delta, using methods with no dynamical structure at all.

| predictor | 6 h | 24 h | 48 h |
|---|---|---|---|
| **ridge regression**, delta_r | **+0.888** | **+0.949** | **+0.948** |
| **ridge regression**, skill | **+0.829** | **+0.965** | **+0.971** |
| cell-line-mean delta, delta_r | +0.174 | +0.131 | +0.124 |
| no-change, by definition | 0.000 | 0.000 | 0.000 |

A single linear map from baseline proteome, perturbation mask and drug
fingerprint to the delta recovers **97% of the response variance** at 48 h. The
signal is not subtle: the response RMS is 1.5× the injected measurement noise.

So the benchmark is fine, and the reading is the first one: at lambda = 0.8 the
neural-ODE architecture, released and improved alike, fails to learn a
perturbation response that ridge regression learns almost perfectly on the same
inputs and the same split.

That looked like a far more consequential finding than the +0.029 AUROC the
optimisation was chasing, and like the same pattern the 2025 benchmarking
literature reports
for transcriptomic perturbation models: Ahlmann-Eltze, Huber & Anders
([Nat Methods 22:1657](https://doi.org/10.1038/s41592-025-02772-6)) found no deep
model consistently beating a linear model or the mean; Wong, Hill & Moccia
([Bioinformatics 41:btaf317](https://doi.org/10.1093/bioinformatics/btaf317))
found ablating scGPT's pretrained weights changed nothing.

**But the comparison as stated was not like-for-like, and that matters — see the
next section.** Ridge optimises the trajectory alone; the ODE variants were
giving it 0.2 weight.

**State this caveat too whenever the ridge number is quoted.** Ridge is *not*
parameter-matched. It fits one coefficient per (input feature, output protein)
per timepoint — about 267,000 coefficients per timepoint at 200 proteins, against
the neural ODE's 92,226 shared across all three. It is a capacity-rich upper
bound on available signal, not a fair architectural rival. The point is not that
ridge is the better model. The point is that the signal is plainly there and the
ODE does not find it.

## The loss weighting was doing the damage

The objection to the section above is immediate: ridge optimises the trajectory
alone, while every ODE variant was trained at lambda = 0.8, which puts only
**0.2** weight on the proteome loss. So the comparison was not like-for-like, and
the flat ladder may say more about the objective than the architecture.

`scripts/lambda_starvation_control.py` sweeps lambda, holding everything else
fixed. Two seeds, same budget, same split.

| lambda | configuration | skill | delta_r |
|---|---|---|---|
| 0.8 *(published)* | released | −2.693 | +0.006 |
| 0.8 *(published)* | all five changes | −3.407 | −0.024 |
| 0.5 | released | −0.176 | +0.186 |
| 0.5 | all five changes | −0.234 | +0.065 |
| **0.0** *(trajectory only)* | released | **+0.036** | **+0.361** |
| **0.0** *(trajectory only)* | all five changes | **+0.424** | **+0.641** |
| — | *ridge, trajectory only, no ODE* | *+0.971* | *+0.948* |

Three things follow, and the first two correct this document's earlier claims.

**The architecture is not incapable — the published objective barely trains the
proteome task.** At lambda = 0 the released architecture reaches delta_r +0.361
and, for the first time anywhere in this work, a *positive* skill score. The
"learns nothing" reading was an artefact of comparing a 0.2-weighted task against
a baseline that optimised it exclusively.

**The five changes do work, once they are not masked.** At lambda = 0 they take
skill from +0.036 to **+0.424** and delta_r from +0.361 to **+0.641**, roughly
doubling the correlation and turning a marginal fit into a substantial one. At
lambda = 0.8 they appeared to do nothing. The improvements were real all along;
the default weighting hid them.

**Ridge still wins, so a real architectural gap remains.** +0.971 / +0.948
against +0.424 / +0.641. The honest statement is "the neural ODE is well behind a
linear map on this task", not "the neural ODE learns nothing".

This also says something about the published model itself. The SI confirms the
reported results used fixed lambda = 0.8. On this simulator that setting leaves
the proteome-dynamics half of a model named for its proteome dynamics close to
untrained, while the drug-efficacy head — which carries 98.3% of the parameters —
gets 0.8 of the gradient. Whether the same holds on the real corpus cannot be
checked without the gated matrix.

## Limits

- **Simulator only.** None of this touches the real corpus, so it says nothing
  about whether ProteinTalks predicts real proteome dynamics. The real matrix is
  approval-gated and the checkpoint's protein index is unpublished.
- **Small budget.** 60 epochs at hidden width 32 against the released config's
  1,000 epochs at width 64. The released architecture may do better with its own
  budget; this ladder shows relative behaviour under a shared one.
- **Two seeds, 112 test conditions.** Differences under ~0.03 AUROC are not
  resolved.
- **The simulator is generous to the ODE**: its ground truth *is* an ODE on a
  protein network, which is exactly the model class the architecture assumes. The
  failure to learn it is therefore harder to explain away, not easier.
- **The ladder is cumulative**, so a change is evaluated on top of the ones above
  it. C scoring worst does not prove C is harmful in isolation.

## Reproducing

```bash
python scripts/verify_ladder_baseline.py        # rung 0 == released model
python scripts/optimize_model.py --n-proteins 200 --hidden 32 --epochs 60 --seeds 0 1
python scripts/delta_learnability_control.py    # is the response learnable at all
python scripts/lambda_starvation_control.py    # is the objective starving Loss1
```
