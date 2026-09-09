# What the paper specifies, and what it leaves open

The ProteinTalks methods section is unusually explicit: Equations 1–18 pin down
the layer types, dimensions, activations, dropout rate, ODE solver, loss
functions and the multi-task weighting rule. Most re-implementations of a
biology foundation model have to guess far more than this.

Below is an audit of what could be implemented directly from the text and what
required a decision. Everything in the "resolved by interpretation" column is
marked with an `# INTERPRETATION` comment at the corresponding line of
`proteintalks/model.py`.

## Fully specified by the paper

| Component | Specification | Source |
|---|---|---|
| Input | 5,585 proteins `P0` concatenated with perturbation channel `D` | Eq. 1 |
| `L1` | Linear, 2 → 32 | Eq. 1 |
| `C1` | Convolutional, 32 → 128 | text after Eq. 1 |
| `L2`, `L3` | Two linear layers, Softplus activation, dropout 0.1 | Eq. 2–3 |
| ODE solver | `rk4`, fixed-step fourth-order Runge–Kutta | Eq. 4 |
| `C2` | Convolutional, 128 → 32 | text after Eq. 4 |
| `L4` | Linear, 32 → proteomics space | Eq. 5 |
| Loss1 | MSE over the 6/24/48 h proteomes | Eq. 6 |
| `C3` | Convolutional over `[P0, P̃6, P̃24, P̃48]` → 128 | Eq. 7 |
| Drug features | 881 fingerprint bits + 55 physicochemical = 935, duplicated for single-drug conditions | Eq. 8 |
| `C4` | Convolutional, (935 × 2) → 128 | Eq. 8 |
| `L5` | Linear + ReLU → 32 | Eq. 9 |
| `L6` | Linear + sigmoid → 1 | Eq. 10 |
| Loss2 | Binary cross-entropy | Eq. 11 |
| λ | 0.8 | Eq. 12 |
| Weight adaptation | cosine similarity of task gradients, step 0.01 | Eq. 13–18 |
| Preprocessing | impute at 0.8 × minimum detected intensity; per-sample min–max | Methods |
| Splits | 0.7/0.2/0.1 → 1070/305/154 conditions | Methods, Setting 1 |

## Resolved by interpretation

**1. Pooling before `L5`.** Eq. 7–10 take `C3` (a convolution over 5,585
proteins) straight into `L5` and then to a scalar. No pooling step is named, but
one is required for the dimensions to close. We use a strided convolution stack
followed by global average pooling over the protein axis. This also keeps the
parameter count independent of the number of proteins, which is consistent with
the paper's emphasis on a small model.

**2. Time rescaling inside the ODE solve.** The paper integrates over
t = 0, 6, 24, 48 in hours. A fixed-step RK4 solve of a learned vector field over
48 unscaled time units diverges in float32. We rescale t by 1/48 before the
solve. This is a reparameterisation absorbed into the learned field, so it does
not change the model class.

**3. Convolution kernel sizes and strides.** The paper names the input and
output channel counts of `C1`–`C4` but no kernel widths. We use kernel 3,
stride 1 for `C1`/`C2` (which must preserve the protein axis for `L4` to decode
per-protein abundances) and kernel 7, stride 4 for `C3`/`C4` (which must reduce
it).

**4. Which parameters the gradient cosine is measured on.** Eq. 13 needs a
parameter set shared by both losses. Loss2 reaches the dynamics module only
through the predicted proteome, so the shared trunk is the whole of module 1.
That is what we use.

**5. The sign convention in Eq. 14–18.** This is the one place where the
equations and the prose disagree, and it matters.

Eq. 15 sets `adjustment_factor = 0.01 × clipped_similarity`. The adjustment is
applied only when the gradients *conflict*, which means the cosine similarity is
negative, which makes `adjustment_factor` negative. Substituted into Eq. 17–18,
a negative factor *increases* `w_Loss1` and *decreases* `w_Loss2` — the exact
opposite of the stated intent, which is: "If the gradient directions conflict,
we prioritize the drug efficacy prediction task by increasing its gradient
weight while decreasing the weights of other tasks."

We follow the prose and use the absolute value. `GradientConflictWeighter(
literal_sign=True)` reproduces the equations as printed, for anyone who wants to
check both. Note also that Eq. 14's clip at 1.0 is a no-op, since a cosine
similarity is already bounded above by 1; it is implemented anyway for fidelity.

## Not specified anywhere in the paper

These have no stated value and are not recoverable from the text. They are
exposed as arguments with the defaults listed:

- optimiser and learning rate (we use Adam, 1e-3)
- batch size (32)
- number of training epochs and the early-stopping rule (300, patience 40 on the
  combined validation objective)
- weight initialisation
- whether the perturbation channel `D` encodes drug target identity, dose, or
  both. The results text says the second module uses "61 targets of 63 drugs",
  which implies target identity is available to the model; Eq. 1's `D` has the
  same dimension as `P0`, so we encode it as a per-protein target-occupancy
  vector.
- how the three timepoints are weighted inside Loss1 (we weight them equally)

## Parameter count

The implementation has **168,226 parameters** at 400 proteins, and the same
count at 5,585 proteins, because every layer is either per-protein
weight-shared or followed by global pooling. For scale, scGPT has ~53 M
parameters and Geneformer ~10–110 M. The paper's claim of "significantly fewer
parameters" is architecturally sound and reproduces directly.
