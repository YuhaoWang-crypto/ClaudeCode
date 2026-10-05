# Applying FFT-PLSR to a new enzyme

A checklist for running this method on a target that is not PylRS. Every step
names the failure it prevents, because most of them are failures that produce a
confident-looking ranking rather than an error.

## 0. Check the method fits the problem

Answer these before writing code. A "no" anywhere means stop.

| Question | Why it matters |
|---|---|
| Do you have **measured** activities for single mutants? | This is a regression model. With no measurements there is nothing to fit — use a PLM/structure tool instead (`enzyme-mutation-ranking`). |
| Are they all on **one parent background**, in **one assay**? | Fitness is a fold-change relative to one parent. Mixing backgrounds silently compares incomparable numbers. |
| Is the **parent itself measured** as the reference? | The parent anchors the scale and is the row every baseline needs. |
| Are you combining mutations at **distinct positions**? | The model recombines; it cannot invent substitutions it has never seen. |
| Do you have ≥ ~10 singles, and a budget to assay ~8–20 combinations? | Below that the model cannot be validated, and the ranking cannot be acted on. |

⚠️ If the mutations are **not** at distinct positions, or you want to explore
positions with no data, this is the wrong tool.

## 1. Encode the data

```python
measured = {
    "WT":    1.00,   # the parent, label must contain no substitution
    "D2N":   3.62,
    "H62Y":  1.97,
    "D2N/H62Y": 6.31,   # include multi-mutants if you have them
}
parent = "MDKKPLNTLISATGLW..."   # the background all labels are relative to
```

Label format is `<wt><1-based pos><mut>`, `/`-joined. `fftplsr.variants` validates
every stated wild-type residue against `parent` and **raises** on a mismatch — the
reference implementation printed a warning and carried on, producing a variant
whose label and sequence disagree. If you hit `VariantError`, your numbering is
off (a signal peptide trimmed, or 0- vs 1-based) — fix it, do not suppress it.

Replicates: average them *before* building `measured`, and keep the spread. A
model fitted to point estimates cannot know that two variants differ by less than
assay noise, and a top-8 list separated by less than that spread is arbitrary —
widen the diversity filter rather than trusting the order.

## 2. Choose sites, and do not over-prune

```python
from fftplsr import design
sites = design.improved_sites(measured, threshold=1.05)
```

Prune **positions** (drop substitutions that do not beat the parent); do not prune
substitutions *within* a position by their solo effect.

✅ The measured reason: on the PylRS data, the best variant in the released panel
(`N7Y/H63L/K67N/V74W`, 2.75× parent) is built from singles ranking 3rd of 3, 1st
of 7, 4th of 5, and 2nd of 2 at their positions. `per_position=3` drops `K67N` and
makes the answer unreachable before any model is fitted. A single mutation's solo
effect is a poor predictor of its value in combination — which is the premise of
the whole exercise.

Space size is `prod(1 + n_subs_at_position)` over positions. Past ~2^20 variants,
cap with `max_order` rather than dropping substitutions: high-order combinations
are the least reliable predictions anyway.

## 3. Run the round and read it in the right order

```python
report = design.design_round(
    parent=parent, measured=measured, sites=sites,
    n_rounds=1,      # 3 for training sets of ~100+
    cv=None,         # LOOCV; switch to an int (10) once n > ~50
    pick=8, max_jaccard=0.6,
)
print(report.summary())
```

Read it like this, and stop early if a step fails:

1. **`report.baselines`.** If nothing beats `mean` (cvR2 ≤ 0), there is no
   learnable signal. Do not order from the ranking; get more or better data.
2. **FFT-PLSR's `cvR2` vs the best baseline.** If a baseline wins, use the
   baseline — it is deterministic, instant, and interpretable. `design_round`
   adds a note when this happens. The FFT encoding's one clear advantage is the
   **cold-start, singles-only round**, where each substitution is seen exactly
   once and one-hot models fall back to the intercept.
3. **`report.notes`.** Automatic flags.
4. **`report.picks`.** The order list. `epistasis` (predicted minus log-additive)
   marks where the model claims to add something over naive stacking; those rows
   are the most informative to assay regardless of rank, because they are what
   tests the model.

## 4. Do not over-read the predictions

- Treat predicted fitness as a **ranking**, not an expected fold-change, unless a
  held-out split has shown the calibration holds. ✅ On the PylRS round-2 split,
  held-out R² was 0.75–0.83 but the model systematically under-predicted the best
  variants (measured 11.1× vs predicted 8.4×).
- Cap `n_components` well below the training-set size. Selection pressure drives
  it upward because more components fit more noise; the paper's round-1 model used
  10 components from 12 points.
- `cvR2` after selecting the best of 566 descriptors is optimistic. Use
  `model.nested_cv_r2` for an honest figure when the decision depends on it.
- If the top candidates are within assay noise of each other, the ordering among
  them carries no information. Diversify and assay breadth instead.

## 5. Close the loop

Assay the picks, add them to `measured`, re-run. ✅ Observed in the paper's data:
the jump from 28 singles-only measurements to 120 singles+doubles is what made the
final winner reachable — doubles are where epistasis first becomes visible, so the
second round is where the model starts earning its keep over an additive baseline.

Keep every round's measurements, including the failures. Variants that came out
*worse* than predicted are the most informative training rows you will get, and
dropping them biases the next round upward.

## Common failure modes

| Symptom | Likely cause |
|---|---|
| `VariantError: states wild-type X but the parent has Y` | residue numbering mismatch — check for a trimmed tag/signal peptide |
| All baselines at cvR2 ≤ 0 | not enough data, or the assay noise exceeds the effect sizes |
| FFT-PLSR cvR2 ≫ held-out R² | descriptor-selection bias; run `nested_cv_r2` |
| Top-8 are near-identical | diversity filter too loose — lower `max_jaccard` |
| In-sample R² ≈ 0.99 | you are reading `train_r2`, not `cv_r2` |
| Scores change between machines | you set `drop_dc=False`; don't |
| The known-best variant is absent from the ranking | over-pruned with `per_position` (see step 2) |
