# trnaplat — the three demos

The platform layer on top of `fftplsr` (FFT-PLSR, *Nat Commun* 2025) and `pylrs`
(the PylRS-libY audit and the MjTyrRS/McTyrRS datasets). Three demos, in the
order they were specified:

| | Demo | Tool | Needs wet lab? |
|---|---|---|---|
| 1 | AzK on a validated PylRS pair, multi-site full-length yield | `demo1_multisite.py` | yes — this module designs the plate and analyses it |
| 1c | libY run once as the out-of-domain negative control | `demo1_libY_control.py` | no |
| 2 | Fixed substrate, compare enzyme frameworks | `demo2_scaffold.py` | **no** — the report is the deliverable |
| 3 | A group of ncAAs, evidence tiers for a new one | `demo3_evidence.py` | no |

Shared: `ncaa_chemotype.py` holds the hand-curated descriptors for 13 TyrRS
ncAAs, libY's 8 training substrates and AzK, in one feature space.

```bash
python3 -m trnaplat.demo1_multisite            # design + power analysis
python3 -m trnaplat.demo1_multisite --ingest yields.csv --baseline-arm unmodified
python3 -m trnaplat.demo1_libY_control --repo /path/to/PylRS-libY
python3 -m trnaplat.demo2_scaffold             # writes DEMO2_REPORT.md
python3 -m trnaplat.demo3_evidence
python3 -m trnaplat.ncaa_chemotype
python3 -m pytest tests/test_trnaplat.py -q    # 34 tests
```

---

## Demo 1 — the transplant is free, and the curve is the readout

✅ **The activity transplant needs no alignment.** The FPFORCOM parent "IFRS" is
wild-type *M. mazei* PylRS with exactly **two** substitutions, **N346S** and
**C348I**, both in the catalytic domain. All **7/7** COM1 activity substitutions
(`D2N/V31I/T56P/R61K/H62Y/T122S/S193R`) sit at positions where wild-type MmPylRS
carries precisely the residue the paper mutated away from — no remapping, no
divergent sites. Contrast `pylrs/mc_scaffold.py`, where the MjTyrRS → McTyrRS
transfer has 7 divergent positions touching 75 of 424 substitutions.

✅ **6 of the 7 lie outside the catalytic domain** (N-terminal tRNA-binding
domain and linker). That is the structural reason this axis is ncAA-agnostic,
and the nearest activity position is **153 residues** from the nearest pocket
position along the chain.

⚠️ "Defined" is not "effective". The evidence that these transfer is the paper's
own 7 synthetases × 6 ncAA types, not this arithmetic.

### Why multi-site, not single-site

The gain compounds as `(p₁/p₀)ⁿ`. A per-site efficiency improvement of
0.55 → 0.72 is **1.31×** at one site — inside assay noise — and **2.94×** at
four. That is the entire argument for making the readout a curve.

### What the curve's shape decides

Fit `log Y(n)/Y(0) = a·n + b·n²`. `b ≈ 0` means sites are independent and the
synthetase is what to engineer; `b < 0` means a shared resource (charged tRNA
pool, EF-Tu, ribosome stalling) is running out and **more synthetase activity
cannot rescue high n** — the tRNA is.

✅ Power, simulated at 15% assay CV over n = 1..4:

| replicates | b = 0 (false positive) | b = −0.02 | b = −0.04 | b = −0.08 |
|---|---|---|---|---|
| 2 | 0.073 | 0.245 | 0.509 | 0.888 |
| 3 | 0.071 | 0.295 | 0.614 | 0.960 |
| 4 | 0.056 | 0.321 | 0.700 | 0.982 |
| 6 | 0.045 | 0.402 | 0.812 | 0.998 |

⚠️ The `F > 4 / b < −0.02` rule is a screening heuristic, not a calibrated test:
even at 6 replicates it mislabels ~1 independent-site curve in 20 as
tRNA-limited, and `b = −0.02` is not reachable at all with this design. Confirm
any positive call by titrating tRNA copy number — the orthogonal experiment, and
the actionable one.

⚠️ Every yield number above is assumed. AzK's charging by the unmodified pair
comes from the client's validated construct, not from anything established here;
fix the exact AzK regioisomer and the tRNA copy number before ordering.

---

## Demo 1c — libY as the negative control

Three separate findings, two of them blockers rather than scores:

❌ **libY structurally cannot score AzK.** Only **18 of 280** features (6.4%) are
substrate-aware — Rosetta `cartesian_ddg`, licence-gated. The other 262 are
identical across all 8 substrate rows of a given variant. And the 18 that do
know the substrate transfer at **AUC 0.487** across substrates — chance — while
the variant-only block reaches 0.704. The one block that knows which substrate
it was asked about is both unobtainable for a new one and worthless for
transfer.

✅ **The categorical test detects the out-of-domain case decisively.** libY's 8
training substrates occupy one chemotype cell (`parent=Tyr`, `attach=para_O`;
confirmed against `3.training/UAA-info.csv` SMILES). All **3/3** of AzK's
categorical levels appear in **zero** libY training rows. A one-hot level no
training row carries has no fitted weight — the model is blind, not
extrapolating.

❌ **Euclidean chemotype distance fails to detect it.** AzK sits **3.70** from
its nearest libY substrate while libY's own substrates span up to **4.32**
internally: **3 of 28** within-training pairs are further apart than AzK is from
the training set. By distance alone a different *parent amino acid* does not look
out of domain.

---

## Demo 2 — scaffold choice, from measured data only

See `DEMO2_REPORT.md`, regenerated from `pylrs/data/` on every run.

✅ Best McTyrRS variant **F+/F⁻ = 3.69** vs best MjTyrRS **2.79**, same paper,
same readout. ✅ Two species converged independently: the best McTyrRS variant is
`Y32G/D158T` in Mj numbering, and `AzPheRS-6`, selected separately for pAzF on
MjTyrRS, carries both. ✅ The substitution that rescued McTyrRS, **Y112F**, maps
to Mj position 108 where the wild type is *already F* — so McTyrRS's advantage
lives in the rest of the scaffold, not at that site.

❌ **The screening model cannot see the mutation that mattered.** `Mc Mut6`
(F+/F⁻ = 1.11) and `Mc Mut6+RM` (3.69) differ only by Y112F, which lies outside
the five scored positions. The screen gives both the **same rank (65/2001)**
despite a **3.3×** measured difference. The hotspot set defines what the model
can see; a substitution outside it is absent, not down-weighted.

✅ Calibration, leave-one-clone-out against 2000 random library members: on
MjTyrRS, 4 of 7 known clones land in the top 10%, median percentile 7.55%
(**7× enrichment**), worst clone at the 96th percentile — a clone that works,
scored as a confident reject.

⚠️ The McTyrRS row of that table is **not** a scaffold comparison: 7 of its 9
pAzF clones are coordinate rewrites of the MjTyrRS clones.

---

## Demo 3 — evidence tiers, and two ideas that did not survive

libY's evaluation architecture, applied to 13 substrates: per-substrate labels,
leave-one-ncAA-out, accuracy per feature block.

✅ Per-block leave-one-ncAA-out (mean over 13 held-out substrates):

| block | mean AUC | range |
|---|---|---|
| onehot | 0.903 | 0.541–1.000 |
| onehot + pssm | 0.902 | 0.534–1.000 |
| pssm | 0.901 | 0.572–1.000 |
| **onehot × chemotype** | **0.664** | 0.092–1.000 |
| chemotype | 0.500 | degenerate in all 13 folds |

✅ The `chemotype` block is constant inside every held-out fold and scores
exactly 0.5 — ncAA descriptors can only enter through an **interaction** with the
variant. ✅ That interaction block reaches **0.664**, within noise of libY's own
honest grouped-CV figure (~0.64) on a different scaffold and different
substrates. Two independent datasets landing in the same place is the most
reliable number in this report.

⚠️ The gap between 0.903 and 0.664 is the decoy artifact: `onehot` wins by
recognising "looks like a selected clone", which is a fact about selection, not
about recognition.

### ❌ Chemotype distance does not predict transfer

Slope **+0.063 AUC per chemotype unit, 95% CI ±0.095, r = +0.36, n = 13** — the
interval spans zero. The tier table refutes itself directly: `pCMF` is the most
distant substrate (4.19, tier D) and scores **0.926**, while `OAY` sits at 1.41
(tier B) and scores **0.541**. Ship the tiers as a **disclosure of how far the
model is reaching**, never as a discount factor on its score.

### ❌ Grouping by product need does not help training

| grouping | within-group AUC | trained on all other substrates |
|---|---|---|
| product need | 0.760 (n=6) | 0.900 |
| chemical class | 0.805 (n=10) | 0.886 |

Focused training wins in only **4 of 16** held-out cases. Both grouping rules
lose, and the product groups lose by more — they are smaller (2 members, 2–8
training clones). On this dataset **more data beats more relevant data**,
decisively.

✅ What the rule is still right for: **defining the product.** A carrier needing
two mutually orthogonal handles needs pAzF + pPRF whatever a model prefers to
train on. Group by product to choose what to build; train on every substrate you
have.

### The categorical gate, carried over from Demo 1c

⚠️ Not validated here either: with `subst_class` excluded as ID-like (13 levels
over 13 molecules — an identifier, not a feature), only **2-NPA** is out of
domain, and it scores 0.947, above the in-domain mean. n = 1 on one side could
not have shown much. The gate remains the right *principle* — it fires 3/3 on
AzK vs libY, where the model genuinely has no parameters — but on these 13
substrates the platform should report *which* levels are unseen as a factual
disclosure and make no accuracy claim from it.

---

## The standing caveat

⚠️ Every AUC in Demo 3 and every calibration in Demo 2 is measured against
**decoys** — unselected library members, not measured rejections. They answer
"does this look like something selection would keep", not "does this recognise
the ncAA". `pylrs/negative_panel.py` designs the score-stratified plate that
would fix this, and `pylrs/mc_scaffold.py` surfaces the only measured non-hits
anywhere in the dataset (sheet 3's Mut4/5/7, for which the paper gives no
numbers). Run that plate before any of these numbers goes to a client.
