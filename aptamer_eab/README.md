# aptamer_eab — ssDNA library design for NPY / PP (PPY) EAB sensors

In-silico prioritisation of ssDNA candidates for electrochemical aptamer-based
(EAB / E-AB) sensors against **neuropeptide Y** and **pancreatic polypeptide**.
Full write-up, including what is already published and what to order:
[`../APTAMER_EAB_REPORT.md`](../APTAMER_EAB_REPORT.md).

**Nothing here has been experimentally tested.** A library is a starting pool
for SELEX / screening, not a binder.

## Run

```bash
pip install numpy pandas viennarna matplotlib
python3 -m aptamer_eab.run_all
```

Runs in ~2 min on CPU. Writes `output/*.csv` and `../figures/aptamer_eab_overview.png`.

## Modules

| file | what it holds |
|---|---|
| `targets.py` | mature NPY / PP / PYY sequences (UniProt), family identity, sliding-window epitope discrimination scan, PP-fold context |
| `known.py` | literature aptamers with citations — ssDNA **4.31** (the one with a published EAB sensor), RNA **DP3**, primer arms, degenerate-oligo specs |
| `eab_filters.py` | EAB switchability metrics + score, the **anchor check**, G4Hunter / `g4_quality`, EAB construct spec |
| `build_library.py` | sub-libraries A–D, S (short), G (G-quadruplex), controls, diversity-aware panel picking |
| `run_all.py` | orchestrates everything, prints the audit trail, writes CSVs and the figure |

## Outputs

| file | contents |
|---|---|
| `output/ORDER_PANEL.csv` | **the deliverable** — 58 curated candidates (22 NPY, 22 PP, 14 target-agnostic), pairwise core edit distance ≥ 6 |
| `output/order_specs.csv` | 4 orderable degenerate oligos (naive N40, NPY maturation, PP retargeting, primers) |
| `output/library_{A,B,C,D,S,G}_full.csv` | every scored candidate per sub-library |
| `output/epitope_scan_{NPY,PP}.csv` | per-window discrimination scores |
| `output/controls.csv` | positive reference + the two negative controls |

## Two things the code deliberately does NOT do

**It does not trust its own score.** `eab_filters.anchor_check()` tests the
scoring function against the only two published E-AB outcomes for this target
(4.31 full 80-nt works at 230 % signal change; its bare 40-nt core gives no
signal). `run_all` aborts rather than emit an order panel if the check fails.
An earlier version of the score ranked the working aptamer *below* its own
scramble control — that is why the check exists. Two anchor points, set after
seeing the sequences, is an anchor, not a validation.

**It does not use the arm-pairing metrics as discriminators.** `arm_arm_bp` and
`arm_core_bp` are reported, never scored: the two fixed arms are identical in
every framed candidate, so 4.31 (9 / 5) and its scramble (9 / 4) look nearly the
same. They exist as a *design constraint* — an orthogonal contact-probability
run puts ~1/3 of the predicted interface on the 5' arm, with strongest contacts
at T14/C15/A16, and this fold independently pairs arm positions 14,15,18,19,20
into the core. Library C dopes those five core partners at 5 % instead of 30 %
so retargeting mutagenesis does not wipe the duplex out (panel retention
9/12 → 11/12).

**It does not score G-quadruplexes with ViennaRNA.** The DNA parameter set has
no G4 term, so a real G4 is scored as coil. The G track is ranked by
`g4_quality` instead, anchored on the thrombin-binding aptamer (G4Hunter 1.13 —
*below* the usual 1.2 threshold, which is why maximising G4Hunter is the wrong
objective and collapses the pool to G-homopolymers).
