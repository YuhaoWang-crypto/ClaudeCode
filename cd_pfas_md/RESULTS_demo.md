# Demo run — proposed dye + CD modification, prefilter output

**Status:** heuristic prefilter, run on CPU. **These are NOT converged MD/FEP
numbers** — they are an order-of-magnitude triage from a transparent scoring
function (`src/prefilter.py`), whose every constant lives in
`config/system.yaml → prefilter:`. The point is to show the *shape* of the
workflow's output and to decide which designs deserve the real GPU FEP run.

## What was proposed and wired in

- **Reporter dye:** **TNS** (6-(p-toluidino)naphthalene-2-sulfonate) — a classic
  β-CD fluorescent probe: dark in water, bright in the cavity, so PFAS
  displacement gives a clean turn-OFF signal. Anionic sulfonate, ΔG anchored to
  the literature TNS·β-CD affinity (~4×10³ M⁻¹, ΔG ≈ −4.9 kcal/mol).
- **Featured modification:** **mono-6-(trimethylammonium)-β-CD** (a permanent
  cationic charge on the primary rim), plus the rest of the library and two
  **fluorophilic (fluorous-lined cavity)** designs as contrast.
- **Analytes:** PFOA, PFOS (modeled as −1 anions).

## Output (reproduce with `python -m cd_pfas_md.src.prefilter`)

```
host                                q  Fφ  ΔG_dye   disp[pfoa]   disp[pfos]    FOM  verdict
-------------------------------------------------------------------------------------------
fluorous_tagged                    +0   Y   -4.90        -0.60        -1.30  +0.95  FAVORABLE
fluorous_cationic                  +1   Y   -5.94        -0.60        -1.30  +0.95  FAVORABLE
beta_cyclodextrin (reference)      +0   ·   -4.90        +1.40        +0.70  -1.05  UNFAVORABLE
mono_6_amino                       +1   ·   -5.94        +1.40        +0.70  -1.05  UNFAVORABLE
mono_6_trimethylammonium           +1   ·   -5.94        +1.40        +0.70  -1.05  UNFAVORABLE
sulfobutylether                    -4   ·   -2.82        +1.40        +0.70  -1.05  UNFAVORABLE
methylated_dimeb                   +0   ·   -4.90        +1.40        +0.70  -1.05  UNFAVORABLE
hepta_6_amino                      +7   ·   -6.98        +1.40        +0.70  -1.05  UNFAVORABLE
```
`ΔG_displace = ΔG_bind(PFAS) − ΔG_bind(dye)`; **negative = PFAS out-competes the
dye = good sensor**. FOM = mean(−ΔG_displace) over the PFAS guests.

## The insight the run surfaces (this is the useful part)

1. **A symmetric rim charge does not buy PFAS selectivity.** TNS and PFAS are
   *both* −1 anions, so a cationic rim stabilizes them almost equally — the
   `disp` columns are **identical** for β-CD, mono-amino, mono-TMA, and
   hepta-amino. Charge tunes **overall affinity / robustness** (note ΔG_dye
   dropping from −4.9 → −6.98 across the +1 → +7 series), **not** the
   displacement that the assay actually reads out. A +7 rim even risks locking
   the dye in too tightly.

2. **Selectivity comes from what PFAS has that the dye doesn't — the fluorous
   tail.** Only the **fluorophilic (fluorous-lined) cavity** flips displacement
   from uphill (+1.4 kcal/mol for PFOA on native β-CD) to downhill (−0.6). That
   is the design lever to pursue.

3. **Combine both levers.** `fluorous_cationic` keeps the favorable displacement
   of the fluorous cavity *and* adds a charge for higher absolute affinity /
   lower detection limit — the natural lead to take into the real FEP screen.

4. **Probe choice matters too.** TNS binds fairly tightly (−4.9); against native
   β-CD it out-competes PFOA, so plain β-CD is a poor host for this pair. Either
   engineer PFAS-selective stabilization (above) or pick a weaker-binding probe.

## Honest caveats

- The fluorophilic bonus (−2.0 kcal/mol) and the electrostatic geometry
  (ε=40, r=4.5 Å, ≤2 contacting rim charges, 0.15 M Debye screening) are
  **assumptions**, not measurements. Different but reasonable choices shift the
  absolute numbers; the *ranking logic* (charge ≠ selectivity; fluorophilic =
  selectivity) is robust to them.
- Real accuracy needs `fep_ti_ddg.py` (relative ΔΔG) anchored by the APR
  calibration in `analyze_apr.py`. The prefilter's job is only to tell you the
  fluorophilic and fluorous-cationic hosts are worth those GPU cycles and the
  symmetric-cationic ones probably are not.

## Suggested next step

Run the rigorous FEP on the two flagged winners against PFOA/PFOS and TNS:

```bash
scripts/run_ddg_screen.sh pfoa   # then pfos, then dye  (needs a GPU)
```
and confirm (or correct) the prefilter's fluorophilic-selectivity prediction.
