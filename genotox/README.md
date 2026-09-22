# genotox — virtual in-vitro genotoxicity assays

**Step 1** — the **umu test** (SOS/umuDC-lacZ, *S. typhimurium*
TA1535/pSK1002).
**Step 2** — a mammalian **GADD45a-GFP reporter line** (p53/Mdm2 delayed
feedback).  It reuses step 1's upstream layer, readout registry and decision
layer unchanged; only the core is new, which was the point of the split.

```bash
pip install numpy scipy matplotlib
python3 -m genotox.run_umu        # step 1: tables + checks + figure
python3 -m genotox.run_p53        # step 2: pulses, cross-endpoint, checks
```

Writes `figures/genotox_umu.png` and `figures/genotox_p53.png`.

## The four seams

```
 chemistry          biology             measurement        decision
 ---------          -------             -----------        --------
 DamageSource  ->  SignalCore     ->    Readout      ->    call_result
 damage.py         core.py              readouts.py        doseresponse.py
        \             |                     |                   |
         `-> DamageFlux (per-channel vector, the only thing that crosses)
```

Each layer is replaceable on its own:

| To add | Change | Leave alone |
|---|---|---|
| a predictive chemistry front end | `DamageSource` subclass | core, readouts, decision |
| endpoint #2 (p53/GADD45a reporter) | new `SignalCore` — *done* | source, decision |
| endpoint #3 (comet, micronucleus) | new readout on a core exposing break density | source |
| a different protocol threshold | a `Protocol` instance | all biology |

Readouts hand the decision layer exactly two generic names, `signal` and
`biomass`; instrument-specific ones (`A410`, `fluorescence_per_cell`) ride
alongside for reporting.  Without that, the decision layer would have to know
that umu means absorbance and a reporter line means fluorescence, and the
layer split would be fictional.

## Why the interface is a vector

`DamageFlux` carries eight lesion channels, not one potency number. Each core
declares its own per-channel weights, so one upstream prediction produces
*different* answers at different endpoints — which is the actual situation:

- `SOSCore.sos_weights` puts **0.0** on the `aneugenic` channel, so a spindle
  poison is structurally incapable of being called SOS-positive. A
  micronucleus core would weight that same channel heavily. With a scalar
  interface, one of those two must be wrong.
- `lesions` is exposed as an observable specifically so a comet readout can
  be calibrated against the same quantity rather than re-deriving damage.

## Step 2: what is different about the mammalian core

Not a relabelled copy — four structural differences:

1. **Delayed negative feedback.** p53 induces Mdm2 through transcription plus
   nuclear import; the delay makes the damaged state *oscillate* (period
   ~5.3 h here) instead of settling. "The p53 level" is not a well-defined
   quantity, and the reporter sees a pulse train.
2. **The reporter integrates.** GFP matures slowly and is stable, so
   fluorescence accumulates in visible steps, one per pulse.
3. **Aneugens enter by a different door.** Spindle interference is not a DNA
   lesion, so it gets its own state variable and reaches p53 through mitotic
   surveillance. Projecting it onto `lesions` would make a future comet
   readout report strand breaks that do not exist.
4. **Cytostasis is endogenous.** p53 -> p21 arrests the cycle, so the assay
   partly causes its own validity-gate failures. In the bacterial core every
   growth effect came from outside.

## Three results that came out of running it

**The reporter integrates promoter activity, not p53.** Tested rather than
asserted, by ranking three predictors of induced fluorescence on the ratio-CV
of each: promoter integral 0.021, p53 integral 0.077, p53 peak 0.508. The
GADD45a promoter is a steep switch, so p53 excursions below its threshold add
area to the p53 integral and almost nothing to the fluorescence — the assay
is structurally blind to low-level sustained p53.

**Growth arrest manufactures induction.** A compound with *zero* DNA lesions
and a demonstrably un-induced SOS promoter (0.51x control — arrest raises
LexA, which deepens repression) still produces an induction ratio of 2.74 in
the umu core, because beta-gal stops being diluted. Only the growth gate
separates that from signal. Both runs keep this as a standing check so it
cannot quietly regress.

**Reporter lines are the wrong instrument for aneugens** — this came out
against expectation. The mitotic-surveillance route is real but weak, and it
is self-limiting: the drug that trips surveillance also stops the cycling
that surveillance requires, so the response peaks mid-range and falls. Genuine
GADD45a engagement reaches only 1.14x while the plate reader shows 2.41x, a
2.1x overstatement that is mostly growth artifact. The bacterial core, by
contrast, shows exactly 1.000x — no route at all. The qualitative
architecture point holds; the quantitative claim "aneugens are
reporter-positive" does not, in this model. That is the argument for reading
the aneugenic channel directly at endpoint #3 (micronucleus), rather than
inferring it from a reporter.

Relatedly, the aneugen's *verdict* is set by the protocol, not the biology:
one unchanged simulation scores INCONCLUSIVE at a 0.80 density gate and
POSITIVE at 0.60.

## What the model does and does not claim

The **topology** is standard and mechanistic in both cores: lesion → ssDNA →
RecA* → LexA autocleavage → umuDC de-repression → lacZ → β-gal → ONPG → A410,
and lesion → ATM/ATR → p53 ⇄ Mdm2 → GADD45a → GFP. LexA self-repression and
the p53-Mdm2 loop supply the feedback that sets each basal set-point.

The **rate constants are illustrative order-of-magnitude values (H), not
fitted to any published dataset.** The p53 parameters were chosen by sweep to
put the damaged-state pulse period near the commonly reported few-hour range
(~5.3 h here) and keep the undamaged state stable — that is a coarse sanity
target, not a fit to a time series. Absolute EC values from this code are not
predictions for any chemical. Curve *shape* and
*relative* behaviour between conditions are what it is currently good for.
Fitting the core to a measured dose series is the next honest step.

The demo compounds are likewise illustrative archetypes; the "-like" names
mark the behaviour being represented, not a claim about those chemicals.

## Two things the model gets right by construction

**Growth normalisation and the validity gate.** A410 scales with total enzyme
in the well, i.e. per-cell enzyme × density, so the induction ratio is formed
from growth-normalised specific activity and every row carries its own growth
factor. `call_result` returns three outcomes, not two: when the threshold is
crossed only in wells that failed the gate, the answer is **INCONCLUSIVE**.
Collapsing that into "negative" is how cytotoxic compounds get mis-scored.

**Lethal vs bacteriostatic toxicity are separate inputs.** They were one
multiplier at first, and that was wrong: suppressing translation also
suppresses the reporter, so every strong cytotoxicant came out non-inducing.
With them split, the model reproduces the opposite artifact too — a
growth-arrested but still-translating culture stops diluting its β-gal and
shows a wildly inflated induction ratio (IR > 25 in the `masked
genotoxicant` series). Both failure modes are real, and the gate is what
separates them from signal.

## Structural checks

Thirteen in total (6 + 7), none of them a curve fit. They assert properties of
the wiring and fail loudly if a layer is connected wrongly.

`run_umu.py` (6): direct-acting compound positive without S9; promutagen call
flips on metabolic activation; aneugen never engages the SOS pathway
(lesions 0, RecA* 0, promoter never induced); the apparent aneugen induction
is a growth-arrest artifact caught by the gate; a pure cytotoxicant is not
called positive; a genotoxicant clearing threshold only in gated wells returns
INCONCLUSIVE rather than NEGATIVE.

`run_p53.py` (7): undamaged state quiet / damaged state pulsing; fluorescence
tracks the promoter integral; direct-acting genotoxicant positive; aneugen
reaches p53 by a route the bacterial core lacks; the reporter overstates the
aneugen response; the aneugen verdict is set by the density gate; a pure
cytotoxicant is not called positive.

Where the failures were fixed matters. Two were fixed in the *assertion*,
because the assertion tested the wrong thing — raw induction ratio where the
claim was about pathway engagement. One was fixed in the *data* (a
"near-silent without S9" promutagen that was 2% active at 50 µM is not
near-silent). Three were fixed in the *model*, and are the two entries above
plus the bounded mitotic-stress variable: modelled as an accumulator it
saturated at the lowest dose tested, so every aneugen returned an identical
induction ratio across three decades and reported an EC extrapolated below
the lowest dose.

## Known limits

- ONPG substrate depletion is not modelled (assumes saturating substrate);
  at very high enzyme this overestimates A410.
- The p53 core has no cell-to-cell variability. Real p53 pulses are
  asynchronous across a population, so a population-average reporter reading
  is smoother than this deterministic single-cell trace implies. Nothing here
  reproduces the single-cell heterogeneity that a real reporter line shows.
- In the SOS core LexA is cleared only by growth dilution, which is why
  arresting the culture deepens repression. Adding proteolysis would damp
  that coupling.
- S9 is a single scalar activation factor. Real S9 composition varies by
  batch and is a major source of false negatives; this is the layer most in
  need of replacement before any absolute claim is made.
- `QSARSource` is a declared seam, deliberately unimplemented — training an
  upstream on assay-outcome labels and then feeding it into a model of that
  same assay is circular. Channels derivable from first principles
  (electrophilicity, tubulin binding, Top2 pharmacophore) are the ones to
  implement first.
