# genotox — virtual in-vitro genotoxicity assays

Step 1 of the plan: the **umu test** (SOS/umuDC-lacZ, *S. typhimurium*
TA1535/pSK1002) as a runnable mechanistic model, built so that steps 2 and 3
plug in without touching what is already here.

```bash
pip install numpy scipy matplotlib
python3 -m genotox.run_umu        # tables + structural checks + figure
```

Writes `figures/genotox_umu.png`.

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
| endpoint #2 (p53/GADD45a reporter) | new `SignalCore` | source, decision |
| endpoint #3 (comet) | new readout on a core exposing break density | source |
| a different protocol threshold | `doseresponse.py` constants | all biology |

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

## What the model does and does not claim

The **topology** is standard and mechanistic: lesion → ssDNA → RecA* → LexA
autocleavage → umuDC de-repression → lacZ → β-gal → ONPG → A410, with LexA
self-repression supplying the feedback that sets the basal set-point.

The **rate constants are illustrative order-of-magnitude values (H), not
fitted to any published induction dataset.** Absolute EC-IR1.5 numbers from
this code are therefore not predictions for any chemical. Curve *shape* and
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

`run_umu.py` asserts five properties of the wiring, none of which is a curve
fit. They fail loudly if a layer is connected wrongly:

1. a direct-acting compound is positive without S9;
2. a promutagen's call flips negative → positive on metabolic activation;
3. an aneugen is negative, because its channel carries weight 0;
4. a pure cytotoxicant is not called positive;
5. a genotoxicant whose signal only clears threshold in gated wells returns
   INCONCLUSIVE rather than NEGATIVE.

## Known limits

- ONPG substrate depletion is not modelled (assumes saturating substrate);
  at very high enzyme this overestimates A410.
- S9 is a single scalar activation factor. Real S9 composition varies by
  batch and is a major source of false negatives; this is the layer most in
  need of replacement before any absolute claim is made.
- `QSARSource` is a declared seam, deliberately unimplemented — training an
  upstream on assay-outcome labels and then feeding it into a model of that
  same assay is circular. Channels derivable from first principles
  (electrophilicity, tubulin binding, Top2 pharmacophore) are the ones to
  implement first.
