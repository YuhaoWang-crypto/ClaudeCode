# genotox — virtual in-vitro genotoxicity assays

**Step 1** — the **umu test** (SOS/umuDC-lacZ, *S. typhimurium*
TA1535/pSK1002).
**Step 2** — a mammalian **GADD45a-GFP reporter line** (p53/Mdm2 delayed
feedback).
**Step 3** — **comet + micronucleus** from one cytogenetic core, sharing the
p53 sub-model with step 2.

Each step reuses the upstream layer, the readout registry and the decision
layer unchanged; only the core is new, which was the point of the split.

```bash
pip install numpy scipy matplotlib
python3 -m genotox.run_umu        # step 1: tables + checks + figure
python3 -m genotox.run_p53        # step 2: pulses, cross-endpoint, checks
python3 -m genotox.run_comet_mn   # step 3: mechanism table + 2 experiments
python3 -m genotox.run_identifiability   # diagnostics: what a fit could recover
```

Writes four figures into `figures/`. 29 structural checks (6 + 7 + 8 + 8).

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
| endpoint #3 (comet, micronucleus) | new `SignalCore` + 2 readouts — *done* | source, decision |
| a different protocol threshold | a `Protocol` instance | all biology |
| ask what a fit could recover | nothing — `identifiability.py` reads any core | everything |

Readouts hand the decision layer exactly two generic names, `signal` and
`biomass`; instrument-specific ones (`A410`, `fluorescence_per_cell`) ride
alongside for reporting.  Without that, the decision layer would have to know
that umu means absorbance and a reporter line means fluorescence, and the
layer split would be fictional.

## Why the interface is a vector

`DamageFlux` carries eight channels, not one potency number. Each core
declares its own per-channel weights, so one upstream prediction produces
*different* answers at different endpoints — which is the actual situation:

- `SOSCore.sos_weights` puts **0.0** on the `aneugenic` channel, so a spindle
  poison is structurally incapable of being called SOS-positive. A
  micronucleus core would weight that same channel heavily. With a scalar
  interface, one of those two must be wrong.
- `lesions` is exposed as an observable specifically so a comet readout can
  be calibrated against the same quantity rather than re-deriving damage.

Two things the channel vector is **not**. It is not eight quantities in one
unit: the DNA-lesion channels are lesion-equivalents per cell per minute,
while `topo` is trapped-complex occupancy and `aneugenic` is spindle
engagement. And a channel left at zero is an assertion, not a default — a
prediction that never examined a channel belongs in `DamageFlux.unknown`,
which travels with the result so a NEGATIVE verdict reports its own coverage
(`"NEGATIVE ...; but 3 channel(s) this core reads had no evidence either
way"`). A channel the core weights at zero raises no caveat, because
ignorance about it cannot mislead that endpoint.

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

## Step 3: two assays, one core

Comet and CBMN are not different biologies — they are the same damaged cells
measured with different instruments at different times. So step 3 adds **one**
core and **two** readouts, and the two assays differ only in exposure length
and readout:

| | comet | micronucleus (CBMN) |
|---|---|---|
| exposure | 4 h | 36 h |
| requires division | no | **yes** |
| measures | break density in situ | what survived into the next mitosis |
| validity gate | relative viability | relative proliferation (CBPI−1) |

The lesion pool had to be split to support this. A single lumped `lesions`
variable drives a transcriptional reporter fine, but the alkaline comet does
**not** see bulky adducts — what migrates is strand breaks and alkali-labile
sites. So the core tracks bulky adducts, alkyl adducts (themselves
alkali-labile), SSB, DSB, acentric fragments and spindle engagement
separately.

Centromere status is carried through the whole chain: acentric fragments give
centromere-negative micronuclei, lagging whole chromosomes give
centromere-positive ones.

## What step 3 buys: naming the mechanism

Step 2's finding was that a transcriptional reporter cannot distinguish a
clastogen from an aneugen — both just raise p53, and most of the apparent
aneugen signal was growth artifact. Two observables fix that. Comet answers
"was the DNA broken"; centromere status answers "was a whole chromosome
lost":

| compound | comet | MN | %MN C+ | mechanism |
|---|---|---|---|---|
| bulky adduct former | POS | POS | 0.0 | clastogenic |
| alkylating agent | POS | POS | 0.0 | clastogenic |
| direct clastogen | POS | POS | 0.0 | clastogenic |
| **aneugen** | **NEG** | **POS** | **98.8** | **aneugenic** |
| **mixed clastogen/aneugen** | **POS** | **POS** | **80.1** | **MIXED — both components** |
| non-genotoxic cytotoxicant | NEG | NEG | — | negative |

The aneugen is comet-negative at every dose because nothing is broken — that
is structural, not tuned: the `aneugenic` channel reaches neither break pool.

## Two more results from running it

**The comet scores excision intermediates, not adducts.** Sweeping the NER
rate at fixed dose, %tail rises from 24 to 58 as repair gets faster, while the
adduct burden falls from 129 to 2.7 lesions/cell. A repair-*proficient* cell
looks more damaged than a repair-deficient one, because the assay is counting
the incisions repair makes, not the lesions it removes. This is the opposite
of the intuition that repair protects, and it is a manipulation of the model,
not a curve fit.

**The micronucleus turnover is a cytotoxicity artifact, not a p53 effect.**
The MN dose-response peaks and falls, which is expected — a cell that does
not divide cannot form a micronucleus. The first version of this attributed
the fall to p53 → p21 arrest. Disabling each cause in turn says otherwise:
the signal falls 84% from peak to top dose normally, 85% with p53 arrest
disabled, and **0%** with lethality disabled (where it rises monotonically to
363/1000). So the turnover is driven almost entirely by cells dying, and p53
arrest contributes ~8%. Reading the high-dose decline as "less genotoxic"
would be exactly backwards — which is what the protocol's cytostasis limit
exists to prevent.

## Diagnostics: what could a fit to these readouts actually recover?

The stated next step has been "fit a core to a measured dose series". The
honest step *before* that is to ask how much such a fit could learn, because
least squares returns a confident number for a parameter the data cannot see.
`identifiability.py` answers it by finite-difference relative sensitivities
(`d log observable / d log parameter`), an SVD of that matrix, and a count of
directions estimable to better than a factor of 1.65 at 5% relative noise —
a count of *directions*, not a condition number, because a condition number
says a model is sloppy without saying what is measurable.

| Readout | parameters | estimable at 5% noise |
|---|---:|---:|
| umu, endpoint (IR + gate) | 17 | **6** |
| umu, time-resolved (4 sampling times) | 17 | **6** |
| GADD45a-GFP, endpoint | 17 | **6** |
| comet alone | 12 | 2 |
| micronucleus alone | 12 | 3 |
| comet + micronucleus | 12 | **5** |
| ...plus an alkylating agent in the design | 12 | **6** |

Three results worth stating separately.

**Some parameters are invisible by construction, not by noise.** The umu
readout is a *ratio*, so anything that scales the reporter linearly cancels
exactly: transcription rate, translation rate, and the instrument gain.
Scaling `k_txn` or `k_tsl` over a 10⁵-fold range moves the induction ratio by
<3×10⁻⁸ (solver tolerance), and a 1000× gain change leaves it identical to 12
decimal places. A fitter handed these would report values for all three.

**The model carries an exact symmetry that was not visible in the equations.**
`(beta_lexA, K_lexA)` scaled *together* by any factor leaves every observable
unchanged (drift ~10⁻⁸ over a 0.5×–4× scan) while basal LexA scales exactly
with the factor. The LexA ODE is homogeneous of degree one under that
rescaling and the promoter only ever sees `L / K_umu`, with `K_umu`
proportional to basal LexA. The concentration scale of LexA is unobservable,
so one of those two parameters should be fixed by convention rather than
fitted.

**More sampling does not fix a symmetry.** Sampling the reporter at four
times instead of one raises every visible singular value by ~1.5× — better
precision — but the estimable count stays at 6, because the flat directions
are symmetries and no sampling schedule removes one. This is worth
distinguishing from the usual demonstration, where extra time points *do*
raise the rank; there the unidentifiability was an artifact of normalising to
the endpoint, which more data genuinely repairs.

**Identifiability is a property of the experiment, not the model.** `k_ber`
(alkyl-adduct excision) has *exactly zero* sensitivity in the cytogenetic
core — not because the model is degenerate but because the probe compound
makes no alkyl adducts. Adding an alkylating agent to the design takes its
sensitivity from 0 to 1.65 and the estimable count from 5 to 6.

All of this is local to one operating point and says nothing about whether
the model is *right* — only about what a fit to this readout could and could
not learn.

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

Twenty-nine in total (6 + 7 + 8 + 8), none of them a curve fit. They assert properties of
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

`run_comet_mn.py` (8): aneugen comet-negative at every dose; aneugen's
micronuclei centromere-positive; clastogen's centromere-negative; the two
mechanisms are told apart; pure cytotoxicant negative in both; comet scores
excision intermediates; MN turnover is a cytotoxicity artifact rather than a
p53-arrest effect — this last one asserts the *attribution*, with all three
decomposition numbers, not just the shape of the curve; and a two-mechanism
compound is reported as MIXED rather than forced into one class.

`run_identifiability.py` (8): the reporter rate constants and the instrument
gain are invisible to a ratio readout; `(beta_lexA, K_lexA)` is an exact
continuous symmetry; each of the three endpoints constrains far fewer
directions than it has parameters; a parameter is identifiable only if the
probe compound exercises it; and extra sampling buys precision rather than
new directions.

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
- The cytogenetic core is deterministic and population-averaged. Real comet
  and MN data are distributions over single cells — %tail is scored per
  nucleoid and micronuclei are counted per binucleate — so the spread, the
  "hedgehog" fraction and the overdispersion that drive real statistical
  power are absent here.
- Micronucleus formation is evaluated as a per-division probability from the
  instantaneous fragment burden, not by tracking individual chromosomes
  through an explicit mitosis.
- Comet %tail saturates by construction (ceiling 92%), which is right, but the
  model has no lysis/electrophoresis step, so slide-to-slide and
  condition-to-condition variation in migration is not represented.
- The identifiability analysis is **local** — a linearisation about one
  operating point, with one probe compound per design. A parameter estimable
  there may not be estimable elsewhere in parameter space, and profile
  likelihoods or a global method would be needed to say more.
- S9 is a single scalar activation factor. Real S9 composition varies by
  batch and is a major source of false negatives; this is the layer most in
  need of replacement before any absolute claim is made.
- `QSARSource` is a declared seam, deliberately unimplemented. The earlier
  wording here ("training on endpoint labels is circular") was too broad and
  has been corrected: training end-to-end on assay outcomes is a perfectly
  good way to *predict* those outcomes. What is not legitimate is relabelling
  such a prediction as a measured lesion flux and feeding it to a mechanistic
  model of the same endpoint — the mechanism then adds no information and the
  apparent agreement is the training signal coming back around. Channels with
  evidence that is not the endpoint itself (electrophilicity, tubulin
  binding, a Top2 pharmacophore) are the ones to implement first.
