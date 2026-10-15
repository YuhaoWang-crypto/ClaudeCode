---
name: state-restoration
description: >-
  Define a cell / tissue / indication state by its perturbation-recovery
  behaviour and ask what returns it to balance. Five steps on omics (or any
  multivariate) data with a reference state and a perturbed or disease state,
  optionally followed over time under an intervention: (1) STEADY STATE (balance
  reference, noise floor, disease displacement, persistence), (2) ORDER
  PARAMETERS (disease coordinate phi, module coordinates, label-free collective
  modes), (3) MINIMAL CORE NETWORK (disrupted modules, couplings, core genes
  replicated across groups), (4) module-level DYNAMICS dy/dt = J(y - yD) + b0 +
  b u (recovery time constants, slow mode = residual, b = where the
  intervention acts), (5) MINIMAL ADJUSTMENT (fewest modules to restore the
  steady state; structural driver set). Use for a new indication, treatment
  trial, challenge study, stimulus time course or knockout; for "does
  treatment restore the balance state", "which components are core", "what
  is the minimal intervention", "is dynamic information better than static".
  Built on netjac (statepipe.py, recovery.py, statmech.py); validated on
  GEO datasets with known positives. Enforces ✅ rigorous / ⚠️ hypothesis
  labels.
---

# State restoration: steady state → order parameters → core network → dynamics → minimal adjustment

Engine: `netjac/statepipe.py`, with helpers in `netjac/recovery.py` (nested M0→M3 tests, driver nodes, greedy restoration) and `netjac/statmech.py` (collective modes, shrinkage couplings, MRA). Data loader: `netjac/data/geo.py`, which reads any GEO series matrix and GPL annotation.

## 0. Before any computation: what kind of state, what data

Answer these four questions first and write the answers down.

1. **Is the state functional or compositional?**
   - *Functional* means how the object answers an input: routing, susceptibility, reserve, a compensated feedback. Static data CANNOT define it; you need a standardised perturbation (stimulus, challenge, treatment).
   - *Compositional* means what the object is made of: genotype, cell type, an integrator such as HbA1c. Static data usually suffice, and adding dynamics adds nothing or makes it worse.

   Evidence:

   | Case | Static | With perturbation |
   |---|---|---|
   | D12 macrophage routing | AUC 0.50 | 1.00 |
   | D12 macrophage genotype | 1.00 | — |
   | D13 influenza susceptibility | 0.59–0.72 | 0.93 |
   | D11 NHANES OGTT | — | ΔC ≈ 0 |

2. **What is the balance reference?** In order of preference:
   - the subject's own unaffected sample (non-lesional tissue, pre-challenge blood, same plate);
   - control tissue from the same study;
   - healthy volunteers from another batch (flag it).
3. **What is the known positive, defined independently of the readout?** Examples: clinical response score, textbook pathway dependence, drug vs placebo, ground-truth simulation. Also define a negative control: placebo, non-targeted module, housekeeping genes.
4. **Which of the five steps can the data support?**

| Data | Steps | Example |
|---|---|---|
| Reference + disease, one time point | 1–3 | CPTAC tumour vs normal (D9) |
| Reference + disease + ≥ 2 follow-up times under an intervention | 1–5 (steady state identifiable) | psoriasis weeks 0–16 (D14/D15) |
| Only 1 follow-up time | 1–4; y* and step 5 flagged unreliable | IBD before/after (D15) |
| Perturbation screen (each module / gene perturbed) | causal J by MRA → causal levers | Perturb-seq (D9–D11), Norman doubles |
| Kinetic model available | exact J, fragility, ground-truth levers | Teusink (D11), network-jacobian skill |

## 1–5. Run

```python
from netjac import statepipe as sp
D = sp.StateData(X=..., features=..., subject=..., kind=..., time=...,   # kind: "ref" | "state"
                 arm={subj: "drug"|"placebo"}, outcome={subj: bool}, group={subj: g})
card = sp.run(D, modules)            # modules: {name: [genes]} or sp.auto_modules(D, sp.steady_state(D))
sp.print_card(card, "my indication")
```

To start a new indication, use the template `assets/run_state_pipeline.py` (GEO or CSV input). It follows the conventions of the worked example `netjac/demos/d15_state_pipeline.py`. Method details and formulas are in `reference/workflow.md`.

What each step returns, and how to read it:

| Step | Output | Read it as |
|---|---|---|
| 1 Steady state | disease genes; displacement in noise units; excluded subjects; control-arm restoration over time | Is the disease state far from balance (≫ 1 noise unit)? Is it itself stable without intervention (control arm ≈ 0)? |
| 2 Order parameters | φ per sample (1 = own disease state, 0 = balance); module coordinates; number of collective modes over the shuffle null; \|r(PC, φ)\| | If PC1 ≈ φ, the state is ONE collective mode, and one high-loading prior gene is a sufficient readout. If not, inter-individual heterogeneity dominates |
| 3 Core network | modules actually disrupted; partial correlations; core genes replicated in every group | The minimal set of components that describes the state. Replication across groups is mandatory |
| 4 Dynamics | J, time constants, slow mode, b0 (natural course), b (intervention), y* or y(T) | Slow mode = what will remain (residual disease). b = mechanism. CV error must beat "no change" |
| 5 Minimal adjustment | greedy module set to 90 %; single-module levers with propagation; structural drivers; applied intervention as a lever | With J from one observational intervention: BOTTLENECKS, not causal drivers. With a measured response library: causal levers (see below) |

## Step 5 needs a response library, not a better formula

The effect direction of lever j in `minimal_adjustment` is J⁻¹e_j. With the MRA
Jacobian J = −R⁻¹ this equals **−R[:, j]**, the measured response to perturbing
module j (verified numerically in `statepipe.causal_J`). So the question is never
"is the formalism right" but "does J⁻¹ reproduce measured perturbation responses".

- **For a skin indication, a default library ships with the skill.**
  `netjac.data.keratinocyte.lever_library()` returns nine measured cytokine
  responses in primary keratinocytes: IFN-α, IFN-γ, **IL-4**, IL-13, IL-17A,
  TNF (GSE36287, 3 donors) plus IL-9, IL-22, IL-26 (GSE305982). Run it from
  `assets/keratinocyte_levers.py`; verify the library with
  `python3 -m netjac.data.keratinocyte`, which re-runs the published checks
  (psoriasis top lever IL-17A; IL-17A +0.175 [0.149, 0.202] and IL-4 +0.103
  [0.074, 0.129] pass their crossed contrast). `VALIDATED` in that module
  lists the approved drug targets per indication, so the crossed check runs
  automatically. **IL-4 is in the default deliberately**: with IL-13 as the
  only type-2 lever the atopic-dermatitis question cannot be answered
  (IL-13 +0.007 n.s.). The sign of the fitted coefficient says what to do —
  positive means the lesion lies along that input, so restoration is blockade.
- **If you have perturbation data** (each module perturbed individually: CRISPR
  screens, cytokine or drug panels hitting different modules, several
  interventions): build R, then use `rank_levers` /
  `minimal_adjustment_from_responses`. Report lever rankings, **and state the
  cell type the library came from** — a library from the wrong cell type scored
  0/4 where the right one scored the clinically validated target.
- **If you have one observational intervention** (a single drug time course):
  you may report bottlenecks and the slow mode, explicitly labelled as such.
  Do NOT report a causal lever ranking. This is the check D15 failed.

Ground truth for the causal route (`demos/d16_causal_lever.py`, 439 held-out
CRISPRi knockdowns, leave-one-out, perturbed gene excluded from module scores):

| Predictor | Names the module that was actually broken |
|---|---|
| causal (measured response library) | **66 %** [61–70 %] |
| naive (largest shift) | 25 % |
| observational (covariance of controls) | 15 % |
| chance / label-permutation null | 7 % / 10 % (p < 0.0005) |

**The right cell type is part of the requirement** (`demos/d17_keratinocyte_lever.py`).
A public cytokine panel on primary keratinocytes (GSE305982: IL-9/IL-13/IL-17A/
IL-22/IL-26/IFN-γ, 3 h and 24 h) answers the psoriasis lever question D15 got
wrong: the top lever at 24 h is **IL-17A** (cos +0.231 vs IFN-γ +0.132), the
clinically validated target, with a positive coefficient, so restoration means
blockade. The same analysis with the FIBROBLAST library misses in all four
settings. A perturbation library from the wrong cell type is not a substitute.
Replicated and extended with IL-4 and TNF in the library (GSE36287, 3 donors,
`demos/d18_il4_and_combinations.py`): psoriasis → IL-17A again (cos +0.220),
and its crossed specificity becomes significant (+0.174 [0.147, 0.201]). The
AD-specific lever is **IL-4, not IL-13** (+0.097 [0.070, 0.125] vs +0.004 n.s.),
so a panel missing IL-4 cannot answer the AD question — **check that the library
contains the levers the clinic actually uses** before concluding anything.
Caveats kept: atopic dermatitis picks IL-22 at 24 h (IL-13 only at 3 h), and
psoriasis flips to IFN-γ under a rank-based metric, so report the metric and the
time point, and prefer the crossed specificity contrast (IL-13→AD +0.181
[0.154, 0.207]; IL-17A→psoriasis +0.026 [−0.009, +0.061] n.s.) over top-1.

Limits found: per-module hit rate ranges 19–90 %, and **no module-level property
predicts it in advance** (|R_j| p = 0.11, collinearity p = 0.69; n = 15 modules).
Errors land on causally nested modules (mitoribosome → complex I), so top-3 is
78 %. For double perturbations a linear J cannot beat additive superposition,
because d_AB = d_A + d_B exactly in linear response: epistasis is outside the
model class (Norman 2019: AUC 0.89 both ways; best-partner ranking 54 % vs 21 %).
Measured combinations put a number on what is lost (D18, keratinocyte
IL-17 + TNF): the measured combination aligns with the psoriasis lesion better
than the additive prediction (+0.053 [0.031, 0.074]). But treat whole-vector
non-additivity with suspicion — it looked like 65 % and was mostly donor noise
(synergy-vector between-donor r 0.10 vs 0.54 for the response itself). Test
synergy per gene with leave-one-donor-out: it was real and concentrated at
IL19, S100A7A, DEFB4A, SERPINB4 and IL23A, i.e. the disease signature genes.

## Validation gates (do not report a step that fails its gate)

1. **Known positive recovered.** Examples: textbook-dependent stimuli separated; drug > placebo; molecular restoration correlates with clinical improvement.
2. **Negative control flat.** Placebo or untreated restoration ≈ 0 at the early visit. A housekeeping module is undefined (not disrupted). Non-targeted modules do not move more than placebo.
3. **Boring baselines beaten or matched.**
   - State definition: the static readout and the clinical early-improvement readout.
   - Dynamics: "no change" and the same-arm mean trajectory.
   - Minimal readouts: random gene sets of the same size.
4. **Transfer.** Leave a whole group out (strain, virus, drug, plate, tissue) and report the per-group numbers.
5. **Published numbers reproduced where they exist.** Verify every PMID before writing it down.

## Validated cases and what passed

| Case | Data | Result |
|---|---|---|
| D11 | Teusink glycolysis (ground truth) | Hidden feedback defect. Static AUC 0.66 (FP 74 %) → one challenge 0.996 (FP 0.4 %). The greedy restoration set contains the true fix ✅ |
| D11 | K562 Perturb-seq / Norman 2019 | Homeostats transfer to unseen families (0.63 vs static 0.06). Single responses predict which double restores the state (0.89 vs 0.55) ✅ |
| D12 | GSE75210 Myd88⁻/⁻ macrophages | Routing: static 0.50 → trajectory 1.00 (replicate r 0.91). Bypass inputs restore NF-κB at an interferon cost ✅ |
| D13 | GSE52428 influenza challenge → GSE17156 HRV/RSV | Static 0.59–0.72 → 0.93 (leave-strain-out). Core = type I IFN (111 genes, null 0). Illness = one collective mode ✅ |
| D14 | Psoriasis, AD vs placebo, IBD | Restoration tracks clinical scores (ρ 0.58–0.66). Dupilumab 0.63 vs placebo −0.04 (paper: 68.8 % vs −10.5 %). Th2 restored, interferon not ✅ |
| D15 | Pipeline end to end on D14 data | See below |

D15 details:
- Steps 1–3 pass on all three datasets.
- Dynamics beats "no change" by 17–31 % in cross-validated error and matches the same-arm mean.
- Slow mode = barrier / epithelium ✅.
- AD direct input b: T cell / DC first, Th2 second, interferon the only positive (specificity ✅, Th2 not first ⚠️).
- Psoriasis lever check FAILED: barrier ranked first, not IL-17. Levers from a single intervention are bottlenecks ⚠️.

## Labelling rules

- ✅ **rigorous**: a measured statistic on real data with its CI or null, an exact computation, or a decimal-exact model.
- ⚠️ **hypothesis**: a mechanism behind a statistic, an extrapolation beyond the sampling window, a lever ranking from observational J, a cross-batch reference, or a virtual cohort.
- Report negative results as findings. Examples: OGTT adds nothing to mortality prediction; the psoriasis lever check failed; IBD static severity predicts response as well as anything.

## Read next

- `reference/workflow.md`: formulas, parameter choices and decision rules for each step.
- `reference/pitfalls.md`: every way this went wrong before it was fixed (regression to the mean, distance vs projection, unstable J, ratio blow-ups, overlapping cohorts, citation errors).
- `reference/validation.md`: the pre-specified checks and their outcomes, case by case.
- The `network-jacobian` skill: kinetic models, fragility, statistical-physics micro → macro, Perturb-seq MRA (causal J for step 5).
