# Composing biological models into longer causal chains

*Why cross-model interfaces break, what to do about it, and a runnable
minimal demo.*

Companion code: [`bioif/`](bioif/) · run `python3 -m bioif.demo` and
`python3 -m bioif.selftest`.

---

## 0. The reframe

The usual framing is "every group has its own data and its own model, so
nothing composes." That is true but not actionable, because it suggests the
fix is standardisation, and standardising all of biology is not a project
anyone finishes.

A sharper framing: **models do not fail to compose at the file format. They
fail at the semantics of the seam.** Moving a tensor is trivial. What breaks
is that each model emits a different *kind* of number, conditioned on a
different context, relative to a different baseline, with no statement of
where it stops being valid. Two models can agree perfectly on `ENSG00000141510`
and still be uncomposable.

That changes the goal. You do not need a universal model. You need a **thin,
typed contract at the seams**, plus the honest admission that a chain of
unvalidated hops produces a hypothesis, not a measurement.

And it changes what "success" means for a long chain. ⚠️ A five-hop chain of
imperfect maps will not give an accurate endpoint prediction, and anyone
promising one is mis-selling. What a well-built chain *can* do, defensibly:

1. propagate uncertainty honestly instead of hiding it,
2. **refuse** when a link is out of domain instead of extrapolating quietly,
3. attribute the endpoint's uncertainty back to individual links, and
4. therefore rank **which single experiment most tightens the chain**.

(4) is the deliverable. A chain is best understood as an experiment-allocation
device, not a prediction device.

---

## 1. Bottlenecks

### 1A. Cross-cutting axes

These apply to *every* pair of models, in every combination. Each is a
distinct failure mode with a distinct detection method and a distinct fix.

| # | Axis | Failure mode at the seam | How it shows up | Mitigation |
|---|---|---|---|---|
| 1 | **Identity** | join on an identifier that is not stable without a build/release/isoform/tautomer | silent row loss, or worse, silent mis-join | namespaced ids with a mandatory `build`; lint at ingest (`check_identity`) |
| 2 | **Quantity semantics** | same name, different thing: `log2FC vs DMSO 6 h` ≠ `log2FC vs isogenic WT` | numbers combine "successfully" and mean nothing | make `reference` part of the type; two quantities are equal only if name+unit+scale+reference match |
| 3 | **Context / conditioning** | model A predicts in K562, model B was fitted in HepG2; joined on gene id | plausible, wrong, unfalsifiable output | context carried on every claim; adapters declare `calibrated_systems` and `requires_context` |
| 4 | **Uncertainty & calibration** | point estimates only; error compounds invisibly across hops | a 5-hop chain reports 3 significant figures | propagate distributions (samples, not moments — the hops are nonlinear); calibrate each model (conformal prediction gives distribution-free intervals cheaply) |
| 5 | **Applicability domain** | no model can say "I don't know", so every model answers everything | confident nonsense far from the training manifold | every adapter returns in-domain / extrapolate / **refuse**; extrapolation widens rather than corrects |
| 6 | **Provenance & versioning** | weights, training set and thresholds not pinned | last quarter's chain is not reproducible; a silent model update changes conclusions | provenance chain on every claim; hash the weights |
| 7 | **Causal direction** | most models are associational; chaining associations does not chain causes | confounded links look predictive in-distribution, then fail under intervention | prefer links fitted on *perturbation* data; mark associational links explicitly; do not present the endpoint as an intervention effect if any link is associational |
| 8 | **Timescale & spatial scale** | binding (µs–ms) → signalling (min) → expression (h) → fitness (days) | a hop implicitly assumes steady state that never obtains | put `time_h` in context; make the steady-state assumption a named, flagged bridge |
| 9 | **Coverage & negative data** | training sets concentrate on well-studied genes, drug-like chemistry, common cell lines; negatives are usually absent or fabricated | headline metrics are driven by the easy majority; the cases you care about are the tail | report performance stratified by data density; treat "no measured negatives" as a domain limit, not a modelling choice |
| 10 | **Evaluation leakage across the seam** | random splits let the same gene family / scaffold / cell line sit on both sides | chain benchmarks look far better than they are | split by gene, scaffold, and cell system; evaluate the *chain*, not each link in isolation |
| 11 | **Operational asymmetry** | a 2-second model chained to a 6-hour GPU model; incompatible licences; one is an API you cannot batch | the chain cannot be run enough times to estimate uncertainty at all | caching + async adapters; treat runtime and licence as interface metadata, not an afterthought |
| 12 | **Reference genome / annotation drift** | GRCh37 vs GRCh38, Ensembl release differences in transcript models | splicing and variant links silently disagree about which exon exists | pin annotation release in `build`; refuse cross-release joins |

### 1B. The layer-crossing matrix

Here the combinations are enumerated. For each crossing: the quantity change,
the *dominant* bottleneck from §1A, and what actually tightens it.

| From → To | Quantity change | Dominant bottleneck | What tightens it |
|---|---|---|---|
| DNA variant → splicing | variant → ΔPSI | (12) annotation drift; (9) exonic variants under-covered | minigene / targeted RNA-seq on the specific variant class |
| DNA variant → expression | variant → log2 expression | (3) tissue/cell context; (7) association vs causation | matched eQTL in the right tissue; allele-specific expression |
| DNA variant → protein structure | variant → ΔΔG, ΔpLDDT | (2) ΔΔG sign/reference conventions differ per tool | one agreed reference state; experimental stability on a subset |
| DNA → chromatin / epigenome | variant → accessibility, methylation | (3) cell-type specificity is the whole signal | matched ATAC/methylation in the same cell type |
| Epigenome → TF occupancy | accessibility → binding | (7) co-occurrence vs causal binding | perturbation (degron on the TF), not correlation |
| Epigenome → transcription | mark → expression | (8) timescale; marks can be consequence, not cause | targeted epigenome editing (dCas9-effector) |
| RNA → protein abundance | log2 RNA → log2 protein | (4)+(7): the single weakest common hop; buffering, autoregulation, complex stoichiometry | matched targeted proteomics on the same samples |
| RNA → pathway activity | expression → pathway score | (2) every signature has a different reference and scaling | fix one signature definition and one baseline; report the baseline |
| Protein abundance → activity | abundance → residual activity | (7) spare capacity, scaffolding, paralogue compensation | knockdown *titration* with a functional readout, not a binary KO |
| Protein structure → function | structure → activity/selectivity | (5) structure prediction confidence ≠ functional confidence | function assays; do not read pLDDT as a functional error bar |
| Small molecule → target affinity | SMILES → pIC50/Ki | (9) chemical-space coverage; (2) assay-format dependence of the "same" IC50 | keep assay annotation on the label; scaffold-split evaluation |
| Affinity → occupancy | IC50 + dose → fraction bound | ✅ arithmetic once dose is present — but (3): **refuses** without a dose | nothing to fix; just never default the dose |
| Occupancy → functional inhibition | fraction bound → residual activity | (8) residence time, resynthesis; (7) occupancy ≠ inhibition | target-engagement assay (CETSA/NanoBRET) plus a pathway readout |
| **Genetic ⇄ chemical perturbation** | KO/KD ⇄ inhibition | (7) the deepest bridge: a degraded protein is not an inhibited one (scaffolding, kinetics) | degron vs inhibitor comparison on the same target; treat agreement as evidence, never as an assumption |
| Pathway state → cell phenotype | activity → growth/viability | (3) cell background dominates; (9) lines ≠ patients | isogenic panels; report per-background, never pooled |
| Compound → morphology/expression signature | structure → profile | (10) batch/plate leakage is severe | batch-aware splits; explicit plate controls |
| Signature → mechanism class | profile → MoA label | (7) similarity is not mechanism | orthogonal perturbation (genetic) matching the chemical hit |
| Cell line → patient | in-vitro effect → clinical effect | (3) the largest context gap in the whole table | ⚠️ do not bridge this computationally; it is a study design question |
| In vitro → in vivo | potency → exposure-adjusted effect | (8) PK/exposure; free fraction | measured PK; free-drug correction |

Two observations worth stating plainly:

- ⚠️ **The RNA→protein and occupancy→function hops are, in most chains,
  where the uncertainty actually lives.** The demo reproduces this: the
  transcript→protein bridge owns roughly half of the endpoint variance while
  the first hop owns almost none.
- ⚠️ **The cell-line → patient crossing is not an interface problem.** No
  contract fixes it. The correct engineering response is a `REFUSE`.

---

## 2. Solution: type the seams, not the world

The design implemented in [`bioif/`](bioif/). It is deliberately small —
about 600 lines — because the leverage is in the contract, not the code.

### 2.1 A six-field claim

Every number that crosses a model boundary carries:

| Field | Kills bottleneck |
|---|---|
| `Entity` — namespace + id + **build** | 1, 12 |
| `Quantity` — name + unit + scale + **reference** | 2 |
| `Context` — system, dose, time, assay, covariates | 3, 8 |
| `Estimate` — a sample distribution, never a point | 4 |
| domain `Verdict` — in / extrapolate / **refuse** | 5 |
| `Provenance` — model, version, calibration set | 6 |

### 2.2 Three kinds of edge, and the distinction that matters

- `COERCION` — algebra (pIC50 ↔ IC50 nM). Lossless, noiseless, **auto-inserted
  by the router**. This is what makes model swapping cheap: a new model that
  reports Ki instead of IC50 needs zero new plumbing.
- `EMPIRICAL` — fitted on data; carries its own residual noise.
- `BRIDGE` — crosses a level of biology on a scientific assumption. The
  assumption is named in the adapter and **flagged onto every downstream
  claim**, so it is still visible on the final number.

A pipeline that does not separate these reports the same confidence for
arithmetic and for a guess.

### 2.3 Convergence currencies

The single most useful primitive. Pick a small number of **shared quantities**
that several modalities are forced through — in the demo,
`residual_target_activity` relative to untreated wild-type.

The convergence condition is strict, and both halves matter:

> the arms must agree on the **entity** and the **quantity**, and whatever
> differed between them — a variant on one side, a compound at a dose on the
> other — must move into **context**, not stay in the identifier.

That is exactly what lets a model fitted on *genetic* perturbation score a
*chemical* one. It is also the point where the genetic⇄chemical bridge
assumption enters, which is why it is flagged rather than buried.

### 2.4 Refusal as a first-class result

An occupancy without a dose is not a number. The adapter declares
`requires_context=('dose_uM',)` and the chain stops with a reason. A pipeline
that cannot refuse will always produce output, which is the same as producing
no information.

### 2.5 Variance attribution → value of information

Re-run the chain with each link's own noise switched off; the drop in endpoint
variance is that link's share. Divide by the cost of the experiment that would
collapse that link, and you have a ranked measurement plan.

### 2.6 Reuse, don't reinvent

| Need | Use |
|---|---|
| variant identity | HGVS + GA4GH VRS |
| gene / transcript / protein | Ensembl (+release), UniProt (+isoform) |
| chemistry | InChIKey, standardised tautomer/salt handling |
| cell system | Cell Ontology, DepMap model ids |
| assay semantics | BAO; units from UO |
| calibrated intervals | conformal prediction (distribution-free, model-agnostic, cheap) |
| provenance | W3C PROV + model cards; hash the weights |
| matrices | AnnData / MuData |
| execution & caching | Nextflow / Snakemake; the contract is orthogonal to the runner |

⚠️ Note what is *not* on this list: a universal schema for all of biology.
The contract governs the seams only.

---

## 3. Ranked plan

Ranked by (expected value × probability it works) ÷ effort. Anything below
the line is worth doing only after the ones above it are real.

| Rank | Option | Effort | Why here |
|---|---|---|---|
| **1** | **Typed seam contract** (entity+quantity+context+provenance) | S | Removes the largest class of silent errors for the least work. Nothing else is trustworthy until identity and reference state are pinned. **Demoed.** |
| **2** | **Refusal + applicability domain per adapter** | S | The cheapest way to stop confident nonsense. Requires only a predicate per model. **Demoed.** |
| **3** | **Distribution propagation + per-model calibration** | S–M | Conformal prediction makes this tractable without retraining anything. Turns a chain from a number into an interval. **Demoed** (calibration itself: not demoed). |
| **4** | **Variance attribution → experiment ranking** | S | The actual deliverable of a long chain. Almost free once (3) exists. **Demoed.** |
| **5** | **Convergence currencies + genetic⇄chemical bridge** | M | Unlocks cross-modality reuse of existing models. High value, real assumption risk — hence flagged, not hidden. **Demoed.** |
| **6** | **Registry + automatic routing** | M | Makes model swapping a registration instead of a rewrite. Pays off only once ≥2 models exist per hop. **Demoed.** |
| **7** | **Leakage-aware chain benchmark** | M–L | Needed before any chain claim is believable. Deliberately not demoed: a benchmark on stub models would be theatre. |
| — | *line* | | |
| 8 | Perturbation-anchored causal layer (fit links on interventions; do-calculus over the chain) | L | Highest ceiling, highest cost. Only meaningful once 1–7 hold. |
| 9 | End-to-end joint training across layers | XL | Needs matched multi-layer data at scale that mostly does not exist. |
| 10 | Full mechanistic ODE coupling between layers | XL | Deepest where it applies; narrow, parameter-hungry. See this repo's `grn_pipeline/` for where that route leads. |

---

## 4. What the demo shows — and what it does not

```bash
python3 -m bioif.demo        # six scenarios
python3 -m bioif.selftest    # 13 interface guarantees
```

⚠️ **Every numeric constant in `bioif/adapters_demo.py` is an illustrative
stub.** No output is a prediction about any real gene, compound or cell line.
What is demonstrated is *interface behaviour*, which is what the ranked plan
above is about.

| Scenario | Demonstrates |
|---|---|
| S0 | a build-less accession is rejected at ingest |
| S1 | 4-hop genetic chain runs; endpoint variance attributed per link |
| S2 | same chain, uncalibrated cell system → extrapolates, interval widens 2.0× **from the declared context mismatch alone** |
| S3 | chemical chain reaches the same endpoint via different adapters, routed automatically, unit algebra inserted for free |
| S3b | both arms land on the same entity *and* quantity; perturbagen moved into context |
| S4 | same chemical chain with the dose omitted → **refuses** |
| S5 | naive point composition returns a number agreeing to ~0.07 — with no interval, no flags, no domain check |
| S6 | ranked experiment plan: variance removed per unit cost |

The S5 contrast is the point of the whole exercise. The naive pipeline and
the typed chain give nearly the same point estimate. Only one of them tells
you the answer spans 1.5 units, rests on three unvalidated cross-level
bridges, and was asked about a cell system no link was ever fitted on.

### Not demoed (deliberately)

- Real models behind any adapter — every one is a stub.
- Conformal calibration (rank 3's other half).
- The leakage-aware benchmark (rank 7). A benchmark over stub models would
  measure nothing.
- Identifier resolution at production scale — `_VARIANT_TO_GENE` is a dict,
  and in a real deployment that component is the most load-bearing piece of
  the system.

### Next demos, in order

1. ~~Replace one adapter with a real model and keep every test green~~ —
   **done, see §5.**
2. ~~Wrap one model in conformal prediction and check the interval's
   empirical coverage~~ — **done, see §6.**
3. ~~Add a second model for the same hop and let the router pick~~ —
   **done, see §7.**
4. Then, and only then, a leakage-aware chain benchmark on a hop with real
   measured truth on both sides.

---

## 5. Swapping in a real source: what broke

```bash
python3 -m bioif.demo_real          # the chain with a real affinity source
python3 -m bioif.real.heterogeneity # the measurement that justifies §2
python3 -m bioif.selftest           # 20/20
```

The stub `affinity-model@0.1-stub` was replaced with measured bioactivity
from the EMBL-EBI ChEMBL REST API (target `CHEMBL2189121`, *GTPase KRas*;
5,000 pChEMBL-bearing records over 2,689 compounds and 206 assays, snapshot
committed under `bioif/real/_snapshot/` with provenance). Everything
downstream of the source is still the illustrative stub set.

**The suite is green at 20/20 — but it did not stay green for free.** Three
of the four changes needed were genuine defects in the original contract that
a stub source structurally could not reveal, because a stub emits exactly one
tidy number of exactly one type with exactly the arity you asked for.

| | Defect the real data exposed | Why stubs hid it | Fix |
|---|---|---|---|
| **A** | A single `PIC50` quantity silently absorbed ChEMBL's `pchembl_value`, which pools IC50, Ki, Kd and EC50 into one column. For the demo compound **11 of 16 records (69%) were not IC50 at all**; across the snapshot 277/5,000 (5.5%). The contract was committing the exact pooling error §1A axis 2 exists to prevent. | the stub only ever emitted a pIC50 | four distinct quantities + `BY_STANDARD_TYPE`; `route(pKd → fitness)` now returns `None` |
| **B** | `Estimate.point(v)` defaulted to **one** sample. A real measurement reports no error bar, so feeding one collapsed the whole Monte Carlo to n=1 and every downstream interval silently became zero-width — the failure the package exists to prevent, inside the package. | stubs always drew N samples | explicit `n`; a measurement is carried as an exact point replicated N times, so any interval downstream is visibly the *adapters'*, not the measurement's |
| **C** | Identity resolution was a dict that always answered. Against real ChEMBL, `"KRAS"` returns **9 candidates spanning 4 target types** (single protein, protein complex, three-way PPI, protein family). | the dict had one key | `resolve_target` refuses unless exactly one candidate survives the organism + target-type constraints |
| **D** | Cheng-Prusoff (pKi → pIC50) had nowhere to declare that it needs `[S]/Km` and a competitive mechanism. | no stub needed it | `requires_covariates`; the bridge refuses on a real record and only fires when someone supplies both on the record |

Defect A is the one worth dwelling on. It is not a coding slip — it is the
contract in §2 being violated by the file that *defines* §2, and it survived
13 passing tests. Only real data with real heterogeneity surfaced it.

### Two real numbers that came out of the swap

**✅ Assay identity moves potency more than readout type does.** Computed over
the snapshot (measured data, not model output):

| comparison | n | median | p90 | max |
|---|---|---|---|---|
| same compound, **different readout type** (IC50 vs Ki vs Kd vs EC50) | 158 compounds | **0.62** | 1.96 | 3.48 |
| same compound, **same readout type, different assay** | 953 pairs | **1.00** | 2.00 | 2.95 |

Units are log10. The rule everyone knows — don't mix IC50 with Ki — is *not*
the dominant term. Holding readout type fixed and changing only the assay is
worse: a median **10× difference in apparent potency for the same compound on
the same target**. That is comparable to the entire uncertainty the stub model
had been assigned (sd 0.60). So assay identity belongs in the context, not in
a column dropped on the way into a model.

**⚠️ The chain saturates, and says so.** Running the same chain once per assay:
1.74 log units of input spread (55× in potency) became **0.290** gene-effect
units of output, because at 0.3 µM every one of those affinities is already
near-saturating. Read as a decision: at this dose, a better affinity number
buys almost nothing, and the measurement worth paying for is further down the
chain. That is the §2.5 value-of-information argument arriving from real data
rather than from a stub tuned to produce it.

### What is still not real

The occupancy → residual-activity → fitness adapters remain stubs, so the
endpoint numbers in `demo_real` are **not** predictions. The spread *between*
them is the result worth reading. `bioif/real/` also adds `epistemics` from
the blueprint's §3 rule 4 — `measured | inferred_association |
calibrated_prediction | mechanistic_hypothesis` — with the invariant that a
chain can only degrade evidence. The real source enters as `measured` and the
endpoint leaves as `mechanistic_hypothesis`, because a bridge sits between
them, and nothing downstream can restore the label.

---

## 6. Conformal prediction: an error bar that was checked

```bash
python3 -m bioif.demo_conformal    # the five experiments below
python3 -m bioif.selftest          # 27/27, ~3 s
```

Every interval elsewhere in this package is **asserted** — a stub declares
`sd=0.55` and the chain believes it. Rank 3 of §3 was the other half of that
promise: replace one assertion with a *checked* interval.

**The hop.** Not an arbitrary one: *may a potency measured in assay A be
reused where assay B is needed?* That is §1A axis 3 in its most concrete
form, it is the question a chain must answer before reusing any number, and
uniquely among the hops in this repo it has real labels — compounds measured
in both assays. Data is the committed ChEMBL snapshot, IC50 records only (so
§5's readout-type confound cannot leak in), predictor is the median shift
between the two assays, and **every split is on compound, never on row**.

### 6.0 A gate that had to exist first

Of 32 candidate assay pairs with ≥10 shared compounds, **13 are not
continuous measurements at all.** They are patent potency bands at one-log
spacing arriving as a float at the band midpoint: in the worst case 73% of
source values are literally the same number, across 4 distinct levels for
n=363. Nothing in the column types says so — only the tie fraction does.
`classify()` refuses them, and `build_transfer()` raises rather than fitting
a regression to an ordinal. That leaves **19 usable pairs, 496 paired
points.**

⚠️ This is a fourth defect class that only real data exposes, and it is not
specific to ChEMBL: any "standard_value" column can carry binned data.

### 6.1 Within a pair, conformal covers and the usual Gaussian does not

α = 0.2 (nominal 0.80), 300 compound-level splits per pair:

| | mean coverage | pairs at/above nominal | typical width |
|---|---|---|---|
| split conformal | **0.830** | **17/18** | 1.0–3.5 |
| Gaussian ±1.28·sd | 0.671 | 3/18 | 0.6–2.0 |

On the one well-powered pair (`CHEMBL5737243→CHEMBL5737244`, the *same*
coupled nucleotide-exchange assay at 2 h vs 20 h, n=226), conformal is
essentially exact at both levels:

| α | nominal | empirical coverage | width |
|---|---|---|---|
| 0.2 | 0.80 | **0.809** | 0.42 |
| 0.1 | 0.90 | **0.910** | 0.66 |

✅ The guarantee holds. ⚠️ It is paid for in width: conformal intervals run
~1.6–2× wider than the Gaussian ones. That is the honest trade — the Gaussian
interval is narrower *and wrong*, under-covering by 10–25 points because
residuals at these sample sizes are skewed and heavy-tailed.

Note also the fitted shift on that pair: **+0.70 log units**. The same assay,
run for 20 h instead of 2 h, reads 5× more potent. "Reuse the number as-is"
is not a neutral default; it is a 5× error.

### 6.2 When the guarantee cannot be bought, the method says so

At α = 0.1, **17 of 18 pairs cannot support an interval at all**: with ~5–7
calibration compounds, `ceil((n+1)·0.9) > n`, so the conformal quantile is
`+∞`. The implementation returns infinity rather than the largest observed
residual. Returning the latter would be claiming a guarantee it does not
have — the same failure mode as an adapter that cannot refuse.

### 6.3 The guarantee does not survive a change of assay pair

Leave-one-assay-pair-out: calibrate on every *other* pair pooled, test on the
held-out one.

| | nominal 0.80 |
|---|---|
| mean coverage across 19 held-out pairs | **0.612** |
| pairs under-covered | **12/19** |
| worst case | **0.002** |

Conformal guarantees coverage *under exchangeability*. Two assay pairs are
not exchangeable — their systematic shifts range from −1.76 to +1.76 log
units — so the guarantee is void, and **the failure is total rather than
graceful**: near-zero coverage, not a slightly optimistic interval.

This is the applicability-domain argument of §2.4 turned into a number. It is
why `ConformalAssayTransfer` **refuses** an uncalibrated assay pair instead of
widening for it: there is no honest width to widen to.

### 6.4 Marginal coverage is met while most classes fail

Pool all pairs into one calibration set, then score coverage per pair:

| | nominal 0.80 |
|---|---|
| pooled, **marginal** coverage over all test points | **0.799** ✅ |
| pooled, pairs under-covered | **11/18**, worst **0.009** |
| Mondrian (per-pair calibration), pairs under-covered | **0/18**, worst **0.804** ✅ |

The pooled interval keeps its promise *on average* and is useless for most
individual pairs. This is the standard marginal-vs-conditional coverage gap,
and the reason it matters here is what repairs it: **the conditioning class
is the assay — the context field the contract already carries.** The
statistics and the interface design land on the same object.

### 6.5 The resulting adapter

`bioif/real/transfer_adapter.py` ships the calibration:

- **kind = EMPIRICAL**, so a `measured` input leaves as
  `calibrated_prediction`. A measurement transferred to another assay is not
  a measurement of that assay — and that demotion is now backed by a coverage
  number rather than by convention.
- **consumes and produces the same quantity**, changing only `context.assay`,
  so routing never inserts it: an assay transfer is a deliberate act and must
  be asked for (tested).
- **refuses** an uncalibrated source assay, **extrapolates** outside the
  calibrated input range, and **refuses** an α its calibration set cannot buy.

Worked: a compound reading pIC50 6.00 in `CHEMBL5737243` gets a guaranteed
80% interval of **[6.47, 6.93]** in `CHEMBL5737244`. In the full chain the
transfer link owns **28.9%** of the endpoint variance — and it is the only
link in that chain whose share rests on a measured coverage figure rather
than an assumed sd.

### What this does and does not establish

✅ Measured, on held-out compounds: within-pair coverage, the Gaussian
baseline's failure, the cross-pair collapse, the marginal/conditional gap, and
the binned-data prevalence.

⚠️ All of it is one target (KRAS) and its assays in one ChEMBL snapshot. It is
not a general claim about conformal prediction, about ChEMBL, or about assay
transfer in other target classes. The downstream occupancy/activity/fitness
adapters remain stubs, so chain endpoints are still hypotheses. And conformal
calibrates an interval — it does not make the underlying predictor correct.

---

## 7. Two models for one hop, and letting the registry choose

```bash
python3 -m bioif.demo_routing     # R1–R5 below
python3 -m bioif.selftest         # 34/34, ~5 s
```

Until now every hop had exactly one implementation, so "which model" was
never a question the interface had to answer. It is the question that matters
most in practice: several groups model the same step, they disagree, and a
pipeline has to pick.

The assay-transfer hop now carries **twelve competing adapters** — four source
assays × three predictors (`identity` y=x, `shift` y=x+b, `linear` y=a+bx),
all fitted on real ChEMBL pairs, all wrapped in the same split-conformal
procedure. Two contract changes make the competition well-posed:

- the transfer **produces a different quantity** (`pIC50` → `pIC50 in the
  reference assay`), so "which assay is this number from?" is a type error
  rather than a footnote, and a raw measurement cannot reach the downstream
  model without an explicit calibrated transfer;
- `Adapter.calibrated_width(alpha)` is a **declared property**, so the
  registry can rank on it. Algebra declares 0; an uncalibrated model declares
  infinity and loses to any model that has been calibrated.

Selection is `in-domain first, then narrowest calibrated interval, then kind`,
evaluated **per claim**. The domain test alone eliminates 11 of 12 (only one
source assay is ever in domain); width settles the rest.

### 7.1 A bug the demo caught, and the headline it reversed

The first version of this section reported that in **14 of 18** assay pairs a
rival model's prediction fell outside the selected model's 80% interval, with
model disagreement running at **74%** of the interval width. **That was
wrong, and it was my bug.**

`ConformalAssayTransfer` emits `predict(x) + a resampled calibration
residual`. Those residuals carry any bias in `predict` with the opposite sign,
so the emitted distribution is *already de-biased* and is **not** centred on
`predict(x)`. Scoring it with a symmetric interval around `predict(x)` — the
textbook |residual| form — punished a model for a bias its own output no
longer had, and comparing raw `predict(x)` values across models measured a
gap that the adapters did not actually exhibit.

The demo surfaced it: `identity` and `shift` were returning **byte-identical
endpoints** while advertising widths that differed 4×. Fixed by scoring on
**two-sided signed** conformal offsets, which are exactly the central
(1−α) region of what the adapter emits.

⚠️ §6's numbers are unaffected — `conformal.py` still uses the symmetric
variant, which is a valid procedure and is what §6 measured. The two forms
now coexist with their roles documented. One real cost: the two-sided form
needs **n_cal ≥ 9** at α = 0.2 against **n_cal ≥ 4** for the symmetric one,
so it supports 7 pairs rather than 18.

### 7.2 Validity does not separate models; width does — and barely

α = 0.2, 300 compound-level splits, 7 pairs with enough calibration data:

| predictor | mean coverage | mean width | pairs valid |
|---|---|---|---|
| identity | 0.815 | 1.31 | 7/7 |
| shift | 0.814 | 1.31 | 7/7 |
| linear | 0.815 | 1.54 | 7/7 |

Two things follow.

**Identity and shift are the same model.** Not similar — identical widths,
identical centres, identical intervals, verified to 1e-9 in the test suite,
while their raw predictors differ by 0.70 log units. Conformalisation absorbs
a constant bias, so "reuse the number as measured" and "reuse it with the
median offset applied" are one conformalised model. Three candidate
predictors, **two** distinct models.

**The flexible model wins only where there is data to fit it.** `linear` is
narrower on the well-powered pair (0.39 vs 0.43, n=226) and wider on all six
small ones (1.54 vs 1.31 on average) — its extra parameter costs more in
variance than it buys in fit. Across the four registered source assays the
winner splits 3 linear / 1 identity≡shift. No model wins everywhere, which is
precisely why this is a routing decision rather than a library-wide default.

### 7.3 What running the alternatives buys

`bioif/ensemble.py` runs **every** in-domain route, not just the chosen one,
and compares the spread between their endpoints against the width the selected
route advertises. Measured on the transfer hop across 7 pairs:

| | |
|---|---|
| median (model spread ÷ selected interval width) | **0.13** |
| worst pair (max spread ÷ width) | **1.22** |
| pairs where some compound's rival falls outside the interval | 4/7 |

So on this data **model choice is a second-order term**: typically it moves
the answer about 13% of the calibrated interval width — far less than the
interval itself, and nothing like the cross-assay failure of §6.3 where
coverage collapsed to 0.00. The first, buggy version of this analysis
overstated it by roughly 6×.

But it is not uniformly small, and *where* it bites is the useful part.
Sweeping a claim across the input range:

| source assay | n | pIC50 in | selected endpoint | sel. width | spread | flag |
|---|---|---|---|---|---|---|
| CHEMBL5737243 | 226 | 4.50 | −0.012 | 0.378 | 0.003 | — |
| CHEMBL5737243 | 226 | 6.00 | −0.396 | 0.488 | 0.037 | — |
| CHEMBL5737243 | 226 | 7.50 | −0.878 | 0.591 | 0.018 | — |
| CHEMBL4368373 | 13 | 4.50 | −0.152 | 0.489 | 0.482 | — |
| CHEMBL4368373 | 13 | 7.00 | −0.854 | 0.694 | 0.716 | **YES** |
| CHEMBL4368373 | 13 | 7.50 | −0.897 | 0.619 | 0.827 | **YES** |

**The well-calibrated source never trips the flag; the thin one trips it
across most of the range.** At pIC50 7.50 from the n=13 source, identity≡shift
say −0.897 [−1.198, −0.579] and linear says −0.070 [−0.375, +0.173] — the two
models disagree about whether the compound does much of anything, and the
envelope (width 1.372) is more than twice the selected interval (0.619).

Model choice bites exactly where there is not enough data to tell the models
apart, which is also where you least want to be guessing. That is why
`RouteComparison.flag()` fires **per claim** rather than per hop: whether the
choice matters depends on where in the input range you are asking.

### What this does and does not establish

✅ Measured: per-predictor coverage and width, the identity≡shift collapse,
the data-dependence of which model wins, and the spread-vs-width ratios.

⚠️ Same scope caveat as §6 — one target, one snapshot, three deliberately
simple predictors. "Model choice is second-order" is a statement about *these*
models on *this* hop; a hop where the candidates embody different mechanistic
assumptions would not be expected to behave this way. And a registry that
ranks on interval width selects for **efficiency among valid models**, which
is not the same as selecting the model that is right. Running the
alternatives prices that gap; it does not close it.

---

## 8. The interface map's 21 pairings, audited and extended

```bash
python3 -m bioif.map21          # the full audit
python3 -m bioif.demo_pairings  # M1–M5 below
```

Working from the uploaded blueprint's map — seven object classes (C compound,
G DNA/variant, E chromatin, R transcript, P protein, F cell phenotype,
D clinical), 21 directed pairings, each graded **A** verifiable under stated
conditions / **B** needs paired experimental calibration / **C** hypothesis
and ranking only.

### 8.1 What is actually buildable

`bioif/map21.py` audits all 21 (plus F→G, which the blueprint's own Demo B
raises and marks `blocked`). The audit keeps two things apart that a
dataset-availability check would conflate: **is there data** versus **are
there labels on BOTH sides**. Sources were probed from this container, not
assumed.

| status | n | pairings |
|---|---|---|
| **built** | 4 | C→P, C→F, C→G, F→G |
| partial | 3 | G→F, P→F, R→P |
| reachable | 1 | G→P |
| blocked | 9 | C→R, C→E, G→E, G→R, E→R, E→P, E→F, R→F, G→D |
| refused | 5 | C→D, R→D, P→D, E→D, F→D |

Probing mattered: SpliceAI's lookup API and Ensembl REST are **not**
reachable here (so G→R stays a stub), while ChEMBL, UniProt/EBI Proteins,
UCSC and the DepMap portal are.

Two things the audit says that the map alone does not:

1. **Every buildable pairing is one where some public source happens to hold
   labels on both sides — 4 of 21.** The binding constraint across the whole
   map is paired measurement, not models and not formats. That is precisely
   why the blueprint grades so much of it B.
2. **R→P is simultaneously the hop contributing the most endpoint variance in
   our chains (~50%, §4 S1) and the hop with no reachable paired labels.** The
   most load-bearing edge is the least verifiable one. Any honest long-chain
   programme should buy that measurement first; no modelling choice
   substitutes for it.

The five `refused` pairings are all level C and all end in D. That is not a
gap in coverage — `REFUSE` is the correct output. F→D in particular is the
largest context gap in the map and is a study-design question, not an
interface one.

### 8.2 Three pairings built on measured labels

Datasets: Tox21 (MoleculeNet release, 7,453 unique skeletons × 12 endpoints)
and the Hansen Ames benchmark (N=6512, 6,499 skeletons), both snapshotted
with provenance. Joined on **InChIKey skeleton** — a declared identity
choice, recorded rather than buried, since joining on full InChIKey loses
most of the overlap to salt and stereo differences.

| node | split | n test | prev | AUROC | AP | cov\|0 | cov\|1 | abstain |
|---|---|---|---|---|---|---|---|---|
| **C→F** SR-p53 | scaffold | 1,933 | 0.074 | 0.718 | 0.286 | 0.920 | **0.832** | 0.535 |
| **C→G** Ames | scaffold | 1,950 | 0.552 | 0.787 | 0.820 | 0.889 | 0.904 | 0.467 |

Both use **label-conditional (Mondrian) conformal**, not marginal: at 6.2%
prevalence a marginal 90% classifier can meet its guarantee by calling
everything inactive. Calibrating per class forbids that, and the price
surfaces as prediction sets of `{0,1}` — "I don't know" — rather than as
silent failure on the minority class. The p53 AUROC of 0.718 sits where the
plan document's own model does (0.733).

⚠️ **The p53 minority class is covered at 0.832 against a nominal 0.90.** Not
a bug: a scaffold split is *designed* to break exchangeability between
calibration and test chemistry, and conformal's guarantee is conditional on
it. Same mechanism as §6.3's cross-assay collapse, and asserted in the suite
so it stays visible.

### 8.3 The blocked edge, graded

The blueprint's Demo B correctly refuses p53 → mutagenicity for want of a
bridge. **2,064 compounds are assayed in both datasets**, so the edge now has
a 2×2:

| endpoint | n | ep+ %Ames+ | ep− %Ames+ | RR | OR | sens | Fisher p |
|---|---|---|---|---|---|---|---|
| SR-p53 | 1,908 | 55.2% | 30.5% | 1.81 | 2.81 | 0.13 | 5.2e-09 |
| SR-ATAD5 (DDR) | 1,979 | 57.7% | 30.8% | 1.87 | 3.06 | 0.09 | 6.2e-08 |
| SR-MMP (cytotox control) | 1,685 | 42.7% | 29.5% | 1.45 | 1.78 | 0.24 | 1.5e-05 |

Both genotoxic-stress reporters beat the general-cytotoxicity control, so the
association is not merely dying cells — the confound the blueprint warns
about is real but does not explain the signal.

The edge moves from `blocked` to **`inferred_association`**, which required a
new edge kind. `ASSOCIATION` is a measured co-occurrence with no model in
between: empirically grounded, not causal, capping evidence at
`inferred_association` — stronger than a prediction, weaker than a
measurement. The blueprint's four-rung ladder has that rung; bioif did not.

And the edge states its own limits: RR 1.81 with **sensitivity 0.13** means it
orders a list and cannot clear one. As a screen it would miss 87% of what it
is screening for.

### 8.4 Does routing through an intermediate beat going direct?

This is the question the whole programme rests on, and for once there are
labels on **both** ends — 1,908 compounds with SR-p53 and Ames. Four
predictors of Ames, same held-out scaffolds, p53 model trained with the test
scaffolds dropped so nothing leaks through the intermediate:

| arm | AUROC | AP | P@50 | P@100 | P@200 |
|---|---|---|---|---|---|
| **direct** C→G | **0.748** | **0.656** | **0.900** | 0.780 | 0.630 |
| chain C→F→G | 0.601 | 0.433 | 0.540 | 0.490 | 0.450 |
| oracle F→G (*measured* p53) | 0.533 | 0.376 | 0.520 | 0.420 | 0.400 |
| augment C+F→G | 0.748 | 0.655 | 0.900 | 0.780 | 0.625 |

The direct one-hop model wins by **+0.223 AP**, and the decomposition says
why:

- chain vs oracle = **−0.056** → the C→F model's error (slightly *negative*:
  a continuous probability ranks better than the binary label it predicts)
- oracle vs direct = **+0.279** → information the intermediate **throws
  away**, which no improvement to the C→F model can recover
- augment vs direct = **−0.000** → what a *measured* p53 label adds on top of
  structure: **nothing**

So the intermediate is, for this endpoint, redundant given the structure.
Routing a 2048-bit input through a 1-bit node destroys information, and the
destruction is irreversible.

⚠️ **Corrected with a seed spread (added later).** The table above is one
canonical split. Repeated over three splits (`permute=True`, which is the only
way `seed` actually moves a scaffold split in this repo — see §9.5):

| | AP, mean [min, max] over 3 splits |
|---|---|
| direct C→G | 0.637 [0.605, 0.656] |
| chain C→F→G | 0.450 [0.433, 0.479] |
| oracle F→G | 0.372 [0.367, 0.376] |
| augment C+F→G | 0.637 [0.611, 0.655] |
| **direct − chain** | **+0.186 [+0.166, +0.223]** |
| **augment − direct** | **+0.001 [−0.003, +0.005]** |

The direction is decisive — every split positive, mean far outside the range —
but **the +0.223 quoted above was the top of the range, i.e. the most
favourable of the three splits.** The honest effect is **+0.186**. The second
result survives exactly as stated: a measured p53 label adds +0.001
[−0.003, +0.005] on top of structure, which straddles zero.

**The conclusion for long chains:** build them where end-to-end labels are
*missing*, and expect a direct model to beat them wherever such labels exist.
A chain's value is reach and interpretability, not accuracy. That is not a
reason to abandon chains — most of this map has no end-to-end labels at all —
but it does mean a chain should never be sold as the more accurate option.

⚠️ Scope: this is about predicting the **Ames label**. The reporter still says
*why*, which a structure model does not, and Ames is not Comet — the endpoint
the blueprint actually wants. The honest use of this result is as a prior on
experiment design: before paying for reporter panels to predict a genotoxicity
endpoint, check whether structure alone already carries that information. On
this evidence, for Ames, it does.

### 8.5 One more thing the evidence ladder enforces

In the typed chain `C→F→G`, the endpoint comes out as
`calibrated_prediction`, not `inferred_association` — because the ASSOCIATION
edge sits *downstream* of a QSAR, and a chain carries its weakest link.
Putting a measured association after a model prediction does not recover the
association's standing. That is the ladder doing its job.

---

## 9. The cheap joint model, and a ranking rule it falsified

```bash
python3 -m bioif.real.copula            # fit + evaluate
python3 -m bioif.real.copula --matched  # the matched-fold comparison
python3 -m bioif.selftest               # 53/53
```

§8.3 left the p53 → mutagenicity edge conditional on a cytotoxicity readout,
implemented by switching between three measured 2×2 tables. That works and
wastes most of what is known: it uses the cytotoxicity call as a binary
stratifier and nothing else. A Gaussian copula over the Tox21 endpoint vector
does the job properly — 12 logistic marginals plus the 66 free correlations
of a latent multivariate probit, giving `P(p53 | fingerprint, MMP = measured)`
in closed form.

### 9.1 What the measurement is worth

✅ Measured on one scaffold split, 1,512 test compounds with both labels
observed (prevalence 0.058), ρ(p53, MMP) = 0.328:

| model | AUROC | AP | abstention |
|---|---|---|---|
| marginal — no MMP | 0.8081 | 0.2940 | 0.610 |
| **conditional — measured MMP** | **0.8327** | **0.3369** | 0.615 |
| gain from the measurement | **+0.0246** | **+0.0429** | +0.005 |

For **66 parameters**. The joint diffusion model of §7 was handed the same
information through inpainting and gained **+0.000** with 1.6M.

Correlations are estimated from **cross-fitted** out-of-fold marginals, which
is not a detail: in-sample marginals explain each label too well on their own,
which drives the latent thresholds to extremes and attenuates the very
correlation the model exists to use.

### 9.2 A leak in my own verification, measured

The first comparison here scored the plain QSAR against the copula on
compounds chosen by the **copula's** scaffold split, while the plain QSAR used
its own (0.4/0.3/0.3 against 0.5/0.2/0.3). **127 of the 1,512 scoring
compounds (8%) sat in the plain QSAR's training fold** — and that alone showed
it at AP 0.509 against an honest 0.294, making it look as though the copula
lost badly.

An 8% leak inflating AP by ~0.21. §1A axis 10, in this repo's own script,
caught before it reached a commit.

### 9.3 The registry's ranking rule was wrong, and this hop proved it

§7 ranks competing models of a hop on **calibrated width** (abstention rate).
On this hop the conditional model is +0.0246 AUROC better and abstains
**0.005 more** — so a width-only rule prefers the *worse* model. Width is a
proxy for sharpness, and **sharpness is not discrimination**: two models can
be equally sharp and differ in ranking quality.

Fixed by letting an adapter declare a measured AUROC together with the
**identity of the fold it was measured on**, and ranking on discrimination
first — but only among candidates declaring the *same* fold. Mixed or missing
fold ids mean discrimination is not comparable and the rule falls back to
width. That is a refusal to guess, and it is motivated directly by §9.2:
comparing AUROCs across splits is how an 8% leak becomes a 0.21 AP swing.

Two further pieces were needed to make it correct rather than merely
plausible:

**The arms had to be made genuinely comparable** — same split *and* same
scoring subset, since a conditional model cannot be scored where there is
nothing to condition on. Refitting the plain p53 QSAR on the copula's split
and scoring it on the same 1,512 compounds gives AUROC 0.8081 / AP 0.2940 /
abstention 0.610: **numerically identical to the copula's own marginal**,
because it is the same estimator on the same rows. That identity is the
result, not an accident — **the +0.0246 is the value of the measurement, not
of a model family.** So the marginal arm is served from the same artefact with
conditioning switched off, rather than from a second copy of the same
coefficients.

**Declared discrimination had to become claim-aware.** A conditional model is
only the better model when the claim actually carries the extra endpoint;
declaring 0.8327 unconditionally would be a lie by omission. So:

| claim | comparable | declared AUROC (marginal / copula) | registry picks |
|---|---|---|---|
| no MMP measured | yes | 0.8081 / 0.8081 | marginal (tie → width) |
| MMP measured | yes | 0.8081 / **0.8327** | **copula** |

The hop is now selected correctly, per claim, for the stated reason.

### 9.4 What this does not do

⚠️ The copula sharpens the **upstream** p53 estimate. It does **not** replace
the stratified 2×2 on F→G, which is conditional for a different reason — the
association is *absent* among cytotoxic compounds (OR 1.13, p=0.77). Both are
needed: one makes the estimate sharper, the other decides whether that
estimate may move the mutagenicity risk at all.

⚠️ It remains an association model over assay endpoints. ρ = 0.328 says these
readouts co-occur; it does not say cytotoxicity causes reporter activation.

⚠️ One split, one seed, one target class. The +0.0246 has no error bar here —
the seed-repeat discipline demanded of the diffusion rerun in §7 applies to
this number too and has not been done.

### 9.5 Two things the diffusion rerun gave back

**The error bar §9.4 said was missing.** The rerun measures the copula's
conditional gain over 3 seeds at **+0.028 [+0.024, +0.033]** AUROC. The
single-seed +0.0246 reported in §9.1 sits at the bottom of that range, so the
number is reproducible and the point estimate here was, if anything, slightly
conservative.

**A bug in this repo's splitter.** `qsar.scaffold_split` accepted a `seed`
and ignored it: the greedy fill always leaves a deficit, so the random branch
was unreachable and every seed returned the identical split. No reported
number is wrong — all were honest measurements on one canonical split, and
§9.4 already said so — but it was a trap, because anyone seeking a seed
spread through that argument would have got zero spread and concluded their
estimate was stable. The split is now documented as canonical by design, with
an explicit `permute=True` for genuine repeats, and a test asserts both
behaviours. A permuted split moves ~58% of the test fold, which is why a
permutation spread covers training and selection noise fully and split
variance only partly.

### 9.6 Propagating the splitter fix, and one headline corrected

`permute` is now plumbed through `qsar.build`, `copula.fit` and
`chain_vs_direct.run`, so a genuine seed repeat is possible from any of them.
Default `False`, so no committed number changes.

Checked rather than assumed: **`conformal.split_by_compound` — the splitter
§6 and §7 rest on — does vary with its seed** (four distinct test folds over
four seeds). So the 300-repeat coverage figures and the routing results are
intact. The broken seed was only ever in the Tox21-side `scaffold_split`, and
every caller of it passed a fixed seed without looping, so no reported number
was a fake average.

Using the fix immediately found one of my own headlines to be optimistic:
§8.4's "+0.223 AP" is the best of three splits against a mean of **+0.186**
(see the correction in §8.4). The direction and the decisiveness hold; the
magnitude was quoted from the luckiest split. That is the same discipline I
asked of the diffusion rerun, applied to my own number, with the same kind of
result — which is an argument for making seed repeats cheap enough that they
are routine rather than a special exercise.
