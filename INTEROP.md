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

1. Replace one adapter with a real model and keep every test green — the
   swap cost is the thing being measured.
2. Wrap one model in conformal prediction and check the interval's empirical
   coverage.
3. Add a second model for the same hop and let the router pick; report how
   much the endpoint moves. Disagreement between two routes is a free,
   experiment-free uncertainty estimate.
4. Then, and only then, a leakage-aware chain benchmark on a hop with real
   measured truth on both sides.
