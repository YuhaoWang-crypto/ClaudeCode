# Interpreting the output and writing the risk call

## From 15-mers to a number

Raw per-15-mer binder counts are not the unit of anything. M5 does five things:

1. **Collapse** overlapping 15-mers into *epitopes*, keyed on (sequence,
   binding core), keeping the best rank and its position.
2. For each epitope, compute the **fraction of the weighted US/EU population
   carrying at least one DR molecule predicted to present it** — Hardy–Weinberg
   phenotypic frequency over the presenting DRB1 set.
3. Apply the **M4 tolerance weight**, so framework cores shared with the human
   proteome stop dominating.
4. Roll epitopes into **positional clusters** (cores within 8 residues join).
5. Emit one headline per sequence.

```
pIRS = 100/L * SUM over foreign epitopes of
               (weighted US/EU presenting fraction) x (tolerance weight)
```

— population-weighted presentable **foreign** epitope content per 100 residues.
`pIRS_no_tolerance_filter` is reported alongside it; the gap between the two is
how much of the raw signal was framework.

Also reported, and much easier to explain to a non-specialist:

```
pop_at_risk = fraction of the weighted US/EU population carrying at least one
              DR molecule predicted to present at least one foreign epitope
```

## pIRS is a relative scale. Never quote it alone

It is interpretable **only against the benchmarks and controls run in the same
batch**. Quote it as a fold-change over the anchor:

> pIRS 0.47 = **1.05× the Protein A Z domain**, 4 non-self epitopes, 39.6 % of
> the weighted US/EU population carry a DR molecule predicted to present at
> least one of them.

The Z domain of rProtein A is the anchor because it is an affinity ligand with
decades of controlled clinical leachate exposure — a real-world tolerated
reference, not a theoretical zero. HSA and human germline VH3-23 are the
tolerised floor; tetanus p30 is the ceiling.

**pIRS is not a predicted ADA incidence.** No in-silico method available today
predicts one. Do not convert it to a percentage of patients, and do not compare
pIRS values across batches run on different panels or thresholds.

## Check the batch before reading the test article

M6 writes the system-suitability verdict. **If a control fails, the run is not
reportable**, however clean the test article looks. In particular the positive
controls must reach a minimum DR breadth at the weak-binder tier
(`positive_control_min_breadth_wb`), or the batch is not sensitive enough for a
negative result to mean anything.

The boundary controls are not pass/fail — they are the measured statement of
where the method is wrong, and they belong in the report:

- **MBP85-99** is a *self* peptide that IS a validated epitope on 10 DR
  molecules. The tolerance filter's core assumption (self ⇒ tolerated) is wrong
  for it, and the run states by how much.
- **CLIP87-101** is a universal DR *ligand* with no positive human T-cell
  record — it tests that binding strength alone is not being read as risk.

## What a flag is worth, and what an absence of a flag is worth

At this pipeline's operating point: sensitivity **0.16**, specificity **0.96**,
and a flagged peptide is a real epitope about **16 %** of the time at 5 %
assumed scan prevalence.

**Peptides below the tier are unflagged, not cleared.** The strong-binder tier
is a high-specificity, low-sensitivity criterion: HA306-318 has positive human
T-cell assays on 25 distinct DR molecules and clears EL %Rank < 1 on only 2 of
the 25 tested here. The tier recovers roughly a quarter of the molecules a
universal epitope is actually presented by. Write that into the report — a
reader who takes "no strong binders" as "no risk" has been misled by the
omission.

## The tolerance filter's boundaries

It is a **screen, not JanusMatrix**. It does not require the human counterpart
peptide to bind the same allele, so it errs toward calling more peptides
tolerised — a conservative-in-the-wrong-direction error worth stating. Every
flagged core is written out with the human protein it matched, so each call is
individually checkable.

## Exposure context

Impurity risk scales with µg delivered per dose. The M8 grid crosses leachate
ppm × dose mg into µg ligand/dose bands. Sub-microgram/dose protein impurity
exposure is the regime in which rProtein A leachate has decades of clinical
use — that is the comparison the grid is for. **Edit `exposure.ligand_mw_kda`,
`leachate_ppm` and `dose_mg` to the actual product**; the defaults are
placeholders and a risk-in-context matrix built on someone else's dose is
worthless.

Risk bands are an **internal triage convention**, not a regulatory
classification. No agency publishes a numeric leachate immunogenicity limit;
the ICH Q6B / EMA expectation is a justified, consistently achieved control
level.

## B-cell layer

Linear B-cell prediction is the weakest model in the pipeline — most real ADA
epitopes are conformational. M7 output prioritises regions for wet-lab work and
flags T/B coincidence; it is never a standalone claim.

## Deliverables

- `report.html` — data-driven, self-contained, every number traced to a file in
  `results/`.
- `report.pptx` — 15 slides via `make_deck.py` → `make_deck.js` (pptxgenjs).
- `check_deck.py` — **run it every time.** It lints slide geometry for
  off-slide shapes, text overflow and overlapping boxes. LibreOffice cannot
  load `.pptx` in a headless container, so this lint stands in for a visual
  render. It does not inspect table cells, so wide tables still deserve a look
  on a real machine.

Escape everything that reaches HTML. Operating-point rule names contain `<`
(e.g. `%Rank < 1`), and emitting them unescaped parses them as tags and
truncates the table silently.

## Writing the risk call

State, in this order: the fold-change over the anchor and what the anchor is;
the number of foreign epitopes and the dominant one with its core, position and
DR breadth; the population fraction at risk; the exposure band for the actual
product; the batch-control verdict; and then the limits — DR-only, prediction
not presentation, relative not absolute, sensitivity 0.16.

Then the wet-lab plan the output scopes: DR competitive binding on the flagged
peptides (days, lowest cost), MAPPs on monocyte-derived dendritic cells from
HLA-typed donors, ex-vivo CD4 proliferation across ~50 donors matched to the
panel. All three are available RUO and all three are scoped by the peptide list
this pipeline produces — which is the practical point of running it.
