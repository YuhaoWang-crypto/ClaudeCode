---
name: hla-dr-immunogenicity
description: >-
  Assess the T-cell immunogenicity risk of a protein impurity or biologic
  against a population-representative HLA-DR panel, and produce a defensible
  risk call rather than a binder count. Use when asked to evaluate the
  immunogenicity of an affinity-chromatography ligand, a leachate or other
  process-related impurity, an HCP, a VHH/nanobody/scFv, or any protein
  sequence; to design an HLA-DR allele panel to a stated population-coverage
  target; to run NetMHCIIpan class-II predictions over a panel via the IEDB
  API; to separate self/tolerised hits from genuinely foreign epitopes; to
  calibrate a decision rule against measured human CD4 T-cell outcomes; or to
  produce the HTML/PPTX risk report. Enforces the discipline that every
  decision rule is measured against labelled data before it is believed.
---

# HLA-DR immunogenicity risk assessment for protein impurities and ligands

A 14-module pipeline that turns a protein sequence into a **ranked, calibrated,
population-weighted immunogenicity risk assessment** — with controls in every
batch and every decision rule scored against labelled human T-cell outcomes.

Built for **process-related impurities** (affinity-ligand leachate, HCP), where
the sequence is often proprietary, the exposure is µg/dose, and the question is
"is this ligand riskier than the ones already in clinical use?" — not "does this
protein contain DR binders?" It runs unchanged on any protein sequence.

## The core problem this solves

A plain NetMHCIIpan scan over a DR panel answers *"how many 15-mers could bind
HLA-DR?"* **That number is not a risk.** It is large for every non-human
protein, it is not comparable between molecules, and for antibody-derived
ligands it is dominated by framework regions that are near-identical to human
germline. Human germline VH3-23 scores in the same range as a camelid VHH on
raw binder count — so the raw count cannot even separate self from foreign.

Six additions turn the count into a decision:

| # | Addition | The failure mode it fixes | Module |
|---|---|---|---|
| 1 | **Panel designed to measured coverage** | The widely used 15-molecule DR reference set reaches only **84.2 %** weighted US/EU DRB1 phenotypic coverage. A greedy build to a stated target reaches 97.3 % with 21 DRB1 + 4 DRB3/4/5. | M2 |
| 2 | **Two orthogonal prediction heads** | EL and BA are both computed. Whether BA should *gate* a call is an empirical question, not an argument — see M11. | M3 |
| 3 | **Self / pre-existing-tolerance filter** | Framework hits inflate every VHH-, scFv- and Fab-derived ligand identically and destroy the ranking. Cores 9/9 or 8/9 identical to a human proteome 9-mer are down-weighted, validated against a shuffled-sequence null in the same run. | M4 |
| 4 | **Population weighting + benchmark anchoring** | "13 strong binders" is uninterpretable. Weight each hit by the fraction of the population carrying the presenting molecule, then express it as a fold-change over a ligand with decades of controlled clinical leachate exposure. | M5, M6 |
| 5 | **B-cell layer and exposure context** | The measured endpoint is an anti-drug **antibody** assay and risk scales with µg/dose. A T-cell-only, dose-free score cannot reach a risk call. | M7, M8 |
| 6 | **Every decision rule measured** | Threshold, gate, tolerance weight, breadth criterion — each was scored against ~9,600 labelled HLA-DR-restricted human CD4 T-cell outcomes from IEDB with a cluster-level bootstrap. | M10–M13 |

## The non-negotiable discipline: rules are measured, not argued

**Any proposed improvement to the decision rule must be scored against the
labelled benchmark before it is adopted.** In this pipeline four plausible
improvements were proposed and *all four were rejected by measurement*:

- **The EL×BA consensus gate did not buy specificity.** It removed 23 true
  positives for 6 false positives and left MCC slightly worse (0.188 → 0.184).
  Turning it off moved the headline from 2.92× the Protein A benchmark to
  **1.05×** — an unvalidated rule had been setting the answer.
- **"3× more promiscuous than universal epitopes" was a string-matching bug.**
  Control epitopes were matched to scanned 15-mers by equality, so a 13-residue
  epitope matched nothing and scored zero breadth. With overlap matching the
  universal epitopes reach 2–6/25 and the test article 6/25 — comparable, not
  above.
- **The tolerance down-weight could not be bounded at all.** IEDB self peptides
  enter mostly through autoimmunity studies, so the labels carry an
  ascertainment bias that the data cannot correct. A sensitivity sweep replaced
  the pretence of an optimum.
- **Panel-wide breadth does not beat the best single-allele rank.** ΔAUC
  +0.0005, 95 % CI [−0.0156, +0.0169]. "How many molecules IEDB tested it on"
  reaches AUC 0.575 by itself — higher than any sequence-derived predictor.

If a run reproduces a claim you cannot trace to a number in `results/`, treat it
as unmeasured. Report negative results as findings. When new data overturns an
earlier claim in the same project, say so explicitly and change the number.

## Run it

```bash
pip install pyyaml biopython pandas numpy matplotlib python-pptx
npm install pptxgenjs                       # deck only

python scripts/m0_fetch_sequences.py        # test article + benchmarks + controls
python scripts/m1_sequence_qc.py
python scripts/m2_panel_design.py           # greedy panel to the coverage target
python scripts/m3_binding_prediction.py     # ~40 min; resumable
python scripts/m4_tolerance_filter.py
python scripts/m5_risk_scoring.py
python scripts/m6_benchmark_calibration.py  # system suitability — batch passes or not
python scripts/m7_bcell_layer.py
python scripts/m8_exposure_context.py
python scripts/m9_deimmunization_scan.py    # optional; only if the ligand can be re-engineered
python scripts/m10_benchmark_fetch.py       # ~10 min; labelled T-cell outcomes
python scripts/m11_threshold_calibration.py # ~1-2 h; resumable
python scripts/m12_tolerance_weight.py
python scripts/m13_promiscuity_vs_bestrank.py  # ~2 h; resumable
python scripts/make_figures.py && python scripts/make_report.py
python scripts/make_deck.py && python scripts/check_deck.py
```

M0–M9 produce the assessment. **M10–M13 are the calibration layer and only need
re-running when a decision rule changes** — their outputs are properties of the
panel and the thresholds, not of the test article. For a new ligand on an
unchanged panel, run M0–M9 and reuse the existing `m1*_*.json`.

`assets/scripts/` holds the whole pipeline and `assets/config/config.yaml` is
the single source of truth for every threshold, weight and population share.
Copy both into the project, drop your de-identified FASTA into
`data/sequences.fasta`, and every module runs unchanged.

The human proteome is fetched once:

```bash
curl -L 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true' \
  | gunzip > data/human_sprot.fasta
```

## The batch: never run the ligand alone

A pIRS number is meaningless in isolation — it is readable only against what
else ran in the same batch. **Every run carries benchmarks and controls, and
the batch is reportable only if they behave.** M6 writes the system-suitability
verdict; if a control fails, the run is not reportable, no matter how clean the
test article looks.

| Role | Members | What it establishes |
|---|---|---|
| Benchmarks | Protein A Z domain, native SpA, protein L, protein G | The clinically qualified leachate anchor — the denominator of the fold-change |
| Clinical anchor | Caplacizumab VHH | The only nanobody in the set with published clinical ADA rates |
| Negative controls | human germline VH3-23, HSA domain 1 | The tolerised floor. VH3-23 scoring like a VHH on raw counts *is* the case for the tolerance filter |
| Positive controls | HA306-318, TT p2, TT p30 | The assay ceiling. HA306-318 has positive human DR T-cell assays on 25 distinct molecules |
| Boundary controls | MBP85-99, CLIP87-101 | Where the method is known to be wrong: a *self* peptide that IS a validated epitope, and a universal DR *ligand* with no T-cell record |

Check every control's role against IEDB before assigning it. PADRE was a
candidate here and was dropped — IEDB holds no positive human DR-restricted
T-cell record for it.

## Reference

Read the relevant file before touching that part of the pipeline:

- `reference/panel-design.md` — the coverage model, the renormalisation that
  most implementations get wrong, DRB3/4/5 handling, hitting a stated target.
- `reference/iedb-api.md` — both IEDB APIs, every quirk that cost time
  (308s, `*` encoding, `offset` without `order`), and the concatenated
  pseudo-protein trick that makes a 5,800-peptide benchmark run feasible.
- `reference/calibration.md` — how to build the labelled benchmark, cluster-
  level bootstrapping, the confounds that must be stratified out, and the four
  measured negative results in full.
- `reference/interpreting.md` — what pIRS is and is not, how to write the risk
  call, the operating point and what a flag is actually worth.

## Limits — state these in every report

In-silico DR screening **ranks and localises** risk. It does not measure it.

- **Prediction, not presentation.** No model of uptake, endosomal proteolysis,
  HLA-DM editing or complex stability.
- **DR only.** DP and DQ contribute to CD4 responses and DQ is implicated in
  several biologic ADA responses. Excluded by the panel specification.
- **pIRS is relative.** Not a predicted ADA incidence. No in-silico method
  available today predicts one.
- **No aggregation or adjuvant effect.** Aggregated impurity is substantially
  more immunogenic than monomer; sequence methods cannot see it.
- **Absolute accuracy numbers are an upper bound.** NetMHCIIpan is trained on
  IEDB data and partial training-set overlap with any IEDB-derived benchmark is
  certain. Rule-vs-rule comparisons are far more robust: leakage inflates both
  arms and largely cancels in the difference.
- **Peptides below the tier are unflagged, not cleared.** At the operating point
  used here, sensitivity is 0.16.

The output is a scoped wet-lab plan, and that is the practical point of running
it: HLA-DR competitive binding on the flagged peptides (days), then MAPPs on
monocyte-derived dendritic cells from HLA-typed donors, then ex-vivo CD4
proliferation across ~50 donors matched to the panel. All three are available
RUO and all three are scoped by the peptide list this pipeline produces.
