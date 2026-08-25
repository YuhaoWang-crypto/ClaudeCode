# HLA-DR immunogenicity assessment of process-related affinity ligands

An in-silico risk-assessment pipeline for **leached affinity-chromatography
ligands** — proprietary VHH/protein-scaffold ligands that enter drug product as
*other* process-related impurities rather than host-cell protein.

The demo runs end to end on a **public stand-in** for a proprietary ligand: the
AAVX affinity ligand (a camelid VHH) from PDB **9DC3**, the ligand of the
CaptureSelect-type AAV affinity resins. Drop your own de-identified FASTA into
`data/sequences.fasta` and every module runs unchanged.

---

## Why this differs from a plain "run NetMHCIIpan over a DR panel"

A single-tool DR scan answers *"how many peptides in this protein could bind
HLA-DR?"*. That number is not a risk, because it is large for every non-human
protein and it is not comparable between ligands. Five things were added so the
output is a decision instead of a count.

| # | Addition | The failure mode it fixes |
|---|---|---|
| 1 | **Panel designed against measured coverage** (M2) | The DR subset of the IEDB class-II reference set is widely used as a "representative" 15-molecule panel. Measured against the IEDB allele-frequency tables it reaches only **84.2 %** weighted US/EU DRB1 phenotypic coverage — it does not meet a 95–98 % requirement. A greedy build to a stated target does. |
| 2 | **Two orthogonal prediction heads** (M3) | Both eluted-ligand and binding-affinity scores are computed for every 15-mer. Whether affinity should *gate* a call was an open question that M11 settled — see below; the gate is off, and both scores are still reported. |
| 3 | **Self / pre-existing-tolerance filter** (M4) | Most DR hits in an antibody-derived ligand sit in *framework* whose 9-mer cores are near-identical to human immunoglobulin V germline. Counting them inflates every VHH-, scFv- and Fab-derived ligand identically and destroys the ranking. Cores 9/9 or 8/9 identical to a human proteome 9-mer are down-weighted, and the cut is validated against a shuffled-sequence null in the same run. |
| 4 | **Population weighting + benchmark calibration with controls** (M5, M6) | "13 strong binders" is uninterpretable. Weighting each hit by the fraction of the US/EU population that carries the presenting molecule, then expressing the result as a fold-change over the Protein A Z-domain — an affinity ligand with decades of controlled clinical leachate exposure — makes it a comparison. Positive and negative controls run in the same batch make the batch *reportable* or not. |
| 5 | **B-cell/ADA layer and exposure context** (M7, M8) | The measured endpoint is an anti-drug **antibody** assay, and impurity risk scales with µg delivered per dose. A T-cell-only, dose-free score cannot reach a risk call. |
| 6 | **The decision rule is measured, not argued** (M10, M11) | Every choice above — which threshold, whether the second head helps, what a flag is worth — was made by reasoning. M10 pulls ~9,600 labelled HLA-DR-restricted human CD4 T-cell outcomes from IEDB; M11 scores each rule against them with a cluster-level bootstrap, and reports sensitivity, specificity and the PPV a flag actually carries at realistic scan prevalence. |

`M9` is optional: an anchor-position deimmunisation scan of the dominant
epitope, for the case where the ligand can be re-engineered.

---

## What the demo run found

The full pipeline was run end to end on the public AAVX VHH. Headline results
(`report.html`, `report.pptx`, and the tables in `results/`):

| | |
|---|---|
| Panel | 21 DRB1 + 4 DRB3/4/5 = **25 DR molecules**, **97.3 %** weighted US/EU DRB1 coverage. The legacy 15-molecule DR reference set reaches **84.2 %** on the same measure. |
| Prediction | **266** eluted-ligand strong-binder calls across the batch. The affinity gate that used to filter them is off: measured against the benchmark it removed 23 true positives and 6 false positives. |
| Tolerance filter | 17.6 % of predicted cores are exact human 9-mers, 4.7 % are one substitution away; the shuffled-sequence null hits 1.1 % at the same cut (**20× enrichment**). |
| Test article | pIRS **0.47** = **1.05× the Protein A Z-domain**; 4 non-self epitopes; **39.6 %** of the weighted US/EU population carry a DR molecule predicted to present at least one of them. |
| Dominant epitope | residues 47–56, core `FVAVQDITA`, peak 15-mer `KEREFVAVQDITASN`, presented by 7 of 25 DR molecules. |
| Accuracy | Best decision rule reaches **ROC AUC 0.66** against measured T-cell outcomes. At the operating point used, sensitivity **0.16**, specificity **0.96**, and a flagged peptide is real about **16 %** of the time at a 5 % assumed scan prevalence. |
| Batch controls | all four system-suitability checks pass. |

### Two claims from the previous revision did not survive measurement

- **The affinity gate was not buying specificity.** Requiring the binding-affinity
  head to agree before calling a strong binder was introduced on a plausible
  argument. Measured against 5,795 labelled outcomes it removed **23 true
  positives and 6 false positives** and left Matthews correlation slightly
  *worse* (0.188 → 0.184). Removing it moved the test article from 2.9× the
  Protein A benchmark to **1.05×** — the gate had been suppressing the
  benchmark's epitopes harder than the test article's, so an unvalidated rule
  was setting the headline.
- **The ligand's core is not ~3× more promiscuous than the universal epitopes.**
  That number came from a bug: control epitopes were matched to scanned 15-mers
  by string equality, so a 13-residue epitope such as HA306-318 matched nothing
  and scored zero breadth. With overlap matching the universal epitopes reach
  2–6 of 25 DR molecules at %Rank < 1 and the ligand's core reaches 6/25 —
  **comparable to**, not above, the most promiscuous epitope with human T-cell
  evidence.

Two findings from the controls are worth reading on their own:

- **The strong-binder tier is a high-specificity, low-sensitivity criterion.**
  HA306-318 has positive human T-cell assays on 25 distinct DR molecules and
  clears EL %Rank < 1 on 2 of the 25 tested here; the tier recovers about a
  quarter of the molecules a universal epitope is actually presented by. That
  matches the benchmark's measured sensitivity of 0.16. Peptides below the tier
  are *unflagged*, not *cleared*.
- **Human germline VH3-23 scores a peak promiscuity of 6/25 and an unfiltered
  pIRS of 0.37, in the same range as the test article.**
  Raw binder counts cannot distinguish a camelid VHH from the human framework
  it resembles. That is the entire case for the tolerance filter, and it is why
  the ranking is unreadable without it.

## Modules

| Module | Script | Output |
|---|---|---|
| M0 | `m0_fetch_sequences.py` | `data/sequences.fasta`, `sequences_metadata.tsv` — every sequence pulled live from RCSB/UniProt with its accession |
| M1 | `m1_sequence_qc.py` | `m1_sequence_qc.tsv` — composition QC, MW/pI, % identity to human germline IGHV3-23, VHH hallmark tetrad |
| M2 | `m2_panel_design.py` | `m2_panel.json`, `m2_panel_alleles.txt` — greedy DR panel + coverage curve |
| M3 | `m3_binding_prediction.py` | `m3_binding_long.tsv` — NetMHCIIpan EL + BA, every 15-mer × every DR molecule |
| M4 | `m4_tolerance_filter.py` | `m4_core_tolerance.tsv` — per-core self/foreign call against Swiss-Prot human |
| M5 | `m5_risk_scoring.py` | `m5_epitopes.tsv`, `m5_clusters.tsv`, `m5_ligand_summary.tsv` — pIRS and population-at-risk |
| M6 | `m6_benchmark_calibration.py` | `m6_calibrated_ranking.tsv`, `m6_system_suitability.json` |
| M7 | `m7_bcell_layer.py` | `m7_bcell_regions.tsv`, `m7_tb_coincidence.tsv` |
| M8 | `m8_exposure_context.py` | `m8_exposure_grid.tsv` — µg ligand/dose bands |
| M9 | `m9_deimmunization_scan.py` | `m9_deimmunization_scan.tsv` |
| M10 | `m10_benchmark_fetch.py` | `m10_benchmark.tsv` — every HLA-DR-restricted human CD4 T-cell assay outcome IEDB holds for the panel, labelled per (peptide, allele) |
| M11 | `m11_threshold_calibration.py` | `m11_calibration.json` — ROC/PR for each decision rule, cluster-bootstrap comparison, calibrated operating point |

Figures: `make_figures.py`. Report: `make_report.py` → `report.html`;
`make_deck.py` (+ `make_deck.js`) → `report.pptx`. `check_deck.py` lints the
deck's geometry for off-slide shapes, text overflow and overlapping boxes —
LibreOffice cannot load `.pptx` in this environment, so that lint stands in for
a visual render. It does not inspect table cells, so wide tables still deserve
a look on a real machine.

## The batch

The ligand is never run alone. Every batch carries controls, and the batch is
only reportable if they behave:

| Sequence | Role | Why it is in the batch |
|---|---|---|
| `AAVX_VHH` (PDB 9DC3) | test article | the ligand under assessment |
| `ProteinA_Z` (PDB 1Q2N) | benchmark | Z domain of rProtein A resins — the calibration anchor |
| `ProteinA_B_native` (P38507) | benchmark | native, non-engineered SpA domain |
| `ProteinL_B1` (Q51918) | benchmark | protein L Ig-binding domain |
| `Caplacizumab_VHH` (PDB 7EOW) | clinical anchor | approved humanised VHH — the only nanobody in the set with published clinical ADA rates |
| `VHH_7D12` (PDB 4KRL) | class comparator | non-humanised camelid VHH background |
| `HumanVH3_23_germline` (P01764) | negative control | human germline VH — tolerised floor |
| `HSA_D1` (P02768) | negative control | human self protein |
| `ProteinG_B1` (P06654) | benchmark | protein G IgG-binding domain — third clinically used affinity ligand |
| `HA_306_318_region` (P03437) | positive control | influenza HA306-318 — **positive human DR-restricted T-cell assays on 25 distinct DR molecules** in IEDB, the best-evidenced promiscuous epitope available |
| `TT_p2_region`, `TT_p30_region` (P04958) | positive controls | tetanus toxin universal T-helper epitopes (6 and 1 DR molecules respectively — weaker than their reputation) |
| `MBP_85_99_region` (P02686) | boundary control | a **self** peptide that IS a validated epitope on 10 DR molecules — measures whether the tolerance filter suppresses real risk |
| `CLIP_87_101_region` (P04233) | boundary control | a universal DR **ligand** with no positive human T-cell record — measures whether binding strength alone is read as risk |

Every control's role was checked against IEDB before being assigned. PADRE
(`AKFVAAWTLKAAA`) was a candidate and was dropped: IEDB holds no positive human
DR-restricted T-cell record for it.

## Running it

```bash
pip install pyyaml biopython pandas numpy matplotlib python-pptx
python scripts/m0_fetch_sequences.py
python scripts/m1_sequence_qc.py
python scripts/m2_panel_design.py
python scripts/m3_binding_prediction.py     # ~40 min against the IEDB cluster
python scripts/m4_tolerance_filter.py
python scripts/m5_risk_scoring.py
python scripts/m6_benchmark_calibration.py
python scripts/m7_bcell_layer.py
python scripts/m8_exposure_context.py
python scripts/m9_deimmunization_scan.py    # optional
python scripts/m10_benchmark_fetch.py       # ~10 min against the IEDB query API
python scripts/m11_threshold_calibration.py # ~1-2 h; resumable, caches as it goes
python scripts/m12_tolerance_weight.py
python scripts/make_figures.py
python scripts/make_report.py
python scripts/make_deck.py
python scripts/check_deck.py
```

`m11_threshold_calibration.py` is the long pole: it scores ~5,800 benchmark
peptides on both prediction heads, and the endpoint's wall-clock cost tracks the
number of requests rather than their size. Peptides are therefore submitted as
concatenated pseudo-proteins — verified against standalone scoring before the
run trusts it — and results are appended to `m11_scores_partial.tsv` as they
land, so an interrupted run resumes instead of restarting.

`m3_binding_prediction.py` is resumable: it reads any existing
`results/m3_binding_long.tsv` and only fetches the (sequence, allele, head)
combinations missing from it, so widening the panel costs only the new allele's
calls rather than a full re-run.

`data/human_sprot.fasta` (UniProt Swiss-Prot *Homo sapiens*, ~20 400 entries)
is fetched once:

```bash
curl -L 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true' \
  | gunzip > data/human_sprot.fasta
```

Everything else is configured in `config/config.yaml` — panel target, %Rank
thresholds, tolerance weights, population weights, exposure grid.

## Methods and their limits

- **Predictors**: NetMHCIIpan (EL + BA heads) served by the IEDB REST API;
  BepiPred-2.0 for linear B-cell propensity. No licensed local binary is
  required.
- **Coverage model**: single-locus Hardy–Weinberg phenotypic frequency over the
  IEDB population-coverage allele-frequency tables, verified against the IEDB
  CLI. DRB3/4/5 carry no frequencies in those tables, so they add
  presentation breadth without changing the DRB1 coverage arithmetic and are
  reported separately, never double-counted.
- **pIRS is a relative scale.** It is interpretable only against the benchmarks
  and controls run in the same batch. It is not a predicted ADA incidence, and
  no in-silico method available today predicts one.
- **The tolerance filter is a screen, not JanusMatrix.** It does not require the
  human counterpart peptide to bind the same allele, so it errs toward calling
  more peptides tolerised. Every flagged core is written out with the human
  protein it matched, so each call is checkable. A 5-of-9 TCR-face variant was
  tested and rejected — it matches the human proteome by chance several times
  per query. The 8/9 whole-core cut used instead is validated against a
  shuffled-sequence null in `m4_filter_validation.json`; if the null hit rate
  is not far below the real one, the run reports the filter as uninformative.
- **Linear B-cell prediction is the weakest model here.** Most real ADA epitopes
  are conformational. M7 output prioritises regions for wet-lab work; it is
  never a standalone claim.
- **Risk bands are an internal triage convention**, not a regulatory
  classification. No agency publishes a numeric leachate immunogenicity limit;
  the ICH Q6B / EMA expectation is a justified, consistently achieved control
  level.

## What this cannot replace

In-silico DR screening ranks and localises risk. It does not measure it. The
confirmatory work this output is designed to scope, in increasing cost:

1. **HLA-DR competitive binding** on the 3–6 flagged peptides against the
   panel's dominant molecules — confirms the predicted binding directly.
2. **MAPPs** (MHC-associated peptide proteomics) on monocyte-derived dendritic
   cells from HLA-typed donors — shows what is actually processed and
   presented, which prediction cannot.
3. **Ex-vivo PBMC / T-cell proliferation** across ~50 HLA-typed donors matched
   to the panel — the closest available surrogate for clinical ADA risk.

All three are available RUO from standard CROs and all three are scoped by the
peptide list this pipeline produces.
