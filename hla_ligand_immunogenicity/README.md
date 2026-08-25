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
| 1 | **Panel designed against measured coverage** (M2) | The DR subset of the IEDB class-II reference set is widely used as a "representative" 15-molecule panel. Measured against the IEDB allele-frequency tables it reaches only **85.3 %** weighted US/EU DRB1 phenotypic coverage — it does not meet a 95–98 % requirement. A greedy build to a stated target does. |
| 2 | **Two orthogonal prediction heads** (M3) | Eluted-ligand (EL) scoring alone over-calls: it rewards motif-like peptides that have poor measured affinity. Requiring EL **and** binding-affinity agreement removes that class of hit and the drop is reported, not hidden. |
| 3 | **Self / pre-existing-tolerance filter** (M4) | Most DR hits in an antibody-derived ligand sit in *framework* whose 9-mer cores are near-identical to human immunoglobulin V germline. Counting them inflates every VHH-, scFv- and Fab-derived ligand identically and destroys the ranking. Cores 9/9 or 8/9 identical to a human proteome 9-mer are down-weighted, and the cut is validated against a shuffled-sequence null in the same run. |
| 4 | **Population weighting + benchmark calibration with controls** (M5, M6) | "13 strong binders" is uninterpretable. Weighting each hit by the fraction of the US/EU population that carries the presenting molecule, then expressing the result as a fold-change over the Protein A Z-domain — an affinity ligand with decades of controlled clinical leachate exposure — makes it a comparison. Positive and negative controls run in the same batch make the batch *reportable* or not. |
| 5 | **B-cell/ADA layer and exposure context** (M7, M8) | The measured endpoint is an anti-drug **antibody** assay, and impurity risk scales with µg delivered per dose. A T-cell-only, dose-free score cannot reach a risk call. |

`M9` is optional: an anchor-position deimmunisation scan of the dominant
epitope, for the case where the ligand can be re-engineered.

---

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

Figures: `make_figures.py`. Report: `make_report.py` (HTML) and
`make_deck.py` (PPTX).

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
| `TT_p2_region`, `TT_p30_region` (P04958) | positive controls | tetanus toxin universal T-helper epitopes — the panel must find them |

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
python scripts/make_figures.py
python scripts/make_report.py
python scripts/make_deck.py
```

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
