# Genome-wide Kinome CRISPR Knockout Screen — Design Package

A complete, execution-ready design for a pooled **CRISPR-Cas9 knockout screen
targeting the entire human kinome** (protein + lipid + selected metabolic
kinases). It covers guide selection, off-target prediction, library
construction, screen scale, sequencing depth, and analysis.

## Contents

| File | What it covers |
|------|----------------|
| [`docs/01_target_kinome.md`](docs/01_target_kinome.md) | Defining the target gene set (~518 protein kinases + extensions) |
| [`docs/02_guide_design.md`](docs/02_guide_design.md) | Optimal guide selection: on-target scoring, filters, controls |
| [`docs/03_offtarget_prediction.md`](docs/03_offtarget_prediction.md) | Off-target enumeration and scoring (CFD/MIT, GuideScan2, Cas-OFFinder) |
| [`docs/04_library_construction.md`](docs/04_library_construction.md) | Oligo pool → cloning → plasmid QC |
| [`docs/05_screen_and_sequencing.md`](docs/05_screen_and_sequencing.md) | Coverage, MOI, timepoints, gDNA input, read depth |
| [`docs/06_analysis.md`](docs/06_analysis.md) | MAGeCK / BAGEL2 / drugZ analysis and QC |
| [`docs/07_two_screen_designs.md`](docs/07_two_screen_designs.md) | Essentiality + chemogenomic/synthetic-lethal in one experiment |
| [`docs/08_drug_arm_sop.md`](docs/08_drug_arm_sop.md) | Drug-arm SOP: dose pre-experiment, dosing, drugZ matrix, expected hits |
| [`docs/09_analysis_commands.md`](docs/09_analysis_commands.md) | Runnable MAGeCK / BAGEL2 / drugZ command cheat-sheet |
| [`docs/10_case_study_A375_BRAFi.md`](docs/10_case_study_A375_BRAFi.md) | Fully worked case: A375 (BRAF^V600E) + vemurafenib, ChEMBL/clinical data |
| [`docs/11_drug_arm_candidate_menu.md`](docs/11_drug_arm_candidate_menu.md) | ChEMBL-verified menu of 7 approved inhibitors: potency, target, matched background |
| [`scripts/coverage_calculator.py`](scripts/coverage_calculator.py) | Runnable scale/depth calculator (cells, gDNA, reads) |
| [`data/kinome_targets.md`](data/kinome_targets.md) | Target-set definition + authoritative source fetch |
| [`templates/sample_sheet.csv`](templates/sample_sheet.csv) | 10-sample sheet layout for the dual-arm design |
| [`docs/zh/`](docs/zh/) | 中文对照版全套流程文档 (Chinese parallel of all docs) |

## Design at a glance

- **Targets:** ~518 protein kinases (Manning classification); extendable to
  ~635 (adds lipid + metabolic) or ~750 (broad "kinome").
- **Guides:** 6 sgRNA/gene (option for 8 for higher confidence) + ~250
  non-targeting + ~50 intergenic-cutting + positive/negative essential
  controls → **~3,600–5,000 guide library**.
- **On-target scoring:** Rule Set 3 (Sanson 2018) / DeepSpCas9, early
  constitutive exons, all-transcript coverage.
- **Off-target:** CFD specificity ≥ 0.2 filter, no perfect/1-mismatch
  exonic off-targets, verified with GuideScan2 + Cas-OFFinder.
- **Scale:** ≥1000× cells/guide maintained throughout, MOI ~0.3, ≥3
  biological replicates.
- **Sequencing:** ~1000 reads/guide (≈5M reads/sample for a 5k library),
  gDNA input sized to preserve ≥1000× genome coverage.

Run the calculator for concrete numbers tuned to your library and cell line:

```bash
python3 scripts/coverage_calculator.py --genes 518 --guides-per-gene 6 --coverage 1000
```
