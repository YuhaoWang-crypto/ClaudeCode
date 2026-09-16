# ADMET & Safety-Panel Coverage Analysis — 30-drug FDA panel

Computational ADMET / drug-likeness profiling of 30 FDA-approved small-molecule
drugs, validated against documented pharmacology, and mapped onto the
secondary-pharmacology safety panels used in drug safety assessment
(Bowes-44 / Safety-77 / SafetyScreen87).

The primary deliverable is the Chinese-language report
`results/admet_analysis_report_zh_v2.pdf` (18 pages, 10 sections, 6 figures,
10 tables).

## Headline results

| Question | Answer |
|---|---|
| Do DMPK / hERG / BBB endpoints recover known pharmacology? | Yes — 19/21 positive controls above 0.5, ROC-AUC 1.00 on hERG, CYP2D6, P-gp, BBB |
| Do Tox21 nuclear-receptor endpoints? | No — 1/5 positive controls. Tamoxifen/ER 0.18, bicalutamide/AR 0.06, rosiglitazone/PPARγ 0.30, anastrozole/aromatase 0.27 |
| How much of the Bowes-44 safety panel can computation cover? | 2/44 targets map nominally (hERG, AR); only 1/44 (hERG, 2.3%) survives positive-control validation |
| Are CYP and P-gp endpoints part of a safety panel? | No — they are DMPK / drug–drug-interaction assays, a different discipline |

## Pipeline

```
scripts/drug_panel.py       30 FDA drugs, SMILES validated against reference formulae
scripts/admet_pipeline.py   standardisation → descriptors → rules → alerts → ADMET-AI
scripts/demo_validation.py  positive/negative control test against documented pharmacology
scripts/make_plots.py       6 figures (PNG + SVG)
scripts/report_data.py      all report tables, computed from the CSVs
scripts/mapping_core.py     data-driven half of the panel coverage mapping
scripts/safety_panel.py     literature half + chapter 6 of the report
scripts/pdfkit.py           ReportLab scaffolding (CJK font, tables, callouts)
scripts/build_report_zh.py  builds the PDF
scripts/export_panel_csv.py exports tables 8 and 9 as CSV
scripts/make_markdown.py    builds results/analysis_report.md
```

Run order:

```bash
python admet_pipeline.py && python demo_validation.py && python make_plots.py \
  && python export_panel_csv.py && python make_markdown.py && python build_report_zh.py
```

## Environment

```bash
uv venv --python 3.11 .venv && . .venv/bin/activate
uv pip install rdkit pandas numpy scikit-learn matplotlib seaborn \
               reportlab joblib pypdf pymupdf admet-ai
```

ADMET-AI downloads its model weights on first use. Everything runs on CPU;
the full pipeline takes a few minutes for 30 molecules. The PDF needs a CJK
font — `/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc` (WenQuanYi Zen Hei).

## Evidence grading in the report

Both source papers (Bowes 2012, Brennan 2024, *Nature Reviews Drug Discovery*)
are paywalled. The per-target Bowes-44 composition used here comes from the
Eurofins SafetyScreen44 product flyer, which the vendor explicitly attributes
to Bowes et al. 2012. **The full Safety-77 target list could not be obtained
and is not reconstructed anywhere in the report** — only its verified kinase
count (20/77, 26%) and qualitative family-enrichment statements are used.
Section 9 of the report states these limits explicitly.

## Scope

A computational screening aid for research triage. Not a regulatory toxicology
assessment, not a clinical safety evaluation, and not a substitute for GLP
in vitro or in vivo safety pharmacology studies.
