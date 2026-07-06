# Biomarker Discovery Workflow (reference)

## 1. BEST framework — define the question first
Map every clinical use to an FDA-NIH BEST category. Same gene/protein can play different roles in different contexts, so anchor on: **which category, versus what comparison group, measurable in what sample?**

| Clinical use | BEST category | Question | Primary engine entry |
|---|---|---|---|
| Onset risk | Susceptibility/Risk | Will they develop disease? | Germline variants (ESM-2 VEP + PROTO) |
| Diagnosis/subtyping | Diagnostic | Present? which subtype? | Geneformer cell-state clustering |
| Enrollment / companion Dx | **Predictive** | Respond to which drug class? | Geneformer response-states + serum + variants |
| Efficacy | Response (PD) | Is the drug working? | Serum PD proteins; Geneformer perturbation |
| Degree of treatment | **Monitoring** | Disease burden trend? | Longitudinal serum + residual cell states |
| Progression | Prognostic | How will it evolve? | Geneformer network drift + clonal evolution |

## 2. Six-step process (indication-agnostic)
1. **Define**: indication + 1–2 BEST categories + explicit comparison groups (responder vs non-responder / progressor vs stable / disease vs healthy).
2. **Ground evidence**: PubMed / bioRxiv / ClinicalTrials / ChEMBL → known biomarkers, competitive baseline, novelty gaps.
3. **Generate candidates** on 3 omic legs: Geneformer (scRNA cell states + in-silico perturbation), ESM-2 VEP (coding variants), PROTO (regulatory/splice variants via Evo2/AlphaGenome/SpliceAI).
4. **Mechanistic filter**: Boltz (structure/binding), EDEN (immunogenicity), ChEMBL (druggability/CDx pairing).
5. **Clinical triage**: three gates — novelty, accessibility (blood/stool measurable?), clinical alignment (endpoints, CDx drug).
6. **Combine** into per-mechanism scorecards (a response probability per drug class + a disease-burden trend), not single markers. Multi-modal cross-corroboration (variant + expression + protein agree) is the core false-positive control.

## 3. Port to a new indication — swap these 6, keep everything else
| # | Swap | Source | File to edit |
|---|---|---|---|
| 1 | Indication + comparison groups | clinical question | panel meta |
| 2 | scRNA/omics dataset | GEO, Broad Single Cell Portal, CELLxGENE, HCA | scRNA_datasets.md |
| 3 | Mechanism/drug classes | drugs in development → Trials/ChEMBL | candidates.yaml `mechanism_classes` |
| 4 | Disease-specific gene module | literature + de novo discovery | `INFLAMMATION_MODULE` in geneformer_denovo_discovery.py |
| 5 | Candidate variant list | GWAS Catalog / Open Targets | proto_variants_*.tsv |
| 6 | Validation cohorts | indication's trials (NCT) | panel `validation_cohort` |

Reused unchanged: candidates.yaml schema, score_candidates.py, both Geneformer pipelines, report layout, the PubMed/Trials/ChEMBL grounding loop.

### Porting examples
| Step | NSCLC immunotherapy | MASH/NASH (liver) |
|---|---|---|
| BEST focus | Predictive/CDx (beyond PD-L1/TMB) | Diagnostic + Prognostic (fibrosis) |
| Comparison | ICI responder vs primary-resistant | progressive fibrosis vs stable |
| scRNA | tumor/immune scRNA (GEO/CELLxGENE) | liver scRNA atlas |
| Disease module | T-cell exhaustion / antigen presentation | stellate-cell activation / fibrosis |
| Variants | HLA / JAK pathway | PNPLA3 / TM6SF2 / HSD17B13 |
| Drug classes | anti-PD1 / anti-CTLA4 / LAG3 | FGF21 / THRβ / GLP1 |

## 4. Worked reference: IBD (autoimmune/inflammation)
Flagship because ~40% primary non-response across 5 mechanism classes (anti-TNF / IL-12·23 / IL-23 / α4β7 / JAK) with no validated drug-selection biomarker. Highest-novelty gaps = IL-23 and JAK response prediction (near-empty literature; genetics supports predictability: IL23R R381Q loss-of-function protective; TYK2 P1104A/A928V/I684S protective). Example candidate panel (`candidates.yaml`) covers OSM/OSMR, serum Olink cluster, TREM1 (discovery-positive but RCT-negative — a triangulation lesson), GIMATS module (Martin 2019 Cell), α4β7/MAdCAM, IL23R, TYK2. Diagnostic subtyping anchor: Smillie 2019 Cell IL13RA2+IL11+ inflammatory fibroblasts (anti-TNF resistance).

## 5. Guardrails
- Attribute PubMed + DOI links for cited articles.
- Discovery ≠ validation; flag grade-1 candidates as hypotheses.
- Record contradictions; stratify by subtype/age/sample; guard batch effects.
- Never fabricate rsIDs/citations; mark unconfirmed `to_validate`.
- Geneformer needs REAL scRNA for real signal (encodes learned co-expression); use pretrained (not fine-tuned) model with matching vocab.
