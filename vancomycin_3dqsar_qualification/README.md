# In silico toxicological qualification of vancomycin B impurities RRT 0.75 & RRT 0.87

3D-QSAR / structure-based **read-across** qualification of two related substances in
Xellia vancomycin, delivered against the sponsor brief:

1. 3D-QSAR / structure-based assessment of the **parent** (Vancomycin B) and the **two daughter** molecules;
2. identification/assessment of any **additional toxophore** present in a daughter but absent from the parent;
3. a documented **read-across justification**;
4. a **regulatory-style report**.

## Test items

| Item | Identity | Formula | Relationship to API |
|---|---|---|---|
| Parent | Vancomycin B (API) | C66H75Cl2N9O24 | — |
| RRT 0.87 | (7S,26R)-Vancomycin B | C66H75Cl2N9O24 | **Diastereomer** (epimerised at C7 & C26); heat/low-water-activity degradant |
| RRT 0.75 | (26R)-DMe-DeLI | C65H73Cl2N9O24 | Residue-1 variant (N-desMe; Leu→Ile) + C26 epimer; fermentation impurity |

## Headline result

**Neither impurity introduces a new toxicophore or a new DNA-reactive (ICH M7) structural alert absent from the parent.**

- **RRT 0.87** is a pure stereoisomer: **all 101 heavy atoms** and the complete functional-group inventory are shared with the API; 2D fingerprints are identical (Tanimoto = 1.00), so the difference is *purely 3D* — hence the 3D/mechanistic treatment.
- **RRT 0.75** shares **99/100 heavy atoms**; its only constitutional change is a secondary (N-methyl) → primary aliphatic amine at residue 1 — a group class already present in the parent (vancosamine) — with unchanged logD (0.42) and near-identical basicity (pKa ≈ 8.1 vs 8.2). The basic-amine / "Red-Man" liability is conserved, not amplified.

Both impurities are qualified by **read-across to the parent**, whose safety is anchored by the
91-day i.v. dog study (NOAEL 75/100 mg/kg/day, ≈2–3× clinical dose), provided each stays within
the Ph. Eur. related-substances limits. Confirm under the company ICH M7 two-(Q)SAR + expert-review workflow.

## Statistical QSAR (12-endpoint Tox21)

An in-house 12-endpoint Tox21 QSAR was trained (class-balanced RandomForest on 2048-bit radius-2
Morgan fingerprints; public MoleculeNet Tox21; **mean 5-fold CV-AUC 0.827**, range 0.72–0.89) and
applied with an explicit applicability-domain (AD) check:

- **Vancomycin is itself a Tox21 training compound** (mol_id TOX25354, nearest-neighbour Tanimoto = 1.00),
  so the QSAR is genuinely *in-domain* for the parent. Its **measured** profile is essentially clean:
  weak **NR-AR** active, **inactive** at every other tested endpoint (incl. genotoxicity SR-ATAD5, SR-p53).
- **RRT 0.87** has a fingerprint identical to the parent → **Δ ≡ 0** at all 12 endpoints (a concrete demonstration
  that a 2D-QSAR is blind to this stereoisomer — hence the 3D read-across carries the weight).
- **RRT 0.75** (genuine external prediction, Tanimoto 0.86): every endpoint probability is within **±0.02**
  of the held-out parent and **none crosses into activity** → no new nuclear-receptor or stress-response liability.

The Tox21 panel covers 12 specific mechanisms only; it does **not** model vancomycin's clinical effects
(nephro-/ototoxicity, Red-Man syndrome) — those are handled mechanistically (toxicophore analysis) and by
read-across to the parent's in vivo package. See report §5.6.

## ICH M7 mutagenicity — two complementary (Q)SAR methodologies

Per ICH M7(R2), bacterial (Ames) mutagenicity was screened by two independent methodologies:

- **Methodology 1 — expert rule-based:** a **Benigni–Bossa** structural-alert rulebase (30 alert classes;
  the public scientific basis of Derek-type systems, as in Toxtree/VEGA).
- **Methodology 2 — statistical:** an **Ames QSAR** (class-balanced RandomForest, Morgan r2) trained on the
  public **Hansen benchmark (N≈6500)**, 5-fold OOF **AUC 0.881**.

| Molecule | Expert rulebase | Statistical Ames | M7 class |
|---|---|---|---|
| Parent | no alert | negative (p=0.41) | **5** |
| RRT 0.87 | no alert | negative (p=0.41) | **5** |
| RRT 0.75 | no alert | negative (p=0.42) | **5** |

Both methodologies **concordant and negative → ICH M7 Class 5** (no structural alert → treated as non-mutagenic;
no dedicated Ames study needed). Caveat: the Ames statistical model is *out-of-domain* for these glycopeptides
(Tanimoto ≈ 0.22), so the domain-independent **expert rulebase is the primary call**, corroborated by the parent's
experimentally inactive Tox21 DNA-damage endpoints (SR-ATAD5, SR-p53). A licensed Derek/Sarah Nexus run should
confirm in the validated workflow. See report §5.7.

## Report formats

- `..._Report.docx` — editable regulatory report (10 sections incl. §5.6 QSAR and §5.7 ICH M7; 7 tables, 8 figures).
- `..._Report.pdf` — 14-page rendered PDF (identical content; produced from `report.html` via headless Chromium,
  since LibreOffice PDF export is unavailable in the build sandbox).
- `report.html` — self-contained HTML twin (embedded figures) used to generate the PDF.

## Deliverables

- `Vancomycin_RRT075_RRT087_3DQSAR_Qualification_Report.docx` — **the regulatory-style report** (main deliverable).
- `figures/` — structures, residue-1 difference, toxicophore-inventory heatmap, 3D/physchem descriptors, similarity metrics.
- `data/` — `analysis_results.json` (physchem + alerts + FG inventory + 3D descriptors + similarity), `mcs_core_alignment.json`, `summary_table.csv`.
- `scripts/` — reproducible pipeline.

## Reproduce

```bash
pip install rdkit matplotlib pandas scikit-learn
python scripts/analyze.py      # physchem, structural alerts, FG inventory, 3D descriptors, O3A/shape
python scripts/render.py       # 2D structures, residue-1 highlight, MCS core 3D alignment
python scripts/figures.py      # summary figures
# QSAR: needs data/tox21.csv.gz (public MoleculeNet Tox21; not committed)
python scripts/qsar.py         # train 12-endpoint Tox21 RF, 5-fold CV, predict, AD
python scripts/qsar_augment.py # measured vancomycin labels + held-out parent predictions
python scripts/qsar_fig.py     # QSAR figure
npm install docx && node scripts/build_report.js   # assemble the DOCX report (reads qsar_*.json)
```

## Methods (summary)

RDKit for structure handling, physicochemical/3D descriptors, PAINS(A/B/C)/Brenk/NIH structural-alert
catalogues, ETKDGv3 conformers + MMFF94, MCS-based core alignment and Open3DAlign/shape metrics;
Inductive Bio (Q)SAR for supporting logD/pKa. Framed by ICH Q3A/Q3B (qualification), ICH M7 (mutagenic
screening), and the OECD (Q)SAR validation principles.

## Honesty / limitations

- This is **3D-structure/shape read-across + mechanistic toxicophore profiling**, not a formal
  CoMFA/CoMSIA predictive 3D-QSAR (no congeneric activity training set exists for these endpoints).
- Small-molecule statistical tox panels (e.g. Tox21) were **deliberately not used as evidence**: a
  ~1435–1449 Da glycopeptide (100+ heavy atoms) is entirely outside their applicability domain; their
  output would convey false precision. logD/pKa predictions were flagged out-of-domain and used only directionally.
- 3D descriptors derive from single low-energy conformers of a large flexible macrocycle → corroborating,
  not decisive; the decisive evidence is the shared skeleton and conserved toxicophore inventory.
- A **screening / hazard-identification and read-across** exercise that supports — does not replace —
  expert toxicological review and any experimental testing a competent authority may require.
