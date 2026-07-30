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

## Deliverables

- `Vancomycin_RRT075_RRT087_3DQSAR_Qualification_Report.docx` — **the regulatory-style report** (main deliverable).
- `figures/` — structures, residue-1 difference, toxicophore-inventory heatmap, 3D/physchem descriptors, similarity metrics.
- `data/` — `analysis_results.json` (physchem + alerts + FG inventory + 3D descriptors + similarity), `mcs_core_alignment.json`, `summary_table.csv`.
- `scripts/` — reproducible pipeline.

## Reproduce

```bash
pip install rdkit matplotlib
python scripts/analyze.py      # physchem, structural alerts, FG inventory, 3D descriptors, O3A/shape
python scripts/render.py       # 2D structures, residue-1 highlight, MCS core 3D alignment
python scripts/figures.py      # summary figures
npm install docx && node scripts/build_report.js   # assemble the DOCX report
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
