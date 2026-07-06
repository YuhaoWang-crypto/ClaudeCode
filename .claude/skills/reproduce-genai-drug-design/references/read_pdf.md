# Reading the paper PDF

`pdftoppm`/`pdftotext` (poppler) are usually missing and apt often can't install
them. Use pypdf for text extraction:

```python
import pypdf                       # pip install pypdf ; may need: pip install --force-reinstall cffi
r = pypdf.PdfReader("paper.pdf")
print(len(r.pages))
for i in range(len(r.pages)):
    print(f"=== PAGE {i+1} ===")
    print(r.pages[i].extract_text())
```

## What to extract (drives the whole reproduction)

- **Target + PDB codes** — wild-type, mutant, and any off-target structure
  (e.g. DprE1 4KW5, mutant 5OEL, CYP2C9 5W0C). Note the native ligand + cofactor
  HET codes for each.
- **Hit / reference compound** and its measured activity.
- **TPP / objectives** — which properties were optimised (activity, ADME, tox).
- **ML model** — type (RF/Bayesian), training-set size, descriptor list, active
  threshold (e.g. pIC50 ≥ 5.75), reported ROC.
- **Candidate SMILES** — almost always in an appendix ("IUPAC and SMILES
  (ChemDraw …)"). These are gold: they let you validate structures and score them.
- **Data availability** — Zenodo/figshare DOI. Get the file download URL from the
  API: `https://zenodo.org/api/records/<id>` → each file's `links.self`
  (the human `/records/<id>/files/<name>?download=1` URL often 404s).
- **Table 3 / Table 4 numbers** — docking scores and MD observables to compare
  against. Transcribe them into the analysis scripts as the reference.
