# GLP1R cross-comparison: ChEMBL × DrugCLIP × humanPPI

Analysis: `analysis/crosscompare_glp1r.py` · Figure: `results/figures/glp1r_crosscompare.png`
· Ranked output: `data/targets/GLP1R/GLP1R_consensus.csv`

## Question
Do DrugCLIP's *predicted* GLP1R binders agree with the *measured* GLP1R actives
in ChEMBL? Are they enriched for known-active chemistry (i.e. is there a
high-confidence consensus set)?

## Method
- ChEMBL actives restricted to the small-molecule subset (MW < 900; n = 693) for a
  like-for-like comparison with DrugCLIP's hits (n = 173).
- Morgan fingerprints (radius 2, 2048 bit). For each DrugCLIP hit, nearest-neighbour
  Tanimoto to any measured active.
- **Enrichment null**: DrugCLIP's predicted hits for an unrelated target, EGFR
  (a kinase; n = 874). A target-specific method should place GLP1R hits closer to
  GLP1R actives than EGFR hits are.

## Result — no consensus, no enrichment (an honest negative)

| metric | value |
|---|---|
| Exact InChIKey overlap (ChEMBL ∩ DrugCLIP) | **0** |
| Connectivity-layer overlap | **0** |
| Median nearest-active Tanimoto — GLP1R hits | 0.167 |
| Median nearest-active Tanimoto — EGFR decoy | 0.165 |
| Mann–Whitney U (GLP1R > decoy), one-sided | **p = 0.20 (n.s.)** |
| Enrichment factor @ Tanimoto ≥ 0.35 | **0.0×** |
| Spearman(DrugCLIP score, active-similarity) | 0.13 (p = 0.1) |

**DrugCLIP's GLP1R predictions are no more similar to measured GLP1R actives than
an unrelated kinase's predictions are.** Not one of the 173 hits is even
scaffold-similar (Tanimoto ≥ 0.35) to a known active.

## Why (Panel B)
The two sets occupy different chemical space by construction:
- **Measured actives**: MW median **607** Da — the modern oral-agonist chemotype
  (danuglipron/orforglipron class), large and specific.
- **DrugCLIP hits**: MW median **257** Da — **fragments** retrieved from
  ZINC/Enamine REAL.

A MW-257 fragment cannot be a close Tanimoto neighbour of a MW-607 drug, so a
fragment-library virtual screen will not recover the known GLP1R chemotype.
DrugCLIP also scored these against a **7TM pocket** (residues TYR250/GLU262/CYS174),
not the extracellular domain we want for a conjugatable targeting ligand.

## What this means for the project
1. **Do not treat the DrugCLIP GLP1R hits as validated binders.** They have no
   measured-pharmacology support and are fragment-sized. Use them only as
   fragment starting points, if at all.
2. **This is the cross-comparison working as intended** — it prevents us from
   chasing 173 unsupported "hits". The method (measured ∩ predicted) is the right
   filter; here it correctly returns ~empty.
3. **Reinforces the peptide-first decision** for GLP1R targeting: the highest-
   affinity, ECD-binding, conjugatable ligand is a GLP-1-derived peptide, not
   these small molecules.
4. **humanPPI** is a separate (protein) modality — 259 predicted partners, 181
   membrane-localized — not fused into this molecular comparison; relevant for
   PPI/peptide-derived targeting, analysed separately.

## When cross-comparison *will* pay off
This same pipeline is valuable where predicted and measured libraries overlap in
size/space — e.g. a target whose DrugCLIP hits are lead-sized, or when comparing
DrugCLIP against a fragment screen. Re-run for any target:
`python analysis/crosscompare_glp1r.py` after mining that target with the
`fetch_drugclip.py` / `fetch_glp1r_ligands.py` scripts.
