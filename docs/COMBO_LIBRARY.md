# Combinatorial lipid library — enumerate → LiON lite → rank (Problem A)

Enumerate: `analysis/enumerate_library.py` · Score: `modal_app/lion.py::screen(folds=1)`
· Ranked output: `data/libraries/combo_v1/combo_v1_ranked.csv` ·
Figure: `results/figures/combo_v1_screen.png`

## What was built
A Su-2026-style combinatorial ionizable-lipid library via **aza-Michael chemistry**
(the real synthetic route): 14 amine heads × {ester, amide} linkers × 16 tails.
Each head N-H is alkylated with an acrylate (ester linker) or acrylamide (amide
linker) tail, so tail count = the head's reactive amine count. Enumerator lives in
`lipidlib/lion_library.py::enumerate_michael_lipids`.

- 448 combos → **323 valid lipids** after RDKit validity + MW 400–1400 + ≥2-tail
  filters. Mix: 139 two-tail, 84 three-tail, 100 four-tail.
- Scored on Modal with the **lite (1-fold, 15-epoch) LiON model** at the KK/HeLa
  in-vitro baseline, in ~6 s.

> These blocks are *inspired by* Su's H1–H14 / T1–T16 (we don't have their exact
> SMILES), but chemically representative of the class.

## The model recovered known structure–activity rules

Two independent literature SAR trends fell out of the ranking — a strong sanity
check that `enumerate → featurize → LiON → rank` produces real chemistry:

1. **Unsaturated tails win** (mean predicted delivery, z):
   di-unsaturated **0.585** ≫ mono-unsat 0.287 > branched 0.162 ≈ saturated 0.146
   > hetero 0.043. All top-15 lipids carry the linoleyl (di-unsaturated) tail.
   → matches Su 2026 ("one or two double bonds dramatically increased expression")
   and the MC3/Moderna lipid design.
2. **3 tails are optimal**: 3-tail mean 0.239 > 2-tail 0.170 > 4-tail 0.108.
   → matches Su 2026's headline "3-tail ionizable lipids (3TILs) superior overall".

Head trend: polyamine cores H6 (1,3-diaminopropane), H4 (aminoethyl-piperazine),
H7 (N-Me-diaminopropane) rank highest. Linker: amide ≈ ester (amide marginally
higher).

## Top candidates (HeLa/in-vitro, KK formulation)

| rank | head | linker | tail | n_tails | MW | pred (z) |
|---|---|---|---|---|---|---|
| 1 | H6 (1,3-DAP) | ester | linoleyl | 4 | 1356 | 0.78 |
| 2 | H6 (1,3-DAP) | amide | linoleyl | 4 | 1352 | 0.78 |
| 3 | **H7 (N-Me-DAP)** | amide | linoleyl | **3** | 1047 | 0.74 |
| 4 | **H4 (AE-piperazine)** | ester | linoleyl | **3** | 1091 | 0.74 |
| 5 | **H4 (AE-piperazine)** | amide | linoleyl | **3** | 1088 | 0.74 |

The **3-tail linoleyl** hits (ranks 3–6: H7/H4) are the practical sweet spot —
they combine the 3-tail optimum with di-unsaturation, at a more synthesizable
MW (~1050–1090) than the 4-tail leads.

## Caveats
- Lite model = single fold, so these are **preliminary rankings**, not validated
  leads; re-score with the full 5-fold ensemble before committing.
- Context is HeLa/in-vitro (KK). For organ-targeted screens change the context in
  `enumerate_library.py` (`Formulation(delivery_target=…, route=…)`) — the model
  conditions on it.

## Reproduce
```bash
python analysis/enumerate_library.py                       # -> 323 lipids + screen input
modal run modal_app/lion.py::screen \
    --smiles-json "$(python -c 'import pandas,json;print(json.dumps(pandas.read_csv("data/libraries/combo_v1/combo_v1.csv").smiles.tolist()))')" \
    --folds 1 --name combo_v1
# then merge screen__combo_v1.csv (volume) with combo_v1_enumerated.csv
```
