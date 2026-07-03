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

## Full 5-fold ensemble re-scoring (more reliable shortlist)

Re-scored all 323 lipids with the **full 5-fold ensemble** (30 epochs/fold, held-out
RMSE 0.78–0.86). The ensemble mean **and its cross-fold std** (a confidence signal)
give a much more trustworthy shortlist. Output:
`data/libraries/combo_v1/combo_v1_shortlist_top50.csv`;
figure `results/figures/combo_v1_ensemble.png`.

**The ensemble reranks meaningfully vs the lite model** (Spearman ρ = 0.78; only
19/50 top-50 overlap) — so the full run was worth it. Two shifts:
- Head/linker matters more than the lite model implied: **H7 (N-Me-1,3-diaminopropane)
  + amide, 3-tail** tops the raw mean across many tails.
- But those raw-top H7+amide picks have **high cross-fold disagreement** (std 0.7–0.8):
  high mean, low confidence.

**Confidence-aware leads** (high mean ≥0.6 AND low std ≤0.4 — all 5 folds agree):

| head | linker | tail | n_tails | MW | ensemble pred |
|---|---|---|---|---|---|
| H7 (N-Me-DAP) | ester | linoleyl | 3 | 1050 | **1.06 ± 0.36** |
| H10 (DMAP-amine) | ester | linoleyl | 2 | 743 | **0.99 ± 0.32** |
| H10 (DMAP-amine) | ester | oleyl | 2 | 747 | **0.88 ± 0.28** |
| H7 (N-Me-DAP) | ester | oleyl | 3 | 1056 | 0.84 ± 0.32 |
| H2 (DEAE-amine) | ester | linoleyl | 2 | 757 | 0.78 ± 0.27 |

These are **ester-linked, unsaturated (linoleyl/oleyl), 2–3 tail, MW 740–1050** —
synthesizable and consistent with proven ionizable-lipid chemistry. Of the top-50:
20 high / 16 medium / 14 low confidence.

SAR trends survive the ensemble: unsaturated ≫ saturated (di-unsat 0.51 > mono 0.46
> branched 0.28 > saturated 0.15 > hetero 0.07); 3-tail > 2-tail ≫ 4-tail
(0.30 / 0.24 / 0.04 — the 4-tail penalty is *larger* in the ensemble, and 4-tail
picks carry the highest uncertainty).

**Recommendation**: prioritise the 5 confidence-aware leads above (not the raw-mean
top, which the folds disagree on). The ensemble was trained preemption-resiliently
via `modal_app/lion.py::train_resilient` (per-fold persist + resume).

## Reproduce
```bash
python analysis/enumerate_library.py                       # -> 323 lipids + screen input
modal run modal_app/lion.py::screen \
    --smiles-json "$(python -c 'import pandas,json;print(json.dumps(pandas.read_csv("data/libraries/combo_v1/combo_v1.csv").smiles.tolist()))')" \
    --folds 1 --name combo_v1
# then merge screen__combo_v1.csv (volume) with combo_v1_enumerated.csv
```
