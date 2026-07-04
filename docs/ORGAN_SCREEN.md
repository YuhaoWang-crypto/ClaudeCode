# Organ-context re-screen of the combinatorial library

Same 323 lipids, 5-fold ensemble, re-scored under four biological contexts by
conditioning the model's context features (delivery target + route + cell/model).
Data: `data/libraries/combo_v1/combo_v1_organ_compare.csv` ·
figure `results/figures/combo_v1_organ.png`.

| context | delivery_target | route | model |
|---|---|---|---|
| HeLa (baseline) | generic_cell | in_vitro | HeLa |
| lung | lung_epithelium | intratracheal | Mouse |
| liver | liver | intravenous | Mouse |
| spleen | spleen | intravenous | Mouse |

## How much does the ranking change?

Overall ordering is **broadly conserved** — a good lipid tends to be good
everywhere — but not identical (Spearman of per-lipid scores):

| | lung | liver | spleen |
|---|---|---|---|
| lung | — | 0.87 | 0.85 |
| liver | | — | 0.93 |
| **HeLa** | 0.76 | 0.80 | **0.60** |

The in-vitro **HeLa ranking transfers worst to spleen (0.60)** — spleen has the
most distinct structure–activity profile.

## The SAR shifts by organ (the interesting part)

1. **Tail-number optimum moves** (panel A): **spleen favours 2-tail** lipids;
   **lung/liver favour 3-tail**; 4-tail is penalised everywhere (worst in liver).
2. **Head preference** (panel B): lung/liver prefer **H7** (N-Me-1,3-diaminopropane);
   **spleen strongly favours the small 2-tail heads H1/H2/H10** (dimethyl/diethyl-
   aminoethyl-amines) — a clear spleen signature.
3. **Tail chemistry flips vs in-vitro** (panel C): in HeLa, di-unsaturated
   (linoleyl) tails dominated; **in the in-vivo organ contexts unsaturated tails are
   disfavoured** and saturated/branched do relatively better (branched peaks for
   spleen). The model learned an in-vitro↔in-vivo tail-chemistry difference.

## Best candidate per organ

`H7 (N-Me-DAP) + amide + C9` (T13, 3-tail, MW 680) is the **top pick in all three
organs** (lung 0.49, liver 0.48, spleen 0.57) — a robust broadly-good lipid.

**Organ-selective** leads (highest *relative* preference for one organ over the
others) are chemically distinct:
- **lung**: H6 + amide + C9/C10, **4-tail** (bigger lipids)
- **liver**: H13 (aromatic xylylene head) + **ester** + C16/C14 (long saturated)
- **spleen**: H1/H10 + short tails (C6/C8), **2-tail** (small lipids)

This mirrors the literature intuition that organ tropism tracks distinct lipid
physicochemistry (Su 2026: linker/head changes redirect organ targeting via the
protein corona), and gives concrete, testable per-organ shortlists.

## Caveat
These are the model **conditioning on context**, and LiON's in-vivo organ data is
sparser than its in-vitro data — the tail-chemistry flip in particular may partly
reflect training-set composition. Treat per-organ picks as **hypotheses to test**,
strongest where the ensemble agrees (low cross-fold std, see COMBO_LIBRARY.md).

## Reproduce
```bash
for c in lung_epithelium_IT:combo_lung liver_IV:combo_liver spleen_IV:combo_spleen; do
  modal run modal_app/lion.py::screen --smiles-json "$(...)" --folds 5 \
      --context "${c%%:*}" --name "${c##*:}"
done
```
