# Connectivity scoring — the math that links any two signatures

Reference: Subramanian et al., *A Next-Generation Connectivity Map: L1000
Platform and the First 1,000,000 Profiles*, **Cell 2017**. Implemented in
`perturbomics/connectivity.py`.

## Why it works cross-modality

Connectivity never looks at *how* a perturbation was made — only at the ranked
gene signature it produced. So a drug's L1000 profile, a CRISPR screen's per-gene
scores, a pseudobulk-DGE signature, and a disease signature are all just ranked
lists, and the same statistic compares any pair. That is the whole reason drug ↔
CRISPR ↔ disease comparison is possible.

## Step 1 — weighted running-sum enrichment score (GSEA ES) ✅ rigorous

`enrichment_score(ranked, tag_set, weight=1.0)`

Walk down the reference signature's genes sorted high→low. Increment a running
sum by a weighted amount at each **hit** (gene in the query tag set), decrement
by a constant at each **miss**. The ES is the signed maximum deviation of that
walk from zero. `weight=0` → classic Kolmogorov–Smirnov; `weight=1` (CMap
default) weights hits by `|score|`. ES ∈ [−1, 1]: near +1 = the tag set sits at
the **top** of the reference, near −1 = at the **bottom**.

## Step 2 — weighted connectivity score (WTCS) ✅ rigorous

`weighted_connectivity_score(query, reference, k=50)`

1. From the **query** signature take its top-`k` **up** and bottom-`k` **down**
   tag sets (`Signature.query_sets`).
2. Score each against the **reference's** ranked list:
   `ES_up`, `ES_down`.
3. Combine:

   ```
   WTCS = (ES_up − ES_down) / 2   if sign(ES_up) ≠ sign(ES_down)
        = 0                        otherwise
   ```

WTCS ∈ [−1, 1]. **`+` = mimic** (query up-tags high & down-tags low in the
reference), **`−` = reversal** (up-tags low, down-tags high — the therapeutic
direction).

**The zero-out rule is a feature, not a bug.** It demands a *coherent* two-tailed
connection. A perturbation that only moves one tail of the query (a **partial**
reverser) scores **0** — which is exactly why single-agent WTCS under-ranks
complementary agents and why `combine.py` exists (see below).

## Step 3 — normalization to NCS (and τ) ✅ rigorous rescaling

Raw WTCS isn't comparable across a heterogeneous reference (different cell lines,
perturbation types). CMap divides each score by the signed mean |score| within
its **cell-line × perturbation-type** class → **NCS**.
`normalized_connectivity(raw_matrix, groups=...)` does this. clue.io then
converts NCS to **τ (tau)**: the percentile of a connection against a large
reference compendium (|τ|=90 ⇒ stronger than 90% of reference connections).
τ needs that external compendium; NCS is what you compute locally.

## Using it

```python
from perturbomics import (weighted_connectivity_score, connectivity_matrix,
                          normalized_connectivity, rank_reversers)

wtcs = weighted_connectivity_score(disease, drug_sig, k=50)   # one pair
M    = connectivity_matrix([disease, *library], k=50)          # all pairs
NCS  = normalized_connectivity(M, groups={s.name: s.meta.get("cell_line","all")
                                          for s in library})
rank_reversers(disease, library)   # ranked table, most-negative first
```

## Interpretation discipline

- ✅ **rigorous:** the ES, the WTCS, the NCS, the coverage fractions — pure
  deterministic statistics on the vectors you supply.
- ⚠️ **hypothesis:** *"WTCS = −0.8 ⇒ this drug treats the disease"*, *"high
  positive connectivity ⇒ same mechanism"*. A similarity number is not a causal
  or therapeutic claim.

## Before you believe a hit

- **Null model.** Random gene sets can score non-trivial WTCS (the demo shows
  decoys reaching ≈−0.35). Compare against a permutation null or the NCS/τ
  percentile; don't threshold raw WTCS blindly.
- **Shared gene universe.** Score only on genes present in *both* signatures,
  mapped to a common ID (Ensembl). Mismatched universes silently deflate scores.
- **`k` sensitivity.** Very small `k` is noisy; very large `k` dilutes the
  signal. 50–150 tags is typical for L1000-scale signatures; report the `k` you
  used.
- **Direction/orientation.** Confirm both signatures encode the same Δ
  convention (perturbation-induced change) before reading a sign.
