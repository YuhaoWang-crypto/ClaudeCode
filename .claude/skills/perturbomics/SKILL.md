---
name: perturbomics
description: >-
  Integrate and combine perturbation-omics data across modalities — drug
  perturbation (Broad Drug Repurposing Hub / Connectivity Map L1000), CRISPR
  perturbation (Broad GPP CRISPick libraries: Brunello/Brie KO, Dolcetto CRISPRi,
  Calabrese CRISPRa; Perturb-seq), single-cell differential expression, and the
  Geneformer foundation model (embeddings + in-silico perturbation) — by
  reducing every perturbation to one common object, a signed ranked gene
  SIGNATURE, then scoring connectivity (CMap WTCS) and mining drug+CRISPR
  COMBINATIONS that jointly reverse a disease state. Use when the task is: build
  a perturbation signature from raw single-cell counts (pseudobulk DGE); compare
  or cluster drug vs CRISPR vs disease signatures; find compounds/genes that
  reverse a signature; propose combination perturbations; screen an approved drug
  for a NEW indication (drug repurposing via the Drug Repurposing Hub —
  launched-status + current-indication novelty); fuse hits with pathway control +
  druggability + clinical status into integrated leads; or pull data from the
  Repurposing Hub, CMap/clue.io, GPP, or Geneformer. Runnable `perturbomics`
  package included. Enforces ✅-rigorous vs ⚠️-hypothesis labeling on every claim.
---

# Perturbation-omics integration & combination analysis

A reusable methodology (and a working `perturbomics` Python package) that turns
four disconnected data worlds — **drug** perturbation, **CRISPR** perturbation,
**single-cell DGE**, and a **foundation-model** — into one comparable space,
so you can ask cross-modality questions like *"which drug + gene knockout
combination reverses this fibrosis signature?"* — and always label what is a
rigorous computation vs. a biological hypothesis.

## The one idea that unifies everything

Every perturbation, no matter the modality, is reduced to the **same object**:

> a **signature** = a signed, ranked vector over genes
> (`+` = induced by the perturbation, `−` = repressed).

| Source | What it natively gives | → coerced to a signature via |
|---|---|---|
| **Drug Repurposing Hub / CMap L1000** | moderated z-score per landmark gene, per compound × cell line × dose | `Signature.from_l1000` |
| **CRISPR screen** (Brunello KO / Dolcetto CRISPRi / Calabrese CRISPRa, or Perturb-seq) | log2FC / MAGeCK score per gene, or DGE per guide | `Signature.from_deseq2` / `from_ranked` |
| **Single-cell perturbation** (sc-best-practices) | pseudobulk DESeq2 Wald stat per gene | `signature_from_pseudobulk` |
| **Geneformer** | cosine embedding shift per gene under in-silico perturbation | `Signature` (see `reference/geneformer.md`) |

Once everything is a `Signature`, **one metric** — the Connectivity-Map
**weighted connectivity score (WTCS)** — compares *any* two of them:
drug↔drug, drug↔CRISPR, perturbation↔disease. `+1` = mimic, `−1` = reversal.
That single fact is what makes combination analysis possible.

## The pipeline

```
raw single-cell counts ──pseudobulk DGE──┐
Drug Hub / CMap L1000 z-scores ──────────┤
CRISPR KO/i/a screen scores ─────────────┼──► Signature (signed ranked genes)
Geneformer in-silico perturbation ───────┘            │
                                                      ▼
                          connectivity (WTCS) ──► reversers / mimics / clusters
                                                      │
                                                      ▼
                          combination analysis ──► drug + CRISPR pairs that
                                                   jointly reverse a disease
```

## Run the package

```bash
pip install numpy pandas scipy          # core; demos need only these
python3 -m perturbomics.demo            # synthetic end-to-end (offline, instant)
python3 -m perturbomics.realdata_ipf    # REAL data: downloads Enrichr libraries
```

`demo.py` plants a disease signature and a mixed drug/CRISPR library, then shows
(1) single-agent reversers by WTCS, (2) a connectivity map, (3) the best
drug+CRISPR **combinations** by disease-coverage, and (4) a rigorous pseudobulk
DGE recovering the planted genes. For real DGE also `pip install pydeseq2`
(auto-detected; falls back to a numpy/scipy path otherwise).

`realdata_ipf.py` runs the SAME pipeline on **genuine public data** — it
downloads three Enrichr libraries (CREEDS disease DE, LINCS L1000 drug, LINCS
CRISPR-KO signatures; ~45 MB, no login), builds a consensus idiopathic-
pulmonary-fibrosis signature, and ranks real drugs + gene-KOs that reverse it
plus the best drug+CRISPR combinations. Reproducible verified run (PYTHONHASHSEED=0):
top drug reversers **trichostatin A** (HDAC inhibitor), **mln4924**,
**canertinib** (pan-ErbB); top CRISPR-KO reversers **CDK13**, **NBAS**,
**KIAA0907**; top cross-modality combination **canertinib + CDK13-KO** (covers
~11% of the IPF signature, +3% from the genetic partner). All hits are ⚠️
hypotheses to validate — the point is the *method* runs on real data. The
Enrichr loaders live in `enrichr.py` (`load_library`, `paired_signatures`,
`crispr_ko_signatures`, `consensus_signature`).

Package layout (`assets/perturbomics/`):

- `signature.py` — the `Signature` class + constructors from every source.
- `connectivity.py` — `enrichment_score` (GSEA ES), `weighted_connectivity_score`
  (WTCS), `connectivity_matrix`, `normalized_connectivity` (NCS).
- `combine.py` — `rank_reversers`, `combination_score`, `best_combinations`
  (with `require_cross_modality` for the drug+CRISPR question).
- `pseudobulk.py` — `pseudobulk` + `signature_from_pseudobulk` (PyDESeq2 or fallback).
- `demo.py` — the runnable end-to-end example above.

## Combination analysis — what it actually computes

Given a **disease** signature and a mixed library:

1. `rank_reversers` — every drug/CRISPR perturbation scored by how strongly it
   *reverses* the disease (most-negative WTCS first), drug and CRISPR hits in
   **one ranked table**.
2. `best_combinations` — all pairs scored on **coverage**: how much of the
   disease signature the two *jointly* reverse, where the combined reversal at
   each gene is the max of the two partners (so covering **different** genes
   adds up, covering the **same** gene doesn't double-count). A
   **complementarity** term rewards pairs that each fix what the other misses.
   `require_cross_modality=True` restricts to **drug + CRISPR** pairs.

The payoff (visible in the demo): a coherent *single* reverser may win WTCS,
but two *partial* agents that individually score ~0 can together cover ~90% of
the disease — combination analysis surfaces exactly what single-agent
connectivity misses.

## Composing with pathway analysis & the drug-discovery pipeline

perturbomics is the transcriptomic **front end** of a larger, multi-scale funnel.
It nominates *what* to perturb; two sibling capabilities then say *where it acts*
and *whether the molecule is real* — all three share one vocabulary (a **target
gene** + a **perturbation magnitude**), so they compose with a thin bridge
(`integrate.py`), no rewrites:

1. **perturbomics** → ranked reversers + drug/CRISPR combinations (`|WTCS|`).
2. **network-biomarker** (pathway/dynamics skill) → is the target an
   **irreducible-core** (M1/M11) or **bistable-switch** (M2/M19) node? what
   **early-warning biomarker** (M4 DNB) tracks a responder?
3. **drug-discovery MCP servers** (ChEMBL/Boltz/ClinicalTrials — the same ones
   network-biomarker's m6–m9 use) → **engageability** (pIC50 + Boltz binding/ADME)
   and **clinical** phase.

`integrated_leads(reversers, network, evidence)` fuses these four orthogonal axes
into one ranked lead table (transparent weights, per-row renormalisation, missing
axes flagged not zeroed). The payoff: a strong-but-promiscuous connectivity hit
gets demoted and a moderate hit with a core+switch, druggable, clinical target
gets promoted — which no single axis would surface. And via network-biomarker's
`m6_integrate` a nominated compound can be pushed *back* through the dynamics to
predict how far it moves the disease network toward its tipping point.

```bash
python3 -m perturbomics.demo_integrate     # offline 4-axis funnel, deterministic
```

Full wiring (which `report()` keys and MCP calls feed each axis) →
`reference/integration.md`.

## Drug repurposing (reuse an approved drug for a new indication)

The Connectivity Map's founding use case, and a thin specialisation of the above:
a repurposing candidate is a drug that (a) **reverses** the disease signature,
(b) is **already launched** (for something else), and (c) is **not already used
for this disease**. `repurpose.py` adds exactly that on top of `rank_reversers`,
using the real **Broad Drug Repurposing Hub** annotations (`clinical_phase`,
`moa`, `target`, `disease_area`, `indication`) for approval status + current use:

```python
from perturbomics import rank_reversers, load_repurposing_hub, screen_repurposing
rev  = rank_reversers(disease_sig, drug_library)
hub  = load_repurposing_hub("repurposing_drugs_20200324.txt")
cands = screen_repurposing(rev, hub, disease_terms=["pulmonary","fibros"],
                           require_launched=True)   # score = |WTCS| × phase × novelty
```

Cross the shortlist against ClinicalTrials to split **de-risked** (already in
trials for the disease) from **novel white-space** (untried) candidates. Real IPF
run (`examples/repurpose_ipf.py`): surfaced **dasatinib** (CML → *already* in an
IPF Ph1 trial, NCT02874989 — de-risked/confirmatory) and **neratinib** (breast
cancer, EGFR inhibitor → no IPF trials — a novel hypothesis that aligns with the
EGFR/canertinib finding). Full method + the real result table →
`reference/repurposing.md`.

## The non-negotiable discipline: honesty labeling

Every result carries one of:

- **✅ rigorous** — a deterministic computation: an enrichment score, a WTCS, a
  coverage fraction, a pseudobulk DESeq2 Wald statistic, an exact z-score read
  from L1000.
- **⚠️ hypothesis** — a biological interpretation that needs orthogonal or
  experimental validation: *"high connectivity ⇒ shared mechanism"*, *"this
  reverser is therapeutic"*, *"this pair is synergistic"*. A WTCS is a
  similarity number; a synergy claim requires an actual combination assay
  (Bliss/Loewe on a dose matrix).

Never blur the two. The package docstrings mark every function this way; carry
it into any report. Report negative/partial results as findings — e.g. random
decoy sets can score as moderate reversers under raw WTCS, which is why
`normalized_connectivity` (NCS) and a null distribution matter before you
believe a hit.

## Common tasks → where to look

- **Get the data** (Repurposing Hub download files & columns, CMap/clue.io
  `.gctx` + `cmapPy`, GPP CRISPick libraries, Geneformer weights, and the
  in-session MCP servers — ChEMBL, PubMed, ClinicalTrials, bioRxiv — that enrich
  a hit) → `reference/data-access.md`.
- **Compute a signature from raw single cells the right way** (pseudobulk,
  avoid the pseudoreplication trap, PyDESeq2 design/contrast) →
  `reference/dge-signatures.md`.
- **The connectivity math** (GSEA ES, WTCS, NCS/τ, thresholds, null models) →
  `reference/connectivity.md`.
- **Geneformer** (embeddings, in-silico perturbation, turning its output into a
  `Signature`) → `reference/geneformer.md`.

## Hard-won gotchas

- **WTCS zeros one-tailed connections by design.** A partial reverser that only
  moves the up-half of a disease scores 0, not weakly-negative, because the
  up/down tags must land at *opposite* ends to count as coherent. Don't read
  that 0 as "no effect" — that's precisely the case `best_combinations` is for.
- **Never test perturbation DGE per-cell.** Cells from one sample aren't
  independent replicates; Wilcoxon-on-cells inflates false positives. Aggregate
  to pseudobulk (sample × condition) first. (`reference/dge-signatures.md`)
- **Match the gene space.** L1000 gives 978 landmark genes (+11k inferred);
  Geneformer tokenizes by Ensembl ID over ~25k genes; CRISPR screens are gene-
  level. Map to a shared identifier (Ensembl) and universe before scoring, or
  connectivity silently compares mismatched sets.
- **Sign conventions bite.** KO/CRISPRi remove a gene's function (loss); CRISPRa
  and cDNA-overexpression add it (gain). A "reversal" between a KO signature and
  a drug signature only means what you think if both are oriented the same way
  (perturbation-induced Δexpression). State the convention.
- **Raw WTCS is not calibrated, and a perfect-looking score can be fragile.**
  Normalize within cell-line × perturbation-type groups
  (`normalized_connectivity`) and run `permutation_pvalue` before calling a hit.
  On real IPF data the top raw reverser scored WTCS −1.000 yet **p ≈ 0.23** (it
  rested on a few genes) — the null caught what the raw score hid. Decoys can
  score moderately; the p-value is the guard.
- **Geneformer's corpus is biased** (~9M brain cells, donors skewed <1yr old,
  non-diseased). In-silico perturbation shifts are hypotheses to confirm in the
  relevant tissue, not ground truth.
