---
name: molecular-toxicology
description: >-
  Structure-based molecular toxicology / side-effect assessment of a small
  molecule (given a SMILES or compound name). Combines ProTox-3.0 (predicted
  acute LD50 + toxicity class, verified working via API) with the session's
  stable structure-based tools (ChEMBL experimental hERG/CYP/ADMET, Inductive
  Bio physchem models, Boltz ADME), and an optional OMIM gene→disease
  annotation layer. Use when asked to predict a compound's toxicity, acute
  lethality, organ/off-target liabilities, or to assess side-effect risk from
  structure; or to run ProTox programmatically. Enforces ✅-computed vs
  ⚠️-inferred honesty labeling and states each prediction's source model.
---

# Structure-based molecular toxicology

Turn a molecular structure (SMILES / name) into a labeled toxicity profile by
combining several independent predictors, and **always** marking what is a real
model output vs. our own inference.

## What the user originally asked, and the honest answer

Two platforms were proposed as toxicology supplements:

- **ProTox-3.0** (tox.charite.de/protox3) — ✅ **usable via API**, ❌ weights are
  **not downloadable** (CC BY-ND, web/API only). Its official sample script is
  404; `assets/protox3_client.py` reproduces the real submission flow and is
  **verified** (aspirin → LD50 250 mg/kg, class 3). Rate limit **250/day/IP**.
- **OMIM geneMap** (omim.org) — ⚠️ has a REST API but needs a personal
  yearly-renewed key (not self-obtainable here; site is Cloudflare-gated → 403),
  and it is **not a structure-based toxicity predictor** — it is a gene↔disease
  catalog. Useful only as a downstream **gene→disorder annotation layer**.

Full evidence and the capability matrix: `reference/platforms.md`.

## The predictors and when each applies

| Tool | Gives you | Rigor | Access |
|---|---|---|---|
| **ProTox-3.0** `predict_acute` | predicted LD50 (mg/kg) + toxicity class 1–6, similarity, accuracy | ✅ model output; ⚠️ *prediction*, not measurement | `assets/protox3_client.py` (250/day) |
| **ProTox-3.0** `get_full_panel` | ~60 more endpoints (organ, Tox21, NRs, MIEs) | ⚠️ best-effort browser scrape | Playwright (preinstalled) |
| **ChEMBL** `get_bioactivity` | **experimental** hERG / CYP / P-gp liabilities | ✅ measured assay data | MCP tool |
| **ChEMBL** `get_admet` | calculated drug-likeness (Lipinski/Veber/QED) | ✅ deterministic from structure | MCP tool |
| **Inductive Bio** `predict_properties` | logD, pKa, … physchem drivers | ✅ model output | MCP tool |
| **Boltz** `start_small_molecule_adme` | structure→ADME | ✅ model output | MCP tool |
| **OMIM** (optional) | target gene → associated inherited disorders | ✅ curated DB | REST API, needs key |

Rule of thumb: **experimental beats predicted.** Use ChEMBL assay data for hERG
/ CYP where it exists; use ProTox for breadth of endpoints and acute lethality
where no assay is available. Never present a prediction as a measurement.

## Standard workflow

```bash
pip install rdkit          # once; Playwright already present for --full
python3 .claude/skills/molecular-toxicology/assets/protox3_client.py "<SMILES>"
```

1. **Acute lethality + breadth** — `predict_acute(smiles)` → LD50 + class.
   Add `--full` (or `get_full_panel`) for the 61-endpoint sweep (⚠️ brittle).
2. **Experimental off-target liabilities** — ChEMBL `target_search` +
   `get_bioactivity` on hERG (CHEMBL240), CYP3A4/2D6/2C9, P-gp (CHEMBL4302).
3. **Drug-likeness & physchem** — ChEMBL `get_admet`; Inductive Bio
   `predict_properties` (logD, pKa); Boltz ADME.
4. **Mechanism / susceptibility (optional)** — map flagged targets → OMIM
   disorders (only once an API key is provided).
5. **Report** — one row per endpoint: value, **source model**, and a
   ✅-computed / ⚠️-inferred tag. Aggregate to a risk summary that never
   over-claims.

## Honesty labeling (non-negotiable)

- ✅ **computed** — a direct output of a named model/DB (e.g. "ProTox LD50 =
  250 mg/kg", "ChEMBL hERG IC50 = 3 µM measured"). Always name the source.
- ⚠️ **inferred** — our own reasoning on top (e.g. "high logP + hERG hit ⇒
  likely cardiotoxic"). Mark it as hypothesis.
- Distinguish **predicted vs measured** every time. ProTox LD50 is a prediction;
  a ChEMBL IC50 may be an experimental measurement. Do not blur them.
- State the constraints when relevant: ProTox 250/day cap and CC BY-ND
  (non-commercial); OMIM needs a key and is not structure-based.

## Files

- `assets/protox3_client.py` — verified ProTox-3.0 client (SMILES → acute tox;
  best-effort full panel). No weights are bundled — none exist to bundle.
- `reference/platforms.md` — full, live-tested capability & access matrix for
  ProTox, OMIM, and the session's structure-based tools.
