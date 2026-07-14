# Platform capability & access matrix (verified 2026-07-14)

Honest assessment of the two platforms the user asked about, plus the
structure-based toxicology tools already wired into this session. Every row was
tested live from this environment, not recalled from memory.

## TL;DR

| Platform | Download weights & run locally? | Programmatic API? | Reachable here? | Structure-based tox? |
|---|---|---|---|---|
| **ProTox-3.0** | ❌ never released (CC BY-ND, web/API only) | ✅ documented POST flow (250/day/IP) | ✅ HTTP 200 | ✅ yes — SMILES in |
| **OMIM geneMap** | ❌ N/A — it is a database, not an ML model | ⚠️ REST API but needs a personal, yearly-renewed key | ❌ HTTP 403 (Cloudflare) | ❌ no — gene→disease, not structure |

## ProTox-3.0 (tox.charite.de/protox3)

- **Weights download: impossible.** ProTox has never published weights or source.
  Licence is **CC BY-ND 4.0, academic / non-commercial**. So "download the weight
  from a cloud platform and run it" cannot be done for ProTox — there is nothing
  to download.
- **API: real and working.** FAQ documents "a simple POST interface", query by
  compound name (via PubChem) or canonical SMILES, CSV output, **250 queries /
  IP / day**, covering 61 endpoints (acute, organ, carcino/mutagen, Tox21
  pathways, nuclear receptors, MIEs, metabolism, targets).
- **Gotcha:** the advertised sample script `protox3_api.py` currently **404s** on
  their server. Our `assets/protox3_client.py` reproduces the real web-form flow
  instead (SMILES → RDKit MOL block → POST `smilesString` to
  `?site=compound_search_similarity`; model checkboxes are client-side display
  filters only).
- **Verified outputs (live):**
  - aspirin `CC(=O)OC1=CC=CC=C1C(=O)O` → **LD50 250 mg/kg, class 3** (matches
    ProTox's published value; 100% similarity = in training set).
  - novel `CCOC(=O)C1=CC=C(NC(=O)C2CC2N)C=C1` → **LD50 400 mg/kg, class 4**,
    similarity 72.6%, accuracy 69.26% (genuine extrapolation, not an echo).
- **Limitation:** the initial POST returns acute oral toxicity immediately; the
  other ~60 endpoints render "Not Calculated" and are filled by in-browser JS.
  `get_full_panel()` drives headless Playwright for those — treat as best-effort.

## OMIM geneMap (omim.org)

- **Category mismatch — call this out to the user.** OMIM is the **Mendelian
  gene ↔ phenotype (inherited-disease) catalog**. Its input is a gene / MIM
  number; it does **not** take a molecular structure and does **not** predict
  compound toxicity. It cannot replace ProTox for structure-based tox.
- **API exists but is gated:** REST (JSON/XML), but requires a **personal API
  key** (register at omim.org/downloads, **renew yearly**) — cannot be
  self-obtained here. The site sits behind Cloudflare; direct fetch returns
  **403** (the `__cf_chl_tk=…` token in the URL the user pasted is the human
  challenge).
- **Where it genuinely helps — an annotation layer, not an engine:** map a
  toxicity *target* (e.g. a ProTox Tox21 / nuclear-receptor hit, or a ChEMBL
  target gene) → OMIM disorders associated with that gene, to reason about
  mechanism / susceptibility. Wire it in only once a key is supplied.

## Already in this session — the stable structure-based backbone

Prefer these for anything programmatic and high-volume; no scraping, no daily
cap. Use ProTox as breadth-of-endpoints supplement, OMIM as gene-disease context.

- **ChEMBL** (`get_admet`, `get_bioactivity`, `get_mechanism`, `target_search`):
  calculated drug-likeness **and experimental** liabilities — hERG (CHEMBL240),
  CYP3A4/2D6/2C9, P-gp — the real assay data ProTox only predicts.
- **Inductive Bio** (`predict_properties`): structure→property models
  (logD, pKa, …).
- **Boltz** (`start_small_molecule_adme`): structure→ADME predictions.

## Recommended composition for a tox assessment of a compound

1. **Acute + endpoint breadth:** `protox3_client.predict_acute(smiles)` (LD50 +
   class, verified) → optionally `get_full_panel` for the 61-endpoint sweep.
2. **Experimental liabilities:** ChEMBL `get_bioactivity` on hERG / CYP targets;
   `get_admet` for drug-likeness.
3. **Physchem drivers:** Inductive Bio `predict_properties` (logD, pKa); Boltz
   ADME.
4. **Mechanistic / susceptibility context (optional, needs key):** map flagged
   targets → OMIM disorders.
5. **Label every claim** ✅ computed vs ⚠️ inferred, and name the source model.
