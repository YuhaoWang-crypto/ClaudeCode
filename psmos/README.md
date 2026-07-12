# PSMOS — Pathway-aware model-organism selection (Evo2-live pilot)

**Question this answers:** *Given a biological signalling pathway, which model
organism is the right one to study it in — and for which purpose?* Not "which
species is most similar to human" (a false premise), but a **layered, provenance-
labelled** recommendation where the sequence-constraint layer is computed by a
genome foundation model (Evo2) rather than asserted.

This is the runnable core behind the Notch / Hippo dashboards. It closes the loop
the dashboards flagged as pending: `Evo2 状态:未接入` → **接入 (live on Modal)**.

## What is actually *computed* vs *curated*

| Layer | Source | Status |
|---|---|---|
| Hard gate — is the pathway even present? | UniProt + Ensembl ortholog search | ✅ **computed** (empirical; yeast/plant confirmed absent) |
| Ortholog CDS (coding DNA) | Ensembl REST canonical transcript | ✅ **computed** (real sequences) |
| Sequence constraint / "naturalness" | **Evo2-7B log-likelihood on Modal H100** | ✅ **computed** (the live layer) |
| Architecture similarity, redundancy (paralogues), tractability/throughput | Notch comparative-genomics priors | ⚠️ **curated** (labelled) |
| Cross-species regulatory grammar (R layer) | AlphaGenome (human/mouse) | ⛔ **not yet wired** (AlphaGenome covers human/mouse only) |

**Honesty boundary (enforced in code):** Evo2 log-likelihood measures sequence
constraint *across the tree of life*, **not** "equivalence to human". So it is
reported as an independent computed axis and its correlation with the curated
fidelity proxy is *reported*, never laundered into the composite silently.

## Pipeline

```
pathways.py         pathway core families + gate + curated priors + ortholog seed
  → orthologs.py    UniProt: real ortholog protein per (family, species)  [gate]
  → cds.py          Ensembl: canonical-transcript CDS (DNA) for each ortholog
  → evo2_modal.py   Modal app: evo2_7b on H100, score_sequences → mean log-LL
  → run_evo2_scoring.py   send CDS → Evo2, cache per-species constraint scores
  → scoring.py      PSMOS G/D/N/R/E/X + hard gate + goal-weighted ranking
  → build_dashboard.py    regenerate notch_dashboard_live.html from computed data
```

## Run

One-time setup (this environment):

```bash
pip install modal python-socks requests
export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt   # trust the agent-proxy CA
# MODAL_TOKEN_ID / MODAL_TOKEN_SECRET already in env
```

End to end (Notch):

```bash
python3 -m psmos.orthologs           # UniProt ortholog search  (gate)
python3 -m psmos.cds                 # Ensembl CDS
python3 -m psmos.run_evo2_scoring Notch   # Evo2 on Modal H100  (live)
python3 -m psmos.build_dashboard     # -> psmos/notch_dashboard_live.html
```

`orthologs`/`cds`/`evo2` results are cached under `psmos/data/` so the dashboard
rebuilds offline.

## Extending to a new pathway

Add a `Pathway` in `pathways.py`: its component families, which are the
gate families, and a per-species ortholog seed table for the gate genes. The
rest of the pipeline is pathway-agnostic. Hippo–YAP (the PSMOS regeneration
dashboard) is the natural next one — same machinery, regeneration-weighted goals.

## Notes / limits

- Evo2 scores DNA (nucleotide alphabet), so the constraint layer uses **CDS**,
  not the UniProt protein. Non-model invertebrates (some without an Ensembl
  gene) resolve to no CDS and are simply left out of the live layer — an
  auditable gap, not an error.
- "Absence of a UniProt hit" ≠ "gene absent". The gate only hard-fails a
  species when the ortholog search finds **nothing** *and* the curated prior
  agrees (the true negative controls); partial misses fall back to curated and
  are flagged as annotation gaps.
