---
name: psmos
description: >-
  Pathway-aware model-organism selection. Given a biological signalling pathway,
  recommend the right model organism(s) — per research purpose (mechanism,
  imaging, genetics, translation, comparative) — using a layered,
  provenance-labelled score whose sequence layer is computed by Evo2 (genome
  foundation model on Modal H100), regulatory layer by AlphaGenome (human/mouse),
  redundancy by Ensembl Compara, and hard gate by an empirical UniProt/Ensembl
  ortholog search. Use when asked which organism to study a pathway in, to add a
  new pathway to the psmos package, to run Evo2 sequence-constraint or
  AlphaGenome regulatory scoring, or to regenerate a model-organism dashboard.
  Enforces computed-vs-curated provenance on every axis and "absence of evidence
  ≠ evidence of absence" on the gate.
---

# PSMOS — pathway-aware model-organism selection

Answers: *given a pathway P, which model organism is right to study it in, and
for what purpose?* Not "which species is most similar to human" (a false
premise) — a **layered, provenance-labelled** recommendation where each axis is
either **computed** (and says by what) or **curated** (and says so). The working
package is `psmos/`.

## The core thesis (what the output proves)

There is **no single best model organism** for a pathway. The right choice
depends on the question, so the output is five role models, not one winner:

| Role | Wants | Notch pick | Hippo pick |
|---|---|---|---|
| **Mechanistic** | pathway complete + low redundancy + tractable | fly | fly |
| **Imaging** | live dynamics, transparent | zebrafish | zebrafish |
| **Genetic** | high-throughput perturbation | fly / worm | fly |
| **Translational** | mammalian, human-like architecture | mouse | mouse |
| **Comparative / negative** | natural rewiring or loss → tests necessity | yeast/plant (loss) | planaria (rewiring) |

The last role is the point: the best organism is often *not* the most similar.

## The layers, and what is computed vs curated

| Layer | Source | Provenance tag |
|---|---|---|
| Hard gate — is the pathway present? | UniProt + Ensembl ortholog search | `computed:uniprot` |
| Ortholog coding DNA (CDS) | Ensembl canonical transcript | `computed:ensembl` |
| **G** sequence constraint / naturalness | **Evo2-7B log-likelihood** (Modal H100) | `computed:evo2` |
| Redundancy / low-copy | **Ensembl Compara** within-species paralogues | `computed:compara` |
| **R** regulatory grammar (human/mouse only) | **AlphaGenome** RNA-seq/ATAC/TF-ChIP | `computed:alphagenome` |
| D domain, N network, E expression, X tractability | comparative-genomics priors | `curated` |

**Complementarity, stated honestly:** Evo2 scores sequence naturalness across
*all* species but cannot judge regulatory equivalence; AlphaGenome predicts
regulatory tracks but only for *human and mouse*. They are complementary, not
redundant. Say this; do not let Evo2 stand in for R or AlphaGenome for G.

## Run

```bash
pip install modal python-socks requests            # Evo2 path
pip install alphagenome --ignore-installed packaging  # R path (separate; see gotcha)
export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt      # trust the agent-proxy CA
export GRPC_DEFAULT_SSL_ROOTS_FILE_PATH=/root/.ccr/ca-bundle.crt  # AlphaGenome gRPC
# MODAL_TOKEN_ID/SECRET and ALPHAGENOME_API_KEY come from the environment

python3 -m psmos.run_all Notch      # UniProt→Ensembl→Compara→AlphaGenome→Evo2→dashboard
python3 -m psmos.run_all Hippo       # the PSMOS six-layer regeneration dashboard
```

Each stage caches under `psmos/data/` so the dashboard rebuilds offline. Single
stages: `python3 -m psmos.{orthologs,cds,paralogs,alphagenome_r,run_evo2_scoring} <Pathway>`.

## The non-negotiable discipline: provenance + honest gate

1. **Every axis is labelled** `computed:<source>` or `curated`. Never silently
   present a curated prior as a measurement. The dashboards carry the legend and
   a `*` on computed cells.
2. **Absence of a hit ≠ gene absent.** The hard gate only *disqualifies* a
   species when the ortholog search finds **nothing** *and* the curated prior
   agrees (the true negative controls — yeast/plant have no Notch). Partial
   misses are annotation gaps → fall back to curated, flagged
   `computed-partial (annotation gap → curated)`. Getting this wrong once
   falsely killed sea anemone (real Notch, missing UniProt CSL annotation).
3. **Report, don't launder.** Evo2 constraint is shown as its own axis and its
   correlation with the curated fidelity proxy is *reported* (Notch: r=+0.73),
   not blended in silently. Where a computed value IS folded into a layer (Evo2
   into Hippo's G; AlphaGenome into mouse's R), the blend is explicit and the
   provenance says so.
4. **Overturned priors are findings.** AlphaGenome measured human↔mouse Hippo
   regulatory conservation at 0.614, below the curated 0.90 — so mouse dropped
   in the ranking. Say it moved; don't hide it.

## Common tasks → where to look

- **Add a new pathway** (e.g. Wnt/β-catenin, cGAS–STING, insulin–mTOR–FOXO) →
  `reference/adding-a-pathway.md` + `assets/pathway_template.py`. One `Pathway`
  entry (families, gate families, per-species ortholog seed) reuses the whole
  pipeline.
- **Infra recipes that actually work behind the proxy** (Modal auth, the Evo2
  image, AlphaGenome gRPC, Ensembl/Compara/liftover) →
  `reference/data-access.md`. These encode several hours of dependency-hell.
- **The layered framework + scoring math + provenance model** →
  `reference/methodology.md`.

## Hard-won gotchas (these were real bugs / blockers)

- **Modal behind the proxy** needs `pip install python-socks` (Modal's grpclib
  proxy support) + `SSL_CERT_FILE=/root/.ccr/ca-bundle.crt`. Without python-socks
  you get `ConnectionError: Could not connect to the Modal server`.
- **Evo2 image**: use base `nvcr.io/nvidia/pytorch:25.04-py3` (Arc's validated
  base). 25.01 gives a `fwd(): incompatible function arguments` at scoring
  (flash-attn/vortex kernel mismatch). Then force `torch.load(weights_only=False)`
  — the checkpoint pickles a `transformer_engine` global that PyTorch ≥2.6
  rejects by default. evo2's real deps are light (`vtx` + torch); NGC provides
  flash-attn/TE prebuilt, so the build is CPU-only.
- **AlphaGenome is API-only** (no open weights → cannot self-host on Modal).
  Needs `ALPHAGENOME_API_KEY`. Its gRPC tunnels through the proxy via
  `GRPC_DEFAULT_SSL_ROOTS_FILE_PATH` + `HTTPS_PROXY`.
- **Mouse assembly**: AlphaGenome mouse = GRCm38/mm10, but Ensembl current is
  GRCm39. Lift mouse coords over with Ensembl `/map/mouse/GRCm39/<region>/GRCm38`.
  Human is GRCh38 both sides (no liftover).
- **protobuf conflict**: `alphagenome` bumps protobuf to 7.x; `modal` needs
  `<7.0`. They don't run in the same step (Evo2 scores are cached), so it's fine
  — but don't try to import both in one process; use separate venvs if you must.
- **Evo2 scores DNA, not protein** — feed CDS (nucleotides), never the UniProt
  protein sequence.
- **WormBase transcript IDs contain real dots** (`R107.8.1`) and a trailing-dot
  artifact from Ensembl; strip only the Ensembl `ENS/FB` version suffix, and
  `rstrip('.')` the artifact.
- **Isoform selection**: pick the ortholog whose length is closest to the human
  reference, not the longest — "longest" grabbed a 1673-aa fusion for fly Su(H)
  instead of the canonical 594-aa.
