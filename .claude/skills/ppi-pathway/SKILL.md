---
name: ppi-pathway
description: >-
  Download protein–protein interaction data from STRING, BioGRID and the Cong-lab
  humanPPI structural predictions, and run omics / pathway analysis on it —
  functional & pathway enrichment (GO, KEGG, Reactome, WikiPathways) of a gene
  list, extraction of the PPI subnetwork / module behind an omics hit list, hub
  & bottleneck ranking, multi-source consensus networks, and InterPro
  domain/family/GO annotation of proteins. Use when a user has a gene list
  (differentially expressed / mutated / hit genes) and wants pathway enrichment,
  network / module analysis, protein-domain annotation, or wants to cache STRING
  / BioGRID / humanPPI data on the cloud for reuse. Online paths (STRING &
  Reactome REST, InterPro) need only `requests`; offline genome-wide path caches
  the bulk files once and uses `networkx`. A built-in demo test
  (`python -m ppi_pathway.demo`) exercises every source end-to-end.
---

# STRING + BioGRID + humanPPI PPI & pathway analysis

Turn a gene list into biology: **pathway enrichment**, a **PPI network module**,
and **protein-domain annotation**, backed by three interaction resources plus
InterPro and Reactome.

- **STRING** — https://string-db.org — scored (probabilistic) functional
  associations, plus a REST API that does functional/pathway enrichment for you.
- **BioGRID** — https://thebiogrid.org — manually curated *physical* and genetic
  interactions with experiment-level evidence.
- **humanPPI** — http://prodata.swmed.edu/humanPPI/ — deep-learning
  (RoseTTAFold2-PPI / AlphaFold) structural predictions of the human interactome
  (Zhang, Humphreys, Pei, … Baker, Cong, *Science* 2025): ~29k confident pairs
  with structural confidence scores. A high-confidence, structure-backed layer.
- **InterPro** — https://www.ebi.ac.uk/interpro/ — protein families/domains/sites
  + GO, to annotate the hubs of your module ("what *is* this protein").
- **Reactome** — https://reactome.org — curated pathway enrichment (note: its
  REST API is Cloudflare-blocked from datacenters; the skill auto-falls back to
  STRING's Reactome-category enrichment — see below).

The code lives in `scripts/ppi_pathway/` and is a small, dependency-light Python
package with a CLI. **The raw data files are never committed** — they are large
(STRING human ≈ 83 MB compressed) and are cached under a git-ignored `ppi_data/`
dir (override with the `PPI_DATA_DIR` env var).

## Decide: online or offline?

| Your situation | Use | Needs |
|---|---|---|
| A gene list (≤ ~2000) → pathway/GO enrichment | **online** STRING API | `requests` only |
| A gene list → small network around it | **online** STRING API | `requests` only |
| Genome-wide work, custom score cutoffs, thousands of genes, no per-call limits, reproducible offline runs | **offline** bulk files | `networkx` + one-time download |
| Physical-interaction (curated) network, experiment evidence | **offline** BioGRID | `networkx` + download |
| High-confidence structure-backed interactions | **offline** humanPPI | `networkx` + download (~14 MB) |
| Annotate a protein: domains, families, GO | **online** InterPro | `requests` only |
| Reactome pathway enrichment | **online** Reactome (→ STRING fallback) | `requests` only |
| Your own gene sets (MSigDB/KEGG `.gmt`), deterministic stats | **offline** `enrich.py` | nothing extra |

Start online — it's instant and needs no download. Drop to offline when the
gene list is large or you need the whole graph.

## Quickstart (from `scripts/`)

```bash
cd .claude/skills/ppi-pathway/scripts
printf 'TP53\nBRCA1\nEGFR\nMYC\nCDK2\nCDK4\nRB1\nMDM2\nATM\nCHEK2\n' > genes.txt

# ONLINE — pathway enrichment (KEGG/Reactome/WikiPathways/GO), FDR-corrected
python -m ppi_pathway.cli enrich   genes.txt --fdr 0.05
python -m ppi_pathway.cli pathways genes.txt --n 20      # pathway categories only

# OFFLINE — cache the bulk data once, then build the module
python -m ppi_pathway.cli download --source all           # STRING + BioGRID + humanPPI
python -m ppi_pathway.cli subnet   genes.txt \
        --min-score 700 --expand 1 --out module            # + hub ranking + Cytoscape files

# ANNOTATE proteins (InterPro domains/families/GO) — accepts UniProt accessions
printf 'P04637\nP00533\n' > accs.txt
python -m ppi_pathway.cli annotate accs.txt

# REACTOME pathway enrichment (auto-falls back to STRING if Cloudflare-blocked)
python -m ppi_pathway.cli reactome genes.txt --n 20
```

## Run the demo test

```bash
cd .claude/skills/ppi-pathway/scripts
python -m ppi_pathway.demo               # exercises STRING, humanPPI, InterPro, Reactome
python -m ppi_pathway.demo --online-only # skip the humanPPI download
python -m ppi_pathway.demo --with-string # also test the STRING∪BioGRID∪humanPPI overlay
```

It runs a fixed DNA-damage/cell-cycle gene set through every source and prints
PASS/FAIL checks — the fastest way to confirm the whole pipeline works.

`enrich` also prints a **PPI-enrichment p-value**: whether the gene set is more
interconnected than random — a quick "is this a real coherent module?" check.

## Library use (inside a notebook / larger analysis)

```python
import sys; sys.path.insert(0, ".claude/skills/ppi-pathway/scripts")
from ppi_pathway import string_api, download, network, subnetwork, enrich

genes = ["TP53", "BRCA1", "EGFR", "MYC", "CDK2"]

# --- online ---
paths  = string_api.top_pathways(genes)              # enriched KEGG/Reactome/WikiPathways
coh    = string_api.ppi_enrichment(genes)            # coherence p-value
net    = string_api.network(genes, required_score=700, add_nodes=10)

# --- offline (after download.string_files() / .biogrid_file()) ---
G   = network.load_string(min_score=700)             # networkx graph, gene symbols
sub = subnetwork.induced(G, genes, expand=1)         # module + 1 neighbour shell
subnetwork.write_cytoscape(sub, "module")            # module.edges.tsv / module.nodes.tsv
hubs = subnetwork.rank_genes(sub)                    # degree/betweenness table

# multi-source consensus: STRING ∪ BioGRID ∪ humanPPI
from ppi_pathway import humanppi, interpro
B   = network.load_biogrid()
H   = humanppi.load(precision=80)                    # after humanppi.download()
ov  = network.overlay(string=G, biogrid=B, humanppi=H)
best = [(u, v) for u, v, d in ov.edges(data=True) if d["support"] == 3]  # all 3 agree

# annotate a hub with InterPro domains/families/GO
doms = interpro.entries_for_protein("P04637")        # TP53

# offline over-representation against your own .gmt gene sets
sets = enrich.read_gmt("h.all.v2023.2.Hs.symbols.gmt")
res  = enrich.enrich(genes, sets, background=my_assay_universe)
```

## The modules

| File | Role |
|---|---|
| `download.py`   | Cache STRING (`links`/`info`/`aliases`) + BioGRID tab3 + humanPPI. Idempotent, size-checked. |
| `string_api.py` | ONLINE STRING REST: `enrichment`, `top_pathways`, `network`, `interaction_partners`, `map_ids`, `ppi_enrichment`. |
| `network.py`    | OFFLINE loaders → `networkx` graphs (gene-symbol nodes); `consensus_graph`, `overlay` (multi-source). |
| `humanppi.py`   | Download + load the Cong-lab structural PPI predictions into a `networkx` graph. |
| `interpro.py`   | ONLINE InterPro REST: `entries_for_protein`, `domains`, `go_terms`, `annotate`. |
| `reactome.py`   | ONLINE Reactome AnalysisService: `top_pathways` (raises `ReactomeBlocked` behind Cloudflare). |
| `subnetwork.py` | `induced` subnetwork, `rank_genes` (degree/betweenness/clustering), `largest_component`, `write_cytoscape`. |
| `enrich.py`     | OFFLINE hypergeometric over-representation + BH-FDR against a `.gmt`. |
| `cli.py`        | `download` / `enrich` / `pathways` / `reactome` / `annotate` / `subnet` subcommands. |
| `demo.py`       | End-to-end integration test across all sources (PASS/FAIL report). |

## Interpreting results (short guide)

- **Enrichment table**: sort by `fdr`. `category` tells you the source (KEGG,
  Reactome Pathways, WikiPathways, Process/GO, COMPARTMENTS, ...). Treat `PMID`
  and `DISEASES` categories as supporting context, not clean pathways.
- **PPI-enrichment p-value** (`ppi_enrichment`): small (< 0.05) ⇒ the genes
  interact more than chance ⇒ a real functional module, so pathway hits are
  trustworthy. Large ⇒ a grab-bag; interpret enrichment cautiously.
- **Hubs vs bottlenecks** (`rank_genes`): high **degree** = hub (many partners,
  often core/essential); high **betweenness** = bottleneck (bridges modules,
  candidate control point / drug target).
- **STRING score** is a confidence, not an affinity: 150 low · 400 medium ·
  700 high · 900 highest. Use ≥ 700 for a clean module.
- **STRING vs BioGRID vs humanPPI**: STRING is broad & predicted (text-mining +
  co-expression); BioGRID is experimentally curated (physical assays); humanPPI
  is structure-based deep-learning prediction. They are complementary evidence
  types — an edge with `support == 3` in `overlay()` is backed by all three and
  is the most defensible.
- **Reactome & Cloudflare**: `reactome.py` works from a normal network but
  raises `ReactomeBlocked` from datacenter IPs (Cloudflare challenge). The CLI
  `reactome` command and the demo auto-fall back to STRING's `RCTM` (Reactome)
  enrichment category, which is *not* blocked. So you always get Reactome
  pathways, just via a different route in the cloud.
- **humanPPI "model weights"**: what the download provides is the *prediction
  table* (~29k scored pairs) — the reusable product for network analysis. The
  underlying RoseTTAFold2-PPI neural-net weights are a separate, multi-GB
  research artifact and are **not** needed (or downloaded) for pathway/network
  work here.

## Notes, licensing, other species

- **No API key needed** for anything here (STRING API + STRING/BioGRID bulk
  *files* are public). Only the BioGRID *webservice* needs a free key — not used.
- **Species**: everything takes `--taxon` / `taxon=`; default `9606` (human).
  Mouse `10090`, rat `10116`, fly `7227`, worm `6239`, yeast `4932` are pre-mapped
  for BioGRID; STRING works for any taxon it publishes.
- **Licensing** (respect for any redistribution): STRING is CC BY 4.0; BioGRID
  is MIT. Cite STRING (Szklarczyk et al.) and BioGRID (Oughtred et al.).
- See `reference/data-access.md` for exact URLs/formats and
  `reference/workflow.md` for the end-to-end omics recipe.
