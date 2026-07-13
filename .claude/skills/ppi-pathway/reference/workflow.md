# End-to-end omics → pathway/network workflow

A recipe for the common case: you ran an omics experiment (RNA-seq, proteomics,
a CRISPR/variant screen), you have a **list of hit genes**, and you want to know
*what pathways they hit* and *how they wire together*.

## 0. Prepare the gene list

One HGNC gene symbol per line in `genes.txt` (blank lines and `#` comments are
ignored). If you have Ensembl/UniProt ids, resolve them first:

```bash
python -m ppi_pathway.string_api  # or:
python - <<'PY'
from ppi_pathway import string_api
for r in string_api.map_ids(["ENSG00000141510","P04637","7157"]):
    print(r["queryItem"], "->", r["preferredName"])
PY
```

Unrecognised symbols simply drop out of enrichment; check with `map_ids` if a
gene you expect is missing.

## 1. Is the set a coherent module? (sanity check first)

```bash
python -m ppi_pathway.cli enrich genes.txt --fdr 0.05
```

Read the `PPI enrichment p-value` line on stderr. If it's small (< 0.05) the
genes interact more than random — enrichment results are meaningful. If it's
large, the list may be noisy; tighten your hit threshold or interpret loosely.

## 2. Pathway enrichment

```bash
# everything (GO, KEGG, Reactome, WikiPathways, COMPARTMENTS, DISEASES, ...)
python -m ppi_pathway.cli enrich   genes.txt --fdr 0.05 --json > enrich.json

# pathway databases only, top 20 by FDR — usually what you report
python -m ppi_pathway.cli pathways genes.txt --n 20 > pathways.tsv
```

Report the pathway rows with lowest FDR; group by `category`. The `PMID` and
`DISEASES` categories are literature/disease context — supporting, not primary.

**Reproducible/offline variant** — bring your own gene sets (download a `.gmt`
from MSigDB: Hallmark, KEGG, Reactome, GO):

```python
from ppi_pathway import enrich
sets = enrich.read_gmt("h.all.v2023.2.Hs.symbols.gmt")
rows = enrich.enrich(hits, sets, background=all_measured_genes, fdr_max=0.05)
```

Passing your real `background` (all genes your assay could detect) makes the
hypergeometric test honest — don't skip it for a publication-grade result.

## 3. Build the network module

```bash
python -m ppi_pathway.cli download --source all           # once per container
python -m ppi_pathway.cli subnet genes.txt \
       --min-score 700 --expand 1 --out module --top 30
```

- `--min-score 700` keeps high-confidence edges only.
- `--expand 1` adds one shell of connector genes — proteins that bridge your
  hits even if they weren't in the list (often the mechanistic link).
- Outputs `module.edges.tsv` + `module.nodes.tsv` → open in Cytoscape/Gephi.
- The printed table ranks genes by degree/betweenness: **hubs** (high degree)
  are the module's core; **bottlenecks** (high betweenness) bridge sub-modules
  and make good intervention/drug-target candidates.

## 4. (optional) Curated-only / consensus network

```python
from ppi_pathway import network
G = network.load_string(min_score=700)        # broad, predicted
B = network.load_biogrid()                    # curated, physical
C = network.consensus_graph(G, B)             # edges tagged in_string / in_biogrid
both = [(u, v) for u, v, d in C.edges(data=True) if d["in_string"] and d["in_biogrid"]]
```

Edges present in **both** resources are the most defensible; use them when you
need a conservative, high-confidence backbone.

## 5. Deeper dynamical analysis (optional, sibling skill)

Once you have an irreducible core module, the **`network-biomarker`** skill in
this repo can take a small pathway/network and compute symmetry cores, bistable-
switch capacity (CRNT deficiency), elementary flux modes, and critical-slowing
early-warning biomarkers. Natural next step after you've isolated a module here.

## Cheat-sheet

| Goal | Command |
|---|---|
| Cache data (STRING+BioGRID+humanPPI) | `python -m ppi_pathway.cli download --source all` |
| Coherence + full enrichment | `python -m ppi_pathway.cli enrich genes.txt` |
| Pathways only | `python -m ppi_pathway.cli pathways genes.txt --n 20` |
| Reactome pathways (→STRING fallback) | `python -m ppi_pathway.cli reactome genes.txt --n 20` |
| InterPro domains/families/GO | `python -m ppi_pathway.cli annotate accs.txt` |
| Network module + hubs | `python -m ppi_pathway.cli subnet genes.txt --min-score 700 --expand 1 --out module` |
| End-to-end demo test | `python -m ppi_pathway.demo` |
| Other species | add `--taxon 10090` (mouse), etc. |

## Multi-source consensus (highest confidence)

After caching all three sources, overlay them and keep edges supported by 2+
(ideally all 3) resources — predicted (STRING) + curated (BioGRID) + structural
(humanPPI):

```python
from ppi_pathway import network, humanppi
ov = network.overlay(string=network.load_string(min_score=700),
                     biogrid=network.load_biogrid(),
                     humanppi=humanppi.load(80))
gold = [(u, v) for u, v, d in ov.edges(data=True) if d["support"] == 3]
```

In a validated run this yields ~8,000 edges backed by all three sources — a
conservative, high-trust human interactome backbone.

## Annotate the module's hubs (mechanism)

Take the top hubs from step 3, resolve to UniProt, and pull InterPro domains:

```python
from ppi_pathway import interpro
for e in interpro.entries_for_protein("P04637"):     # TP53
    print(e["accession"], e["type"], e["name"])
```

`humanppi.load()` conveniently stores each gene's UniProt accession on the edge
(`up1`/`up2`), so you can map hub symbols → accessions without another service.
