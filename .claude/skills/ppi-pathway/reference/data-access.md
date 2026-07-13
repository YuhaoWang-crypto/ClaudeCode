# Data access reference — STRING & BioGRID

Exact download URLs, file formats, and identifiers. All URLs below are public
and need **no API key**. Human = taxon `9606`.

## STRING (v12.0)

Bulk downloads: `https://stringdb-downloads.org/download/<file>` (the
`string-db.org/cgi/download` page is the human-facing index of these).

| File | What it is | Human size (.gz) |
|---|---|---|
| `protein.links.v12.0/9606.protein.links.v12.0.txt.gz` | scored interactions: `protein1 protein2 combined_score` (space-separated, score 0–1000) | ≈ 83 MB |
| `protein.info.v12.0/9606.protein.info.v12.0.txt.gz` | `string_protein_id  preferred_name  size  annotation` — the gene-symbol map | ≈ 2 MB |
| `protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz` | every external id → STRING id (Ensembl, UniProt, gene symbol, …) | ≈ 20 MB |

Extra bulk files you can add if needed: `protein.links.detailed` (per-channel
subscores: neighborhood, fusion, coexpression, experiments, database,
textmining), `protein.physical.links` (physical subnetwork only),
`protein.enrichment.terms` (term → protein memberships for offline enrichment).

Node ids look like `9606.ENSP00000269305`; `network.py` relabels them to the
`preferred_name` gene symbol from `protein.info`.

### STRING REST API (online, no download)

Base: `https://string-db.org/api/<format>/<endpoint>` (format = `tsv`/`json`/…).
Wrapped in `string_api.py`. Key endpoints:

- `enrichment` — functional/pathway enrichment (FDR-corrected by STRING).
- `network`, `interaction_partners` — interactions / neighbours.
- `get_string_ids` — identifier resolution.
- `ppi_enrichment` — whole-set connectivity p-value.

Etiquette: send `caller_identity`, keep to a few requests/second, prefer POST
for long id lists. Practical ceiling ~2000 identifiers/call — beyond that use
the offline bulk files.

## BioGRID (Latest Release, tab3)

Per-organism archive (one zip, one file per organism inside):

```
https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip
```

`download.py` fetches the zip and extracts the member whose name contains the
organism (e.g. `Homo_sapiens`). Other archives on the same path if you need them:
`BIOGRID-ALL-LATEST.tab3.zip` (everything), `BIOGRID-MV-Physical-LATEST` (multi-
validated physical only).

### tab3 columns used (0-based)

| idx | column | used for |
|---|---|---|
| 7 / 8 | Official Symbol Interactor A / B | node names |
| 11 | Experimental System | edge `systems` |
| 12 | Experimental System Type | `physical` vs `genetic` filter |
| 17 | Throughput | (available; low/high) |
| 35 / 36 | Organism ID Interactor A / B | same-organism filter |

`load_biogrid()` keeps **physical, same-organism** interactions by default.

### BioGRID webservice (optional, needs a free key)

`https://webservice.thebiogrid.org/interactions?...&accesskey=KEY` — only needed
for programmatic slices without downloading files. Not used by this skill; get a
key at https://webservice.thebiogrid.org/ if you want it.

## Caching & storage

- Default cache dir: `./ppi_data` (git-ignored). Override with `PPI_DATA_DIR`.
- Downloads are size-checked against `Content-Length`; a complete file is
  skipped on re-run (`--force` to override).
- On an ephemeral cloud container the cache is wiped when the container is
  reclaimed — just re-run `download`. It is intentionally **not** committed to
  the repo (too large, and always re-fetchable).

## Licensing / citation

- **STRING**: CC BY 4.0. Cite Szklarczyk et al., *Nucleic Acids Research*.
- **BioGRID**: MIT license. Cite Oughtred et al., *Protein Science / NAR*.
