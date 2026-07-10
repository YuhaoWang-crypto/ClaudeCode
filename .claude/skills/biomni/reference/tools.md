# Biomni tool reference

218 tool APIs across 21 modules (from `A1(...).module2api`). Import any tool as
`from biomni.tool.<module> import <fn>`.

## Module breakdown (tool counts)

| Module | # | Module | # |
|---|---|---|---|
| database | 40 | pharmacology | 23 |
| genomics | 19 | molecular_biology | 18 |
| microbiology | 12 | physiology | 11 |
| immunology | 10 | bioimaging | 10 |
| genetics | 9 | literature | 8 |
| synthetic_biology | 8 | bioengineering | 7 |
| pathology | 7 | systems_biology | 7 |
| biochemistry | 6 | cancer_biology | 6 |
| cell_biology | 5 | biophysics | 3 |
| glycoengineering | 3 | lab_automation | 3 |
| support_tools | 3 | | |

## No-key tools (structured args → public API, no LLM)

Use these for demos, tests, and deterministic lookups.

**Literature** (`biomni.tool.literature`)
- `query_pubmed(query, max_papers=10)` — NCBI E-utilities
- `query_arxiv(query, max_papers=10)` — arXiv API
- `extract_pdf_content(url)`, `extract_url_content(url)`

**Database** (`biomni.tool.database`) — LLM-free subset:
- `query_alphafold(uniprot_id, endpoint='prediction', ...)` — EBI AlphaFold
- `query_pdb_identifiers(identifiers, return_type='entry', ...)` — RCSB PDB
- `blast_sequence(...)` — NCBI BLAST (slow; submits a job)
- `get_genes_near_ccre(...)`, `region_to_ccre_screen(...)` — SCREEN/ENCODE cCREs
- `get_hpo_names(hpo_terms, data_lake_path)` — needs local data-lake file

## Key-required tools (LLM in the loop)

- **The autonomous agent:** `agent.go("<free-text biomedical task>")`,
  `agent.configure(...)`.
- **Natural-language DB wrappers** — the `query_*(prompt="...")` form calls an
  LLM to build the API request. These include: `query_uniprot`, `query_pubchem`,
  `query_ensembl`, `query_kegg`, `query_reactome`, `query_stringdb`, `query_chembl`,
  `query_clinvar`, `query_gwas_catalog`, `query_opentarget`, `query_clinicaltrials`,
  `query_dbsnp`, `query_gnomad`, `query_cbioportal`, `query_geo`, `query_encode`,
  `query_uniprot`, and most other `query_<db>` functions in `database.py`.
  Without a key they raise *"Could not resolve authentication method."*

## Full database tool list (`biomni.tool.database`, 40)

blast_sequence, get_genes_near_ccre, get_hpo_names, query_alphafold,
query_cbioportal, query_chembl, query_clinicaltrials, query_clinvar,
query_dailymed, query_dbsnp, query_emdb, query_encode, query_ensembl,
query_geo, query_gnomad, query_gtopdb, query_gwas_catalog, query_interpro,
query_iucn, query_jaspar, query_kegg, query_monarch, query_mpd, query_openfda,
query_opentarget, query_paleobiology, query_pdb, query_pdb_identifiers,
query_pride, query_pubchem, query_quickgo, query_reactome, query_regulomedb,
query_remap, query_stringdb, query_synapse, query_ucsc, query_unichem,
query_uniprot, query_worms, region_to_ccre_screen.

## Overlap with this session's MCP servers

Several biomni tools duplicate MCP servers already connected in Claude Code
(ChEMBL, Clinical_Trials, PubMed, bioRxiv). For a quick lookup, prefer the MCP
tool (no key, no install). Use biomni when you want the **agent** to chain many
tools autonomously, or a database not covered by an MCP server (AlphaFold, PDB,
Ensembl, GWAS Catalog, ...).
