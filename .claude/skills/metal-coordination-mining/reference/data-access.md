# Data access — building the structure set

`assets/fetch_structures.py` wraps steps 1–2; this file documents the endpoints,
their quirks, and what a run costs.

## 1. Define the sequence space (InterPro / Pfam / SSF)

Start from the domain that defines the fold, not from a single seed sequence — a
seed reintroduces the sequence bias the method exists to escape. The published
run started from the Pfam clan of the cupin superfold (clan CL0029) and kept the
member families that are plausibly enzymes.

Practical route: browse the clan or superfamily on InterPro, export the member
family table, and mark which families to include. Then:

```bash
python assets/fetch_structures.py accessions \
  --families PF13640,PF03171 --families-file more_families.txt \
  --outdir acc/ --min-length 150
```

Endpoint used (paginated, `next` cursor):

```
https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/InterPro/IPR005123/?page_size=200
                                                        .../PFAM/PF13640/?page_size=200
                                                        .../ssf/SSF51182/?page_size=200
```

Prefix routing is automatic: `IPR…` → InterPro, `PF…` → PFAM, `SSF…` → ssf.

Behaviour to expect:
- **204** means no content — treated as end of pages.
- **408/429/5xx** are transient; the script backs off and retries.
- Large families take a while at 200/page with a 0.5 s politeness delay; use
  `--max-pages` for a pilot before committing. Re-runs skip families already on
  disk unless `--force`.

Two filters are applied when the per-family tables are merged, both from the
published pipeline:
- keep only `in_alphafold == true` — no model, nothing to mine;
- keep `Length > 150 aa` — fragment/partial-sequence cutoff;
- deduplicate by accession (one protein can belong to several families; the
  `families` column keeps the full provenance).

## 2. Retrieve AlphaFold models

```bash
python assets/fetch_structures.py download \
  --accessions acc/accessions_af2_len_gt_150.csv --outdir structs/ --workers 8
```

Per accession: `GET https://alphafold.ebi.ac.uk/api/prediction/<ACC>` → take
`pdbUrl` → download. Going through the API rather than guessing
`AF-<ACC>-F1-model_v4.pdb` means the current model version is always used and
missing entries are detected cleanly.

- Existing non-empty files are skipped, so an interrupted run resumes.
- Failures are logged to `download_failures.txt` with a reason and never abort
  the batch.
- Keep `--workers` modest (8–12). AFDB throttles; more workers mostly buys 429s.
- Files must be named `<accession>.pdb` — `mcmine` uses the stem as the protein
  identifier throughout, and every downstream table joins on it.

**Budget**: an AF2 monomer PDB is ~0.2–1 MB. 500k models ≈ 200–400 GB and many
hours — the download, not the mining, is the bottleneck. Options: mine in
batches and delete each batch after its sites CSV is written (the CSV is the
durable artifact); restrict to fewer families; or pilot with `--limit`.
`mcmine` also reads `.pdb.gz` / `.cif` / `.cif.gz`, so keeping models compressed
is a cheap ~4× disk saving.

## 3. Annotate the hits (UniProt)

The mining step returns accessions; function/name text is what removes non-enzyme
members of the fold.

- ID mapping / retrieval: <https://www.uniprot.org/id-mapping> — upload
  `hit_accessions_*.txt`, retrieve FASTA or TSV.
- REST equivalent: `https://rest.uniprot.org/uniprotkb/search?query=accession:(A OR B …)&format=tsv&fields=accession,protein_name,organism_name,length,ec,lineage`.
- Published exclusion keywords (case-insensitive, on the protein name/function):
  `transcription`, `regulator`, `AraC`, `globin`, `AlkB`, `glutelin`, `TehB`,
  `tet`, `fragment`, `chemotaxis`, `helix-turn-helix`, `tellurite`, `adenosyl`,
  `SAM`, `ferredoxin`. Tune per fold, and record what the filter removed —
  keyword filtering is a blunt instrument and belongs in the ⚠️ column.

## 4. Alternatives to AFDB

- **ESM Metagenomic Atlas** — metagenomic space AFDB does not cover; same mining
  logic applies (download bulk, name files by accession).
- **PDB / mmCIF experimental structures** — use to sanity-check a motif against
  *holo* structures where the metal and ligands are actually present; `mcmine`
  parses mmCIF, ignores HETATM residues, and uses only the first model.
- **Locally predicted structures** (AF2/3, ESMFold, Boltz) — fine, as long as
  files are named by accession. Check pLDDT in the anchor region before trusting
  a hit.

## 5. Provenance to record for every run

Database releases move. Record: InterPro release, AFDB model version, the family
list, the exact motif spec (`mcmine` copies it into every `summary_*.json`), and
the date. Without those, a hit count is not reproducible.
