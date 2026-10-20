# Data access recipes

Everything below was exercised end-to-end from a sandboxed container. Where an
endpoint was blocked, the working fallback is given rather than the ideal one.

---

## Foldseek web server — a whole-AFDB sweep with no downloads

Base: `https://search.foldseek.com/api`

```bash
# submit
curl -s -X POST https://search.foldseek.com/api/ticket \
  -F q=@query.pdb -F mode=3diaa \
  -F 'database[]=afdb-proteome' -F 'database[]=afdb-swissprot'
# -> {"id":"...","status":"PENDING"}

# poll until COMPLETE
curl -s https://search.foldseek.com/api/ticket/<ID>

# fetch
curl -s https://search.foldseek.com/api/result/<ID>/0 -o hits.json
```

`scripts/fsweb.py` wraps all three.

### Databases (from `GET /api/databases`)

| path | what it is | taxonomy | notes |
|---|---|---|---|
| `afdb50` | whole AFDB clustered at 50% id | yes | broadest; `--max-accept 1000` |
| `afdb-proteome` | AFDB reference proteomes | yes | **best taxonomic spread per hit slot** |
| `afdb-swissprot` | AFDB over Swiss-Prot | yes | **best annotations** |
| `pdb100` | experimental structures | yes | for finding solved relatives |
| `cath50` | CATH domains | yes | fold-level context |
| `bfvd_v3`, `BFVD_folddisco` | viral protein structures | yes | if viruses are in scope |
| `mgnify_esm30`, `gmgcl_id` | metagenome-predicted | no | environmental breadth, no taxonomy |

### Fields returned per hit

`target, taxId, taxName, seqId, eval, prob, score, alnLength, qStartPos,
qEndPos, dbStartPos, dbEndPos, qLen, dbLen, qAln, dbAln, tSeq, tCa`

`tCa` is a flat comma-separated `x,y,z` list of length `3 * dbLen`, in register
with `tSeq`. **This is the useful part**: it means every hit arrives with its
full Calpha trace, so the chemistry filter and the catalytic-shell geometry
check run with no structure downloads at all.

Not returned: TM-score, lddt, per-residue pLDDT. `fsweb.tm_score()` recomputes
TM from the coordinates; pLDDT is genuinely unavailable in web mode and must
be reported missing.

### The cap, and what it does to your conclusions

`afdb50` accepts at most 1000 hits and fills them with the query's **nearest**
relatives. Measured: human cGAS vs `afdb50` returned 1000 hits spanning 382
species, overwhelmingly molluscs and cnidarians (102 from *Dreissena
polymorpha* alone) — and reached no plant or fungal protein at all. The same
query against `afdb-proteome` + `afdb-swissprot` returned 190 hits that
included bacterial CD-NTases.

**A taxon missing from a capped search is not a taxon lacking the fold.**
Say so in writing.

Note also that `afdb50` is clustered: a specific accession may be absent
because a cluster representative stands in for it. Three proteins verified
present in `afdb-proteome` (O64642, A0A0P0VE25, A0A0P0XSC2) were all absent
from `afdb50` for this reason.

---

## AlphaFold DB

### Metadata API — works, and carries the sequence

```bash
curl -s https://alphafold.ebi.ac.uk/api/prediction/O64642
```

Gives `uniprotAccession`, `organismScientificName`, `taxId`, `sequence`,
`globalMetricValue` (mean pLDDT), `latestVersion`, and `pdbUrl`/`cifUrl`/`bcifUrl`.

### Per-file downloads — frequently blocked

`https://alphafold.ebi.ac.uk/files/AF-<ACC>-F1-model_v6.pdb` returned **403**
for every accession in this environment, including URLs taken verbatim from
the API. The body is a Google frontend error page, not a proxy denial.

Fallbacks, in order of preference:

1. **Foldseek `tCa`** (above) — no download needed at all.
2. **EBI FTP proteome archives** — a different host, which worked:
   ```bash
   https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/
   # e.g. UP000002311_559292_YEAST_v6.tar
   ```
   48 model-organism reference proteomes are available. Each tar holds both
   `.cif.gz` and `.pdb.gz` per entry; extract one format:
   ```bash
   tar xf YEAST.tar -C dir --wildcards '*.pdb.gz'
   ```
3. **RCSB** for anything experimentally solved (`files.rcsb.org/download/`).
4. **3D-Beacons** (`www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/uniprot/summary/<ACC>.json`)
   enumerates every model provider for an accession — useful when AFDB itself
   is unreachable but another provider is not.

### Long transfers get cut — always resume

Measured throughput: ~0.9 MB/s single-stream, ~2.5-3 MB/s across four parallel
streams, with `curl: (56) Recv failure: Connection reset by peer` on several.
Never trust one attempt:

```bash
want=$(curl -sI "$url" | grep -i '^content-length' | awk '{print $2}' | tr -d '\r')
for i in $(seq 1 40); do
  have=$(stat -c %s "$out" 2>/dev/null || echo 0)
  [ "$have" -ge "$want" ] && break
  curl -sSL -C - --max-time 1800 -o "$out" "$url"
done
```

Verify against `Content-Length`; a truncated tar fails late and confusingly.

Budget disk: a proteome tar is 0.2-5 GB (METJA 173 MB, YEAST 976 MB, ARATH
3.7 GB, HUMAN 4.9 GB). Process one at a time — download, extract, search,
delete — to keep peak usage bounded.

---

## Local Foldseek

```bash
curl -sSL -o foldseek.tar.gz https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xzf foldseek.tar.gz           # -sse2 build if AVX2 is unavailable
./foldseek/bin/foldseek easy-search query.pdb targets_dir hits.tsv tmp \
  --format-output query,target,qstart,qend,tstart,tend,qaln,taln,evalue,prob,alntmscore,qtmscore,ttmscore,lddt,fident,alnlen,qcov,tcov,qlen,tlen
```

Gzipped inputs are read directly. Whole archaeal proteome (1,773 models)
searched in **12 s** on 8 threads — local mode is cheap once the data is down,
and it is the only mode that gives side chains and pLDDT.

---

## UniProt

```bash
# one entry, selected fields
curl -s 'https://rest.uniprot.org/uniprotkb/O64642.json?fields=accession,protein_name,gene_names,cc_function,xref_pfam'

# sequence
curl -s https://rest.uniprot.org/uniprotkb/O64642.fasta

# by domain across a taxon
curl -s 'https://rest.uniprot.org/uniprotkb/search?query=xref:pfam-PF01909+AND+taxonomy_id:4751&size=500&fields=accession,organism_name,protein_name'
```

### Taxonomy lineage, in batches

```bash
curl -s 'https://rest.uniprot.org/taxonomy/search?query=(tax_id:3702)OR(tax_id:4896)&size=500&fields=id,scientific_name,lineage'
```

Each result's `lineage` is an ordered list with `rank` annotations; match
against kingdom-level markers (`Viridiplantae`, `Metazoa`, `Fungi`,
`Discoba`, `Alveolata`, `Amoebozoa`, `Bacteria`, `Archaea`). Batch ~40 ids per
request and cache to disk — `scripts/annotate.py` does both.

Relevant Pfam families for the polymerase-beta-like nucleotidyltransferases:
**PF01909** (`NTP_transf_2`, the core domain — 24,600+ proteins),
**PF03828** (`PAP_assoc`, the associated domain), **PF18765**.

---

## RCSB

```bash
# full-text search
curl -s https://search.rcsb.org/rcsbsearch/v2/query -H 'Content-Type: application/json' \
  -d '{"query":{"type":"terminal","service":"full_text",
       "parameters":{"value":"terminal uridylyltransferase URT1 Arabidopsis"}},
       "return_type":"entry",
       "request_options":{"paginate":{"start":0,"rows":10}}}'

curl -sSL -o 6L8K.pdb https://files.rcsb.org/download/6L8K.pdb
```

**Always grep the header for engineered mutations before trusting a catalytic
residue:**

```bash
grep -E '^SEQADV|MUTATION' 6L8K.pdb
# SEQADV 6L8K ALA A 547 UNP O64642 ASP 547 ENGINEERED MUTATION
```

Both AtURT1 structures are D547A inactive mutants. Finding this *after*
building a pipeline on them wastes a day.
