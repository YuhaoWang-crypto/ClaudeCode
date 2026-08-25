# The IEDB APIs: what works, and every quirk that costs time

Two entirely separate services. Neither is documented well enough to use
without the notes below.

## 1. Prediction — `tools_api` (NetMHCIIpan, BepiPred)

```
https://tools-cluster-interface.iedb.org/tools_api/mhcii/    # class II binding
https://tools-cluster-interface.iedb.org/tools_api/bcell/    # BepiPred-2.0
```

Form-encoded POST. No key, no local licensed binary required.

**Use `https`.** The http endpoint answers a POST with a **308 redirect**, and
`urllib` does not follow 308 for POST — the body is silently dropped and you
get a confusing parse error rather than a redirect error.

**Batch by allele block, not by sequence.** The endpoint takes parallel
comma-separated `allele` and `length` lists plus a *multi-record* FASTA, so one
call covers every sequence for a block of alleles. `seq_num` in the response is
the 1-based FASTA record index — map it back to your IDs yourself.

```python
_post(endpoint, {
    "method": "netmhciipan_el",          # or netmhciipan_ba
    "sequence_text": "".join(f">{i}\n{seqs[i]}\n" for i in ids),
    "allele": ",".join(alleles),
    "length": ",".join(["15"] * len(alleles)),
})
```

`alleles_per_request: 5` is a reasonable block size. Retry with backoff; the
cluster is shared and transient failures are normal.

### Wall-clock scales with the number of requests

Twenty separate small records took 89 s; one 60-aa record took 45 s under load.
Whatever the mechanism (it was not isolated), **cost tracks request count far
more than payload size**. For a 5,800-peptide benchmark, submitting one request
per peptide is not feasible.

**Solution: concatenated pseudo-proteins.** Join the peptides into one long
record separated by `(length-1)` glycines, submit in ~3,000-residue chunks, and
map result positions back to peptides via the layout you recorded when
building the chunk. The spacer length guarantees no scored 15-mer window can
span two peptides.

```python
SPACER_AA = "G"
CHUNK_AA  = 3000
```

**Always verify the trick before trusting a run built on it.**
`verify_concatenation()` scores a handful of peptides standalone and
concatenated and requires exact agreement. It is cheap and it is the difference
between an optimisation and a silent corruption.

### Make long runs resumable

Append each result to a partial TSV *as it lands*, and read that file back on
start to skip completed work. A 1–2 h run that writes its cache only at the end
loses everything to one timeout. `results/m11_scores_partial.tsv` and
`results/m13_panel_scores.tsv` are these caches — regenerable, so they are
gitignored, while the analysis outputs derived from them are tracked.

## 2. Labelled data — `query-api` (PostgREST)

```
https://query-api.iedb.org/tcell_search
```

Plain GET, PostgREST filter syntax. Two rules that are not in any documentation
and both fail with unhelpful errors:

**Do not percent-encode `*` in an `eq.` filter.** `HLA-DRB1*01:01` goes in raw.
Encoding it as `%2A` is rejected. The space in `structure_type=eq.Linear
peptide` *does* need encoding, so the query is a mix — which is exactly why it
looks like it should be uniformly encoded.

**`offset` requires `order`.** The API refuses an offset without an order, to
keep paging stable. Add `&order=tcell_id`.

```
?select=...&mhc_class=eq.II&host_organism_iri=eq.NCBITaxon:9606
&structure_type=eq.Linear%20peptide
&mhc_restriction=eq.HLA-DRB1*01:01
&order=tcell_id&limit=1000&offset=0
```

**Surface the HTTP response body on error.** This API explains what is wrong
with a query in the body of a 4xx; swallowing it costs hours. The paging rule
above was found exactly that way. Also: do not retry a 4xx — a malformed query
will not fix itself.

```python
except urllib.error.HTTPError as e:
    body = e.read().decode()[:400]
    last = f"{e} - {body}"
    if e.code < 500:
        break
```

## Fetching the human proteome

Once, ~20,400 reviewed entries:

```bash
curl -L 'https://rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta&compressed=true' \
  | gunzip > data/human_sprot.fasta
```

## Fetching sequences by motif, not by offset

`m0_fetch_sequences.py` pulls every sequence live from RCSB/UniProt and locates
each region **by searching for its motif**, never by a hard-coded residue
offset. Numbering conventions differ between PDB entries, UniProt isoforms and
the literature; an offset that is right today breaks silently on the next
release, and a control peptide that has quietly shifted by two residues is
close to undetectable downstream.
