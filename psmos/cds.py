"""
cds.py — resolve REAL coding DNA sequences (CDS) for pathway orthologs.

Evo2 is a *genome* foundation model: it scores nucleotide sequences, not
proteins. So the fidelity/constraint layer must be computed on coding DNA, not
on the UniProt protein. This module maps each ortholog to its species' Ensembl
canonical transcript and pulls the CDS nucleotides.

Route (all real network calls to rest.ensembl.org, which serves vertebrates and
metazoa on one endpoint):
    symbol --xrefs/symbol--> gene id --lookup(expand)--> canonical transcript
           --sequence(type=cds)--> ACGT string

Species not on the main Ensembl REST endpoint (e.g. some non-model inverts) or
the negative controls simply resolve to nothing — which is the correct,
auditable outcome, not an error.
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict

ENSEMBL = "https://rest.ensembl.org"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "data", "cds")

# species_key -> Ensembl production name
ENSEMBL_SPECIES = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "chicken": "gallus_gallus",
    "frog": "xenopus_tropicalis",
    "zebrafish": "danio_rerio",
    "fly": "drosophila_melanogaster",
    "worm": "caenorhabditis_elegans",
    "tunicate": "ciona_intestinalis",
    "sea_anemone": "nematostella_vectensis",
    "yeast": "saccharomyces_cerevisiae",
    "plant": "arabidopsis_thaliana",
    "naked_mole_rat": "heterocephalus_glaber_female",
    # axolotl (Ambystoma mexicanum) and planaria (Schmidtea mediterranea) are
    # not on the main Ensembl REST endpoint -> honest CDS/Evo2 gaps.
}


@dataclass
class CDS:
    family: str
    species_key: str
    ensembl_species: str
    gene_id: str | None
    transcript_id: str | None
    nt_length: int
    sequence: str          # ACGT...
    found: bool
    note: str = ""

    def slim(self) -> dict:
        d = asdict(self)
        d["sequence"] = f"{self.sequence[:48]}...({self.nt_length} nt)" if self.sequence else ""
        return d


def _get(url: str, retries: int = 4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None      # symbol/gene not in this species — expected
            last = e
            time.sleep(2 ** i)
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"Ensembl request failed: {last}")


def resolve_cds(family: str, species_key: str, symbols: list[str]) -> CDS:
    esp = ENSEMBL_SPECIES.get(species_key)
    if not esp:
        return CDS(family, species_key, "", None, None, 0, "", False, "no Ensembl species map")

    gene_id = None
    for sym in symbols:
        if not sym:
            continue
        url = f"{ENSEMBL}/xrefs/symbol/{esp}/{urllib.parse.quote(sym)}?content-type=application/json"
        data = _get(url)
        if data:
            genes = [x["id"] for x in data if x.get("type") == "gene"]
            if genes:
                gene_id = genes[0]
                break
    if not gene_id:
        return CDS(family, species_key, esp, None, None, 0, "", False, "no Ensembl gene for symbols")

    look = _get(f"{ENSEMBL}/lookup/id/{gene_id}?expand=1&content-type=application/json")
    if not look:
        return CDS(family, species_key, esp, gene_id, None, 0, "", False, "lookup failed")
    tid = look.get("canonical_transcript")
    if not tid:
        # fall back to the longest transcript
        trs = look.get("Transcript", [])
        if not trs:
            return CDS(family, species_key, esp, gene_id, None, 0, "", False, "no transcripts")
        tid = max(trs, key=lambda t: t.get("length", 0))["id"]
    # Ensembl sometimes returns worm canonical ids with a trailing-dot artifact
    # (e.g. "R107.8.1."). Strip that first.
    tid = tid.rstrip(".")
    # Strip only the Ensembl version suffix (e.g. ENST....1, FBtr....2). Do NOT
    # split on every dot: WormBase transcript names contain structural dots
    # (lin-12 = R107.8.1), so a blanket split truncated them to "R107".
    if re.match(r"^(ENS|FB)", tid):
        tid = tid.split(".")[0]

    seq = _get(f"{ENSEMBL}/sequence/id/{tid}?type=cds&content-type=application/json")
    if not seq or "seq" not in seq:
        return CDS(family, species_key, esp, gene_id, tid, 0, "", False, "no CDS")
    s = seq["seq"].upper()
    return CDS(family, species_key, esp, gene_id, tid, len(s), s, True)


def fetch_pathway_cds(pathway, orthologs, use_cache: bool = True) -> list[CDS]:
    """For each found ortholog, resolve its CDS. Symbols tried in fallback order:
    the species-specific seed symbols (+ the UniProt-resolved gene name)."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache = os.path.join(CACHE_DIR, f"{pathway.key}_gates_cds.json")
    if use_cache and os.path.exists(cache):
        with open(cache) as fh:
            return [CDS(**c) for c in json.load(fh)]

    seed = pathway.ortholog_seed
    out: list[CDS] = []
    for o in orthologs:
        symbols = list(seed.get(o.family, {}).get(o.species_key, []))
        if o.gene and o.gene not in symbols:
            symbols.append(o.gene)
        c = resolve_cds(o.family, o.species_key, symbols)
        out.append(c)
        print(f"  [{o.family:16s}] {o.species_key:12s} -> "
              f"{c.transcript_id or '—':20s} {c.nt_length:6d} nt  "
              f"{'ok' if c.found else c.note}")
    with open(cache, "w") as fh:
        json.dump([asdict(c) for c in out], fh, indent=1)
    return out


if __name__ == "__main__":
    import sys
    from psmos.pathways import get_pathway
    from psmos.orthologs import fetch_pathway_gates
    pk = sys.argv[1] if len(sys.argv) > 1 else "Notch"
    pw = get_pathway(pk)
    orths = fetch_pathway_gates(pw, use_cache=True)
    cdss = fetch_pathway_cds(pw, orths, use_cache=False)
    print(f"\nresolved CDS for {sum(c.found for c in cdss)}/{len(cdss)} orthologs")
