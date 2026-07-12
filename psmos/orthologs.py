"""
orthologs.py — pull REAL ortholog protein sequences from UniProt.

For each (family, species) in a pathway's curated ortholog seed table, query the
UniProt REST API for that species' gene and return the canonical protein
sequence. Reviewed (Swiss-Prot) entries are preferred; otherwise the longest
UniProtKB hit is taken (and flagged). Absence of any hit is itself a signal —
it is exactly what should happen for the negative-control species and it drives
the hard gate.

Everything here is real network data: the accession chosen and the sequence
returned are logged so the provenance of every downstream Evo2 score is
auditable. Nothing is asserted.

Requires only `requests` and network egress to rest.uniprot.org (verified
reachable from this environment). Set SSL_CERT_FILE to the agent-proxy CA.
"""

import json
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "data", "orthologs")


@dataclass
class Ortholog:
    family: str
    species_key: str
    taxon: int
    query_symbol: str
    accession: str | None
    reviewed: bool
    protein_name: str
    gene: str
    length: int
    sequence: str
    found: bool
    note: str = ""

    def slim(self) -> dict:
        d = asdict(self)
        d["sequence"] = self.sequence[:60] + ("..." if len(self.sequence) > 60 else "")
        return d


def _get(url: str, retries: int = 4) -> dict:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode())
        except Exception as e:  # network hiccup: exponential backoff
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"UniProt request failed after {retries} tries: {last}")


def fetch_one(family: str, species_key: str, taxon: int, symbol: str,
              ref_len: int | None = None) -> Ortholog:
    """Query UniProt for one species' ortholog of one family.

    Selection: prefer reviewed (Swiss-Prot); then, if a reference length is
    known (the human ortholog of this family), prefer the entry whose length is
    closest to it — orthologues are near-equal length, so this rejects spurious
    fusion/mis-annotated isoforms that a naive "take the longest" rule grabs.
    With no reference length yet (the human query itself), fall back to longest.
    """
    q = f"(gene:{symbol}) AND (organism_id:{taxon})"
    params = {
        "query": q,
        "format": "json",
        "size": "25",
        "fields": "accession,reviewed,protein_name,gene_names,length,sequence",
    }
    url = f"{UNIPROT}?{urllib.parse.urlencode(params)}"
    data = _get(url)
    results = data.get("results", [])

    if not results:
        return Ortholog(family, species_key, taxon, symbol, None, False, "", "", 0, "",
                        found=False, note="no UniProtKB hit")

    def rank(entry):
        reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
        length = entry.get("sequence", {}).get("length", 0)
        if ref_len:
            return (1 if reviewed else 0, -abs(length - ref_len))
        return (1 if reviewed else 0, length)

    best = max(results, key=rank)
    reviewed = best.get("entryType", "").startswith("UniProtKB reviewed")
    seq = best.get("sequence", {}).get("value", "")
    genes = best.get("genes", [])
    gene = ""
    if genes:
        gene = genes[0].get("geneName", {}).get("value", "") or symbol
    pname = (best.get("proteinDescription", {})
             .get("recommendedName", {}).get("fullName", {}).get("value", ""))
    if not pname:
        sub = best.get("proteinDescription", {}).get("submissionNames", [])
        if sub:
            pname = sub[0].get("fullName", {}).get("value", "")

    note = "" if reviewed else "unreviewed (TrEMBL) — longest hit"
    return Ortholog(
        family, species_key, taxon, symbol,
        best.get("primaryAccession"), reviewed, pname, gene,
        best.get("sequence", {}).get("length", len(seq)), seq,
        found=True, note=note,
    )


def fetch_pathway_gates(pathway, species_keys: list[str] | None = None,
                        use_cache: bool = True) -> list[Ortholog]:
    """Fetch ortholog sequences for every GATE family across the panel."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache = os.path.join(CACHE_DIR, f"{pathway.key}_gates.json")
    if use_cache and os.path.exists(cache):
        with open(cache) as fh:
            raw = json.load(fh)
        return [Ortholog(**o) for o in raw]

    out: list[Ortholog] = []
    for fam_key, per_species in pathway.ortholog_seed.items():
        # Resolve the human reference first to get this family's reference length.
        ref = pathway.reference_species
        ref_len = None
        if ref in per_species:
            sp = pathway_species(pathway, ref)
            ho = fetch_one(fam_key, ref, sp.taxon, per_species[ref][0])
            if ho.found:
                ref_len = ho.length

        targets = species_keys or list(per_species.keys())
        for sk in targets:
            if sk not in per_species:
                continue
            # Symbol list is a fallback order: take the first symbol that hits.
            o = None
            for sym in per_species[sk]:
                sp = pathway_species(pathway, sk)
                cand = fetch_one(fam_key, sk, sp.taxon, sym,
                                 ref_len=ref_len if sk != ref else None)
                if cand.found:
                    o = cand
                    break
                o = cand  # keep last (not-found) for reporting
            out.append(o)
            print(f"  [{fam_key:16s}] {sk:12s} {o.query_symbol:8s} -> "
                  f"{o.accession or '—':10s} len={o.length:5d} "
                  f"{'reviewed' if o.reviewed else o.note}")
    with open(cache, "w") as fh:
        json.dump([asdict(o) for o in out], fh, indent=1)
    return out


def pathway_species(pathway, species_key):
    from .pathways import SPECIES
    return SPECIES[species_key]


if __name__ == "__main__":
    import sys
    from psmos.pathways import get_pathway
    pk = sys.argv[1] if len(sys.argv) > 1 else "Notch"
    orths = fetch_pathway_gates(get_pathway(pk), use_cache=False)
    print(f"\nfetched {sum(o.found for o in orths)}/{len(orths)} ortholog sequences")
