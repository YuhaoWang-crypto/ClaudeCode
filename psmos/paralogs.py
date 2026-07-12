"""
paralogs.py — compute WITHIN-species paralogue counts via Ensembl Compara.

The "low-redundancy / deconvolution" axis (a single knockout gives a clean
phenotype only if the gene isn't buffered by paralogues) was a curated prior.
This module makes it computed: for each pathway component family, ask Ensembl
Compara how many within-species paralogues the family's gene has, per species.

    /homology/symbol/{species}/{gene}?type=paralogues  -> within-species paralogues

A species with many extra copies of the pathway's genes (e.g. teleost
whole-genome-duplication fish) scores LOW redundancy-simplicity; a single-copy
pathway (fly, basal metazoans) scores HIGH. The per-species mean paralogue count
across families becomes the redundancy signal, normalised to a 0..1 simplicity
score (fewer paralogues -> simpler -> higher).

Real network calls to rest.ensembl.org (Compara). Species not on Ensembl resolve
to no data (auditable gap), same as the CDS step.
"""

import json
import os
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict

from .cds import ENSEMBL_SPECIES

ENSEMBL = "https://rest.ensembl.org"
CACHE_DIR = os.path.join(os.path.dirname(__file__), "data", "paralogs")


@dataclass
class ParalogCount:
    family: str
    species_key: str
    gene: str
    n_paralogues: int      # extra within-species copies (0 = single copy)
    found: bool
    note: str = ""


def _get(url, retries=4):
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=40) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (400, 404):
                return None
            last = e
            time.sleep(2 ** i)
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"Ensembl Compara request failed: {last}")


def count_one(family, species_key, symbols):
    esp = ENSEMBL_SPECIES.get(species_key)
    if not esp:
        return ParalogCount(family, species_key, symbols[0] if symbols else "", 0,
                            False, "no Ensembl species map")
    for sym in symbols:
        url = (f"{ENSEMBL}/homology/symbol/{esp}/{urllib.parse.quote(sym)}"
               f"?type=paralogues&format=condensed&content-type=application/json")
        data = _get(url)
        if data and data.get("data") and data["data"][0].get("homologies") is not None:
            hom = data["data"][0]["homologies"]
            # count distinct within-species paralogue targets
            n = len({h.get("id") for h in hom if h.get("id")})
            return ParalogCount(family, species_key, sym, n, True)
    return ParalogCount(family, species_key, symbols[0] if symbols else "", 0,
                        False, "no Compara paralogues record")


def fetch_pathway_paralogs(pathway, use_cache=True):
    """Count paralogues for every seeded family across the panel's species."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    cache = os.path.join(CACHE_DIR, f"{pathway.key}_paralogs.json")
    if use_cache and os.path.exists(cache):
        return [ParalogCount(**c) for c in json.load(open(cache))]

    out = []
    for fam_key, per_species in pathway.ortholog_seed.items():
        for sk, syms in per_species.items():
            c = count_one(fam_key, sk, syms)
            out.append(c)
            print(f"  [{fam_key:12s}] {sk:15s} {c.gene:8s} paralogues={c.n_paralogues:2d} "
                  f"{'ok' if c.found else c.note}")
    json.dump([asdict(c) for c in out], open(cache, "w"), indent=1)
    return out


def species_simplicity(pathway, use_cache=True):
    """Per-species simplicity (0..1): fewer paralogues across the gate families
    => simpler => higher. Normalised across the species that have Compara data."""
    counts = fetch_pathway_paralogs(pathway, use_cache=use_cache)
    per_sp = {}
    for c in counts:
        if c.found:
            per_sp.setdefault(c.species_key, []).append(c.n_paralogues)
    mean_par = {sk: sum(v) / len(v) for sk, v in per_sp.items()}
    if not mean_par:
        return {}, {}
    lo, hi = min(mean_par.values()), max(mean_par.values())
    simp = {sk: (1.0 if hi == lo else 1 - (m - lo) / (hi - lo)) for sk, m in mean_par.items()}
    return mean_par, simp


if __name__ == "__main__":
    import sys
    from psmos.pathways import get_pathway
    pk = sys.argv[1] if len(sys.argv) > 1 else "Notch"
    pw = get_pathway(pk)
    mp, simp = species_simplicity(pw, use_cache=False)
    print(f"\n{pk}: computed within-species mean paralogue counts (gate families):")
    for sk in sorted(mp, key=mp.get):
        print(f"  {sk:15s} mean_paralogues={mp[sk]:.2f}  simplicity={simp[sk]:.2f}")
