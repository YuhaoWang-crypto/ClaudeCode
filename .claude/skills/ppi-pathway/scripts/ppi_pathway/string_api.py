"""
string_api.py — thin, dependency-light wrapper over the STRING REST API.

Only needs ``requests``. Use this when you have a short-to-medium gene list
(an omics hit list) and want pathway/functional enrichment or a small PPI
network *without* downloading the multi-hundred-MB bulk files.

Endpoints wrapped
-----------------
* ``enrichment``            — functional/pathway enrichment (GO, KEGG, Reactome,
                              WikiPathways, Pfam, ...), with FDR correction.
* ``network``               — scored interactions among the input genes.
* ``interaction_partners``  — nearest STRING neighbours of the input genes.
* ``get_string_ids``        — map arbitrary identifiers to STRING preferred names.
* ``ppi_enrichment``        — is the network more connected than random?

STRING asks callers to identify themselves (``caller_identity``) and to rate-
limit; for >1 request add a short delay. See https://string-db.org/help/api/ .
"""
from __future__ import annotations

import csv
import io
import time
from typing import Iterable

import requests

from . import HUMAN_TAXON

BASE = "https://string-db.org/api"
CALLER = "ppi_pathway_skill"


def _post(endpoint: str, fmt: str, params: dict) -> requests.Response:
    url = f"{BASE}/{fmt}/{endpoint}"
    params = {"caller_identity": CALLER, **params}
    resp = requests.post(url, data=params, timeout=120)
    resp.raise_for_status()
    return resp


def _tsv_to_dicts(text: str) -> list[dict]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return [dict(row) for row in reader]


def _ids_param(genes: Iterable[str]) -> str:
    # STRING accepts carriage-return-separated identifiers in a form field.
    return "\r".join(str(g).strip() for g in genes if str(g).strip())


def enrichment(genes: Iterable[str], taxon: int = HUMAN_TAXON,
               background: Iterable[str] | None = None) -> list[dict]:
    """Functional/pathway enrichment for ``genes``.

    Returns a list of rows (one per enriched term) with keys including
    ``category`` (KEGG / Process / Reactome / WikiPathways / ...), ``term``,
    ``description``, ``p_value``, ``fdr``, ``number_of_genes`` and the matching
    ``inputGenes`` / ``preferredNames``. Already FDR-corrected by STRING.
    """
    params = {"identifiers": _ids_param(genes), "species": taxon}
    if background is not None:
        params["background_string_identifiers"] = _ids_param(background)
    rows = _tsv_to_dicts(_post("enrichment", "tsv", params).text)
    for r in rows:
        for k in ("p_value", "fdr"):
            if k in r:
                try:
                    r[k] = float(r[k])
                except ValueError:
                    pass
    return rows


def network(genes: Iterable[str], taxon: int = HUMAN_TAXON,
            required_score: int = 400, add_nodes: int = 0) -> list[dict]:
    """Scored interactions among ``genes`` (optionally padded with the top
    ``add_nodes`` connected partners). ``required_score`` is 0-1000
    (400 = medium, 700 = high confidence)."""
    params = {
        "identifiers": _ids_param(genes),
        "species": taxon,
        "required_score": required_score,
        "add_nodes": add_nodes,
    }
    return _tsv_to_dicts(_post("network", "tsv", params).text)


def interaction_partners(genes: Iterable[str], taxon: int = HUMAN_TAXON,
                         required_score: int = 400, limit: int = 20) -> list[dict]:
    """Nearest STRING neighbours of each input gene (up to ``limit`` each)."""
    params = {
        "identifiers": _ids_param(genes),
        "species": taxon,
        "required_score": required_score,
        "limit": limit,
    }
    return _tsv_to_dicts(_post("interaction_partners", "tsv", params).text)


def map_ids(genes: Iterable[str], taxon: int = HUMAN_TAXON,
            limit: int = 1) -> list[dict]:
    """Resolve arbitrary identifiers to STRING preferred names / STRING ids.

    Useful to sanity-check that your gene symbols are recognised before
    running enrichment. ``limit`` caps matches per query identifier."""
    params = {
        "identifiers": _ids_param(genes),
        "species": taxon,
        "limit": limit,
        "echo_query": 1,
    }
    return _tsv_to_dicts(_post("get_string_ids", "tsv", params).text)


def ppi_enrichment(genes: Iterable[str], taxon: int = HUMAN_TAXON,
                   required_score: int = 400) -> dict:
    """Whole-network PPI enrichment p-value: are these genes more densely
    connected than a random set of the same size? A small p-value means the
    gene set is functionally coherent (a real module, not noise)."""
    params = {
        "identifiers": _ids_param(genes),
        "species": taxon,
        "required_score": required_score,
    }
    rows = _tsv_to_dicts(_post("ppi_enrichment", "tsv", params).text)
    return rows[0] if rows else {}


def top_pathways(genes: Iterable[str], taxon: int = HUMAN_TAXON,
                 fdr_max: float = 0.05,
                 categories: tuple[str, ...] = ("KEGG", "Reactome Pathways",
                                                "WikiPathways"),
                 n: int = 25) -> list[dict]:
    """Convenience: enrichment filtered to pathway categories, FDR-sorted."""
    rows = enrichment(genes, taxon=taxon)
    keep = [r for r in rows
            if r.get("category") in categories and float(r.get("fdr", 1)) <= fdr_max]
    keep.sort(key=lambda r: float(r.get("fdr", 1)))
    return keep[:n]


if __name__ == "__main__":  # tiny manual smoke test
    import json
    demo = ["TP53", "BRCA1", "EGFR", "MYC", "CDK2", "CDK4", "RB1", "MDM2"]
    print("PPI enrichment:", json.dumps(ppi_enrichment(demo), indent=2))
    for row in top_pathways(demo, n=8):
        print(f"  {row['category']:20s} {row['fdr']:.2e}  {row['description']}")
        time.sleep(0.2)
