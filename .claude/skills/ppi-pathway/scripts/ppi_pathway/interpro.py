"""
interpro.py — protein family / domain annotation via the InterPro REST API
(https://www.ebi.ac.uk/interpro/, https://www.ebi.ac.uk/interpro/api/).

InterPro integrates member databases (Pfam, PANTHER, PROSITE, SMART, CDD, ...)
into non-redundant protein families, domains and functional sites, each with GO
term mappings. Use it to annotate the hubs of your network module — "what *is*
this protein / what domains does it carry / what molecular function" — which
turns a gene list into mechanistic hypotheses.

Only needs ``requests``. Query by UniProt accession (e.g. P04637).
"""
from __future__ import annotations

import time
from typing import Iterable

import requests

BASE = "https://www.ebi.ac.uk/interpro/api"


def _get(path: str, params: dict | None = None) -> dict:
    r = requests.get(f"{BASE}/{path}", params=params or {},
                     headers={"Accept": "application/json"}, timeout=60)
    r.raise_for_status()
    return r.json()


def entries_for_protein(uniprot_acc: str, page_size: int = 100) -> list[dict]:
    """All InterPro entries (families/domains/sites) annotating a protein.

    Returns a list of dicts with ``accession`` (IPRxxxxxx), ``name``, ``type``
    (family/domain/homologous_superfamily/...), ``member_databases`` and
    ``go_terms``.
    """
    data = _get(f"entry/interpro/protein/reviewed/{uniprot_acc}",
                {"page_size": page_size})
    out = []
    for item in data.get("results", []):
        md = item.get("metadata", {})
        out.append({
            "accession": md.get("accession"),
            "name": md.get("name"),
            "type": md.get("type"),
            "member_databases": md.get("member_databases"),
            "go_terms": [{"id": g.get("identifier"), "name": g.get("name"),
                          "category": g.get("category", {}).get("name")}
                         for g in (md.get("go_terms") or [])],
        })
    return out


def domains(uniprot_acc: str) -> list[dict]:
    """Just the domain-type entries for a protein (structural/functional units)."""
    return [e for e in entries_for_protein(uniprot_acc)
            if e.get("type") == "domain"]


def go_terms(uniprot_acc: str) -> list[dict]:
    """Deduplicated GO terms implied by a protein's InterPro entries."""
    seen, out = set(), []
    for e in entries_for_protein(uniprot_acc):
        for g in e["go_terms"]:
            if g["id"] and g["id"] not in seen:
                seen.add(g["id"])
                out.append(g)
    return out


def annotate(uniprot_accs: Iterable[str], sleep: float = 0.34) -> dict[str, list[dict]]:
    """Annotate several proteins. ``sleep`` throttles to be polite to EBI."""
    result = {}
    for acc in uniprot_accs:
        try:
            result[acc] = entries_for_protein(acc)
        except Exception as e:
            result[acc] = [{"error": str(e)}]
        time.sleep(sleep)
    return result


if __name__ == "__main__":
    for e in entries_for_protein("P04637"):
        gos = ", ".join(g["id"] for g in e["go_terms"][:3])
        print(f"{e['accession']:12s} {e['type']:22s} {e['name']}  [{gos}]")
