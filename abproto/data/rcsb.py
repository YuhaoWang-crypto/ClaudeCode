"""Fetch PDB entry + polymer-entity metadata from the RCSB GraphQL batch API.

One POST returns dates, resolution, and per-entity descriptions + sequences for
many entries at once. No API key required.
"""

from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, field
from typing import Optional

GRAPHQL_URL = "https://data.rcsb.org/graphql"

_QUERY = """
{ entries(entry_ids: %s) {
    rcsb_id
    rcsb_accession_info { deposit_date initial_release_date }
    rcsb_entry_info { resolution_combined }
    polymer_entities {
      rcsb_id
      rcsb_polymer_entity { pdbx_description }
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
    }
} }
"""


@dataclass
class Entity:
    entity_id: str
    description: str
    sequence: str
    polymer_type: str


@dataclass
class Entry:
    pdb_id: str
    deposit_date: str            # "YYYY-MM-DD" or ""
    release_date: str
    resolution: Optional[float]
    entities: list[Entity] = field(default_factory=list)


def _post(ids: list[str], timeout: float) -> dict:
    query = _QUERY % json.dumps(ids)
    body = json.dumps({"query": query}).encode()
    req = urllib.request.Request(
        GRAPHQL_URL, data=body, headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:   # noqa: S310 (trusted host)
        return json.loads(r.read().decode())


def fetch_entries(pdb_ids: list[str], chunk: int = 50, timeout: float = 60.0) -> list[Entry]:
    """Return parsed entries for ``pdb_ids`` (batched; unknown ids are dropped)."""
    out: list[Entry] = []
    ids = [i.upper() for i in pdb_ids]
    for i in range(0, len(ids), chunk):
        data = _post(ids[i : i + chunk], timeout)
        for e in (data.get("data", {}) or {}).get("entries", []) or []:
            if e is None:
                continue
            res = (e.get("rcsb_entry_info") or {}).get("resolution_combined")
            ents = []
            for pe in e.get("polymer_entities") or []:
                poly = pe.get("entity_poly") or {}
                ents.append(Entity(
                    entity_id=pe.get("rcsb_id", ""),
                    description=(pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "",
                    sequence=poly.get("pdbx_seq_one_letter_code_can") or "",
                    polymer_type=poly.get("rcsb_entity_polymer_type") or "",
                ))
            out.append(Entry(
                pdb_id=e.get("rcsb_id", ""),
                deposit_date=((e.get("rcsb_accession_info") or {}).get("deposit_date") or "")[:10],
                release_date=((e.get("rcsb_accession_info") or {}).get("initial_release_date") or "")[:10],
                resolution=float(res[0]) if res else None,
                entities=ents,
            ))
    return out


def search_antibody_antigen(max_rows: int = 200, timeout: float = 60.0) -> list[str]:
    """Discover candidate Ab-Ag complexes via the RCSB Search API.

    Full-text 'Fab' restricted to protein-protein assemblies with >= 3 protein
    entities (two antibody chains + >=1 antigen). This is a recall-oriented
    seed; classify.to_complexes does the real filtering. Returns PDB ids.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
        "query": {
            "type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "full_text",
                 "parameters": {"value": "Fab antibody"}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein",
                                "operator": "greater_or_equal", "value": 3}},
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": max_rows},
                            "results_content_type": ["experimental"]},
    }
    req = urllib.request.Request(
        url + "?json=" + urllib.parse.quote(json.dumps(query)),
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
            data = json.loads(r.read().decode())
    except Exception:
        return []
    return [hit["identifier"] for hit in data.get("result_set", [])]


# urllib.parse referenced in search_antibody_antigen
import urllib.parse  # noqa: E402
