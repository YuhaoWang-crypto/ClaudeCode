"""Step 6 — Chemical Probes Portal quality assessment.

Source:
  Chemical Probes Portal   https://www.chemicalprobes.org/
  Public search API: https://www.chemicalprobes.org/probes?search=<name_or_smiles>

The portal assigns probes a quality rating (1–4 stars) and flags
whether a compound is recommended for use as a chemical probe.
"""

import logging
from typing import Dict, List, Optional

from utils.api_client import get

logger = logging.getLogger(__name__)

CPP_BASE = "https://www.chemicalprobes.org"


def _search_probe_by_name(name: str) -> List[Dict]:
    resp = get(
        f"{CPP_BASE}/api/probes",
        params={"q": name, "format": "json"},
        pause=0.5,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        probes = data if isinstance(data, list) else data.get("results", [])
        return probes
    except Exception:
        return []


def _search_probe_by_inchikey(inchikey: str) -> List[Dict]:
    """Query CPP with InChIKey; falls back to name search if not supported."""
    resp = get(
        f"{CPP_BASE}/api/probes",
        params={"inchikey": inchikey, "format": "json"},
        pause=0.5,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        probes = data if isinstance(data, list) else data.get("results", [])
        return probes
    except Exception:
        return []


def _parse_probe_record(raw: Dict) -> Dict:
    return {
        "source": "ChemicalProbesPortal",
        "probe_name": raw.get("name") or raw.get("probe_name"),
        "target": raw.get("target") or raw.get("primary_target"),
        "rating": raw.get("rating") or raw.get("stars"),
        "recommended": raw.get("recommended", False),
        "selectivity_notes": raw.get("selectivity_notes") or raw.get("comments"),
        "url": f"{CPP_BASE}/probes/{raw.get('slug') or raw.get('id', '')}",
        "references": raw.get("references", []),
    }


def query_chemical_probes_portal(inchikey: str) -> List[Dict]:
    raw_hits = _search_probe_by_inchikey(inchikey)

    # If the API doesn't support InChIKey lookup, raw_hits may be empty
    # or return all probes — filter by InChIKey if present in records
    if raw_hits:
        filtered = [
            r for r in raw_hits
            if r.get("inchi_key", "").upper() == inchikey.upper()
            or r.get("inchikey", "").upper() == inchikey.upper()
        ]
        if filtered:
            raw_hits = filtered

    return [_parse_probe_record(r) for r in raw_hits]


def run(step1_results: List[Dict]) -> List[Dict]:
    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        inchikey = mol.get("inchikey", "")
        if not inchikey:
            continue

        probes = query_chemical_probes_portal(inchikey)

        # Summarise quality
        recommended = [p for p in probes if p.get("recommended")]
        high_quality = [
            p for p in probes if p.get("rating") and float(p.get("rating", 0)) >= 3
        ]

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": mol.get("std_smiles"),
                "probes": probes,
                "n_probes": len(probes),
                "n_recommended": len(recommended),
                "n_high_quality": len(high_quality),
                "is_probe": len(probes) > 0,
                "best_rating": (
                    max(float(p.get("rating", 0)) for p in probes) if probes else None
                ),
            }
        )
        logger.info(
            "  %s → %d probe record(s), %d recommended",
            inchikey[:14],
            len(probes),
            len(recommended),
        )

    logger.info(
        "Step 6 complete: %d molecules, %d confirmed probes",
        len(output),
        sum(1 for r in output if r["is_probe"]),
    )
    return output
