"""Step 2 — Experimental target retrieval.

Sources:
  - ChEMBL REST API   https://www.ebi.ac.uk/chembl/api/data/
  - BindingDB REST    https://bindingdb.org/axis2/services/BDBService
  - PubChem BioAssay  https://pubchem.ncbi.nlm.nih.gov/rest/pug/
  - GtoPdb REST       https://www.guidetopharmacology.org/services/
"""

import logging
from typing import Dict, List, Optional

from utils.api_client import get, post

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
BINDINGDB_BASE = "https://bindingdb.org/axis2/services/BDBService"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
GTOPDB_BASE = "https://www.guidetopharmacology.org/services"


# ── ChEMBL ──────────────────────────────────────────────────────────────────

def _chembl_by_inchikey(inchikey: str) -> Optional[str]:
    """Return ChEMBL molecule ID for an InChIKey."""
    resp = get(
        f"{CHEMBL_BASE}/molecule/{inchikey}.json",
        pause=0.3,
    )
    if resp:
        return resp.json().get("molecule_chembl_id")
    return None


def _chembl_activities(chembl_id: str, limit: int = 100) -> List[Dict]:
    resp = get(
        f"{CHEMBL_BASE}/activity.json",
        params={
            "molecule_chembl_id": chembl_id,
            "limit": limit,
            "format": "json",
        },
        pause=0.3,
    )
    if not resp:
        return []
    data = resp.json()
    activities = data.get("activities", [])
    targets = []
    for act in activities:
        if not act.get("target_chembl_id"):
            continue
        targets.append(
            {
                "source": "ChEMBL",
                "target_id": act["target_chembl_id"],
                "target_name": act.get("target_pref_name"),
                "target_type": act.get("target_organism"),
                "assay_type": act.get("assay_type"),
                "standard_type": act.get("standard_type"),
                "standard_value": act.get("standard_value"),
                "standard_units": act.get("standard_units"),
                "pchembl_value": act.get("pchembl_value"),
                "doc_id": act.get("document_chembl_id"),
            }
        )
    return targets


def query_chembl(inchikey: str) -> List[Dict]:
    chembl_id = _chembl_by_inchikey(inchikey)
    if not chembl_id:
        return []
    return _chembl_activities(chembl_id)


# ── BindingDB ────────────────────────────────────────────────────────────────

def query_bindingdb(inchikey: str) -> List[Dict]:
    resp = get(
        f"{BINDINGDB_BASE}/getLigandsByInchikey",
        params={"inchikey": inchikey, "response": "json"},
        pause=0.5,
    )
    if not resp:
        return []
    try:
        data = resp.json()
    except Exception:
        return []

    entries = data.get("getLigandsByInchikeyResponse", {}).get("affinities", [])
    if isinstance(entries, dict):
        entries = [entries]

    results = []
    for e in entries:
        results.append(
            {
                "source": "BindingDB",
                "target_id": e.get("UniProt_Id"),
                "target_name": e.get("Target_Name"),
                "target_type": e.get("Target_Source_Organism"),
                "assay_type": "binding",
                "standard_type": e.get("Affinity_Type"),
                "standard_value": e.get("Affinity"),
                "standard_units": "nM",
                "pchembl_value": None,
                "doc_id": e.get("Publication"),
            }
        )
    return results


# ── PubChem BioAssay ─────────────────────────────────────────────────────────

def query_pubchem(inchikey: str) -> List[Dict]:
    # Resolve CID from InChIKey
    resp = get(
        f"{PUBCHEM_BASE}/compound/inchikey/{inchikey}/cids/JSON",
        pause=0.3,
    )
    if not resp:
        return []
    cids = resp.json().get("IdentifierList", {}).get("CID", [])
    if not cids:
        return []
    cid = cids[0]

    # Get bioassay activity
    resp2 = get(
        f"{PUBCHEM_BASE}/compound/cid/{cid}/assaysummary/JSON",
        pause=0.5,
    )
    if not resp2:
        return []
    try:
        data = resp2.json()
    except Exception:
        return []

    table = data.get("Table", {})
    columns = [c.get("name", "") for c in table.get("Columns", {}).get("Column", [])]
    rows = table.get("Row", [])

    results = []
    for row in rows[:100]:
        cells = row.get("Cell", [])
        rec = dict(zip(columns, cells))
        if rec.get("Activity Outcome") not in {"Active", "active"}:
            continue
        results.append(
            {
                "source": "PubChem_BioAssay",
                "target_id": rec.get("Target GI"),
                "target_name": rec.get("Target Name"),
                "target_type": rec.get("Target Taxonomy"),
                "assay_type": rec.get("Assay Type"),
                "standard_type": "Activity",
                "standard_value": rec.get("Activity Value [uM]"),
                "standard_units": "uM",
                "pchembl_value": None,
                "doc_id": f"AID{rec.get('AID', '')}",
            }
        )
    return results


# ── GtoPdb ───────────────────────────────────────────────────────────────────

def query_gtopdb(inchikey: str) -> List[Dict]:
    resp = get(
        f"{GTOPDB_BASE}/ligands",
        params={"inchiKey": inchikey, "type": "Approved"},
        pause=0.5,
    )
    if not resp:
        return []
    try:
        ligands = resp.json()
    except Exception:
        return []
    if not ligands:
        return []

    ligand_id = ligands[0].get("ligandId")
    if not ligand_id:
        return []

    resp2 = get(
        f"{GTOPDB_BASE}/ligands/{ligand_id}/interactions",
        pause=0.5,
    )
    if not resp2:
        return []
    try:
        interactions = resp2.json()
    except Exception:
        return []

    results = []
    for ia in interactions:
        results.append(
            {
                "source": "GtoPdb",
                "target_id": ia.get("targetId"),
                "target_name": ia.get("targetName"),
                "target_type": ia.get("targetSpecies"),
                "assay_type": ia.get("type"),
                "standard_type": ia.get("affinityParameter"),
                "standard_value": ia.get("affinity"),
                "standard_units": ia.get("affinityUnits"),
                "pchembl_value": None,
                "doc_id": str(ia.get("primaryReference", "")),
            }
        )
    return results


# ── Aggregator ───────────────────────────────────────────────────────────────

def run(step1_results: List[Dict]) -> List[Dict]:
    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        inchikey = mol.get("inchikey", "")
        if not inchikey:
            continue

        targets: List[Dict] = []
        targets += query_chembl(inchikey)
        targets += query_bindingdb(inchikey)
        targets += query_pubchem(inchikey)
        targets += query_gtopdb(inchikey)

        # Deduplicate by (source, target_id)
        seen = set()
        deduped = []
        for t in targets:
            key = (t["source"], str(t.get("target_id", "")))
            if key not in seen:
                seen.add(key)
                deduped.append(t)

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": mol.get("std_smiles"),
                "targets": deduped,
                "n_targets": len(deduped),
                "sources": list({t["source"] for t in deduped}),
            }
        )
        logger.info(
            "  %s → %d experimental targets", inchikey[:14], len(deduped)
        )

    logger.info(
        "Step 2 complete: %d molecules queried across ChEMBL/BindingDB/PubChem/GtoPdb",
        len(output),
    )
    return output
