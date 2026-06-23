"""Step 5 — Pathway and target-disease context.

Sources:
  - Reactome REST API   https://reactome.org/ContentService/
  - Open Targets GraphQL https://api.platform.opentargets.org/api/v4/graphql
"""

import logging
from typing import Dict, List, Optional

from utils.api_client import get, post

logger = logging.getLogger(__name__)

REACTOME_BASE = "https://reactome.org/ContentService"
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"


# ── Reactome ─────────────────────────────────────────────────────────────────

def _reactome_pathways_for_uniprot(uniprot: str) -> List[Dict]:
    """Retrieve pathways containing a given protein."""
    resp = get(
        f"{REACTOME_BASE}/data/pathways/low/entity/{uniprot}/allForms",
        params={"species": "9606"},  # Homo sapiens NCBI taxon ID
        pause=0.5,
    )
    if not resp:
        return []
    try:
        pathways = resp.json()
        return [
            {
                "reactome_id": p.get("stId"),
                "name": p.get("displayName"),
                "is_inferred": p.get("isInferred", False),
            }
            for p in (pathways if isinstance(pathways, list) else [])
        ]
    except Exception as e:
        logger.debug("Reactome parse error for %s: %s", uniprot, e)
        return []


# ── Open Targets ─────────────────────────────────────────────────────────────

_OT_TARGET_QUERY = """
query TargetAssociation($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    approvedName
    biotype
    associatedDiseases(page: { index: 0, size: 20 }) {
      rows {
        disease {
          id
          name
          therapeuticAreas { id name }
        }
        score
        datatypeScores { id score }
      }
    }
    pathways {
      pathway
      pathwayId
      topLevelTerm
    }
  }
}
"""

_OT_CHEMBL_QUERY = """
query DrugTargets($chemblId: String!) {
  drug(chemblId: $chemblId) {
    id
    name
    mechanismsOfAction {
      rows {
        mechanismOfAction
        targets { id approvedSymbol }
      }
    }
    linkedDiseases {
      rows { id name }
    }
  }
}
"""


def _ot_target_by_ensembl(ensembl_id: str) -> Optional[Dict]:
    resp = post(
        OT_GRAPHQL,
        json={"query": _OT_TARGET_QUERY, "variables": {"ensemblId": ensembl_id}},
        headers={"Content-Type": "application/json"},
        pause=0.5,
    )
    if not resp:
        return None
    try:
        return resp.json().get("data", {}).get("target")
    except Exception:
        return None


def _ot_drug_by_chembl(chembl_id: str) -> Optional[Dict]:
    resp = post(
        OT_GRAPHQL,
        json={"query": _OT_CHEMBL_QUERY, "variables": {"chemblId": chembl_id}},
        headers={"Content-Type": "application/json"},
        pause=0.5,
    )
    if not resp:
        return None
    try:
        return resp.json().get("data", {}).get("drug")
    except Exception:
        return None


def _uniprot_to_ensembl(uniprot: str) -> Optional[str]:
    """Map UniProt accession → Ensembl gene ID via UniProt REST."""
    resp = get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot}.json",
        pause=0.3,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "Ensembl":
                for prop in xref.get("properties", []):
                    if prop.get("key") == "GeneId":
                        return prop.get("value")
    except Exception:
        pass
    return None


# ── Aggregator ───────────────────────────────────────────────────────────────

def run(
    step1_results: List[Dict],
    step2_results: List[Dict],
    step3_results: List[Dict],
) -> List[Dict]:
    """Combine Reactome pathways and Open Targets disease associations."""

    # Index step2 by inchikey for quick lookup
    step2_idx = {r["inchikey"]: r for r in step2_results}
    step3_idx = {r["inchikey"]: r for r in step3_results}

    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        inchikey = mol.get("inchikey", "")

        s2 = step2_idx.get(inchikey, {})
        s3 = step3_idx.get(inchikey, {})

        # Collect UniProt IDs from experimental target data
        uniprots = set()
        for t in s2.get("targets", []):
            uid = t.get("uniprot") or t.get("target_id", "")
            if uid and uid.startswith("P") or (uid and len(uid) in {6, 10}):
                uniprots.add(uid)

        # Query Reactome for each UniProt
        reactome_pathways = []
        for uid in list(uniprots)[:10]:  # cap to avoid excessive calls
            for pw in _reactome_pathways_for_uniprot(uid):
                pw["uniprot"] = uid
                reactome_pathways.append(pw)

        # Deduplicate pathways
        seen_pw = set()
        deduped_pw = []
        for pw in reactome_pathways:
            rid = pw["reactome_id"]
            if rid not in seen_pw:
                seen_pw.add(rid)
                deduped_pw.append(pw)

        # Open Targets via ChEMBL ID (if molecule is a known drug)
        ot_drug_data: Optional[Dict] = None
        dc_info = s3.get("drugcentral", {})
        # Try to get a ChEMBL ID for this molecule from ChEMBL step
        chembl_id = None
        for t in s2.get("targets", []):
            if t.get("source") == "ChEMBL":
                pass  # we don't store the mol-level ChEMBL ID here directly

        # Use Open Targets target lookup for known targets
        ot_disease_assocs: List[Dict] = []
        for uid in list(uniprots)[:5]:
            ensembl = _uniprot_to_ensembl(uid)
            if ensembl:
                ot_data = _ot_target_by_ensembl(ensembl)
                if ot_data:
                    for row in ot_data.get("associatedDiseases", {}).get("rows", []):
                        disease = row.get("disease", {})
                        ot_disease_assocs.append(
                            {
                                "uniprot": uid,
                                "ensembl": ensembl,
                                "target_symbol": ot_data.get("approvedSymbol"),
                                "disease_id": disease.get("id"),
                                "disease_name": disease.get("name"),
                                "association_score": round(row.get("score", 0), 4),
                                "therapeutic_areas": [
                                    ta["name"]
                                    for ta in disease.get("therapeuticAreas", [])
                                ],
                            }
                        )

        ot_disease_assocs.sort(key=lambda x: x["association_score"], reverse=True)

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": mol.get("std_smiles"),
                "reactome_pathways": deduped_pw,
                "n_reactome_pathways": len(deduped_pw),
                "ot_disease_associations": ot_disease_assocs[:50],
                "n_ot_diseases": len(ot_disease_assocs),
                "top_diseases": [
                    d["disease_name"] for d in ot_disease_assocs[:5]
                ],
            }
        )
        logger.info(
            "  %s → Reactome: %d pathways | OT: %d disease assocs",
            inchikey[:14],
            len(deduped_pw),
            len(ot_disease_assocs),
        )

    logger.info(
        "Step 5 complete: %d molecules with pathway/disease context", len(output)
    )
    return output
