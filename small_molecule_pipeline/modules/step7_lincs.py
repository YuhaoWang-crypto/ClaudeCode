"""Step 7 — Downstream functional signatures via LINCS / CMap.

Sources:
  - LINCS L1000 API  (iLINCS)   http://www.ilincs.org/api/
  - CMap Query API               https://clue.io/api

iLINCS provides a public REST API for querying perturbation signatures.
CMap requires an API key (set CLUE_API_KEY env var).
"""

import logging
import os
from typing import Dict, List, Optional

from utils.api_client import get, post

logger = logging.getLogger(__name__)

ILINCS_BASE = "https://www.ilincs.org/api"
CMAP_BASE = "https://api.clue.io/api"


# ── iLINCS ───────────────────────────────────────────────────────────────────

def _ilincs_pert_id(smiles: str) -> Optional[str]:
    """Search for a perturbagen by SMILES in iLINCS."""
    resp = get(
        f"{ILINCS_BASE}/ilincsR/findCompounds",
        params={"smiles": smiles},
        pause=1.0,
        timeout=60,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        compounds = data.get("data", [])
        if compounds:
            return compounds[0].get("pertid") or compounds[0].get("compound_id")
    except Exception:
        pass
    return None


def _ilincs_signature(pert_id: str) -> List[Dict]:
    """Retrieve gene expression signature for a LINCS perturbagen."""
    resp = get(
        f"{ILINCS_BASE}/SignatureFiles",
        params={"pertid": pert_id, "cellline": "MCF7"},  # default cell line
        pause=1.0,
        timeout=60,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        sigs = data.get("data", [])
        return [
            {
                "source": "LINCS_iLINCS",
                "signature_id": s.get("signatureid"),
                "cell_line": s.get("cellline"),
                "time_point": s.get("timepoint"),
                "dose": s.get("concentration"),
                "dose_unit": s.get("concentration_unit"),
                "n_genes_up": s.get("up_genes"),
                "n_genes_dn": s.get("dn_genes"),
            }
            for s in sigs
        ]
    except Exception:
        return []


def _ilincs_similar_compounds(pert_id: str, top_n: int = 10) -> List[Dict]:
    """Find compounds with similar transcriptional signatures."""
    resp = get(
        f"{ILINCS_BASE}/ilincsR/findSimilarCompounds",
        params={"pertid": pert_id, "topN": top_n},
        pause=1.0,
        timeout=60,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        return [
            {
                "pert_id": c.get("pertid"),
                "name": c.get("compoundname"),
                "similarity_score": c.get("concordance"),
            }
            for c in data.get("data", [])
        ]
    except Exception:
        return []


# ── CMap (Clue.io) ───────────────────────────────────────────────────────────

def _cmap_available() -> bool:
    return bool(os.environ.get("CLUE_API_KEY", ""))


def _cmap_query(smiles: str) -> Optional[Dict]:
    """Query CMap for transcriptional connectivity."""
    api_key = os.environ.get("CLUE_API_KEY", "")
    if not api_key:
        return None

    # Step 1: find pert_id
    resp = get(
        f"{CMAP_BASE}/perts",
        params={"filter": f'{{"where":{{"canonical_smiles":"{smiles}"}}}}'},
        headers={"user_key": api_key},
        pause=1.0,
    )
    if not resp:
        return None
    try:
        perts = resp.json()
        if not perts:
            return None
        pert_id = perts[0].get("pert_id")
    except Exception:
        return None

    # Step 2: fetch signature metadata
    resp2 = get(
        f"{CMAP_BASE}/sigs",
        params={
            "filter": f'{{"where":{{"pert_id":"{pert_id}"}},"limit":20}}',
        },
        headers={"user_key": api_key},
        pause=1.0,
    )
    if not resp2:
        return None
    try:
        sigs = resp2.json()
        return {
            "pert_id": pert_id,
            "signatures": [
                {
                    "sig_id": s.get("sig_id"),
                    "cell_id": s.get("cell_id"),
                    "pert_dose": s.get("pert_dose"),
                    "pert_time": s.get("pert_time"),
                    "distil_cc_q75": s.get("distil_cc_q75"),
                }
                for s in sigs
            ],
        }
    except Exception:
        return None


# ── Aggregator ───────────────────────────────────────────────────────────────

def run(step1_results: List[Dict]) -> List[Dict]:
    cmap_on = _cmap_available()
    if cmap_on:
        logger.info("CMap enabled (CLUE_API_KEY set)")
    else:
        logger.info("CMap disabled (set CLUE_API_KEY to enable)")

    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        smiles = mol.get("std_smiles", "")
        inchikey = mol.get("inchikey", "")
        if not smiles:
            continue

        # iLINCS
        pert_id = _ilincs_pert_id(smiles)
        sigs: List[Dict] = []
        similar_compounds: List[Dict] = []
        if pert_id:
            sigs = _ilincs_signature(pert_id)
            similar_compounds = _ilincs_similar_compounds(pert_id)

        # CMap
        cmap_result: Optional[Dict] = None
        if cmap_on:
            cmap_result = _cmap_query(smiles)

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": smiles,
                "ilincs_pert_id": pert_id,
                "lincs_signatures": sigs,
                "n_lincs_sigs": len(sigs),
                "similar_compounds_lincs": similar_compounds,
                "cmap": cmap_result,
                "has_lincs_data": pert_id is not None,
                "has_cmap_data": cmap_result is not None,
            }
        )
        logger.info(
            "  %s → LINCS: %s (%d sigs) | CMap: %s",
            inchikey[:14],
            pert_id or "not found",
            len(sigs),
            "found" if cmap_result else "not found",
        )

    logger.info(
        "Step 7 complete: %d molecules, %d with LINCS signatures",
        len(output),
        sum(1 for r in output if r["has_lincs_data"]),
    )
    return output
