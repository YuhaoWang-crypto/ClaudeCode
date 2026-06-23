"""Step 3 — Known drug mechanisms.

Sources:
  - DrugCentral REST API   https://drugcentral.org/api/
  - Repurposing Hub        downloadable TSV (Broad Institute)
    Local path configured via REPURPOSING_HUB_PATH env var, falls back to
    fetching from GitHub mirror.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

from utils.api_client import get

logger = logging.getLogger(__name__)

DRUGCENTRAL_BASE = "https://drugcentral.org/api/v1"

# Public mirror of the Repurposing Hub drug-target table
REPURPOSING_HUB_URL = (
    "https://raw.githubusercontent.com/broadinstitute/repurposing/"
    "master/drug-target-annotations.csv"
)


# ── DrugCentral ──────────────────────────────────────────────────────────────

def _drugcentral_drug_id(inchikey: str) -> Optional[int]:
    resp = get(
        f"{DRUGCENTRAL_BASE}/drug",
        params={"inchikey": inchikey},
        pause=0.5,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        drugs = data.get("data", [])
        if drugs:
            return drugs[0].get("id")
    except Exception:
        pass
    return None


def _drugcentral_indications(drug_id: int) -> List[Dict]:
    resp = get(
        f"{DRUGCENTRAL_BASE}/drug/{drug_id}/indications",
        pause=0.5,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        return [
            {
                "source": "DrugCentral",
                "record_type": "indication",
                "disease": ind.get("concept_name"),
                "umls_cui": ind.get("concept_id"),
                "snomed": ind.get("snomed_id"),
            }
            for ind in data.get("data", [])
        ]
    except Exception:
        return []


def _drugcentral_targets(drug_id: int) -> List[Dict]:
    resp = get(
        f"{DRUGCENTRAL_BASE}/drug/{drug_id}/targets",
        pause=0.5,
    )
    if not resp:
        return []
    try:
        data = resp.json()
        return [
            {
                "source": "DrugCentral",
                "record_type": "target",
                "target_name": t.get("gene"),
                "uniprot": t.get("accession"),
                "act_type": t.get("act_type"),
                "act_source": t.get("act_source"),
                "moa": t.get("moa"),
            }
            for t in data.get("data", [])
        ]
    except Exception:
        return []


def query_drugcentral(inchikey: str) -> Dict:
    drug_id = _drugcentral_drug_id(inchikey)
    if not drug_id:
        return {"drug_id": None, "indications": [], "targets": []}
    return {
        "drug_id": drug_id,
        "indications": _drugcentral_indications(drug_id),
        "targets": _drugcentral_targets(drug_id),
    }


# ── Repurposing Hub ───────────────────────────────────────────────────────────

_REPURPOSING_DF: Optional[pd.DataFrame] = None


def _load_repurposing_hub() -> Optional[pd.DataFrame]:
    global _REPURPOSING_DF
    if _REPURPOSING_DF is not None:
        return _REPURPOSING_DF

    local_path = os.environ.get("REPURPOSING_HUB_PATH", "")
    if local_path and Path(local_path).exists():
        logger.info("Loading Repurposing Hub from local file: %s", local_path)
        try:
            _REPURPOSING_DF = pd.read_csv(local_path, sep="\t", low_memory=False)
            return _REPURPOSING_DF
        except Exception as e:
            logger.warning("Failed to read local Repurposing Hub: %s", e)

    # Fall back to downloading the public CSV
    logger.info("Fetching Repurposing Hub from remote …")
    resp = get(REPURPOSING_HUB_URL, pause=1.0, timeout=60)
    if not resp:
        logger.warning("Could not fetch Repurposing Hub data")
        return None
    try:
        import io
        _REPURPOSING_DF = pd.read_csv(io.StringIO(resp.text), low_memory=False)
        return _REPURPOSING_DF
    except Exception as e:
        logger.warning("Failed to parse Repurposing Hub data: %s", e)
        return None


def query_repurposing_hub(inchikey: str) -> List[Dict]:
    df = _load_repurposing_hub()
    if df is None:
        return []

    # The Hub uses InChIKey in column 'InChIKey' (or similar variants)
    key_col = next(
        (c for c in df.columns if "inchi" in c.lower()), None
    )
    if not key_col:
        logger.warning("No InChIKey column found in Repurposing Hub table")
        return []

    matches = df[df[key_col].str.upper() == inchikey.upper()]
    results = []
    for _, row in matches.iterrows():
        results.append(
            {
                "source": "RepurposingHub",
                "record_type": "drug_annotation",
                "drug_name": row.get("pert_iname") or row.get("Name") or "",
                "moa": row.get("moa") or row.get("MOA") or "",
                "target": row.get("target") or row.get("Target") or "",
                "clinical_phase": row.get("clinical_phase") or "",
                "disease_area": row.get("disease_area") or "",
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

        dc = query_drugcentral(inchikey)
        rh = query_repurposing_hub(inchikey)

        rec = {
            "inchikey": inchikey,
            "std_smiles": mol.get("std_smiles"),
            "drugcentral": dc,
            "repurposing_hub": rh,
            "is_known_drug": dc["drug_id"] is not None or len(rh) > 0,
            "n_indications": len(dc.get("indications", [])),
            "n_dc_targets": len(dc.get("targets", [])),
            "n_rh_annotations": len(rh),
        }
        output.append(rec)
        logger.info(
            "  %s → DrugCentral: %s | RepurposingHub annotations: %d",
            inchikey[:14],
            "found" if dc["drug_id"] else "not found",
            len(rh),
        )

    logger.info(
        "Step 3 complete: %d molecules, %d known drugs",
        len(output),
        sum(1 for r in output if r["is_known_drug"]),
    )
    return output
