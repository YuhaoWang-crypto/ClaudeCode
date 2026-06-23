"""Step 8 — ADMET, toxicity, and assay-interference risk.

Sources:
  - ADMETlab 2.0   https://admetmesh.scbdd.com/  (POST SMILES, no key required)
  - CompTox (EPA)  https://comptox.epa.gov/dashboard/api/
    Set COMPTOX_API_KEY env var for authenticated access.
  - ProTox-3.0     https://tox.charite.de/protox3/  (POST endpoint)
"""

import logging
import os
from typing import Dict, List, Optional

from utils.api_client import get, post

logger = logging.getLogger(__name__)

ADMETLAB_BASE = "https://admetmesh.scbdd.com"
COMPTOX_BASE = "https://comptoxdashboard.epa.gov/api/chemical/detail"
PROTOX_BASE = "https://tox.charite.de/protox3/prediction/result"


# ── ADMETlab 2.0 ─────────────────────────────────────────────────────────────

def query_admetlab(smiles: str) -> Optional[Dict]:
    """Submit SMILES to ADMETlab and return predicted ADMET properties."""
    resp = post(
        f"{ADMETLAB_BASE}/service/predict",
        json={"smiles": smiles},
        headers={"Content-Type": "application/json"},
        pause=2.0,
        timeout=120,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        result = data.get("result") or data
        return {
            "source": "ADMETlab2",
            # Absorption
            "Caco2_permeability": result.get("Caco2"),
            "HIA": result.get("HIA"),
            "oral_bioavailability": result.get("F20") or result.get("Bioavailability"),
            "pgp_substrate": result.get("Pgp-substrate"),
            "pgp_inhibitor": result.get("Pgp-inhibitor"),
            # Distribution
            "VD": result.get("VD"),
            "fu": result.get("Fu"),
            "BBB": result.get("BBB"),
            # Metabolism
            "CYP1A2_inhibitor": result.get("CYP1A2-inhibitor"),
            "CYP2C9_inhibitor": result.get("CYP2C9-inhibitor"),
            "CYP2C19_inhibitor": result.get("CYP2C19-inhibitor"),
            "CYP2D6_inhibitor": result.get("CYP2D6-inhibitor"),
            "CYP3A4_inhibitor": result.get("CYP3A4-inhibitor"),
            # Excretion
            "t_half": result.get("T12") or result.get("Half-life"),
            "clearance": result.get("CL"),
            # Toxicity
            "hERG_inhibitor": result.get("hERG"),
            "DILI": result.get("DILI"),
            "ames_mutagenicity": result.get("Ames"),
            "FDAMDD": result.get("FDAMDD"),
            "skin_sensitization": result.get("SkinSen"),
            "carcinogenicity": result.get("Carcinogenicity"),
            "rat_oral_ld50": result.get("Rat-oral-LD50"),
        }
    except Exception as e:
        logger.warning("ADMETlab parse error for SMILES %s…: %s", smiles[:20], e)
        return None


# ── CompTox (EPA) ─────────────────────────────────────────────────────────────

def query_comptox(inchikey: str) -> Optional[Dict]:
    """Query EPA CompTox Dashboard for hazard data."""
    api_key = os.environ.get("COMPTOX_API_KEY", "")
    headers = {}
    if api_key:
        headers["x-api-key"] = api_key

    resp = get(
        f"https://comptox.epa.gov/dashboard-api/ccdapp1/chemical-files/inchikey/{inchikey}",
        headers=headers,
        pause=0.5,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        return {
            "source": "CompTox",
            "dtxsid": data.get("dtxsid"),
            "casrn": data.get("casrn"),
            "preferred_name": data.get("preferredName"),
            "mol_weight": data.get("molWeight"),
            "is_qsar_ready": data.get("isQsarReady"),
            "has_experimental_data": data.get("hasExpData"),
            "exposure_flags": data.get("exposureFlags", []),
            "toxicity_data": {
                "ToxCast_active": data.get("toxcastActive"),
                "tox21_active": data.get("tox21Active"),
                "in_vitro_active": data.get("inVitroActive"),
            },
        }
    except Exception as e:
        logger.debug("CompTox error for %s: %s", inchikey, e)
        return None


# ── ProTox-3.0 ────────────────────────────────────────────────────────────────

def query_protox(smiles: str) -> Optional[Dict]:
    """Predict toxicity endpoints via ProTox-3.0."""
    resp = post(
        "https://tox.charite.de/protox3/api/predict",
        json={"smiles": smiles},
        headers={"Content-Type": "application/json"},
        pause=2.0,
        timeout=120,
    )
    if not resp:
        return None
    try:
        data = resp.json()
        preds = data.get("predictions") or data
        return {
            "source": "ProTox3",
            "ld50_predicted_mg_kg": preds.get("LD50"),
            "toxicity_class": preds.get("Class"),
            "hepatotoxicity": preds.get("Hepatotoxicity"),
            "immunotoxicity": preds.get("Immunotoxicity"),
            "neurotoxicity": preds.get("Neurotoxicity"),
            "carcinogenicity": preds.get("Carcinogenicity"),
            "mutagenicity": preds.get("Mutagenicity"),
            "cytotoxicity": preds.get("Cytotoxicity"),
            "NR_AR": preds.get("NR-AR"),
            "NR_ER": preds.get("NR-ER"),
            "SR_p53": preds.get("SR-p53"),
        }
    except Exception as e:
        logger.warning("ProTox error for SMILES %s…: %s", smiles[:20], e)
        return None


# ── Interference / assay-flag summary ────────────────────────────────────────

def _flag_interference(mol_result: Dict, admet: Optional[Dict]) -> Dict:
    """Compile a risk summary from upstream alert data and ADMET predictions."""
    flags = []

    # From step1 structural alerts
    alert_flags = mol_result.get("alert_flags", [])
    if "PAINS" in alert_flags:
        flags.append("PAINS_alert")
    if "Brenk" in alert_flags:
        flags.append("Brenk_alert")
    if "NIH" in alert_flags:
        flags.append("NIH_alert")

    if admet:
        if admet.get("hERG_inhibitor") in {"Positive", True, 1, "1"}:
            flags.append("hERG_risk")
        if admet.get("ames_mutagenicity") in {"Positive", True, 1, "1"}:
            flags.append("AMES_mutagen")
        if admet.get("DILI") in {"Positive", True, 1, "1"}:
            flags.append("DILI_risk")
        if admet.get("carcinogenicity") in {"Positive", True, 1, "1"}:
            flags.append("carcinogen")

    return {
        "interference_flags": flags,
        "n_flags": len(flags),
        "risk_level": (
            "high" if len(flags) >= 3
            else "medium" if len(flags) >= 1
            else "low"
        ),
    }


# ── Aggregator ───────────────────────────────────────────────────────────────

def run(step1_results: List[Dict]) -> List[Dict]:
    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        smiles = mol.get("std_smiles", "")
        inchikey = mol.get("inchikey", "")
        if not smiles:
            continue

        admet = query_admetlab(smiles)
        comptox = query_comptox(inchikey)
        protox = query_protox(smiles)
        interference = _flag_interference(mol, admet)

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": smiles,
                "admetlab": admet,
                "comptox": comptox,
                "protox": protox,
                "interference": interference,
                "risk_level": interference["risk_level"],
                "n_tox_flags": interference["n_flags"],
            }
        )
        logger.info(
            "  %s → ADMET: %s | CompTox: %s | ProTox LD50: %s | Risk: %s",
            inchikey[:14],
            "ok" if admet else "failed",
            "ok" if comptox else "not found",
            protox.get("ld50_predicted_mg_kg", "N/A") if protox else "failed",
            interference["risk_level"],
        )

    logger.info(
        "Step 8 complete: %d molecules | high-risk: %d | medium: %d | low: %d",
        len(output),
        sum(1 for r in output if r["risk_level"] == "high"),
        sum(1 for r in output if r["risk_level"] == "medium"),
        sum(1 for r in output if r["risk_level"] == "low"),
    )
    return output
