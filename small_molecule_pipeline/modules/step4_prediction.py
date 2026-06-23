"""Step 4 — Predicted target identification.

Sources:
  - SwissTargetPrediction  http://www.swisstargetprediction.ch/
    (reverse-engineered public HTTP endpoint; no formal public API key needed)
  - DrugCLIP               local embedding-based model (optional)
    Set DRUGCLIP_MODEL_PATH env var to enable.

SwissTargetPrediction returns a ranked list of human targets with
probability scores estimated from chemical similarity to known ligands.
"""

import logging
import os
from typing import Dict, List, Optional

from utils.api_client import get, post

logger = logging.getLogger(__name__)

STP_BASE = "http://www.swisstargetprediction.ch"
STP_ORGANISM = "Homo sapiens"


# ── SwissTargetPrediction ────────────────────────────────────────────────────

def query_swisstargetprediction(smiles: str, organism: str = STP_ORGANISM) -> List[Dict]:
    """Submit SMILES to STP and parse returned predictions."""
    # STP uses a form POST then redirects to a result page
    resp = post(
        f"{STP_BASE}/target_predictions",
        data={"smiles": smiles, "organism": organism},
        pause=1.0,
        timeout=60,
    )
    if not resp:
        return []

    # STP returns JSON when Accept: application/json is set; fall back to HTML parse
    try:
        resp_json = resp.json()
    except Exception:
        resp_json = None

    if resp_json:
        predictions = resp_json if isinstance(resp_json, list) else resp_json.get("predictions", [])
        return [
            {
                "source": "SwissTargetPrediction",
                "target_name": p.get("common_name") or p.get("target"),
                "uniprot": p.get("uniprot_id"),
                "probability": p.get("probability"),
                "known_actives": p.get("known_actives"),
                "target_class": p.get("target_class"),
                "organism": organism,
            }
            for p in predictions
            if p.get("probability", 0) > 0.1
        ]

    # Minimal HTML fallback: look for JSON embedded in the page
    import re
    text = resp.text
    match = re.search(r"var\s+predictions\s*=\s*(\[.*?\]);", text, re.DOTALL)
    if match:
        import json
        try:
            raw = json.loads(match.group(1))
            return [
                {
                    "source": "SwissTargetPrediction",
                    "target_name": r.get("common_name") or r.get("target"),
                    "uniprot": r.get("uniprot_id"),
                    "probability": r.get("probability"),
                    "known_actives": r.get("known_actives"),
                    "target_class": r.get("target_class"),
                    "organism": organism,
                }
                for r in raw
                if r.get("probability", 0) > 0.1
            ]
        except Exception as e:
            logger.debug("STP HTML parse failed: %s", e)

    logger.warning("SwissTargetPrediction: could not parse response for SMILES %s…", smiles[:30])
    return []


# ── DrugCLIP (optional local model) ─────────────────────────────────────────

def _drugclip_available() -> bool:
    path = os.environ.get("DRUGCLIP_MODEL_PATH", "")
    if not path:
        return False
    try:
        import torch  # noqa: F401
        return True
    except ImportError:
        logger.warning("DrugCLIP requested but torch not installed")
        return False


def query_drugclip(smiles: str) -> List[Dict]:
    """
    Embed SMILES with DrugCLIP and return nearest protein targets.

    Expects:
      DRUGCLIP_MODEL_PATH   – path to saved DrugCLIP checkpoint
      DRUGCLIP_TARGET_DB    – path to pre-embedded target database (.pt)
    """
    if not _drugclip_available():
        return []

    import torch

    model_path = os.environ["DRUGCLIP_MODEL_PATH"]
    target_db_path = os.environ.get("DRUGCLIP_TARGET_DB", "")

    try:
        from drugclip import DrugCLIP  # type: ignore
    except ImportError:
        logger.warning("drugclip package not installed; skipping DrugCLIP")
        return []

    try:
        model = DrugCLIP.load_from_checkpoint(model_path)
        model.eval()

        drug_emb = model.encode_smiles([smiles])  # (1, D)

        if not target_db_path:
            logger.warning("DRUGCLIP_TARGET_DB not set; skipping target lookup")
            return []

        db = torch.load(target_db_path)  # {"embeddings": Tensor, "targets": [...]}
        target_embs = db["embeddings"]  # (N, D)
        targets_meta = db["targets"]

        sims = torch.nn.functional.cosine_similarity(drug_emb, target_embs)
        top_k = min(20, len(targets_meta))
        top_idx = sims.topk(top_k).indices.tolist()

        return [
            {
                "source": "DrugCLIP",
                "target_name": targets_meta[i].get("name"),
                "uniprot": targets_meta[i].get("uniprot"),
                "similarity": round(float(sims[i]), 4),
                "target_class": targets_meta[i].get("class"),
                "organism": targets_meta[i].get("organism"),
            }
            for i in top_idx
            if float(sims[i]) > 0.5
        ]
    except Exception as e:
        logger.error("DrugCLIP inference failed: %s", e)
        return []


# ── Aggregator ───────────────────────────────────────────────────────────────

def run(step1_results: List[Dict]) -> List[Dict]:
    drugclip_on = _drugclip_available()
    if drugclip_on:
        logger.info("DrugCLIP enabled (model: %s)", os.environ.get("DRUGCLIP_MODEL_PATH"))
    else:
        logger.info("DrugCLIP disabled (set DRUGCLIP_MODEL_PATH to enable)")

    output = []
    for mol in step1_results:
        if not mol.get("passed"):
            continue
        smiles = mol.get("std_smiles", "")
        inchikey = mol.get("inchikey", "")
        if not smiles:
            continue

        stp_preds = query_swisstargetprediction(smiles)
        dc_preds = query_drugclip(smiles) if drugclip_on else []

        all_preds = stp_preds + dc_preds
        # Sort by probability / similarity descending
        all_preds.sort(
            key=lambda x: float(x.get("probability") or x.get("similarity") or 0),
            reverse=True,
        )

        output.append(
            {
                "inchikey": inchikey,
                "std_smiles": smiles,
                "predicted_targets": all_preds,
                "n_stp": len(stp_preds),
                "n_drugclip": len(dc_preds),
                "top_target": all_preds[0].get("target_name") if all_preds else None,
            }
        )
        logger.info(
            "  %s → STP: %d | DrugCLIP: %d predicted targets",
            inchikey[:14],
            len(stp_preds),
            len(dc_preds),
        )

    logger.info(
        "Step 4 complete: %d molecules with predicted targets", len(output)
    )
    return output
