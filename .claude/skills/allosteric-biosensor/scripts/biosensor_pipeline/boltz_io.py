"""
boltz_io.py -- build Boltz-2.1 structure-and-binding payloads and parse results.

The actual paid calls go through the Boltz remote MCP tools
(``mcp__Boltz_API__boltz_estimate_structure_and_binding`` /
``...boltz_start_structure_and_binding``).  This module only *constructs* the
JSON ``input`` payloads (so they are reviewable and reproducible) and parses
the returned confidence / affinity metrics.  No secrets, no network here.

Predictions we use per construct:
  * holo : chimera protein (chain A) + analyte ligand (chain L), with
           ligand_protein_binding on chain L  -> fold + ligand binding
  * apo  : chimera protein only (chain A)                          -> fold of the OFF-state scaffold
  * ctrl : native receptor (chain A) + analyte ligand (chain L)    -> positive control that the binder+ligand pair is modellable
"""

from __future__ import annotations


def holo_payload(chimera_seq: str, ligand_smiles: str, *, num_samples: int = 1,
                 msa_empty: bool = True) -> dict:
    """Chimera + ligand, binding scored on the ligand chain."""
    prot = {"type": "protein", "chain_ids": ["A"], "value": chimera_seq}
    if msa_empty:
        prot["msa"] = {"type": "empty"}   # de-novo/engineered chain: single-sequence mode
    return {
        "input": {
            "entities": [
                prot,
                {"type": "ligand_smiles", "chain_ids": ["L"], "value": ligand_smiles},
            ],
            "binding": {"type": "ligand_protein_binding", "binder_chain_id": "L"},
            "num_samples": num_samples,
        }
    }


def apo_payload(chimera_seq: str, *, num_samples: int = 1, msa_empty: bool = True) -> dict:
    """Chimera protein only (OFF-state scaffold fold)."""
    prot = {"type": "protein", "chain_ids": ["A"], "value": chimera_seq}
    if msa_empty:
        prot["msa"] = {"type": "empty"}
    return {"input": {"entities": [prot], "num_samples": num_samples}}


def control_payload(receptor_seq: str, ligand_smiles: str, *, num_samples: int = 1,
                    msa_empty: bool = True) -> dict:
    """Native receptor + ligand -- positive control for the binder/ligand pair."""
    return holo_payload(receptor_seq, ligand_smiles, num_samples=num_samples, msa_empty=msa_empty)


# ---------------------------------------------------------------------------
# result parsing
# ---------------------------------------------------------------------------
def parse_confidence(result: dict) -> dict:
    """Pull the confidence + affinity fields out of a Boltz result blob.

    Boltz returns per-model confidence (confidence_score / ptm / iptm /
    complex_plddt) and, when binding is requested, affinity metrics
    (affinity_pred_value, affinity_probability_binary).  Field names vary
    slightly across API versions, so we search defensively.
    """
    out = {}
    blob = result or {}

    def dig(d, keys):
        for k in keys:
            if isinstance(d, dict) and k in d and d[k] is not None:
                return d[k]
        return None

    # confidence can live at top level, under 'confidence', or per-model
    cand = [blob, blob.get("confidence", {}) if isinstance(blob, dict) else {}]
    models = blob.get("models") or blob.get("predictions") or []
    if models and isinstance(models, list):
        cand.append(models[0])
        cand.append(models[0].get("confidence", {}) if isinstance(models[0], dict) else {})

    for c in cand:
        if not isinstance(c, dict):
            continue
        out.setdefault("confidence_score", dig(c, ["confidence_score", "confidence", "aggregate_score"]))
        out.setdefault("ptm", dig(c, ["ptm", "pTM"]))
        out.setdefault("iptm", dig(c, ["iptm", "ipTM"]))
        out.setdefault("complex_plddt", dig(c, ["complex_plddt", "plddt", "mean_plddt"]))
        out.setdefault("affinity_pred_value", dig(c, ["affinity_pred_value", "affinity", "affinity_value"]))
        out.setdefault("affinity_probability_binary",
                       dig(c, ["affinity_probability_binary", "affinity_probability", "binding_probability"]))

    return {k: v for k, v in out.items() if v is not None}
