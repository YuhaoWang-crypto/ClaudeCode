"""Run provenance: what data, what predictor, what code, what evidence level.

Adopted from the audit discipline in the uploaded Neoantigen_Vaccine_Package,
which does this better than the first version of this pipeline did. A ranked
CSV with no record of which predictor produced it, which reference proteome it
was tiled against, or which weights were applied is not auditable, and a
personalized therapy pipeline that is not auditable is not usable.

Evidence levels attach to the *presentation* numbers, which are the ones a
reader is most likely to over-trust:

  E3  measured on this patient's own tumor (mass-spec immunopeptidomics)
  E2  a real, published predictor run on real sequence (NetMHCpan, MHCflurry)
  E1  a real predictor run on synthetic or fixture sequence
  E0  a placeholder/proxy -- no biological meaning, refuse to interpret

This package only ever produces E2 (or E1 if you feed it a fixture). It has no
E0 mode: there is no hash-based stand-in predictor, because a number that looks
like an affinity but is a hash gets quoted as an affinity eventually.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional, Sequence

from . import __version__


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: str) -> Optional[str]:
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def sha256_frame(df) -> str:
    """Content hash of a DataFrame, stable across column order."""
    import pandas as pd
    if df is None or not len(df):
        return hashlib.sha256(b"").hexdigest()
    payload = df.reindex(sorted(df.columns), axis=1).to_csv(index=False).encode()
    return hashlib.sha256(payload).hexdigest()


def git_revision() -> Optional[str]:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                             text=True, timeout=10,
                             cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        return out.stdout.strip() or None
    except Exception:                       # noqa: BLE001 - provenance is best effort
        return None


def capabilities() -> Dict[str, object]:
    """What is actually installed. Availability is not permission: a licensed
    binary on PATH says nothing about whether you may use it commercially."""
    import importlib.util
    import shutil
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "mhcflurry": bool(importlib.util.find_spec("mhcflurry")),
        "netMHCpan_binary": shutil.which("netMHCpan"),
        "netMHCIIpan_binary": shutil.which("netMHCIIpan"),
        "iedb_cloud_rest": True,
        "note": ("A binary on PATH establishes neither model-weight availability "
                 "nor a commercial licence."),
    }


def manifest(res: Dict[str, object], cfg, proteome_path: Optional[str] = None,
             input_path: Optional[str] = None,
             evidence_level: str = "E2") -> Dict[str, object]:
    """One JSON that lets someone else re-derive this run, or reject it."""
    preds = res.get("predictions")
    sel = res.get("selected")
    con = res.get("construct") or {}
    p = cfg.patient
    return {
        "generated_at": now(),
        "package_version": __version__,
        "git_revision": git_revision(),
        "research_use_only": True,
        "not_a_reproduction_of": ("intismeran autogene / mRNA-4157 / V940 "
                                  "proprietary scoring or construct rules"),
        "patient": {
            "patient_id": p.patient_id,
            "hla_class1": list(p.hla_class1),
            "hla_class2": list(p.hla_class2),
            "hla_source": "supplied by caller -- NOT predicted by this package",
            "tumor_purity": p.tumor_purity,
            "mhc1_lengths": list(p.mhc1_lengths),
        },
        "inputs": {
            "variant_table_sha256": sha256_frame(res.get("variants")),
            "variant_input_path": input_path,
            "variant_input_sha256": sha256_file(input_path) if input_path else None,
            "reference_proteome_sha256": sha256_file(proteome_path) if proteome_path else None,
            "iedb_reference_epitopes": res.get("n_iedb_reference_epitopes"),
        },
        "presentation": {
            "predictor": cfg.mhc1_method,
            "backend": "IEDB tools-cluster REST",
            "evidence_level": evidence_level,
            "peptide_allele_predictions": int(len(preds)) if preds is not None else 0,
            "predictions_sha256": sha256_frame(preds),
        },
        "scoring": {
            "weights": cfg.weight_dict(),
            "gates": cfg.gates.__dict__,
            "selection_rules": cfg.selection.__dict__,
        },
        "outputs": {
            "n_selected": int(len(sel)) if sel is not None else 0,
            "selected_sha256": sha256_frame(sel),
            "construct_cds_nt": len(con.get("cds", "")) if con else 0,
            "construct_qc_pass": (con.get("qc", {}) or {}).get("pass"),
            "junction_binders_flagged": (
                int(con["junction_scan"]["flagged"].sum())
                if con and hasattr(con.get("junction_scan"), "empty")
                and not con["junction_scan"].empty else None),
        },
        "capabilities": capabilities(),
        "limitations": [
            "Every immunogenicity number here is an unverified prediction.",
            "HLA typing is an input; a wrong allele invalidates the ranking silently.",
            "The composite score is not validated and does not beat the binding "
            "predictor on either shipped benchmark.",
            "The CDS is a research-use cassette draft: no UTRs, cap, poly(A), "
            "modified nucleosides, LNP formulation or release specifications.",
            "Frameshift / neo-ORF and fusion neoantigens need transcript-level "
            "annotation and are reported as skipped rather than invented.",
        ],
    }


def write_manifest(res: Dict[str, object], cfg, outdir: str, **kw) -> str:
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, "audit_manifest.json")
    with open(path, "w") as fh:
        json.dump(manifest(res, cfg, **kw), fh, indent=2, default=str)
    return path
