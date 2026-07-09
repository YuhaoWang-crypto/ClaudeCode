"""Integration layer: unify the `immunogenicity-multimodel` skill backends
(IEDB NetMHCpan MHC-I cloud + BioLib DTU/ImmunoGeNN MHC-II) with this package's
own DTU adapters, into one normalized EpitopeRecord table.

Why a bridge: the skill runs the two backends that this package could not reach
locally — IEDB cloud NetMHCpan (MHC-I, no license, no token) and DTU/ImmunoGeNN
(MHC-II population immunogenicity). Its results come back as a pandas DataFrame
with columns (id, source, model, allele, length, peptide, position,
affinity_nM, percentile_rank, score, is_binder). We map each row onto this
package's `EpitopeRecord` so the skill's output flows through the same
`aggregate.consensus` / report machinery as the native adapters.
"""

from __future__ import annotations

import os
import sys
from typing import List, Optional, Sequence

from .models import EpitopeRecord
from .runner import ModelResult


# --------------------------------------------------------------------------- #
# Locate and import the installed skill kernel
# --------------------------------------------------------------------------- #
def _load_kernel(skill_dir: Optional[str] = None):
    """Import the skill's kernel.py. Searches common install locations."""
    candidates = [skill_dir] if skill_dir else []
    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(os.path.dirname(here))          # .../ClaudeCode
    candidates += [
        os.path.join(repo, ".claude", "skills", "immunogenicity-multimodel"),
        os.path.expanduser("~/.claude/skills/immunogenicity-multimodel"),
    ]
    for c in candidates:
        if c and os.path.isfile(os.path.join(c, "kernel.py")):
            if c not in sys.path:
                sys.path.insert(0, c)
            import kernel  # noqa
            return kernel
    raise ImportError(
        "immunogenicity-multimodel kernel.py not found. Install the skill under "
        ".claude/skills/immunogenicity-multimodel/ or pass skill_dir=...")


# --------------------------------------------------------------------------- #
# DataFrame(skill schema) -> EpitopeRecord list
# --------------------------------------------------------------------------- #
_CATEGORY = {"immunogenn": "mhc_ii", "netmhcpan_el": "mhc_i"}

# Map the skill's binder vocabulary onto this package's (SB/WB) convention so
# skill rows flow through aggregate.consensus unchanged.
_BIND_MAP = {"strong": "SB", "weak": "WB", "immunogenic": "SB",
             "none": "", "": ""}


def _df_to_records(df, model_name: str, category: str) -> List[EpitopeRecord]:
    import pandas as pd
    out: List[EpitopeRecord] = []
    for _, r in df.iterrows():
        rank = r.get("percentile_rank")
        rank = float(rank) if pd.notna(rank) else None
        pos = r.get("position")
        raw_level = str(r.get("is_binder", "")).strip()
        out.append(EpitopeRecord(
            model=model_name, category=category,
            seq_id=str(r.get("id", "")), peptide=str(r.get("peptide", "")),
            position=int(pos) + 1 if pd.notna(pos) else None,   # -> 1-indexed
            allele=r.get("allele"),
            score=float(r["score"]) if pd.notna(r.get("score")) else None,
            rank=rank,
            affinity_nM=float(r["affinity_nM"]) if pd.notna(r.get("affinity_nM")) else None,
            bind_level=_BIND_MAP.get(raw_level, raw_level),
        ))
    return out


# --------------------------------------------------------------------------- #
# Public entry points
# --------------------------------------------------------------------------- #
def run_iedb_mhci(fasta_path: str, alleles: Sequence[str], lengths=(9,),
                  skill_dir: Optional[str] = None) -> ModelResult:
    """MHC-I via IEDB cloud NetMHCpan (no token). Returns a ModelResult so it
    plugs into aggregate/report exactly like a native adapter."""
    import time
    k = _load_kernel(skill_dir)
    seqs = k.read_fasta(fasta_path)
    t0 = time.time()
    try:
        df = k.run_iedb_mhci(seqs, list(alleles), tuple(lengths))
    except Exception as exc:
        return ModelResult("NetMHCpan(IEDB)", "mhc_i", "error",
                           message=str(exc)[:300])
    recs = _df_to_records(df, "NetMHCpan(IEDB)", "mhc_i")
    return ModelResult("NetMHCpan(IEDB)", "mhc_i", "ok", records=recs,
                       message=f"{len(recs)} peptide×allele (IEDB cloud)",
                       command="IEDB tools_api/mhci netmhcpan_el",
                       duration_s=time.time() - t0)


def run_immunogenn(fasta_path: str, skill_dir: Optional[str] = None,
                   timeout: int = 1200) -> ModelResult:
    """MHC-II population immunogenicity via BioLib DTU/ImmunoGeNN (anonymous).

    NOTE: anonymous BioLib jobs do not complete in a non-interactive sandbox
    (they submit but never leave 'in_progress'); set BIOLIB_TOKEN to enable.
    Returns a ModelResult; status 'error' with a clear message on non-completion.
    """
    import time
    k = _load_kernel(skill_dir)
    t0 = time.time()
    try:
        _, df = k.run_immunogenn(fasta_path, mode="screen", timeout=timeout)
    except Exception as exc:
        return ModelResult("ImmunoGeNN", "mhc_ii", "error",
                           message=str(exc)[:300])
    if df is None:
        return ModelResult("ImmunoGeNN", "mhc_ii", "error",
                           message="job did not complete (set BIOLIB_TOKEN; "
                                   "anonymous cloud jobs stall in this sandbox)",
                           duration_s=time.time() - t0)
    recs = _df_to_records(df, "ImmunoGeNN", "mhc_ii")
    return ModelResult("ImmunoGeNN", "mhc_ii", "ok", records=recs,
                       message=f"{len(recs)} 15mers (BioLib)",
                       duration_s=time.time() - t0)


def run_integrated(fasta_path: str, alleles_i: Sequence[str], lengths_i=(9,),
                   use=("iedb", "immunogenn"),
                   skill_dir: Optional[str] = None) -> List[ModelResult]:
    """Run the skill backends and return ModelResults ready for aggregate.*.

    Combine with native adapters via `runner.run_models(...)` and concatenate
    the result lists before calling `aggregate.consensus`.
    """
    results: List[ModelResult] = []
    if "iedb" in use:
        results.append(run_iedb_mhci(fasta_path, alleles_i, lengths_i, skill_dir))
    if "immunogenn" in use:
        results.append(run_immunogenn(fasta_path, skill_dir))
    return results
