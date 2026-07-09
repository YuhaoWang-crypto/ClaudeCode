"""Parallel orchestrator: push one protein through many models at once."""

from __future__ import annotations

import os
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import List, Optional

from . import fasta
from .models import EpitopeRecord, ModelSpec


@dataclass
class ModelResult:
    model: str
    category: str
    status: str                     # ok | not_installed | error | needs_structure
    records: List[EpitopeRecord] = field(default_factory=list)
    message: str = ""
    command: str = ""
    duration_s: float = 0.0


def _run_biolib(spec: ModelSpec, fasta_path: str, structure: Optional[str],
                workdir: str, timeout: int, t0: float) -> ModelResult:
    """Dispatch a model to the BioLib cloud via pybiolib."""
    from . import biolib_backend as bb

    if spec.needs_structure and not structure:
        return ModelResult(spec.name, spec.category, "needs_structure",
                           message=(f"{spec.app_uri} needs a PDB structure; "
                                    "pass --structure to enable"))
    model_wd = os.path.join(workdir, spec.name.replace("/", "_") + "_biolib")
    os.makedirs(model_wd, exist_ok=True)
    args = spec.build_args(spec, fasta_path, structure, model_wd)
    try:
        outdir = bb.run_app(spec.app_uri, args, model_wd, timeout=timeout)
    except bb.BiolibNotAvailable as exc:
        return ModelResult(spec.name, spec.category, "not_installed",
                           message=str(exc), command=f"{spec.app_uri} {args}")
    except Exception as exc:
        return ModelResult(spec.name, spec.category, "error",
                           message=f"biolib run failed: {str(exc)[:300]}",
                           command=f"{spec.app_uri} {args}",
                           duration_s=time.time() - t0)
    try:
        recs = spec.parse_files(spec, outdir, [])
    except Exception as exc:
        return ModelResult(spec.name, spec.category, "error",
                           message=f"parse failed: {exc}",
                           command=f"{spec.app_uri} {args}",
                           duration_s=time.time() - t0)
    return ModelResult(spec.name, spec.category, "ok", records=recs,
                       message=f"{len(recs)} records (biolib cloud)",
                       command=f"{spec.app_uri} {args}",
                       duration_s=time.time() - t0)


def _run_one(spec: ModelSpec, fasta_path: str, alleles: List[str],
             lengths, workdir: str, binary_override: Optional[str],
             timeout: int, structure: Optional[str] = None) -> ModelResult:
    t0 = time.time()

    if spec.backend == "biolib":
        return _run_biolib(spec, fasta_path, structure, workdir, timeout, t0)

    if spec.needs_structure:
        return ModelResult(spec.name, spec.category, "needs_structure",
                           message=("requires a 3D structure (PDB), not a bare "
                                    "FASTA; supply --structure to enable"))

    binary = spec.resolve_binary(binary_override)
    if not binary:
        return ModelResult(
            spec.name, spec.category, "not_installed",
            message=f"binary '{binary_override or spec.binary}' not on PATH. "
                    f"{spec.install_hint} ({spec.homepage})")

    if spec.needs_alleles and not alleles:
        return ModelResult(spec.name, spec.category, "error",
                           message="model needs HLA alleles but none supplied")

    model_wd = os.path.join(workdir, spec.name.replace("/", "_"))
    os.makedirs(model_wd, exist_ok=True)
    cmd = spec.build_command(spec, fasta_path,
                             alleles if spec.needs_alleles else [],
                             lengths, model_wd)
    cmd[0] = binary  # honor resolved absolute path / override
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout, cwd=model_wd)
    except subprocess.TimeoutExpired:
        return ModelResult(spec.name, spec.category, "error",
                           message=f"timed out after {timeout}s",
                           command=" ".join(cmd), duration_s=time.time() - t0)
    except OSError as exc:
        return ModelResult(spec.name, spec.category, "error",
                           message=str(exc), command=" ".join(cmd))

    if proc.returncode != 0 and not proc.stdout.strip():
        return ModelResult(
            spec.name, spec.category, "error",
            message=f"exit {proc.returncode}: {proc.stderr.strip()[:400]}",
            command=" ".join(cmd), duration_s=time.time() - t0)

    try:
        recs = spec.parse(spec, proc.stdout, [])
    except Exception as exc:  # parsing must never crash the whole run
        return ModelResult(spec.name, spec.category, "error",
                           message=f"parse failed: {exc}",
                           command=" ".join(cmd), duration_s=time.time() - t0)

    return ModelResult(spec.name, spec.category, "ok", records=recs,
                       message=f"{len(recs)} records",
                       command=" ".join(cmd), duration_s=time.time() - t0)


def run_models(specs: List[ModelSpec], fasta_path: str, alleles_i: List[str],
               alleles_ii: List[str], lengths_i, lengths_ii,
               workdir: str, binary_overrides: dict,
               max_workers: int = 8, timeout: int = 3600,
               structure: Optional[str] = None) -> List[ModelResult]:
    """Run every model concurrently and collect ModelResults.

    MHC-I and MHC-II models receive their respective allele lists.
    """
    fasta_path = os.path.abspath(fasta_path)
    if structure:
        structure = os.path.abspath(structure)
    records = fasta.parse_fasta(fasta_path)
    warnings = fasta.validate(records)
    if warnings:
        print("[fasta] warnings:")
        for w in warnings:
            print(f"  - {w}")

    os.makedirs(workdir, exist_ok=True)
    results: List[ModelResult] = []
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futs = {}
        for spec in specs:
            if spec.category == "mhc_ii":
                alleles, lengths = alleles_ii, lengths_ii
            elif spec.category == "mhc_i":
                alleles, lengths = alleles_i, lengths_i
            else:
                alleles, lengths = [], None
            override = binary_overrides.get(spec.binary)
            fut = pool.submit(_run_one, spec, fasta_path, alleles, lengths,
                              workdir, override, timeout, structure)
            futs[fut] = spec.name
        for fut in as_completed(futs):
            results.append(fut.result())

    # stable ordering by category then name for reproducible reports
    order = {"mhc_i": 0, "mhc_ii": 1, "bcell_linear": 2,
             "bcell_conf": 3, "aux": 4}
    results.sort(key=lambda r: (order.get(r.category, 9), r.model))
    return results
