"""Shared helpers: config loading, logging, unit conversions, platform choice."""
from __future__ import annotations

import logging
import math
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

# Gas constant in kcal / (mol K)
R_KCAL = 1.987204259e-3

REPO_ROOT = Path(__file__).resolve().parents[1]


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=os.environ.get("CDPFAS_LOGLEVEL", "INFO"),
        format="%(asctime)s | %(levelname)-7s | %(name)s | %(message)s",
    )
    return logging.getLogger(name)


log = get_logger("cd_pfas.utils")


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.is_absolute():
        path = REPO_ROOT / path
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    log.info("loaded config: %s", path)
    return cfg


def resolve(path: str | Path) -> Path:
    """Resolve a possibly-repo-relative path to an absolute Path."""
    p = Path(path)
    return p if p.is_absolute() else (REPO_ROOT / p)


def require_tool(binary: str) -> str:
    """Fail early with a clear message if an external binary is missing."""
    found = shutil.which(binary)
    if not found:
        raise RuntimeError(
            f"required external tool '{binary}' not found on PATH. "
            "Did you `conda activate cd-pfas-md`?  See environment.yml."
        )
    return found


def run(cmd: list[str], cwd: Path | None = None, check: bool = True) -> subprocess.CompletedProcess:
    log.info("$ %s", " ".join(str(c) for c in cmd))
    return subprocess.run(
        [str(c) for c in cmd], cwd=cwd, check=check,
        capture_output=True, text=True,
    )


def ka_to_dg(ka_m_inv: float, temperature_k: float, standard_conc_m: float = 1.0) -> float:
    """Convert an association constant Ka (1/M) to ΔG_bind (kcal/mol) at 1 M standard state.

    ΔG = -RT ln(Ka * C°)   with C° = 1 M   =>   ΔG = -RT ln(Ka[1/M])
    """
    ka_m_inv = float(ka_m_inv)  # YAML may hand us '3.9e3' as a string
    return -R_KCAL * temperature_k * math.log(ka_m_inv * standard_conc_m)


def dg_to_ka(dg_kcal: float, temperature_k: float) -> float:
    return math.exp(-dg_kcal / (R_KCAL * temperature_k))


@dataclass
class Guest:
    key: str
    name: str
    smiles: str | None
    mol2: str | None
    net_charge: int
    charge_method: str
    fluorine_policy: str
    exp_dg: float | None
    structure: str | None = None   # path to a 3D structure file (sdf/pdb/mol2)

    @classmethod
    def from_config(cls, key: str, gcfg: dict, temperature_k: float) -> "Guest":
        exp = gcfg.get("experimental", {}) or {}
        exp_dg = exp.get("dG_kcal_per_mol")
        if exp_dg is None and exp.get("Ka_M_inv"):
            exp_dg = ka_to_dg(exp["Ka_M_inv"], temperature_k)
        return cls(
            key=key,
            name=gcfg.get("name", key),
            smiles=gcfg.get("smiles"),
            mol2=gcfg.get("mol2"),
            net_charge=int(gcfg.get("net_charge", 0)),
            charge_method=gcfg.get("charge_method", "bcc"),
            fluorine_policy=gcfg.get("fluorine_policy", "standard"),
            exp_dg=exp_dg,
            structure=gcfg.get("structure"),
        )


def pick_platform(requested: str):
    """Return an OpenMM Platform, honoring config but degrading gracefully."""
    from openmm import Platform

    order = (
        [requested] if requested and requested != "auto"
        else ["CUDA", "OpenCL", "CPU"]
    )
    for name in order:
        try:
            plat = Platform.getPlatformByName(name)
            log.info("using OpenMM platform: %s", name)
            return plat
        except Exception:  # noqa: BLE001
            continue
    raise RuntimeError(f"no usable OpenMM platform among {order}")
