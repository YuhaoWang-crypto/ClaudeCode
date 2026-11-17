"""Real Folder adapter around Boltz-2 (MIT, commercial-friendly).

Boltz runs as a CLI: you write a YAML describing the chains, call
``boltz predict`, and read back predicted structures (mmCIF) plus a confidence
JSON (plddt / pae / pde / iptm). This adapter wires that up. It is NOT exercised
by the smoke test (needs the model + weights); the integration points are
marked clearly and the import/CLI presence is checked lazily.

Install (see https://github.com/jwohlwend/boltz):
    pip install boltz -U        # pulls torch; first run downloads weights

The pipeline only ever calls ``fold`` (unconditional refold for self-consistency).
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from ..types import ANTIGEN, HEAVY, LIGHT, Confidence, Design, DesignSpec, Pose


class BoltzFolder:
    def __init__(self, boltz_bin: str = "boltz", use_msa_server: bool = True,
                 workdir: Optional[str] = None):
        self.boltz_bin = boltz_bin
        self.use_msa_server = use_msa_server
        self.workdir = Path(workdir) if workdir else None

    def _check(self) -> None:
        if shutil.which(self.boltz_bin) is None:
            raise RuntimeError(
                f"'{self.boltz_bin}' not found on PATH. Install with `pip install boltz -U` "
                "or use MockFolder for a dry run."
            )

    def _write_yaml(self, design: Design, spec: DesignSpec, path: Path) -> None:
        """Boltz input schema: one entry per polymer chain."""
        import yaml  # local import; pyyaml is a light dep

        doc = {
            "version": 1,
            "sequences": [
                {"protein": {"id": HEAVY, "sequence": design.heavy_seq}},
                {"protein": {"id": LIGHT, "sequence": design.light_seq}},
                {"protein": {"id": ANTIGEN, "sequence": spec.antigen_seq}},
            ],
            # Ask for the confidence/affinity heads.
            "properties": [{"affinity": {"binder": HEAVY}}],
        }
        path.write_text(yaml.safe_dump(doc, sort_keys=False))

    def _parse_confidence(self, out_dir: Path) -> Confidence:
        # Boltz writes confidence_*.json alongside the predicted structure.
        hits = list(out_dir.rglob("confidence_*.json"))
        if not hits:
            raise RuntimeError(f"no confidence json under {out_dir}")
        data = json.loads(hits[0].read_text())
        return Confidence(
            plddt=float(data.get("complex_plddt", data.get("plddt", 0.0))),
            pae=float(data.get("complex_pae", data.get("pae", 0.0))),
            pde=float(data.get("complex_pde", data.get("pde", 0.0))),
            iptm=float(data.get("iptm", data.get("complex_iptm", 0.0))),
        )

    def _parse_ca(self, out_dir: Path) -> dict[str, np.ndarray]:
        """Parse C-alpha coords per chain from the first predicted mmCIF."""
        cifs = list(out_dir.rglob("*_model_0.cif")) or list(out_dir.rglob("*.cif"))
        if not cifs:
            raise RuntimeError(f"no predicted mmCIF under {out_dir}")
        return _read_cif_ca(cifs[0])

    def fold(self, design: Design, spec: DesignSpec, seed: Optional[int] = None) -> Pose:
        self._check()
        base = self.workdir or Path(tempfile.mkdtemp(prefix="boltz_"))
        base.mkdir(parents=True, exist_ok=True)
        yaml_path = base / f"{design.design_id}.yaml"
        self._write_yaml(design, spec, yaml_path)

        cmd = [self.boltz_bin, "predict", str(yaml_path),
               "--out_dir", str(base), "--output_format", "mmcif"]
        if self.use_msa_server:
            cmd.append("--use_msa_server")
        if seed is not None:
            cmd += ["--seed", str(seed)]
        subprocess.run(cmd, check=True)

        out_dir = base / f"boltz_results_{design.design_id}"
        return Pose(backbone=self._parse_ca(out_dir), confidence=self._parse_confidence(out_dir))


def _read_cif_ca(cif_path: Path) -> dict[str, np.ndarray]:
    """Very small mmCIF reader: pull CA atoms grouped by auth chain id.

    Good enough for backbone RMSD; not a general parser.
    """
    chains: dict[str, list[list[float]]] = {}
    lines = cif_path.read_text().splitlines()
    # Locate the atom_site loop and its column order.
    cols: list[str] = []
    in_loop = False
    for line in lines:
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1])
            in_loop = True
            continue
        if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
            f = line.split()
            if len(f) < len(cols):
                continue
            rec = dict(zip(cols, f))
            if rec.get("label_atom_id") != "CA":
                continue
            chain = rec.get("auth_asym_id", rec.get("label_asym_id", "?"))
            chains.setdefault(chain, []).append(
                [float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"])]
            )
        elif in_loop and line.startswith("#"):
            break
    return {c: np.asarray(v, dtype=float) for c, v in chains.items()}
