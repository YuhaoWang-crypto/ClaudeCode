"""I/O helpers: reading molecule lists and writing per-step results."""

import json
import logging
from pathlib import Path
from typing import List

import pandas as pd

logger = logging.getLogger(__name__)


def read_smiles_list(path: str) -> List[str]:
    """Read SMILES from a file (one per line, or first column of CSV/TSV)."""
    p = Path(path)
    suffix = p.suffix.lower()

    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(p, sep=sep, header=None)
        smiles = df.iloc[:, 0].dropna().astype(str).tolist()
    else:
        smiles = [
            line.strip()
            for line in p.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]

    logger.info("Read %d SMILES from %s", len(smiles), path)
    return smiles


def save_results(results: List[dict], output_dir: str, step_name: str) -> Path:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    json_path = out / f"{step_name}.json"
    csv_path = out / f"{step_name}.csv"

    with open(json_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    if results:
        df = pd.json_normalize(results)
        df.to_csv(csv_path, index=False)

    logger.info("Saved %d records → %s", len(results), json_path)
    return json_path


def load_step_results(output_dir: str, step_name: str) -> List[dict]:
    p = Path(output_dir) / f"{step_name}.json"
    if not p.exists():
        return []
    with open(p) as f:
        return json.load(f)
