"""Small shared helpers: logging, file staging, timing."""
from __future__ import annotations

import json
import logging
import os
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator


def get_logger(name: str = "mdscreen") -> logging.Logger:
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(
            logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                              datefmt="%H:%M:%S")
        )
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
    return logger


log = get_logger()


@contextmanager
def timer(label: str) -> Iterator[None]:
    start = time.perf_counter()
    log.info("▶ %s ...", label)
    try:
        yield
    finally:
        log.info("✔ %s (%.1fs)", label, time.perf_counter() - start)


def ensure_dir(path: os.PathLike | str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(path: os.PathLike | str, obj: dict) -> None:
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_json_default)


def _json_default(o):
    # numpy scalars / arrays -> native python
    if hasattr(o, "tolist"):
        return o.tolist()
    if hasattr(o, "item"):
        return o.item()
    return str(o)


def read_bytes(path: os.PathLike | str) -> bytes:
    with open(path, "rb") as fh:
        return fh.read()


def is_smiles(s: str) -> bool:
    """Cheap heuristic: does this string look like a SMILES rather than a path?"""
    if os.path.exists(s):
        return False
    if any(s.lower().endswith(ext) for ext in (".sdf", ".mol", ".mol2", ".pdb")):
        return False
    # SMILES contain no path separators and no whitespace
    return "/" not in s and "\\" not in s and " " not in s.strip()
