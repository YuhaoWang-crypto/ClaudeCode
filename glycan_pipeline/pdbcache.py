"""Fetch and cache structures from RCSB.

Cached under ``figures/_pdb_cache/`` so repeat runs are offline, mirroring the
``figures/_sbml_cache/`` convention already used by ``grn_pipeline``.
"""

from __future__ import annotations

import os
import urllib.request

CACHE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "figures", "_pdb_cache",
)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
RCSB_ASSEMBLY_URL = "https://files.rcsb.org/download/{pdb_id}-assembly{n}.cif"


def fetch(pdb_id: str, assembly: int | None = None,
          cache_dir: str | None = None, timeout: int = 90) -> str:
    """Return a local path to a structure, downloading it once if needed.

    ``assembly=1`` fetches the biological assembly rather than the asymmetric
    unit. That matters for glycan work: in an oligomeric glycoprotein the
    asymmetric unit often holds a single protomer, so every glycan contributed
    by a neighbouring subunit is simply missing from the coordinates and any
    glycan-glycan shielding term is silently undercounted.
    """
    pdb_id = pdb_id.upper().strip()
    cache_dir = cache_dir or CACHE_DIR
    os.makedirs(cache_dir, exist_ok=True)
    suffix = "" if assembly is None else f"-assembly{assembly}"
    path = os.path.join(cache_dir, f"{pdb_id}{suffix}.cif")
    if os.path.exists(path) and os.path.getsize(path) > 10_000:
        return path

    if assembly is None:
        url = RCSB_URL.format(pdb_id=pdb_id)
    else:
        url = RCSB_ASSEMBLY_URL.format(pdb_id=pdb_id, n=int(assembly))
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        if resp.status != 200:
            raise RuntimeError(f"RCSB returned HTTP {resp.status} for {pdb_id}")
        data = resp.read()
    if len(data) < 10_000:
        raise RuntimeError(f"{pdb_id}: suspiciously small response ({len(data)} B)")
    with open(path, "wb") as fh:
        fh.write(data)
    return path
