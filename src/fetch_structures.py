"""
fetch_structures.py
-------------------
Resolve a target's conformer pair to local PDB files.

Resolution order for each structure:
  1. data/raw/<PDB>.pdb if already present (cache).
  2. An explicit offline file from targets.yaml `offline_files` (e.g. ProDy-bundled
     AdK), copied into data/raw/. This lets the AdK demo run with NO network.
  3. Network download, trying several mirrors in order. NOTE: in this sandbox all
     of these hosts are blocked by the egress allowlist (HTTP 403). The moment you
     add `files.rcsb.org` (and optionally `www.rcsb.org`) to the environment's
     network egress settings, this path starts working with no code change.

Mirrors tried (first that returns a real PDB wins):
  https://files.rcsb.org/download/<ID>.pdb
  https://www.ebi.ac.uk/pdbe/entry-files/download/pdb<idlower>.ent
  https://files.pdbj.org/pub/pdb/data/structures/divided/pdb/<mid>/pdb<idlower>.ent.gz
"""
from __future__ import annotations
import gzip
import os
import shutil
import urllib.request
import urllib.error

RAW_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "raw")


def _prody_datafile(name: str) -> str:
    import prody
    base = os.path.dirname(prody.__file__)
    path = os.path.join(base, "tests", "datafiles", name)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"ProDy bundled file not found: {path}")
    return path


def _mirrors(pdb_id: str):
    i = pdb_id.lower()
    return [
        f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
        f"https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{i}.ent",
        f"https://files.pdbj.org/pub/pdb/data/structures/divided/pdb/{i[1:3]}/pdb{i}.ent.gz",
    ]


def _download(url: str, dest: str) -> bool:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "allostery-pipeline/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            data = r.read()
        if url.endswith(".gz"):
            data = gzip.decompress(data)
        if b"ATOM" not in data[:200000]:
            return False
        with open(dest, "wb") as fh:
            fh.write(data)
        return True
    except (urllib.error.URLError, urllib.error.HTTPError, OSError, ValueError):
        return False


def resolve_structure(pdb_id: str, offline_file: str | None = None) -> str:
    """Return a local path to the PDB for `pdb_id`, fetching/copying as needed."""
    os.makedirs(RAW_DIR, exist_ok=True)
    dest = os.path.join(RAW_DIR, f"{pdb_id.upper()}.pdb")
    if os.path.isfile(dest) and os.path.getsize(dest) > 0:
        return dest

    # 2. offline bundled file
    if offline_file:
        if offline_file.startswith("prody:"):
            src = _prody_datafile(offline_file.split(":", 1)[1])
        else:
            src = offline_file
        if os.path.isfile(src):
            shutil.copyfile(src, dest)
            return dest

    # 3. network mirrors
    for url in _mirrors(pdb_id):
        if _download(url, dest):
            return dest

    raise RuntimeError(
        f"Could not obtain {pdb_id}. All structural-data hosts are blocked by the "
        f"network egress allowlist in this environment. Add 'files.rcsb.org' to the "
        f"environment network settings, or place {pdb_id.upper()}.pdb in {RAW_DIR}/."
    )
