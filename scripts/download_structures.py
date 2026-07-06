#!/usr/bin/env python3
"""Download PDB coordinate files and canonical sequences for every target system.

Reads data/targets.json, fetches each PDB from RCSB into data/pdb/, and writes a
combined FASTA (one record per PDB chain) into data/fasta/. Idempotent: existing
files are skipped unless --force is given.
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGETS = ROOT / "data" / "targets.json"
PDB_DIR = ROOT / "data" / "pdb"
FASTA_DIR = ROOT / "data" / "fasta"

PDB_URL = "https://files.rcsb.org/download/{pdb}.pdb"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb}"


def fetch(url: str, retries: int = 4) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            wait = 2 ** (attempt + 1)
            print(f"  ! {url} failed ({e}); retry in {wait}s", file=sys.stderr)
            time.sleep(wait)
    raise RuntimeError(f"failed after {retries} attempts: {url} ({last})")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    args = ap.parse_args()

    PDB_DIR.mkdir(parents=True, exist_ok=True)
    FASTA_DIR.mkdir(parents=True, exist_ok=True)

    data = json.loads(TARGETS.read_text())
    pdbs = []
    for sysrec in data["systems"]:
        for st in sysrec["structures"]:
            pdbs.append((sysrec["id"], st["pdb"].upper()))

    # de-duplicate while preserving order
    seen = set()
    ordered = [(s, p) for s, p in pdbs if not (p in seen or seen.add(p))]

    print(f"Downloading {len(ordered)} unique PDB entries...")
    fasta_records = []
    for sysid, pdb in ordered:
        pdb_path = PDB_DIR / f"{pdb}.pdb"
        if pdb_path.exists() and not args.force:
            print(f"  = {pdb} (cached, {sysid})")
        else:
            print(f"  + {pdb} ({sysid})")
            pdb_path.write_bytes(fetch(PDB_URL.format(pdb=pdb)))
        fasta_records.append(fetch(FASTA_URL.format(pdb=pdb)).decode())

    combined = FASTA_DIR / "all_targets.fasta"
    combined.write_text("\n".join(r.strip() for r in fasta_records) + "\n")
    print(f"\nWrote {combined.relative_to(ROOT)} ({len(fasta_records)} entries)")
    print(f"PDB files in {PDB_DIR.relative_to(ROOT)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
