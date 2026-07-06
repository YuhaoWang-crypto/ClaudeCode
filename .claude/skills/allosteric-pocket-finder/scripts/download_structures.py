#!/usr/bin/env python3
"""Download PDB coordinate files (+ FASTA sequences) from RCSB.

Two modes:
  --pdb 1V4S 1V4T ...            explicit list of PDB IDs
  --manifest targets.json        systems manifest (see assets/targets.json schema)
Defaults: --out data/ (creates data/pdb/ and data/fasta/). Idempotent + retrying.
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

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


def pdbs_from_manifest(path: Path):
    data = json.loads(path.read_text())
    out = []
    for sysrec in data.get("systems", []):
        for st in sysrec.get("structures", []):
            out.append((sysrec.get("id", "?"), st["pdb"].upper()))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", nargs="*", default=[], help="explicit PDB IDs")
    ap.add_argument("--manifest", type=Path, help="systems manifest JSON")
    ap.add_argument("--out", type=Path, default=Path("data"), help="output root (default: data/)")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    args = ap.parse_args()

    entries = [("cli", p.upper()) for p in args.pdb]
    if args.manifest:
        entries += pdbs_from_manifest(args.manifest)
    if not entries:
        ap.error("give --pdb <IDs> and/or --manifest <file>")

    pdb_dir = args.out / "pdb"
    fasta_dir = args.out / "fasta"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    seen = set()
    ordered = [(s, p) for s, p in entries if not (p in seen or seen.add(p))]
    print(f"Downloading {len(ordered)} unique PDB entries -> {args.out}/")
    fasta_records = []
    for tag, pdb in ordered:
        dest = pdb_dir / f"{pdb}.pdb"
        if dest.exists() and not args.force:
            print(f"  = {pdb} (cached, {tag})")
        else:
            print(f"  + {pdb} ({tag})")
            dest.write_bytes(fetch(PDB_URL.format(pdb=pdb)))
        fasta_records.append(fetch(FASTA_URL.format(pdb=pdb)).decode())

    combined = fasta_dir / "all_targets.fasta"
    combined.write_text("\n".join(r.strip() for r in fasta_records) + "\n")
    print(f"\nWrote {combined} ({len(fasta_records)} sequence entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
