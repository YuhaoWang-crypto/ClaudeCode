"""Download the full ArchMap atlas files (expression matrices + trained models).

All 17 atlases together are ~103 GB (PBMC counts alone ~50 GB), so pick what
you need:

    python download_atlas.py --list
    python download_atlas.py HLCA Heart --kinds reference_data model
    python download_atlas.py --all --kinds reference_data --out /data/archmap

Kinds: reference_data (data.h5ad, the model-feature matrix ArchMap maps onto),
reference_counts (data_only_count.h5ad, raw counts), model (model.pt /
model_params.pt / attr.pkl / var_names.csv).  Downloads resume with curl -C -.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from archmap_api import ArchMap, remote_size
from build_db import file_kind


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("atlases", nargs="*", help="atlas names or ids (case-insensitive)")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--kinds", nargs="+", default=["reference_data", "model"],
                    choices=["reference_data", "reference_counts", "model"])
    ap.add_argument("--out", default="archmap_data")
    args = ap.parse_args()

    api = ArchMap()
    atlases = api.atlases()
    if args.list:
        for a in atlases:
            print(f"{a['_id']}  {a['name']:<36} {a.get('numberOfCells') or 0:>10,} cells  {a['species']}")
        return
    want = {w.lower() for w in args.atlases}
    chosen = [a for a in atlases if args.all or a["name"].lower() in want or a["_id"] in want]
    missing = want - {a["name"].lower() for a in chosen} - {a["_id"] for a in chosen}
    if missing:
        raise SystemExit(f"unknown atlas: {', '.join(sorted(missing))} (see --list)")

    for a in chosen:
        for f in api.atlas_files(a["_id"]):
            kind = file_kind(f["fileName"])
            if not (kind in args.kinds or (kind.startswith("model") and "model" in args.kinds)):
                continue
            dest = Path(args.out) / a["name"].replace("/", "_").replace(" ", "_") / f["fileName"]
            dest.parent.mkdir(parents=True, exist_ok=True)
            size = remote_size(f["presignedUrl"])
            if dest.exists() and dest.stat().st_size == size:
                print(f"ok   {dest}")
                continue
            print(f"get  {dest}  ({(size or 0) / 1e9:.2f} GB)")
            subprocess.run(["curl", "-fL", "--retry", "5", "-C", "-", "-o", str(dest),
                            f["presignedUrl"]], check=True)


if __name__ == "__main__":
    main()
