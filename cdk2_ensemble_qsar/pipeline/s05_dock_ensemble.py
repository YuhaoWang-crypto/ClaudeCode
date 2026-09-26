"""Stage 5 - dock every ligand into every conformational slice.

This is the expensive stage (n_ligands x n_slices docking runs), so it is
checkpointed: each finished job is appended to a JSONL ledger and re-runs skip
whatever is already in it. Killing and restarting the script is safe.

Poses are kept, not just scores. The whole point of the benchmark is to test
whether the *interaction pattern* carries signal that the scalar score does
not, and that pattern can only be read off a pose.

Docking is deterministic given --seed, so a re-run reproduces the ledger.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
STRUCT = ROOT / "structures"
RESULTS = ROOT / "results"
POSES = RESULTS / "poses"
SMINA = ROOT / "bin" / "smina"

EXHAUSTIVENESS = 8
SEED = 42
N_WORKERS = 4
LEDGER = RESULTS / "docking_ledger.jsonl"


def dock_one(job: dict) -> dict:
    """Dock one ligand into one receptor; return score and pose SDF text."""
    with tempfile.TemporaryDirectory() as tmp:
        out_sdf = Path(tmp) / "pose.sdf"
        cmd = [
            str(SMINA),
            "-r", job["receptor_pdbqt"],
            "-l", job["ligand_pdbqt"],
            "-o", str(out_sdf),
            "--center_x", str(job["cx"]),
            "--center_y", str(job["cy"]),
            "--center_z", str(job["cz"]),
            "--size_x", str(job["sx"]),
            "--size_y", str(job["sy"]),
            "--size_z", str(job["sz"]),
            "--exhaustiveness", str(EXHAUSTIVENESS),
            "--num_modes", "1",
            "--seed", str(SEED),
            "--cpu", "1",
        ]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        dt = time.time() - t0

        rec = {
            "ligand_id": job["ligand_id"],
            "receptor": job["receptor"],
            "seconds": round(dt, 1),
        }
        if proc.returncode != 0 or not out_sdf.exists():
            rec["status"] = "failed"
            rec["stderr"] = proc.stderr[-200:]
            return rec

        text = out_sdf.read_text()
        score = None
        for line in proc.stdout.splitlines():
            parts = line.split()
            if len(parts) >= 2 and parts[0] == "1":
                try:
                    score = float(parts[1])
                    break
                except ValueError:
                    pass
        rec["status"] = "ok"
        rec["score"] = score
        rec["sdf"] = text
        return rec


def load_done() -> set[tuple[str, str]]:
    done: set[tuple[str, str]] = set()
    if LEDGER.exists():
        with open(LEDGER) as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue  # truncated final line from a hard kill
                done.add((r["ligand_id"], r["receptor"]))
    return done


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    ap.add_argument("--limit", type=int, default=0, help="debug: cap job count")
    ap.add_argument(
        "--max-ligands",
        type=int,
        default=0,
        help="dock only the first N ligands (validation-scale run). Ligands are\ntaken in order, and every slice of a ligand is queued together, so the\nresult is a complete ensemble for a prefix of the screening set rather\nthan a ragged matrix.",
    )
    args = ap.parse_args()

    POSES.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((STRUCT / "ensemble_manifest.json").read_text())
    box = manifest["box"]
    ligands = pd.read_csv(DATA / "screening_set.csv")
    if args.max_ligands:
        ligands = ligands.head(args.max_ligands)
        print(f"validation-scale run: first {len(ligands)} ligands only")

    done = load_done()
    print(f"ledger already holds {len(done)} finished jobs")

    jobs = []
    for _, lig in ligands.iterrows():
        for sl in manifest["slices"]:
            key = (lig["ligand_id"], sl["pdb_id"])
            if key in done:
                continue
            jobs.append(
                {
                    "ligand_id": lig["ligand_id"],
                    "ligand_pdbqt": str(DATA / "ligands" / lig["pdbqt"]),
                    "receptor": sl["pdb_id"],
                    "receptor_pdbqt": str(STRUCT / sl["receptor_pdbqt"]),
                    "cx": box["center_x"], "cy": box["center_y"], "cz": box["center_z"],
                    "sx": box["size_x"], "sy": box["size_y"], "sz": box["size_z"],
                }
            )
    if args.limit:
        jobs = jobs[: args.limit]

    total = len(ligands) * len(manifest["slices"])
    print(
        f"{len(jobs)} jobs to run ({total} total, {len(done)} cached) "
        f"on {args.workers} workers"
    )
    if not jobs:
        print("nothing to do")
        return

    # Poses are appended per receptor as results arrive, in PLAIN text.
    #
    # An earlier version appended to gzip. That silently corrupts across a
    # restart: a killed process leaves an unterminated gzip member with no
    # trailer, the next run opens in append mode and starts a fresh member
    # after it, and the whole file then fails to decompress a few records in.
    # It cost a full re-dock to discover, because the ledger said 600/600 while
    # only 5 poses per receptor were readable. Plain text append has no such
    # failure mode; s06 reads .sdf or .sdf.gz, and compressing afterwards is a
    # separate, idempotent step.
    pose_handles = {
        sl["pdb_id"]: open(POSES / f"{sl['pdb_id']}.sdf", "a")
        for sl in manifest["slices"]
    }

    t0 = time.time()
    n_ok = n_fail = 0
    with open(LEDGER, "a") as ledger, ProcessPoolExecutor(
        max_workers=args.workers
    ) as ex:
        futures = {ex.submit(dock_one, j): j for j in jobs}
        for i, fut in enumerate(as_completed(futures), 1):
            rec = fut.result()
            sdf = rec.pop("sdf", None)
            if rec["status"] == "ok" and sdf:
                fh = pose_handles[rec["receptor"]]
                # Tag the record so poses stay traceable after concatenation.
                fh.write(sdf.replace("\n$$$$", f"\n> <LIGAND_ID>\n{rec['ligand_id']}\n\n$$$$", 1))
                fh.flush()
                n_ok += 1
            else:
                n_fail += 1
            ledger.write(json.dumps(rec) + "\n")
            ledger.flush()

            if i % 50 == 0 or i == len(jobs):
                elapsed = time.time() - t0
                rate = i / elapsed
                eta = (len(jobs) - i) / rate if rate else 0
                print(
                    f"  {i}/{len(jobs)}  ok={n_ok} fail={n_fail}  "
                    f"{rate * 60:.1f} jobs/min  eta {eta / 60:.0f} min",
                    flush=True,
                )

    for fh in pose_handles.values():
        fh.close()

    # --- flatten the ledger into a score matrix -----------------------------
    rows = [json.loads(line) for line in open(LEDGER)]
    df = pd.DataFrame([r for r in rows if r.get("status") == "ok"])
    wide = df.pivot_table(index="ligand_id", columns="receptor", values="score")
    wide.to_csv(RESULTS / "docking_scores.csv")
    print(f"\nwrote {RESULTS / 'docking_scores.csv'}  shape={wide.shape}")
    print(f"failures: {sum(1 for r in rows if r.get('status') != 'ok')}")
    print(f"mean {df['seconds'].mean():.1f}s per dock; wall {(time.time() - t0) / 60:.0f} min")


if __name__ == "__main__":
    main()
