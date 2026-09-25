"""Stage 3a - cross-docking control: does the ensemble premise hold?

Every ATP-site co-crystal ligand is docked into every slice receptor, all in
the one superposed frame. Two things fall out:

* the diagonal is the classic redocking control - can the setup put a ligand
  back where crystallography found it (RMSD < 2 A)?
* the off-diagonal is the actual justification for ensemble docking. If any
  receptor reproduced any ligand's pose equally well, conformational slices
  would be redundant and a single structure would do. A spread of
  cross-docking RMSDs is the evidence that the slices carry distinct
  information.

RMSD is symmetry-aware and template-free. The reference is the *docking input*
PDBQT converted back to SDF rather than the deposited PDB: it carries the same
crystallographic coordinates but has been through the same Open Babel bond
perception and hydrogen addition as the pose, so the two molecules are
guaranteed to share a topology and rdMolAlign.CalcRMS can find the best
symmetry mapping without re-aligning.

An earlier version used the RCSB ideal-geometry SDF as a bond-order template
and returned n/a for every single ligand. Open Babel adds polar hydrogens when
it writes the PDBQT, and RDKit silently keeps explicit hydrogens when a file is
read with sanitize=False, so a 23-heavy-atom template was being matched against
a 26-atom probe. Index-order RMSD is wrong here too: on LS1 it reads 1.23 A
where the symmetry-corrected value is 0.74 A, because the sulfonamide-phenyl
flip is a genuine equivalence rather than a displacement.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "structures"
RESULTS = ROOT / "results"
WORK = RESULTS / "crossdock"
SMINA = ROOT / "bin" / "smina"

EXHAUSTIVENESS = 16
SEED = 42
RMSD_PASS = 2.0
N_WORKERS = 4


def load_heavy(sdf: Path) -> Chem.Mol | None:
    """Read an SDF, sanitize, and drop hydrogens.

    Sanitizing before RemoveHs matters: RDKit cannot strip hydrogens from an
    unsanitized molecule, and silently leaving them in breaks any downstream
    atom matching.
    """
    mol = Chem.MolFromMolFile(str(sdf), removeHs=False, sanitize=False)
    if mol is None:
        return None
    if Chem.SanitizeMol(mol, catchErrors=True) != Chem.SanitizeFlags.SANITIZE_NONE:
        return None
    try:
        return Chem.RemoveHs(mol)
    except Exception:
        return None


def dock_one(job: dict) -> dict:
    """Dock one ligand into one receptor. Runs in a worker process."""
    out_pdbqt = WORK / f"{job['ligand_slice']}_{job['ccd']}_into_{job['receptor']}.pdbqt"
    cmd = [
        str(SMINA),
        "-r", str(STRUCT / job["receptor_pdbqt"]),
        "-l", str(job["ligand_pdbqt"]),
        "-o", str(out_pdbqt),
        "--center_x", str(job["box"]["center_x"]),
        "--center_y", str(job["box"]["center_y"]),
        "--center_z", str(job["box"]["center_z"]),
        "--size_x", str(job["box"]["size_x"]),
        "--size_y", str(job["box"]["size_y"]),
        "--size_z", str(job["box"]["size_z"]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", "5",
        "--seed", str(SEED),
        "--cpu", "1",
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    if proc.returncode != 0:
        return {**_key(job), "status": "smina_failed", "stderr": proc.stderr[-300:]}

    score = None
    for line in proc.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                score = float(parts[1])
                break
            except ValueError:
                pass

    top_sdf = out_pdbqt.with_name(out_pdbqt.stem + "_top.sdf")
    subprocess.run(
        ["obabel", str(out_pdbqt), "-O", str(top_sdf), "-f", "1", "-l", "1"],
        check=True,
        capture_output=True,
    )
    return {
        **_key(job),
        "status": "ok",
        "score_kcal_mol": score,
        "seconds": round(dt, 1),
        "top_sdf": str(top_sdf),
    }


def _key(job: dict) -> dict:
    return {
        "ligand_slice": job["ligand_slice"],
        "ccd": job["ccd"],
        "receptor": job["receptor"],
        "is_self_dock": job["ligand_slice"] == job["receptor"],
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--reuse",
        action="store_true",
        help="recompute RMSD from cached poses instead of re-docking",
    )
    args = ap.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((STRUCT / "ensemble_manifest.json").read_text())
    box = manifest["box"]
    slices = manifest["slices"]

    # Ligands: the ATP-site co-crystal ligands, one per chemotype.
    ligands = []
    seen_ccd: set[str] = set()
    for entry in slices:
        if not entry["ligand_in_atp_site"]:
            continue
        ccd = entry["native_ligand"]
        if ccd in seen_ccd:
            print(f"  {entry['pdb_id']}: {ccd} already represented, skipping duplicate")
            continue
        seen_ccd.add(ccd)
        lig_pdb = STRUCT / entry["native_ligand_pdb"]
        lig_pdbqt = WORK / f"{entry['pdb_id']}_{ccd}.pdbqt"
        subprocess.run(
            ["obabel", str(lig_pdb), "-O", str(lig_pdbqt), "-h"],
            check=True,
            capture_output=True,
        )
        # Round-trip the docking input back to SDF: this, not the deposited PDB,
        # is the RMSD reference (see the module docstring).
        ref_sdf = WORK / f"{entry['pdb_id']}_{ccd}_reference.sdf"
        subprocess.run(
            ["obabel", str(lig_pdbqt), "-O", str(ref_sdf)],
            check=True,
            capture_output=True,
        )
        ligands.append(
            {
                "slice": entry["pdb_id"],
                "ccd": ccd,
                "pdb": lig_pdb,
                "pdbqt": lig_pdbqt,
                "ref_sdf": ref_sdf,
            }
        )
    print(f"ligands: {[(l['slice'], l['ccd']) for l in ligands]}")
    print(f"receptors: {[s['pdb_id'] for s in slices]}")

    jobs = [
        {
            "ligand_slice": lig["slice"],
            "ccd": lig["ccd"],
            "ligand_pdbqt": lig["pdbqt"],
            "receptor": rec["pdb_id"],
            "receptor_pdbqt": rec["receptor_pdbqt"],
            "box": box,
        }
        for lig in ligands
        for rec in slices
    ]

    cache = RESULTS / "crossdock_control.json"
    if args.reuse and cache.exists():
        # Docking is the expensive part and is deterministic given --seed, so a
        # re-run that only changes the RMSD calculation reuses the saved poses.
        prev = {(r["ccd"], r["receptor"]): r for r in json.loads(cache.read_text())["results"]}
        results = []
        for job in jobs:
            r = dict(prev.get((job["ccd"], job["receptor"]), {}))
            r.update(_key(job))
            r.setdefault("status", "missing")
            r["top_sdf"] = str(
                WORK
                / f"{job['ligand_slice']}_{job['ccd']}_into_{job['receptor']}_top.sdf"
            )
            if not Path(r["top_sdf"]).exists():
                r["status"] = "missing"
            results.append(r)
        n_reused = sum(1 for r in results if r["status"] == "ok")
        print(f"\nreusing {n_reused}/{len(jobs)} cached poses (no docking run)")
    else:
        print(f"\nrunning {len(jobs)} cross-docking jobs on {N_WORKERS} workers ...")
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
            results = list(ex.map(dock_one, jobs))
        print(f"docking wall time: {time.time() - t0:.0f}s")

    # --- symmetry-aware RMSD against each ligand's crystallographic pose ----
    by_ccd = {lig["ccd"]: lig for lig in ligands}
    for res in results:
        if res["status"] != "ok":
            res["rmsd_A"] = None
            continue
        lig = by_ccd[res["ccd"]]
        ref = load_heavy(lig["ref_sdf"])
        pose = load_heavy(Path(res["top_sdf"]))
        if ref is None or pose is None:
            res["rmsd_A"] = None
            res["rmsd_note"] = "molecule unreadable"
            continue
        if Chem.MolToSmiles(ref) != Chem.MolToSmiles(pose):
            res["rmsd_A"] = None
            res["rmsd_note"] = "reference and pose topologies disagree"
            continue
        try:
            res["rmsd_A"] = round(float(rdMolAlign.CalcRMS(pose, ref)), 2)
        except Exception as exc:
            res["rmsd_A"] = None
            res["rmsd_note"] = f"CalcRMS failed: {exc}"

    # --- report -------------------------------------------------------------
    rec_ids = [s["pdb_id"] for s in slices]
    lig_ids = [l["ccd"] for l in ligands]
    print("\ncross-docking RMSD (A)   rows = ligand, cols = receptor")
    print("        " + "".join(f"{r:>9s}" for r in rec_ids))
    for lig in ligands:
        cells = []
        for rid in rec_ids:
            m = next(
                (
                    r
                    for r in results
                    if r["ccd"] == lig["ccd"] and r["receptor"] == rid
                ),
                None,
            )
            v = m.get("rmsd_A") if m else None
            mark = "*" if m and m["is_self_dock"] else " "
            cells.append(f"{v:>8.2f}{mark}" if v is not None else f"{'n/a':>8s} ")
        print(f"{lig['ccd']:>6s}  " + "".join(cells))
    print("  (* = self-dock, i.e. the classic redocking control)")

    print("\ncross-docking score (kcal/mol)")
    print("        " + "".join(f"{r:>9s}" for r in rec_ids))
    for lig in ligands:
        cells = []
        for rid in rec_ids:
            m = next(
                (r for r in results if r["ccd"] == lig["ccd"] and r["receptor"] == rid),
                None,
            )
            v = m.get("score_kcal_mol") if m else None
            cells.append(f"{v:>9.2f}" if v is not None else f"{'n/a':>9s}")
        print(f"{lig['ccd']:>6s}  " + "".join(cells))

    self_rmsd = [
        r["rmsd_A"] for r in results if r["is_self_dock"] and r.get("rmsd_A") is not None
    ]
    cross_rmsd = [
        r["rmsd_A"]
        for r in results
        if not r["is_self_dock"] and r.get("rmsd_A") is not None
    ]
    timings = [r["seconds"] for r in results if r.get("seconds")]

    summary = {
        "exhaustiveness": EXHAUSTIVENESS,
        "rmsd_pass_threshold_A": RMSD_PASS,
        "n_self_dock": len(self_rmsd),
        "n_self_dock_pass": sum(1 for v in self_rmsd if v < RMSD_PASS),
        "self_dock_rmsd_median": round(float(np.median(self_rmsd)), 2) if self_rmsd else None,
        "cross_dock_rmsd_median": round(float(np.median(cross_rmsd)), 2) if cross_rmsd else None,
        "cross_dock_rmsd_frac_under_2A": round(
            float(np.mean([v < RMSD_PASS for v in cross_rmsd])), 3
        )
        if cross_rmsd
        else None,
        "mean_seconds_per_dock": round(float(np.mean(timings)), 1) if timings else None,
        "results": [
            {k: v for k, v in r.items() if k != "top_sdf"} for r in results
        ],
    }
    (RESULTS / "crossdock_control.json").write_text(json.dumps(summary, indent=2))
    print(
        f"\nself-dock pass {summary['n_self_dock_pass']}/{summary['n_self_dock']} "
        f"(median {summary['self_dock_rmsd_median']} A); "
        f"cross-dock median {summary['cross_dock_rmsd_median']} A, "
        f"{summary['cross_dock_rmsd_frac_under_2A']} under {RMSD_PASS} A"
    )
    print(f"mean {summary['mean_seconds_per_dock']}s per dock at exh={EXHAUSTIVENESS}")


if __name__ == "__main__":
    main()
