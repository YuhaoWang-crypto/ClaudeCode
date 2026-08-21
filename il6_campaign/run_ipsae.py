"""Run the reference ipSAE implementation (Dunbrack lab, ipsae.py v4) over a directory of
Boltz design outputs and collect ipSAE_min per design.

Each design directory must contain pae.npz + predicted_structure.cif (+ metrics.json),
which is exactly what the Boltz protein-design archive unpacks to.

ipSAE is asymmetric; the paper's ranking metric is ipSAE_min, the smaller of the two
chain-direction values. pDockQ, pDockQ2 and LIS come along for free and are recorded.

Usage:
    python run_ipsae.py --ipsae ipsae.py --root pae/ --out ipsae_scores.csv \
        [--pae-cutoff 10] [--dist-cutoff 15] [--workers 8]
"""

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


def parse_ipsae_table(path: Path):
    """Parse the *_<pae>_<dist>.txt table written by ipsae.py."""
    rows = []
    with path.open() as fh:
        header = []
        for line in fh:  # the table is preceded by a blank line
            header = line.split()
            if header:
                break
        for line in fh:
            parts = line.split()
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def score_one(args):
    ipsae_script, design_dir, pae_cut, dist_cut = args
    design_dir = Path(design_dir)
    pae = design_dir / "pae.npz"
    cif = design_dir / "predicted_structure.cif"
    if not (pae.exists() and cif.exists()):
        return {"design_id": design_dir.name, "error": "missing inputs"}

    tag = f"{int(pae_cut):02d}_{int(dist_cut):02d}"
    out = design_dir / f"predicted_structure_{tag}.txt"
    if not out.exists():
        proc = subprocess.run(
            [sys.executable, str(ipsae_script), str(pae), str(cif), str(pae_cut), str(dist_cut)],
            capture_output=True, text=True, timeout=600,
        )
        if not out.exists():
            return {"design_id": design_dir.name, "error": proc.stderr.strip()[:200] or "no output"}

    rows = parse_ipsae_table(out)
    asym = [r for r in rows if r.get("Type") == "asym"]
    if len(asym) < 2:
        return {"design_id": design_dir.name, "error": "no asym rows"}

    vals = [float(r["ipSAE"]) for r in asym]
    rec = {
        "design_id": design_dir.name,
        "ipsae_min": min(vals),
        "ipsae_max": max(vals),
        "ipsae_binder_to_target": next(
            (float(r["ipSAE"]) for r in asym if r["Chn1"] == "X"), float("nan")),
        "ipsae_target_to_binder": next(
            (float(r["ipSAE"]) for r in asym if r["Chn1"] != "X"), float("nan")),
        "pdockq": float(asym[0]["pDockQ"]),
        "pdockq2": float(asym[0]["pDockQ2"]),
        "lis": max(float(r["LIS"]) for r in asym),
        "n0res_min": min(int(r["n0res"]) for r in asym),
    }
    metrics = design_dir / "metrics.json"
    if metrics.exists():
        m = json.loads(metrics.read_text())
        rec.update(iptm=m.get("iptm"), min_interaction_pae=m.get("min_interaction_pae"),
                   structure_confidence=m.get("structure_confidence"),
                   binding_confidence=m.get("binding_confidence"),
                   complex_plddt=m.get("complex_plddt"),
                   helix_fraction=m.get("helix_fraction"),
                   sheet_fraction=m.get("sheet_fraction"))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ipsae", required=True, type=Path, help="path to ipsae.py")
    ap.add_argument("--root", required=True, type=Path, help="dir of <run>/<design_id>/ dirs")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--pae-cutoff", type=float, default=10)
    ap.add_argument("--dist-cutoff", type=float, default=15)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    jobs, runs = [], {}
    for run_dir in sorted(p for p in args.root.iterdir() if p.is_dir()):
        for d in sorted(p for p in run_dir.iterdir() if p.is_dir()):
            jobs.append((args.ipsae, d, args.pae_cutoff, args.dist_cutoff))
            runs[d.name] = run_dir.name
    print(f"scoring {len(jobs)} designs with ipSAE (pae<{args.pae_cutoff}, dist<{args.dist_cutoff})")

    with ProcessPoolExecutor(args.workers) as ex:
        recs = list(ex.map(score_one, jobs))

    ok = [r for r in recs if "error" not in r]
    bad = [r for r in recs if "error" in r]
    for r in ok:
        r["run"] = runs[r["design_id"]]
    print(f"scored {len(ok)}, failed {len(bad)}")
    for r in bad[:5]:
        print("  FAIL", r["design_id"], r["error"])

    fields = ["run", "design_id", "ipsae_min", "ipsae_max", "ipsae_binder_to_target",
              "ipsae_target_to_binder", "pdockq", "pdockq2", "lis", "n0res_min",
              "iptm", "min_interaction_pae", "structure_confidence", "binding_confidence",
              "complex_plddt", "helix_fraction", "sheet_fraction"]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(ok, key=lambda r: (r["run"], -r["ipsae_min"])))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
