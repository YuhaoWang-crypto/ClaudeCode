"""Interface metrics + epitope recall for every design in a directory tree of Boltz
design outputs (each subdir holding predicted_structure.cif).

Usage:
    python batch_interface.py --root pae/ --out interface_all.csv [--workers 8]
"""

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from analyze_interface import analyze, load  # noqa: E402

# 1P9M chain-B numbering; the design models number the target entity from 1 with a
# three-residue N-terminal extension, so 1P9M number = model number - 2.
OFFSET = -2
SITE_I = [30, 33, 54, 61, 66, 69, 73, 74, 75, 78, 172, 175, 178, 179, 180, 182, 183]
SITE_II = [19, 24, 27, 28, 30, 31, 34, 110, 111, 113, 114, 117, 118, 121, 124, 125, 128]


def one(d):
    d = Path(d)
    cif = d / "predicted_structure.cif"
    try:
        r = analyze(load(cif), "B", "X")
    except Exception as e:  # keep the batch going, record the failure
        return {"design_id": d.name, "error": str(e)[:120]}
    iface = {k + OFFSET for k in r["target_interface_residues"]}
    n = max(1, len(iface))
    return {
        "design_id": d.name,
        "binder_len": r["binder_length"],
        "sequence": r["binder_sequence"],
        "bsa_A2": r["bsa_total_A2"],
        "n_target_iface": r["n_target_interface"],
        "atom_pairs_4A": r["n_atom_pairs_within_4A"],
        "hbonds": len(r["hbond_like_contacts"]),
        "salt_bridges": len(r["salt_bridges"]),
        "mean_kd_hydropathy": r["binder_mean_kd_hydropathy"],
        "siteI_recall": round(len(iface & set(SITE_I)) / len(SITE_I), 3),
        "siteI_footprint_frac": round(len(iface & set(SITE_I)) / n, 3),
        "siteII_recall": round(len(iface & set(SITE_II)) / len(SITE_II), 3),
        "siteII_footprint_frac": round(len(iface & set(SITE_II)) / n, 3),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    jobs, runs = [], {}
    for run_dir in sorted(p for p in args.root.iterdir() if p.is_dir()):
        for d in sorted(p for p in run_dir.iterdir() if p.is_dir()):
            jobs.append(str(d))
            runs[d.name] = run_dir.name
    print(f"analysing {len(jobs)} interfaces")

    with ProcessPoolExecutor(args.workers) as ex:
        recs = list(ex.map(one, jobs, chunksize=4))

    ok = [r for r in recs if "error" not in r]
    for r in ok:
        r["run"] = runs[r["design_id"]]
    print(f"done {len(ok)}, failed {len(recs) - len(ok)}")

    fields = ["run", "design_id", "binder_len", "bsa_A2", "n_target_iface", "atom_pairs_4A",
              "hbonds", "salt_bridges", "mean_kd_hydropathy", "siteI_recall",
              "siteI_footprint_frac", "siteII_recall", "siteII_footprint_frac", "sequence"]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(ok, key=lambda r: (r["run"], r["design_id"])))
    print("wrote", args.out)


if __name__ == "__main__":
    main()
