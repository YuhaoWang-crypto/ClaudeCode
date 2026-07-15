#!/usr/bin/env python3
"""
Step 4 -- DeePMD vs MACE vs NequIP bake-off. Reads each trainer's native output,
normalises to validation energy-RMSE (meV/atom) and force-RMSE (meV/A), prints a
ranking table, and writes a comparison bar chart.

    # after training each model on the SAME data (same r_max/rcut):
    python compare_potentials.py \
        --deepmd ../data/deepmd_run/lcurve.out \
        --mace   mace_demo_run/train.log \
        --nequip nequip_runs/desal_nequip/metrics_epoch.csv \
        --out ../06_dft_surrogate_study/results/bakeoff.png

    python compare_potentials.py --demo   # parser self-check on sample data

Notes on normalisation
----------------------
Each code reports errors differently; we convert everything to
(meV/atom for energy, meV/A for forces) on the VALIDATION set:
  * MACE   -- PerAtomRMSE error table is already meV/atom and meV/A.
  * DeePMD -- lcurve.out rmse_e_val is total-energy eV; pass --natoms to get
              per-atom, and eV->meV. rmse_f_val is eV/A -> meV/A.
  * NequIP -- metrics_epoch.csv validation columns in eV and eV/A -> meV.
Always train the three with identical r_max/rcut and data for a fair comparison.
"""
import argparse
import csv
import os
import re


def parse_mace_log(path, natoms=None):
    """Parse the 'Error-table on TRAIN and VALID' block MACE prints at the end.
    Returns validation (e_rmse_meV_atom, f_rmse_meV_A)."""
    e = f = None
    with open(path) as fh:
        lines = fh.readlines()
    for i, ln in enumerate(lines):
        if "valid" in ln.lower() and "|" in ln:
            cols = [c.strip() for c in ln.split("|")]
            # columns: '', config_type, RMSE E/meV/atom, RMSE F/meV/A, relF%
            nums = [c for c in cols if re.match(r"^[\d.]+$", c)]
            if len(nums) >= 2:
                e, f = float(nums[0]), float(nums[1])
    return e, f


def parse_deepmd_lcurve(path, natoms=None):
    """Parse DeePMD lcurve.out; return validation (e_rmse_meV_atom, f_meV_A).
    Energy is total eV -> needs --natoms for per-atom."""
    header, last = None, None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                header = ln.lstrip("#").split()
            elif ln.strip():
                last = ln.split()
    if not header or not last:
        return None, None
    idx = {name: k for k, name in enumerate(header)}

    def get(key):
        return float(last[idx[key]]) if key in idx else None

    e_ev = get("rmse_e_val") or get("rmse_e")
    f_ev = get("rmse_f_val") or get("rmse_f")
    e = (e_ev * 1000.0 / natoms) if (e_ev is not None and natoms) else (
        e_ev * 1000.0 if e_ev is not None else None)
    f = f_ev * 1000.0 if f_ev is not None else None
    return e, f


def parse_nequip_metrics(path, natoms=None):
    """Parse NequIP metrics_epoch.csv; return validation (e_meV_atom, f_meV_A)."""
    rows = list(csv.DictReader(open(path)))
    if not rows:
        return None, None
    last = rows[-1]

    def find(*subs):
        for k in last:
            kl = k.lower()
            if all(s in kl for s in subs):
                try:
                    return float(last[k])
                except ValueError:
                    pass
        return None

    e = find("valid", "e", "rmse")
    f = find("valid", "f", "rmse")
    # NequIP energy metric may be per-atom already ('e/N'); heuristics:
    if e is not None:
        e = e * 1000.0 / (natoms if (natoms and "n" not in "e/n") else 1) \
            if e < 1 else e * 1000.0
    if f is not None:
        f = f * 1000.0
    return e, f


def make_table(results):
    print(f"\n{'model':10s} {'E-RMSE (meV/atom)':>18s} {'F-RMSE (meV/A)':>16s}")
    print("-" * 46)
    ranked = sorted([r for r in results if r[2] is not None],
                    key=lambda r: r[2])  # rank by force RMSE
    for name, e, f in results:
        es = f"{e:.1f}" if e is not None else "n/a"
        fs = f"{f:.1f}" if f is not None else "n/a"
        print(f"{name:10s} {es:>18s} {fs:>16s}")
    if ranked:
        print(f"\nBest force accuracy: {ranked[0][0]} "
              f"({ranked[0][2]:.1f} meV/A)")
    return ranked


def plot(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    names = [r[0] for r in results]
    e = [r[1] if r[1] is not None else 0 for r in results]
    f = [r[2] if r[2] is not None else 0 for r in results]
    x = np.arange(len(names))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    a1.bar(x, e, color="#48a"); a1.set_xticks(x); a1.set_xticklabels(names)
    a1.set_ylabel("validation E-RMSE (meV/atom)"); a1.set_title("Energy")
    a2.bar(x, f, color="#c74"); a2.set_xticks(x); a2.set_xticklabels(names)
    a2.set_ylabel("validation F-RMSE (meV/Å)"); a2.set_title("Forces")
    fig.suptitle("MLP bake-off (same data, same cutoff) — lower is better")
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")


def _demo():
    """Parser self-check on the real MACE demo log + synthetic DeePMD/NequIP."""
    here = os.path.dirname(__file__)
    scratch = os.environ.get("TMPDIR", "/tmp")
    # synthetic deepmd lcurve
    dp = os.path.join(scratch, "lcurve_demo.out")
    with open(dp, "w") as fh:
        fh.write("# step rmse_e_val rmse_f_val rmse_e_trn rmse_f_trn lr\n")
        fh.write("1000 0.020 0.150 0.018 0.140 1e-3\n")
        fh.write("50000 0.0021 0.045 0.0019 0.040 3e-5\n")
    # synthetic nequip metrics
    nq = os.path.join(scratch, "metrics_demo.csv")
    with open(nq, "w") as fh:
        fh.write("epoch,validation_e_rmse,validation_f_rmse\n")
        fh.write("1,0.5,2.0\n100,0.030,0.055\n")

    mace_log = os.path.join(here, "mace_demo_run", "train.log")
    results = []
    if os.path.exists(mace_log):
        e, f = parse_mace_log(mace_log)
        results.append(("MACE", e, f))
    e, f = parse_deepmd_lcurve(dp, natoms=38)
    results.append(("DeePMD", e, f))
    e, f = parse_nequip_metrics(nq, natoms=38)
    results.append(("NequIP", e, f))
    make_table(results)
    print("\n(demo: DeePMD/NequIP are synthetic sample files; MACE is the real "
          "demo-run log if present)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--deepmd"); ap.add_argument("--mace")
    ap.add_argument("--nequip")
    ap.add_argument("--natoms", type=int, default=None,
                    help="atoms/frame (for per-atom energy normalisation)")
    ap.add_argument("--out", default="bakeoff.png")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    if args.demo:
        _demo(); return

    results = []
    if args.deepmd:
        results.append(("DeePMD", *parse_deepmd_lcurve(args.deepmd, args.natoms)))
    if args.mace:
        results.append(("MACE", *parse_mace_log(args.mace, args.natoms)))
    if args.nequip:
        results.append(("NequIP", *parse_nequip_metrics(args.nequip, args.natoms)))
    if not results:
        raise SystemExit("Provide at least one of --deepmd/--mace/--nequip "
                         "(or --demo)")
    make_table(results)
    plot(results, args.out)


if __name__ == "__main__":
    main()
