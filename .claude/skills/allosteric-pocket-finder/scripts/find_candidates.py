#!/usr/bin/env python3
"""Find distal (allosteric-candidate) conformational hotspots from an active/inactive
structure pair by robust Ca superposition, and propose His-pair/His-triplet metal sites.

Method
------
1. Match chain CA atoms between the two states by residue number.
2. Robust core superposition: iterative outlier rejection locks the reference frame onto
   the rigid domain so hinge/mobile regions stand out instead of being averaged away.
3. Per-residue Ca displacement + distance from the catalytic/active site (ligand centroid
   or named catalytic residues).
4. Residues that are BOTH high-displacement AND distal -> clustered into contiguous
   regions, ranked by (mean displacement x length). Target 5-20 regions.
5. Geometric scan flags residue pairs/triples whose Ca-Ca spacing is compatible with an
   engineered His-pair / His-triplet metal site.

Usage
-----
Single pair (CLI):
  find_candidates.py --bound 1V4S --unbound 1V4T --chain A --active-ligand GLC --name GCK
  find_candidates.py --bound 1T49 --unbound 1SUG --active-residues 215 --name PTP1B
Batch (JSON of {name: {bound,unbound,chain,active_site:{ligand|residues},note}}):
  find_candidates.py --pairs pairs.json
Common:
  --pdb-dir data/pdb  --out results  --top 12
"""
import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

DISPLACEMENT_MIN = 2.0    # A; residue considered "moved"
DISTAL_MIN = 15.0         # A; residue considered distal from active site
CORE_REJECT_FRAC = 0.30   # fraction of worst residues dropped each superposition iteration
CORE_ITERS = 5
HIS_PAIR_MIN, HIS_PAIR_MAX = 4.5, 12.0  # Ca-Ca window plausible for an engineered His site


def load_ca(pdb_path, chain_id):
    st = PDBParser(QUIET=True).get_structure(pdb_path.stem, pdb_path)
    model = next(iter(st))
    if chain_id not in model:
        chain_id = next(iter(model)).id
    out = {}
    for res in model[chain_id]:
        if res.id[0] != " " or "CA" not in res:
            continue
        ca = res["CA"]
        if ca.is_disordered() and "A" in ca.disordered_get_id_list():
            ca.disordered_select("A")
        out[res.id[1]] = (ca.get_coord().astype(float), res.get_resname())
    return out


def ligand_centroid(pdb_path, resname):
    st = PDBParser(QUIET=True).get_structure(pdb_path.stem, pdb_path)
    coords = [a.get_coord() for a in st.get_atoms()
              if a.get_parent().get_resname() == resname and a.element != "H"]
    if not coords:
        raise ValueError(f"ligand {resname} not found in {pdb_path.name}")
    return np.mean(coords, axis=0).astype(float)


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def robust_superpose(mob, ref):
    mask = np.ones(len(mob), dtype=bool)
    R, t = kabsch(mob, ref)
    for _ in range(CORE_ITERS):
        R, t = kabsch(mob[mask], ref[mask])
        dev = np.linalg.norm((mob @ R.T + t) - ref, axis=1)
        keep_n = max(3, int(round(len(mob) * (1 - CORE_REJECT_FRAC))))
        mask = dev <= np.sort(dev)[keep_n - 1]
    return R, t, mask


def cluster(resnums):
    runs, cur = [], []
    for r in sorted(resnums):
        if cur and r - cur[-1] > 2:
            runs.append(cur); cur = []
        cur.append(r)
    if cur:
        runs.append(cur)
    return runs


def analyze(cfg, pdb_dir):
    bp = pdb_dir / f"{cfg['bound'].upper()}.pdb"
    up = pdb_dir / f"{cfg['unbound'].upper()}.pdb"
    chain = cfg.get("chain", "A")
    b, u = load_ca(bp, chain), load_ca(up, chain)
    common = sorted(set(b) & set(u))
    if len(common) < 4:
        raise ValueError(f"only {len(common)} residues match between "
                         f"{cfg['unbound']} and {cfg['bound']} chain {chain}")
    P = np.array([u[r][0] for r in common])
    Q = np.array([b[r][0] for r in common])
    R, t, core = robust_superpose(P, Q)
    disp = np.linalg.norm((P @ R.T + t) - Q, axis=1)

    a = cfg["active_site"]
    if "ligand" in a:
        asite = ligand_centroid(bp, a["ligand"])
    else:
        asite = np.mean([b[r][0] for r in a["residues"] if r in b], axis=0)
    dist_as = np.linalg.norm(Q - asite, axis=1)

    rows = [{"resnum": r, "resname": b[r][1],
             "displacement": round(float(disp[i]), 2),
             "dist_active_site": round(float(dist_as[i]), 2),
             "in_core": bool(core[i])} for i, r in enumerate(common)]
    dmap = {row["resnum"]: row for row in rows}
    cand = [r["resnum"] for r in rows
            if r["displacement"] >= DISPLACEMENT_MIN and r["dist_active_site"] >= DISTAL_MIN]

    coord = {r: b[r][0] for r in common}
    regions = []
    for run in cluster(cand):
        d = [dmap[r]["displacement"] for r in run]
        da = [dmap[r]["dist_active_site"] for r in run]
        near = [r for r in coord if run[0] - 3 <= r <= run[-1] + 3]
        pairs = []
        for i in range(len(near)):
            for j in range(i + 1, len(near)):
                if near[j] - near[i] < 2:
                    continue
                dd = float(np.linalg.norm(coord[near[i]] - coord[near[j]]))
                if HIS_PAIR_MIN <= dd <= HIS_PAIR_MAX:
                    pairs.append({"i": near[i], "j": near[j], "ca_ca": round(dd, 2)})
        pairs.sort(key=lambda p: p["ca_ca"])
        regions.append({"start": run[0], "end": run[-1], "length": len(run),
                        "mean_displacement": round(float(np.mean(d)), 2),
                        "max_displacement": round(float(np.max(d)), 2),
                        "mean_dist_active_site": round(float(np.mean(da)), 2),
                        "residues": run, "his_pair_candidates": pairs[:6]})
    regions.sort(key=lambda x: x["mean_displacement"] * x["length"], reverse=True)

    stats = {"n_common": len(common), "n_core": int(core.sum()),
             "core_rmsd": round(float(np.sqrt(np.mean(disp[core] ** 2))), 2),
             "overall_rmsd": round(float(np.sqrt(np.mean(disp ** 2))), 2),
             "max_displacement": round(float(disp.max()), 2)}
    return rows, regions, stats


def write_csv(path, rows):
    lines = ["resnum,resname,displacement,dist_active_site,in_core"]
    lines += [f"{r['resnum']},{r['resname']},{r['displacement']},"
              f"{r['dist_active_site']},{int(r['in_core'])}" for r in rows]
    path.write_text("\n".join(lines) + "\n")


def build_pairs(args):
    if args.pairs:
        return json.loads(Path(args.pairs).read_text())
    if not (args.bound and args.unbound):
        raise SystemExit("give --pairs FILE, or --bound and --unbound (+ active site)")
    if args.active_ligand:
        asite = {"ligand": args.active_ligand}
    elif args.active_residues:
        asite = {"residues": [int(x) for x in args.active_residues]}
    else:
        raise SystemExit("define the active site with --active-ligand or --active-residues")
    name = args.name or f"{args.unbound}_{args.bound}"
    return {name: {"bound": args.bound, "unbound": args.unbound,
                   "chain": args.chain, "active_site": asite, "note": args.note or ""}}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", help="JSON of named pair configs (batch mode)")
    ap.add_argument("--bound"); ap.add_argument("--unbound")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--active-ligand", help="HETATM resname defining the active site")
    ap.add_argument("--active-residues", nargs="*", help="catalytic residue numbers")
    ap.add_argument("--name")
    ap.add_argument("--note", default="")
    ap.add_argument("--pdb-dir", type=Path, default=Path("data/pdb"))
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--top", type=int, default=12)
    args = ap.parse_args()

    pairs = build_pairs(args)
    args.out.mkdir(parents=True, exist_ok=True)
    summary = {}
    for name, cfg in pairs.items():
        rows, regions, stats = analyze(cfg, args.pdb_dir)
        top = regions[: args.top]
        write_csv(args.out / f"{name}_candidates.csv", rows)
        (args.out / f"{name}_regions.json").write_text(
            json.dumps({"system": name, "config": cfg, "stats": stats, "regions": top}, indent=2))
        summary[name] = {"stats": stats, "n_regions": len(regions)}
        print(f"\n=== {name}  ({cfg['unbound']} -> {cfg['bound']}) ===")
        if cfg.get("note"):
            print(f"  {cfg['note']}")
        print(f"  matched {stats['n_common']} | core {stats['n_core']} | "
              f"core-RMSD {stats['core_rmsd']} A | overall-RMSD {stats['overall_rmsd']} A | "
              f"max disp {stats['max_displacement']} A")
        print(f"  distal candidate regions (top {len(top)} of {len(regions)}):")
        for k, reg in enumerate(top, 1):
            hp = ", ".join(f"{p['i']}/{p['j']}({p['ca_ca']}A)" for p in reg["his_pair_candidates"][:3])
            print(f"   {k:2d}. res {reg['start']}-{reg['end']} (len {reg['length']}) "
                  f"meanDisp {reg['mean_displacement']} A, {reg['mean_dist_active_site']} A "
                  f"from active site | His-pairs: {hp or '-'}")
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote CSV/JSON + summary.json to {args.out}/")


if __name__ == "__main__":
    main()
