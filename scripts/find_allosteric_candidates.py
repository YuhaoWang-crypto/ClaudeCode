#!/usr/bin/env python3
"""Find distal (allosteric-candidate) conformational hotspots from an active/inactive
structure pair by robust Ca superposition.

Method
------
1. Parse both structures, keep chain A CA atoms, match residues by residue number.
2. Superpose iteratively on a stable "core": fit -> reject residues whose deviation is
   in the top fraction -> refit, repeat. This locks the reference frame onto the rigid
   domain so that hinge/mobile regions stand out instead of being averaged away.
3. Per-residue Ca displacement between states is computed in the core frame.
4. Distance of each residue from the catalytic/active site is measured (centroid of the
   defining ligand in the bound structure, or a named catalytic residue).
5. Residues that are BOTH high-displacement AND distal from the active site are the
   allosteric candidates. Contiguous residues are clustered into regions and the top
   regions (target: 5-20) are reported.
6. Within candidate regions a geometric scan flags residue pairs/triples whose Ca-Ca
   spacing is compatible with engineering a His-pair / His-triplet metal site.

Output: results/<SYS>_candidates.csv (per-residue), results/<SYS>_regions.json (regions
+ His-site suggestions), and a printed summary.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent.parent
PDB_DIR = ROOT / "data" / "pdb"
RESULTS = ROOT / "results"

# ---- analysis configuration: which pairs to compare and how to define the active site
# active_site: {"ligand": <resname in bound structure>} centroid of that HETATM group,
#              or {"residues": [resnum,...]} centroid of those CA in the bound structure.
PAIRS = {
    "GCK": {
        "bound": "1V4S", "unbound": "1V4T", "chain": "A",
        "active_site": {"ligand": "GLC"},  # glucose in the catalytic cleft
        "note": "large open<->closed domain motion; look for hinge + activator-site loops",
    },
    "PTP1B": {
        "bound": "1T49", "unbound": "1SUG", "chain": "A",
        "active_site": {"residues": [215]},  # catalytic Cys215
        "note": "distal allosteric inhibitor site ~20 A from Cys215; WPD loop 177-185",
    },
    "AdK": {
        "bound": "1AKE", "unbound": "4AKE", "chain": "A",
        "active_site": {"ligand": "AP5"},  # Ap5A two-substrate mimic
        "note": "LID + NMP domain closure; hinge/counterweight = de novo Zn-switch target",
    },
}

DISPLACEMENT_MIN = 2.0   # A; residue considered "moved" above this
DISTAL_MIN = 15.0        # A; residue considered distal from active site above this
CORE_REJECT_FRAC = 0.30  # fraction of worst residues dropped each superposition iteration
CORE_ITERS = 5
HIS_PAIR_MIN, HIS_PAIR_MAX = 4.5, 12.0  # Ca-Ca window plausible for an engineered His site


def load_ca(pdb_id, chain_id):
    """Return {resnum: (CA xyz, resname)} for one chain, first model, altloc A/blank."""
    st = PDBParser(QUIET=True).get_structure(pdb_id, PDB_DIR / f"{pdb_id}.pdb")
    model = next(iter(st))
    out = {}
    if chain_id not in model:
        chain_id = next(iter(model)).id
    for res in model[chain_id]:
        if res.id[0] != " ":  # skip HETATM/water
            continue
        if "CA" not in res:
            continue
        ca = res["CA"]
        if ca.is_disordered():
            ca.disordered_select("A") if "A" in ca.disordered_get_id_list() else None
        out[res.id[1]] = (ca.get_coord().astype(float), res.get_resname())
    return out


def ligand_centroid(pdb_id, resname):
    st = PDBParser(QUIET=True).get_structure(pdb_id, PDB_DIR / f"{pdb_id}.pdb")
    coords = [a.get_coord() for a in st.get_atoms()
              if a.get_parent().get_resname() == resname and a.element != "H"]
    if not coords:
        raise ValueError(f"ligand {resname} not found in {pdb_id}")
    return np.mean(coords, axis=0).astype(float)


def kabsch(P, Q):
    """Rotation+translation mapping P onto Q (both Nx3). Returns R, t."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    return R, Qc - R @ Pc


def robust_superpose(mob, ref):
    """Iterative core superposition. mob,ref are Nx3 aligned by index. Returns R,t,core_mask."""
    mask = np.ones(len(mob), dtype=bool)
    R, t = kabsch(mob, ref)
    for _ in range(CORE_ITERS):
        R, t = kabsch(mob[mask], ref[mask])
        dev = np.linalg.norm((mob @ R.T + t) - ref, axis=1)
        keep_n = max(3, int(round(len(mob) * (1 - CORE_REJECT_FRAC))))
        thresh = np.sort(dev)[keep_n - 1]
        mask = dev <= thresh
    return R, t, mask


def cluster(resnums):
    """Group sorted residue numbers into contiguous runs (gap<=1)."""
    runs, cur = [], []
    for r in sorted(resnums):
        if cur and r - cur[-1] > 2:
            runs.append(cur); cur = []
        cur.append(r)
    if cur:
        runs.append(cur)
    return runs


def analyze(sysid, cfg):
    bound, unbound, chain = cfg["bound"], cfg["unbound"], cfg["chain"]
    b = load_ca(bound, chain)
    u = load_ca(unbound, chain)
    common = sorted(set(b) & set(u))
    P = np.array([u[r][0] for r in common])   # move unbound onto bound frame
    Q = np.array([b[r][0] for r in common])
    R, t, core_mask = robust_superpose(P, Q)
    Pf = P @ R.T + t
    disp = np.linalg.norm(Pf - Q, axis=1)

    # active-site reference in bound frame (bound coords are the reference frame)
    a = cfg["active_site"]
    if "ligand" in a:
        asite = ligand_centroid(bound, a["ligand"])
    else:
        asite = np.mean([b[r][0] for r in a["residues"] if r in b], axis=0)
    dist_as = np.linalg.norm(Q - asite, axis=1)

    rows = []
    for i, r in enumerate(common):
        rows.append({
            "resnum": r, "resname": b[r][1],
            "displacement": round(float(disp[i]), 2),
            "dist_active_site": round(float(dist_as[i]), 2),
            "in_core": bool(core_mask[i]),
        })

    # candidates: moved AND distal
    cand = [row for row in rows
            if row["displacement"] >= DISPLACEMENT_MIN
            and row["dist_active_site"] >= DISTAL_MIN]
    cand_nums = [row["resnum"] for row in cand]
    dmap = {row["resnum"]: row for row in rows}

    regions = []
    for run in cluster(cand_nums):
        disps = [dmap[r]["displacement"] for r in run]
        dists = [dmap[r]["dist_active_site"] for r in run]
        regions.append({
            "start": run[0], "end": run[-1], "length": len(run),
            "mean_displacement": round(float(np.mean(disps)), 2),
            "max_displacement": round(float(np.max(disps)), 2),
            "mean_dist_active_site": round(float(np.mean(dists)), 2),
            "residues": run,
        })
    regions.sort(key=lambda x: (x["mean_displacement"] * x["length"]), reverse=True)

    # His-pair/triplet scan within top regions (use bound-frame coords)
    coord = {r: b[r][0] for r in common}
    for reg in regions:
        near = [r for r in coord
                if reg["start"] - 3 <= r <= reg["end"] + 3]
        pairs = []
        for i in range(len(near)):
            for j in range(i + 1, len(near)):
                ri, rj = near[i], near[j]
                if rj - ri < 2:  # avoid i,i+1 (can't both face a metal)
                    continue
                d = float(np.linalg.norm(coord[ri] - coord[rj]))
                if HIS_PAIR_MIN <= d <= HIS_PAIR_MAX:
                    pairs.append({"i": ri, "j": rj, "ca_ca": round(d, 2)})
        pairs.sort(key=lambda p: p["ca_ca"])
        reg["his_pair_candidates"] = pairs[:6]

    return rows, regions, {
        "n_common": len(common), "n_core": int(core_mask.sum()),
        "core_rmsd": round(float(np.sqrt(np.mean(disp[core_mask] ** 2))), 2),
        "overall_rmsd": round(float(np.sqrt(np.mean(disp ** 2))), 2),
        "max_displacement": round(float(disp.max()), 2),
    }


def write_csv(path, rows):
    lines = ["resnum,resname,displacement,dist_active_site,in_core"]
    for r in rows:
        lines.append(f"{r['resnum']},{r['resname']},{r['displacement']},"
                     f"{r['dist_active_site']},{int(r['in_core'])}")
    path.write_text("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("systems", nargs="*", default=list(PAIRS),
                    help="subset of systems to run (default: all monomeric pairs)")
    ap.add_argument("--top", type=int, default=12, help="max regions to report")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)

    summary = {}
    for sysid in (args.systems or list(PAIRS)):
        cfg = PAIRS[sysid]
        rows, regions, stats = analyze(sysid, cfg)
        top = regions[: args.top]
        write_csv(RESULTS / f"{sysid}_candidates.csv", rows)
        (RESULTS / f"{sysid}_regions.json").write_text(
            json.dumps({"system": sysid, "config": cfg, "stats": stats,
                        "regions": top}, indent=2))
        summary[sysid] = {"stats": stats, "n_regions": len(regions)}

        print(f"\n=== {sysid}  ({cfg['unbound']} -> {cfg['bound']}) ===")
        print(f"  {cfg['note']}")
        print(f"  matched {stats['n_common']} residues | core {stats['n_core']} "
              f"| core-RMSD {stats['core_rmsd']} A | overall-RMSD {stats['overall_rmsd']} A "
              f"| max disp {stats['max_displacement']} A")
        print(f"  distal candidate regions (top {len(top)} of {len(regions)}):")
        for k, reg in enumerate(top, 1):
            hp = reg["his_pair_candidates"]
            hp_str = ", ".join(f"{p['i']}/{p['j']}({p['ca_ca']}A)" for p in hp[:3])
            print(f"   {k:2d}. res {reg['start']}-{reg['end']} "
                  f"(len {reg['length']}) meanDisp {reg['mean_displacement']} A, "
                  f"{reg['mean_dist_active_site']} A from active site | "
                  f"His-pairs: {hp_str or '-'}")

    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nWrote per-system CSV/JSON + summary.json to {RESULTS.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
