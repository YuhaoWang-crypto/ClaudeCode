#!/usr/bin/env python3
"""Is the clash problem intrinsic, or a fixable exit-vector problem?

The RGD site of alpha-v/beta-6 sits in a cleft between the alpha-v beta-propeller
and the beta-6 betaI domain. If almost no direction out of that cleft is open,
then NO bulky scaffold can present 4-Delta to the receptor and the whole grafting
concept fails on geometry alone. If a wide cone is open, the clashes seen with a
rigidly superposed graft are an EXIT-VECTOR problem, i.e. a junction-design
problem, which is fixable.

This script measures, in the 4UM9 bound frame:
  1. the escape cone: the fraction of directions out of the bound RGD that a
     ~6 A-radius bulky body could travel along without hitting the receptor,
  2. the natural exit vectors, i.e. where the pro-TGF-beta1 LAP chain itself
     leaves the site (these directions are viable by construction - a whole
     domain hangs off them),
  3. for each candidate, the angle between its scaffold direction and the
     nearest open direction.
"""
import glob
import json
import os
import sys

import mdtraj as md
import numpy as np

MODELS = sys.argv[1] if len(sys.argv) > 1 else "models"
OUT = "results/07_exit_vector.json"
FOURD = "FETLRGDLRILSILRHQNLLKEL"
MOTIF_IDX = list(range(4, 11))
LAP_MOTIF_RESSEQ = [241, 242, 243, 244, 245, 246, 247]
BB = ("N", "CA", "C")
BODY_RADIUS = 6.0      # a domain-sized probe
PROBE_LEN = 25.0       # how far out it has to clear

ref = md.load("data/4UM9.pdb")
ref = ref.atom_slice(ref.topology.select("protein"))
rtop = ref.topology
rxyz = ref.xyz[0] * 10.0
rch = {ch.index: [r for r in ch.residues if r.is_protein] for ch in rtop.chains}
LAP_CI = max((i for i, rs in rch.items() if len(rs) <= 15), key=lambda i: len(rch[i]))
lap_res = {r.resSeq: r for r in rch[LAP_CI]}
lap_cent = rxyz[[a.index for a in rtop.chain(LAP_CI).atoms]].mean(axis=0)
big = [i for i, rs in rch.items() if len(rs) > 400]
RECEPTOR_CIS = sorted(big, key=lambda ci: float(np.linalg.norm(
    rxyz[[a.index for a in rtop.chain(ci).atoms]] - lap_cent, axis=1).min()))[:2]
rec_atoms = [a.index for ci in RECEPTOR_CIS for a in rtop.chain(ci).atoms
             if a.element.symbol != "H"]
REC = rxyz[rec_atoms]


def atom(res, name, xyz=rxyz):
    a = next((x for x in res.atoms if x.name == name), None)
    return None if a is None else xyz[a.index]


def bb_coords(res_list, xyz):
    out = []
    for r in res_list:
        for nm in BB:
            a = next((x for x in r.atoms if x.name == nm), None)
            if a is None:
                return None
            out.append(xyz[a.index])
    return np.array(out)


LAP_TARGET = bb_coords([lap_res[s] for s in LAP_MOTIF_RESSEQ], rxyz)
RGD_CENT = np.array([atom(lap_res[s], "CA") for s in (241, 242, 243)]).mean(axis=0)
print(f"bound RGD centroid in 4UM9 frame: {RGD_CENT.round(1)}")

# ------------------------------------------------------------------ 1. escape cone
def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


DIRS = fibonacci_sphere(2000)
rel = REC - RGD_CENT
dist = np.linalg.norm(rel, axis=1)
near = rel[dist < PROBE_LEN + BODY_RADIUS]
free_mask = np.zeros(len(DIRS), dtype=bool)
for k, d in enumerate(DIRS):
    proj = near @ d
    sel = (proj > 2.0) & (proj < PROBE_LEN)       # ahead of the site
    if not sel.any():
        free_mask[k] = True
        continue
    perp = np.linalg.norm(near[sel] - np.outer(proj[sel], d), axis=1)
    free_mask[k] = bool((perp > BODY_RADIUS).all())

free_frac = float(free_mask.mean())
print(f"escape cone: {free_frac*100:.1f}% of directions let a {BODY_RADIUS:.0f} A-radius "
      f"body clear {PROBE_LEN:.0f} A  ({int(free_mask.sum())}/{len(DIRS)} directions)")

# ------------------------------------------------------------------ 2. natural exit vectors
nat = {}
for label, resseq in (("N_side", 238), ("C_side", 249)):
    r = lap_res.get(resseq)
    if r is None:
        continue
    c = atom(r, "CA")
    if c is None:
        continue
    v = c - RGD_CENT
    v = v / np.linalg.norm(v)
    ang_free = np.degrees(np.arccos(np.clip(DIRS[free_mask] @ v, -1, 1)))
    nat[label] = {
        "lap_residue": resseq,
        "unit_vector": [round(float(x), 3) for x in v],
        "is_itself_open": bool(free_mask[np.argmax(DIRS @ v)]),
        "angle_to_nearest_open_direction_deg": round(float(ang_free.min()), 1),
    }
    print(f"natural {label} exit (LAP {resseq}): nearest open direction "
          f"{nat[label]['angle_to_nearest_open_direction_deg']} deg away")

FREE_DIRS = DIRS[free_mask]

# ------------------------------------------------------------------ 3. per candidate
def kabsch(mobile, target):
    mc, tc = mobile.mean(0), target.mean(0)
    m, t = mobile - mc, target - tc
    u, s, vt = np.linalg.svd(m.T @ t)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return r, mc, tc


def evaluate(cid, payload):
    rows = []
    for f in sorted(glob.glob(os.path.join(MODELS, cid, "*.cif"))):
        t = md.load(f)
        t = t.atom_slice(t.topology.select("protein"))
        top = t.topology
        xyz = t.xyz[0] * 10.0
        hit = None
        for ch in top.chains:
            rs = [r for r in ch.residues if r.is_protein]
            p = "".join(r.code or "X" for r in rs).find(payload)
            if p >= 0:
                hit = (rs, p)
                break
        if hit is None:
            continue
        rs, p0 = hit
        pay = rs[p0:p0 + len(payload)]
        idx = [i for i in MOTIF_IDX if i < len(payload)]
        mob = bb_coords([pay[i] for i in idx], xyz)
        if mob is None or len(mob) != len(LAP_TARGET):
            continue
        r, mc, tc = kabsch(mob, LAP_TARGET)
        moved = (xyz - mc) @ r + tc

        pay_set = {x.index for x in pay}
        # the folded scaffold: everything that is not payload and not a flexible
        # G/S/P junction residue contiguous with it
        flex = set()
        pos = [k for k, x in enumerate(rs) if x.index in pay_set]
        if pos:
            k = min(pos) - 1
            while k >= 0 and (rs[k].code or "X") in "GSP":
                flex.add(rs[k].index)
                k -= 1
            k = max(pos) + 1
            while k < len(rs) and (rs[k].code or "X") in "GSP":
                flex.add(rs[k].index)
                k += 1
        core_ca = [a.index for a in top.atoms if a.name == "CA"
                   and a.residue.index not in pay_set and a.residue.index not in flex]
        if not core_ca:
            continue
        body = moved[core_ca].mean(axis=0)
        v = body - RGD_CENT
        dist_body = float(np.linalg.norm(v))
        v = v / dist_body
        ang = float(np.degrees(np.arccos(np.clip(FREE_DIRS @ v, -1, 1))).min())
        rows.append({
            "file": os.path.basename(f),
            "scaffold_direction_unit": [round(float(x), 3) for x in v],
            "scaffold_centroid_distance_A": round(dist_body, 1),
            "angle_to_nearest_open_direction_deg": round(ang, 1),
        })
    return rows


if __name__ == "__main__":
    cand = json.load(open("results/03_candidates.json"))
    by_id = {c["id"]: c for c in cand["candidates"] + cand["controls"]}
    out = {}
    for cid, c in by_id.items():
        payload = c.get("payload")
        if payload is None:
            if cid in ("REF_4delta_free", "BM_CH3_Cterm_fusion"):
                payload = FOURD
            else:
                continue
        rows = evaluate(cid, payload)
        if not rows:
            continue
        ang = [r["angle_to_nearest_open_direction_deg"] for r in rows]
        dd = [r["scaffold_centroid_distance_A"] for r in rows]
        out[cid] = {
            "n_samples": len(rows),
            "angle_to_open_deg": {"mean": round(float(np.mean(ang)), 1),
                                  "min": round(float(np.min(ang)), 1)},
            "scaffold_centroid_distance_A": {"mean": round(float(np.mean(dd)), 1)},
            "per_sample": rows,
        }
        print(f"{cid:<22} scaffold {np.mean(dd):>5.1f} A from RGD, "
              f"{np.mean(ang):>5.1f} deg off the nearest open direction "
              f"(best {np.min(ang):.1f})")

    json.dump({
        "escape_cone": {
            "body_radius_A": BODY_RADIUS, "clearance_required_A": PROBE_LEN,
            "free_direction_fraction": round(free_frac, 4),
            "n_directions_sampled": len(DIRS),
            "interpretation": ("fraction of directions out of the receptor-bound RGD along "
                               "which a domain-sized body can escape the alpha-v/beta-6 cleft"),
        },
        "natural_exit_vectors": nat,
        "candidates": out,
    }, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")
