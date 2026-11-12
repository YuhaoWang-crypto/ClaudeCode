#!/usr/bin/env python3
"""Task 3b: docking-free test of whether a real integrin can actually engage the
grafted peptide.

Rationale. 4-Delta carries RGDLRIL. The alpha-v/beta-6 co-crystal 4UM9 contains
the pro-TGF-beta1 LAP peptide GRGDLGRL bound in the alpha-v/beta-6 cleft, i.e.
the same RGDLxxL specificity motif in its RECEPTOR-BOUND conformation. So the
bound-state geometry of the motif is known experimentally.

For each predicted candidate we therefore:
  1. least-squares superpose the candidate's RGDLRIL backbone onto the LAP
     RGDLGRL backbone in the 4UM9 complex (fit RMSD = how bound-competent the
     grafted motif's conformation is),
  2. apply that transform to the whole candidate, placing the antibody module
     exactly where it would sit if its RGD engaged the integrin,
  3. count heavy-atom clashes and contacts between the antibody scaffold and the
     alpha-v/beta-6 headpiece.

Zero clashes = the integrin can reach this graft. Many clashes = the motif can
only be presented in a pose where the antibody collides with the receptor, which
no amount of affinity can rescue.
"""
import glob
import json
import os
import sys

import mdtraj as md
import numpy as np

MODELS = sys.argv[1] if len(sys.argv) > 1 else "models"
OUT = "results/06_integrin_steric.json"

FOURD = "FETLRGDLRILSILRHQNLLKEL"
MOTIF_IDX = list(range(4, 11))        # R G D L R I L  (payload 0-based)
LAP_MOTIF_RESSEQ = [241, 242, 243, 244, 245, 246, 247]   # R G D L G R L in 4UM9
BB = ("N", "CA", "C")
CLASH_A = 3.0
CONTACT_A = 4.5

# ------------------------------------------------------------------ reference complex
ref = md.load("data/4UM9.pdb")
ref = ref.atom_slice(ref.topology.select("protein"))
rtop = ref.topology
rxyz = ref.xyz[0] * 10.0

chains = {}
for ch in rtop.chains:
    rs = [r for r in ch.residues if r.is_protein]
    if rs:
        chains[ch.index] = rs
lap_chains = [i for i, rs in chains.items() if len(rs) <= 15]
big_chains = [i for i, rs in chains.items() if len(rs) > 400]

def chain_atoms(ci):
    return [a.index for a in rtop.chain(ci).atoms if a.element.symbol != "H"]

# pick the LAP copy with the most complete motif, then the two large chains
# closest to it -> one alpha-v/beta-6/LAP functional unit
best = None
for lc in lap_chains:
    have = {r.resSeq: r for r in chains[lc]}
    if not all(s in have for s in LAP_MOTIF_RESSEQ):
        continue
    lap_at = chain_atoms(lc)
    cent = rxyz[lap_at].mean(axis=0)
    dists = []
    for bc in big_chains:
        at = chain_atoms(bc)
        dists.append((float(np.linalg.norm(rxyz[at] - cent, axis=1).min()), bc))
    dists.sort()
    near = [bc for _, bc in dists[:2]]
    score = sum(d for d, _ in dists[:2])
    if best is None or score < best[0]:
        best = (score, lc, near)

_, LAP_CI, RECEPTOR_CIS = best
print(f"reference unit: LAP chain idx {LAP_CI}, receptor chain idxs {RECEPTOR_CIS}")
print(f"  receptor sizes: {[len(chains[c]) for c in RECEPTOR_CIS]}")

lap_res = {r.resSeq: r for r in chains[LAP_CI]}
lap_bb = []
for s in LAP_MOTIF_RESSEQ:
    r = lap_res[s]
    for nm in BB:
        a = next((x for x in r.atoms if x.name == nm), None)
        if a is None:
            raise SystemExit(f"LAP residue {s} missing {nm}")
        lap_bb.append(a.index)
LAP_TARGET = rxyz[lap_bb]                         # (21, 3)
print(f"LAP motif sequence: {''.join(lap_res[s].code for s in LAP_MOTIF_RESSEQ)}")

receptor_atoms = [i for ci in RECEPTOR_CIS for i in chain_atoms(ci)]
RECEPTOR_XYZ = rxyz[receptor_atoms]
# the receptor's own contacts with LAP, as the scale for "a real interface"
d_ref = np.linalg.norm(rxyz[chain_atoms(LAP_CI)][:, None, :] - RECEPTOR_XYZ[None, :, :], axis=2)
print(f"LAP/receptor reference: {int((d_ref < CONTACT_A).sum())} contacts, "
      f"{int((d_ref < CLASH_A).sum())} pairs < {CLASH_A} A")


def kabsch(mobile, target):
    mc, tc = mobile.mean(0), target.mean(0)
    m, t = mobile - mc, target - tc
    u, s, vt = np.linalg.svd(m.T @ t)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    rms = float(np.sqrt(((m @ r - t) ** 2).sum(1).mean()))
    return r, mc, tc, rms


def evaluate(cid, payload):
    files = sorted(glob.glob(os.path.join(MODELS, cid, "*.cif")))
    rows = []
    for f in files:
        t = md.load(f)
        t = t.atom_slice(t.topology.select("protein"))
        top = t.topology
        xyz = t.xyz[0] * 10.0
        # locate payload
        hit = None
        for ch in top.chains:
            rs = [r for r in ch.residues if r.is_protein]
            s = "".join(r.code or "X" for r in rs)
            p = s.find(payload)
            if p >= 0:
                hit = (rs, p)
                break
        if hit is None:
            continue
        rs, p0 = hit
        pay_res = rs[p0:p0 + len(payload)]
        idx = [i for i in MOTIF_IDX if i < len(payload)]
        mob_bb = []
        ok = True
        for i in idx:
            for nm in BB:
                a = next((x for x in pay_res[i].atoms if x.name == nm), None)
                if a is None:
                    ok = False
                    break
                mob_bb.append(a.index)
            if not ok:
                break
        if not ok or len(mob_bb) != len(LAP_TARGET):
            continue

        r, mc, tc, fit_rmsd = kabsch(xyz[mob_bb], LAP_TARGET)
        moved = (xyz - mc) @ r + tc

        pay_set = {x.index for x in pay_res}
        # Gly/Ser junction residues flanking the payload are flexible and will
        # rearrange, so a second, fairer count excludes them and scores only the
        # folded Ig core against the receptor.
        flex = set()
        for r in rs:
            if r.index in pay_set:
                continue
            # walk outward from the payload through any contiguous G/S/P stretch
            pass
        chain_rs = rs
        pay_positions = [k for k, r in enumerate(chain_rs) if r.index in pay_set]
        if pay_positions:
            lo_p, hi_p = min(pay_positions), max(pay_positions)
            k = lo_p - 1
            while k >= 0 and (chain_rs[k].code or "X") in "GSP":
                flex.add(chain_rs[k].index)
                k -= 1
            k = hi_p + 1
            while k < len(chain_rs) and (chain_rs[k].code or "X") in "GSP":
                flex.add(chain_rs[k].index)
                k += 1

        scaffold = [a.index for a in top.atoms
                    if a.residue.index not in pay_set and a.element.symbol != "H"]
        core = [a.index for a in top.atoms
                if a.residue.index not in pay_set and a.residue.index not in flex
                and a.element.symbol != "H"]
        if scaffold:
            d = np.linalg.norm(moved[scaffold][:, None, :] - RECEPTOR_XYZ[None, :, :], axis=2)
            clashes = int((d < CLASH_A).sum())
            contacts = int((d < CONTACT_A).sum())
            nearest = float(d.min())
            clash_res = len({top.atom(scaffold[i]).residue.index
                             for i in np.where((d < CLASH_A).any(axis=1))[0]})
            dc = np.linalg.norm(moved[core][:, None, :] - RECEPTOR_XYZ[None, :, :], axis=2)
            core_clashes = int((dc < CLASH_A).sum())
            core_clash_res = len({top.atom(core[i]).residue.index
                                  for i in np.where((dc < CLASH_A).any(axis=1))[0]})
            core_nearest = float(dc.min())
        else:
            clashes = contacts = clash_res = core_clashes = core_clash_res = 0
            nearest = core_nearest = float("nan")

        rows.append({
            "file": os.path.basename(f),
            "motif_fit_rmsd_A": round(fit_rmsd, 3),
            "scaffold_receptor_clashes": clashes,
            "scaffold_receptor_clash_residues": clash_res,
            "scaffold_receptor_contacts": contacts,
            "nearest_scaffold_receptor_A": round(nearest, 2),
            "core_receptor_clashes": core_clashes,
            "core_receptor_clash_residues": core_clash_res,
            "nearest_core_receptor_A": round(core_nearest, 2),
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
            print(f"{cid}: no models")
            continue
        fit = [r["motif_fit_rmsd_A"] for r in rows]
        cl = [r["scaffold_receptor_clashes"] for r in rows]
        ct = [r["scaffold_receptor_contacts"] for r in rows]
        cc = [r["core_receptor_clashes"] for r in rows]
        out[cid] = {
            "n_samples": len(rows),
            "core_receptor_clashes": {"mean": round(float(np.mean(cc)), 1),
                                      "min": int(np.min(cc)), "max": int(np.max(cc))},
            "n_samples_core_clash_free": int(sum(1 for x in cc if x == 0)),
            "motif_fit_rmsd_A": {"mean": round(float(np.mean(fit)), 3),
                                 "min": round(float(np.min(fit)), 3),
                                 "sd": round(float(np.std(fit)), 3)},
            "scaffold_receptor_clashes": {"mean": round(float(np.mean(cl)), 1),
                                          "min": int(np.min(cl)), "max": int(np.max(cl))},
            "scaffold_receptor_contacts": {"mean": round(float(np.mean(ct)), 1)},
            "n_samples_clash_free": int(sum(1 for x in cl if x == 0)),
            "per_sample": rows,
        }
        print(f"{cid:<22} motif-fit {np.mean(fit):>5.2f} (best {np.min(fit):>4.2f})  "
              f"all-clash {np.mean(cl):>7.1f}  core-clash {np.mean(cc):>7.1f}  "
              f"core-free {out[cid]['n_samples_core_clash_free']}/{len(rows)}")

    json.dump({
        "reference": {
            "complex": "4UM9 alpha-v/beta-6 + pro-TGF-beta1 LAP peptide",
            "lap_motif": "RGDLGRL (resSeq 241-247)",
            "candidate_motif": "RGDLRIL (4-Delta payload indices 4-10)",
            "clash_cutoff_A": CLASH_A, "contact_cutoff_A": CONTACT_A,
            "free_vs_bound_motif_rmsd_A": {
                "value": "1.77-2.43 (mean 2.08) over 15 NMR models",
                "meaning": ("the experimental FREE-state motif (6R2X CgA43-49) already sits "
                            "~2.1 A from the receptor-BOUND pose, because free 4-Delta has its "
                            "RGD inside the helix (DSSP SHHHHH) whereas the bound LAP motif has "
                            "RGD in a turn with the helix starting at +4 (DSSP STTHHH). "
                            "alpha-v/beta-6 engagement therefore requires the RGD to unwind "
                            "locally out of the helix. A candidate motif_fit_rmsd near 2.1 A is "
                            "behaving like the free peptide; this is the calibration baseline, "
                            "NOT zero."),
            },
            "core_clash_definition": ("clashes counted after excluding the payload and the "
                                      "flexible Gly/Ser/Pro junction residues flanking it, "
                                      "since those will rearrange; this scores the folded Ig "
                                      "core against the receptor"),
        },
        "results": out,
    }, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")
