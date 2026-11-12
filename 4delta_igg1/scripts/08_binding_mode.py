#!/usr/bin/env python3
"""Tier 2: did the predicted complex reproduce real alpha-v/beta-6 RGD chemistry?

Integrin RGD recognition has three experimentally defined hallmarks:
  1. the ligand Asp carboxylate directly coordinates the beta-subunit MIDAS
     divalent metal (~2.0-2.3 A),
  2. the ligand Arg guanidinium salt-bridges alpha-v Asp218 in the beta-propeller,
  3. for alpha-v/beta-6 specifically the RGDLxxL helix packs against the beta-6
     specificity-determining loop.

Rather than trusting ipTM, this script measures all three, and also superposes on
the RECEPTOR (not the ligand) to compare the predicted ligand pose with the LAP
crystal pose in 4UM9. A method-validation run (the LAP peptide itself) and a
specificity control (RGD->RGE) bound the scale.
"""
import glob
import json
import os
import sys

import mdtraj as md
import numpy as np

MODELS = sys.argv[1] if len(sys.argv) > 1 else "models"
OUT = "results/08_binding_mode.json"
ALPHA_V_ASP = 218            # alpha-v Asp218, canonical Arg-binding carboxylate

# ------------------------------------------------------------------ crystal reference
ref = md.load("data/4UM9.pdb")
rprot = ref.atom_slice(ref.topology.select("protein"))
rtop = rprot.topology
rxyz = rprot.xyz[0] * 10.0
rchains = {ch.index: [r for r in ch.residues if r.is_protein] for ch in rtop.chains}
LAP_CI = max((i for i, rs in rchains.items() if len(rs) <= 15),
             key=lambda i: len(rchains[i]))
lap_res = {r.resSeq: r for r in rchains[LAP_CI]}
# the receptor copy nearest that LAP peptide
lap_cent = rxyz[[a.index for a in rtop.chain(LAP_CI).atoms]].mean(axis=0)
big = [i for i, rs in rchains.items() if len(rs) > 400]
near = sorted(big, key=lambda ci: float(np.linalg.norm(
    rxyz[[a.index for a in rtop.chain(ci).atoms]] - lap_cent, axis=1).min()))[:2]
# alpha-v is the larger of the two
ALPHA_CI, BETA_CI = sorted(near, key=lambda ci: -len(rchains[ci]))
print(f"crystal reference: LAP chain {LAP_CI}, alpha-v chain {ALPHA_CI} "
      f"({len(rchains[ALPHA_CI])} aa), beta-6 chain {BETA_CI} ({len(rchains[BETA_CI])} aa)")


def res_atom(residue, name, xyz):
    a = next((x for x in residue.atoms if x.name == name), None)
    return None if a is None else xyz[a.index]


# crystal pose of the LAP motif, in receptor frame
LAP_MOTIF = [241, 242, 243, 244, 245, 246, 247]
lap_ca_cryst = np.array([res_atom(lap_res[s], "CA", rxyz) for s in LAP_MOTIF])
# receptor CA frame used for superposition
ref_rec_ca = {}
for ci, tag in ((ALPHA_CI, "V"), (BETA_CI, "W")):
    for r in rchains[ci]:
        c = res_atom(r, "CA", rxyz)
        if c is not None:
            ref_rec_ca[(tag, r.resSeq)] = c

# crystal hallmark distances, as the yardstick
av_asp = next((r for r in rchains[ALPHA_CI] if r.resSeq == ALPHA_V_ASP), None)
cryst = {}
if av_asp is not None:
    arg = lap_res[241]
    cz = res_atom(arg, "CZ", rxyz)
    best = None
    for nm in ("OD1", "OD2"):
        o = res_atom(av_asp, nm, rxyz)
        if o is not None and cz is not None:
            d = float(np.linalg.norm(cz - o))
            best = d if best is None else min(best, d)
    cryst["arg_to_alphaV_Asp218_A"] = round(best, 2) if best else None
    print(f"  crystal: LAP Arg241 CZ -> alpha-v Asp{ALPHA_V_ASP} carboxylate = "
          f"{cryst.get('arg_to_alphaV_Asp218_A')} A")

# crystal MIDAS metal: the metal in the original file closest to the LAP Asp
full = md.load("data/4UM9.pdb")
ftop = full.topology
fxyz = full.xyz[0] * 10.0
metals = [a.index for a in ftop.atoms
          if a.element.symbol in ("Mg", "Ca", "Mn") and a.residue.name in ("MG", "CA", "MN")]
asp_o = [res_atom(lap_res[243], nm, rxyz) for nm in ("OD1", "OD2")]
asp_o = [x for x in asp_o if x is not None]
if metals and asp_o:
    dm = min(float(np.linalg.norm(fxyz[m] - o)) for m in metals for o in asp_o)
    cryst["asp_to_MIDAS_metal_A"] = round(dm, 2)
    print(f"  crystal: LAP Asp243 carboxylate -> nearest metal = {dm:.2f} A")


def kabsch_apply(mobile_pairs, target_pairs, xyz_all):
    m = np.array(mobile_pairs)
    t = np.array(target_pairs)
    mc, tc = m.mean(0), t.mean(0)
    u, s, vt = np.linalg.svd((m - mc).T @ (t - tc))
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    rms = float(np.sqrt((((m - mc) @ r - (t - tc)) ** 2).sum(1).mean()))
    return (xyz_all - mc) @ r + tc, rms


def analyse(cid, ligand_seq, motif_offset=4):
    """motif_offset = index of the Arg of R-G-D inside ligand_seq"""
    files = sorted(glob.glob(os.path.join(MODELS, cid, "*.cif")))
    rows = []
    for f in files:
        t = md.load(f)
        top = t.topology
        xyz = t.xyz[0] * 10.0
        # identify chains by sequence length / content
        lig_ch = rec_chs = None
        chs = {}
        for ch in top.chains:
            rs = [r for r in ch.residues if r.is_protein]
            if rs:
                chs[ch.index] = rs
        for ci, rs in chs.items():
            if "".join(r.code or "X" for r in rs) == ligand_seq:
                lig_ch = ci
                break
        if lig_ch is None:
            continue
        rec_chs = [ci for ci in chs if ci != lig_ch and len(chs[ci]) > 200]
        if len(rec_chs) < 2:
            continue
        alpha_ci, beta_ci = sorted(rec_chs, key=lambda ci: -len(chs[ci]))

        lig = chs[lig_ch]
        i_r, i_g, i_d = motif_offset, motif_offset + 1, motif_offset + 2
        arg_cz = res_atom(lig[i_r], "CZ", xyz)
        acid = lig[i_d]
        acid_o = [res_atom(acid, nm, xyz) for nm in ("OD1", "OD2", "OE1", "OE2")]
        acid_o = [x for x in acid_o if x is not None]

        # 1. acidic carboxylate -> MIDAS metal
        met = [a.index for a in top.atoms
               if a.element.symbol in ("Mg", "Mn", "Ca")
               and a.residue.name in ("MG", "MN", "CA")]
        d_metal = (min(float(np.linalg.norm(xyz[m] - o)) for m in met for o in acid_o)
                   if met and acid_o else None)

        # 2. Arg -> alpha-v Asp218
        av_res = next((r for r in chs[alpha_ci] if r.resSeq == ALPHA_V_ASP), None)
        d_arg = None
        if av_res is not None and arg_cz is not None:
            cands = [res_atom(av_res, nm, xyz) for nm in ("OD1", "OD2")]
            cands = [x for x in cands if x is not None]
            if cands:
                d_arg = min(float(np.linalg.norm(arg_cz - c)) for c in cands)

        # 3. receptor-frame superposition -> is the ligand in the crystal pose?
        mob, tgt = [], []
        for ci, tag in ((alpha_ci, "V"), (beta_ci, "W")):
            for r in chs[ci]:
                key = (tag, r.resSeq)
                if key in ref_rec_ca:
                    c = res_atom(r, "CA", xyz)
                    if c is not None:
                        mob.append(c)
                        tgt.append(ref_rec_ca[key])
        pose_rmsd = rec_rmsd = None
        if len(mob) > 200:
            moved, rec_rmsd = kabsch_apply(mob, tgt, xyz)
            mca = []
            for k in range(7):
                j = motif_offset + k
                if j < len(lig):
                    a = next((x for x in lig[j].atoms if x.name == "CA"), None)
                    if a is not None:
                        mca.append(moved[a.index])
            if len(mca) == 7:
                pose_rmsd = float(np.sqrt(
                    ((np.array(mca) - lap_ca_cryst) ** 2).sum(1).mean()))

        # ligand/receptor interface size
        lig_at = [a.index for a in top.atoms
                  if a.residue.index in {r.index for r in lig} and a.element.symbol != "H"]
        rec_at = [a.index for a in top.atoms
                  if a.residue.index in {r.index for ci in (alpha_ci, beta_ci) for r in chs[ci]}
                  and a.element.symbol != "H"]
        dd = np.linalg.norm(xyz[lig_at][:, None, :] - xyz[rec_at][None, :, :], axis=2)
        rows.append({
            "file": os.path.basename(f),
            "acid_to_MIDAS_metal_A": round(d_metal, 2) if d_metal else None,
            "arg_to_alphaV_Asp218_A": round(d_arg, 2) if d_arg else None,
            "ligand_pose_rmsd_vs_LAP_crystal_A": round(pose_rmsd, 2) if pose_rmsd else None,
            "receptor_superposition_rmsd_A": round(rec_rmsd, 2) if rec_rmsd else None,
            "ligand_receptor_contacts_4.5A": int((dd < 4.5).sum()),
        })
    return rows


if __name__ == "__main__":
    spec = json.load(open("results/04_job_ids.json"))["tier2_integrin_complexes"]
    out = {"crystal_reference": cryst, "results": {}}
    for name, meta in spec.items():
        if name.startswith("_"):
            continue
        lig = meta["ligand"]
        off = lig.find("RGD")
        if off < 0:
            off = lig.find("RGE")
        rows = analyse(name, lig, motif_offset=off)
        if not rows:
            print(f"{name}: no models yet")
            continue
        def best(key, lo=True):
            v = [r[key] for r in rows if r.get(key) is not None]
            if not v:
                return None
            return round(float(min(v) if lo else max(v)), 2)
        out["results"][name] = {
            "ligand": lig, "purpose": meta["purpose"], "n_samples": len(rows),
            "best_acid_to_MIDAS_metal_A": best("acid_to_MIDAS_metal_A"),
            "best_arg_to_alphaV_Asp218_A": best("arg_to_alphaV_Asp218_A"),
            "best_pose_rmsd_vs_crystal_A": best("ligand_pose_rmsd_vs_LAP_crystal_A"),
            "max_contacts_4.5A": best("ligand_receptor_contacts_4.5A", lo=False),
            "per_sample": rows,
        }
        r = out["results"][name]
        print(f"{name:<22} Asp->metal {str(r['best_acid_to_MIDAS_metal_A']):>6}  "
              f"Arg->avD218 {str(r['best_arg_to_alphaV_Asp218_A']):>6}  "
              f"pose-RMSD {str(r['best_pose_rmsd_vs_crystal_A']):>6}  "
              f"contacts {r['max_contacts_4.5A']}")
    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")
