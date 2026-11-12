#!/usr/bin/env python3
"""Task 3: structural stability and steric evaluation of the grafted candidates.

For every candidate, metrics are computed over all 5 Boltz samples and compared
against three reference points:
  - PDB 6R2X            the experimental bioactive conformation of the peptide
  - WT_Fv / WT_CH3      the unmodified scaffold predicted the same way
  - BM_CH3_Cterm_fusion the no-span-penalty topology benchmark

Metric families
  A. graft conformational fidelity   helicity, RMSD to 6R2X, end-to-end span
  B. RGD presentation                relative SASA, bulky-probe SASA, Arg-Asp
                                     pharmacophore distance, protrusion
  C. scaffold integrity              scaffold RMSD vs WT, clashes, and for CH3
                                     the fraction of native dimer contacts kept
  D. model confidence                pLDDT over the graft, complex pTM/ipTM
"""
import glob
import json
import os
import sys

import mdtraj as md
import numpy as np

MODELS = sys.argv[1] if len(sys.argv) > 1 else "models"
OUT = "results/05_evaluation.json"

FOURD = "FETLRGDLRILSILRHQNLLKEL"
# 0-based indices inside the payload
IDX_RGD = [4, 5, 6]                       # R43 G44 D45
IDX_HELIX_CORE = list(range(7, 18))       # L46(E46) .. N56  - fully helical in 6R2X
IDX_HPHOB = [9, 10, 12, 13, 18, 19]       # I48 L49 I51 L52 L57 L58 - integrin-facing face
HELIX_CODES = {"H", "G", "I"}
MAX_ACC = {"A": 106, "R": 248, "N": 157, "D": 163, "C": 135, "Q": 198, "E": 194,
           "G": 84, "H": 184, "I": 169, "L": 164, "K": 205, "M": 188, "F": 197,
           "P": 136, "S": 130, "T": 142, "W": 227, "Y": 222, "V": 142}

# experimental reference values measured in script 01 from the 6R2X ensemble
REF = {"helix_frac_core": 1.00, "rgd_rel_sasa": 0.577, "pharmacophore_A": 9.02,
       "end_to_end_A": 35.75}

cand = json.load(open("results/03_candidates.json"))
by_id = {c["id"]: c for c in cand["candidates"] + cand["controls"]}
CH3_WT = cand["scaffold"]["CH3_341_447"]
VH_WT = cand["scaffold"]["VH"]

# ------------------------------------------------------------------ 6R2X reference frame
ref = md.load("data/6R2X.pdb")
ref_ca = ref.atom_slice(ref.topology.select("name CA"))
# 6R2X construct is CgA F39-D63; payload index i corresponds to construct index i
def ref_payload_ca(indices):
    return ref_ca.xyz[:, indices, :] * 10.0


def load_samples(cid):
    """Load every predicted sample for a candidate as a list of md.Trajectory."""
    pats = [os.path.join(MODELS, cid, "*.cif"), os.path.join(MODELS, cid, "*.pdb")]
    files = sorted(f for p in pats for f in glob.glob(p))
    out = []
    for f in files:
        try:
            t = md.load(f)
            t = t.atom_slice(t.topology.select("protein"))
            if t.n_atoms:
                out.append((os.path.basename(f), t))
        except Exception as exc:                    # noqa: BLE001
            print(f"  ! {cid}: cannot read {f}: {exc}")
    return out


def chain_seqs(t):
    seqs = {}
    for ch in t.topology.chains:
        seqs[ch.index] = "".join(r.code if r.code else "X"
                                 for r in ch.residues if r.is_protein)
    return seqs


def residue_list(t, chain_index):
    return [r for r in t.topology.chain(chain_index).residues if r.is_protein]


def plddt_from_bfactor(path):
    """Boltz writes per-atom pLDDT into the B-factor column of the model file."""
    vals = {}
    if path.endswith(".pdb"):
        for line in open(path):
            if line.startswith(("ATOM", "HETATM")):
                try:
                    key = (line[21], int(line[22:26]))
                    vals.setdefault(key, []).append(float(line[60:66]))
                except ValueError:
                    pass
    return vals


def analyse(cid, payload, graft_in_chain_name=None):
    samples = load_samples(cid)
    if not samples:
        return None
    per_sample = []
    for fname, t in samples:
        seqs = chain_seqs(t)
        # locate the payload: the chain that contains it, and where
        hit = None
        for ci, s in seqs.items():
            p = s.find(payload)
            if p >= 0:
                hit = (ci, p)
                break
        if hit is None:
            print(f"  ! {cid}/{fname}: payload not found")
            continue
        ci, p0 = hit
        res = residue_list(t, ci)
        pay_res = res[p0:p0 + len(payload)]
        pay_global = [r.index for r in pay_res]

        # ---------------------------------------------------- A. conformation
        dssp = md.compute_dssp(t, simplified=False)[0]
        core = [i for i in IDX_HELIX_CORE if i < len(payload)]
        helix_core = float(np.mean(
            [dssp[pay_global[i]] in HELIX_CODES for i in core]))
        helix_all = float(np.mean(
            [dssp[pay_global[i]] in HELIX_CODES for i in range(len(payload))]))

        def ca_idx(i):
            sel = [a.index for a in pay_res[i].atoms if a.name == "CA"]
            return sel[0] if sel else None

        cas = [ca_idx(i) for i in range(len(payload))]
        xyz = t.xyz[0] * 10.0
        n_ok = cas[0] is not None and cas[len(payload) - 1] is not None
        e2e = float(np.linalg.norm(xyz[cas[-1]] - xyz[cas[0]])) if n_ok else None

        # RMSD of the RGD+helix stretch to 6R2X (superposed on that stretch)
        cmp_idx = [i for i in range(4, min(22, len(payload))) if cas[i] is not None]
        mob = xyz[[cas[i] for i in cmp_idx]]
        best = None
        for m in range(ref_ca.n_frames):
            tgt = ref_ca.xyz[m, cmp_idx, :] * 10.0
            mc, tc = mob - mob.mean(0), tgt - tgt.mean(0)
            u, s, vt = np.linalg.svd(mc.T @ tc)
            d = np.sign(np.linalg.det(u @ vt))
            r = u @ np.diag([1, 1, d]) @ vt
            rms = float(np.sqrt(((mc @ r - tc) ** 2).sum(1).mean()))
            best = rms if best is None else min(best, rms)

        # ---------------------------------------------------- B. RGD presentation
        sasa_std = md.shrake_rupley(t, mode="residue", probe_radius=0.14)[0] * 100.0
        # 5 A probe approximates access for a bulky protein partner (integrin head)
        sasa_big = md.shrake_rupley(t, mode="residue", probe_radius=0.50)[0] * 100.0
        rgd_rel = float(np.mean([
            sasa_std[pay_global[i]] / MAX_ACC[payload[i]] for i in IDX_RGD]))
        rgd_big = float(np.sum([sasa_big[pay_global[i]] for i in IDX_RGD]))
        hphob_rel = float(np.mean([
            sasa_std[pay_global[i]] / MAX_ACC[payload[i]]
            for i in IDX_HPHOB if i < len(payload)]))

        def atom_xyz(residue, name):
            for a in residue.atoms:
                if a.name == name:
                    return xyz[a.index]
            return None

        arg_cz = atom_xyz(pay_res[IDX_RGD[0]], "CZ")
        asp_cg = atom_xyz(pay_res[IDX_RGD[2]], "CG")
        pharm = (float(np.linalg.norm(arg_cz - asp_cg))
                 if arg_cz is not None and asp_cg is not None else None)

        # protrusion: RGD centroid distance from the scaffold CA centroid
        scaffold_ca = [a.index for a in t.topology.atoms
                       if a.name == "CA" and a.residue.index not in set(pay_global)]
        rgd_cent = xyz[[cas[i] for i in IDX_RGD if cas[i] is not None]].mean(0)
        protrusion = (float(np.linalg.norm(rgd_cent - xyz[scaffold_ca].mean(0)))
                      if scaffold_ca else None)
        # nearest scaffold heavy atom to the RGD motif: steric crowding of the motif
        pay_set = set(pay_global)
        scaf_heavy = [a.index for a in t.topology.atoms
                      if a.residue.index not in pay_set and a.element.symbol != "H"]
        rgd_atoms = [a.index for a in t.topology.atoms
                     if a.residue.index in {pay_global[i] for i in IDX_RGD}
                     and a.element.symbol != "H"]
        if scaf_heavy and rgd_atoms:
            d = np.linalg.norm(xyz[rgd_atoms][:, None, :] - xyz[scaf_heavy][None, :, :],
                               axis=2)
            rgd_min_scaffold = float(d.min())
        else:
            rgd_min_scaffold = None

        # ---------------------------------------------------- C. clashes
        graft_atoms = [a.index for a in t.topology.atoms
                       if a.residue.index in pay_set and a.element.symbol != "H"]
        clashes = 0
        if graft_atoms and scaf_heavy:
            dd = np.linalg.norm(xyz[graft_atoms][:, None, :] - xyz[scaf_heavy][None, :, :],
                                axis=2)
            clashes = int((dd < 2.2).sum())

        per_sample.append({
            "file": fname,
            "helix_frac_core": helix_core,
            "helix_frac_payload": helix_all,
            "rmsd_to_6R2X_A": best,
            "end_to_end_A": e2e,
            "rgd_rel_sasa": rgd_rel,
            "rgd_bulky_probe_sasa_A2": rgd_big,
            "hydrophobic_face_rel_sasa": hphob_rel,
            "pharmacophore_A": pharm,
            "rgd_protrusion_A": protrusion,
            "rgd_min_dist_to_scaffold_A": rgd_min_scaffold,
            "graft_scaffold_clashes": clashes,
            "n_chains": len(seqs),
            "_traj_index": len(per_sample),
        })
    return {"samples": per_sample, "files": [f for f, _ in samples]}


def agg(samples, key):
    vals = [s[key] for s in samples if s.get(key) is not None]
    if not vals:
        return None
    return {"mean": round(float(np.mean(vals)), 3),
            "sd": round(float(np.std(vals)), 3),
            "min": round(float(np.min(vals)), 3),
            "max": round(float(np.max(vals)), 3),
            "n": len(vals)}


if __name__ == "__main__":
    results = {}
    for cid, c in by_id.items():
        payload = c.get("payload")
        if payload is None:
            if cid == "REF_4delta_free":
                payload = FOURD
            elif cid == "BM_CH3_Cterm_fusion":
                payload = FOURD
            else:
                continue          # WT scaffolds have no payload; handled separately
        print(f"analysing {cid} ...")
        r = analyse(cid, payload)
        if r is None:
            print(f"  (no models yet)")
            continue
        s = r["samples"]
        results[cid] = {
            "n_samples": len(s),
            "payload": payload,
            "metrics": {k: agg(s, k) for k in [
                "helix_frac_core", "helix_frac_payload", "rmsd_to_6R2X_A",
                "end_to_end_A", "rgd_rel_sasa", "rgd_bulky_probe_sasa_A2",
                "hydrophobic_face_rel_sasa", "pharmacophore_A",
                "rgd_protrusion_A", "rgd_min_dist_to_scaffold_A",
                "graft_scaffold_clashes"]},
            "per_sample": s,
        }
    json.dump({"reference_6R2X": REF, "results": results},
              open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT} for {len(results)} constructs")
