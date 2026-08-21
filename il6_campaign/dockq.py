"""Self-consistency DockQ between a design model and an independent prediction.

DockQ (Basu & Wallner 2016) combines three terms:
    fnat   fraction of the reference interface residue contacts recovered
    LRMSD  binder RMSD after superposing on the target
    iRMSD  backbone RMSD of interface residues after superposing on them

    DockQ = (fnat + 1/(1+(LRMSD/8.5)^2) + 1/(1+(iRMSD/1.5)^2)) / 3

Used here as a *gate*, not a ranking term: it answers "did the independent co-fold
reproduce the pose the design was built for?" — the paper keeps sc-DockQ for exactly
this purpose and reports that adding it to the ranking score leaves discrimination
unchanged.

Residue numbering differs between the two files, so both sides are mapped onto a
common frame with --ref-offset / --mod-offset (added to each file's target numbering).

Usage:
    python dockq.py design.cif cofold.cif --ref-target B --ref-binder X \
        --mod-target A --mod-binder X --ref-offset -2 --mod-offset 1
"""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

BACKBONE = ("N", "CA", "C", "O")


def load(path: Path):
    p = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return p.get_structure("s", str(path))[0]


def residues(chain):
    return [r for r in chain if r.id[0] == " "]


def heavy(res):
    return [a for a in res if a.element != "H"]


def contacts(model, tchain, bchain, toff, boff=0, cut=5.0):
    """Set of (target_resnum, binder_resnum) pairs with any heavy-atom pair < cut."""
    tres, bres = residues(model[tchain]), residues(model[bchain])
    tc = [np.array([a.coord for a in heavy(r)]) for r in tres]
    bc = [np.array([a.coord for a in heavy(r)]) for r in bres]
    out = set()
    for i, ta in enumerate(tc):
        if not len(ta):
            continue
        for j, ba in enumerate(bc):
            if not len(ba):
                continue
            if np.min(np.linalg.norm(ta[:, None] - ba[None], axis=-1)) < cut:
                out.add((tres[i].id[1] + toff, bres[j].id[1] + boff))
    return out


def common_atoms(ref_chain, mod_chain, ref_off, mod_off, resnums=None, names=BACKBONE):
    """Matched coordinate arrays for residues present in both chains."""
    rmap = {r.id[1] + ref_off: r for r in residues(ref_chain)}
    mmap = {r.id[1] + mod_off: r for r in residues(mod_chain)}
    keys = sorted(set(rmap) & set(mmap)) if resnums is None else \
        sorted(set(rmap) & set(mmap) & set(resnums))
    a, b = [], []
    for k in keys:
        for n in names:
            if n in rmap[k] and n in mmap[k]:
                a.append(rmap[k][n].coord)
                b.append(mmap[k][n].coord)
    return np.array(a), np.array(b)


def kabsch(P, Q):
    """Rotation+translation mapping Q onto P; returns (R, t, rmsd)."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (Q - qc).T @ (P - pc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = pc - R @ qc
    rmsd = float(np.sqrt((((R @ Q.T).T + t - P) ** 2).sum(1).mean()))
    return R, t, rmsd


def dockq(ref, mod, ref_t, ref_b, mod_t, mod_b, ref_off, mod_off):
    ref_contacts = contacts(ref, ref_t, ref_b, ref_off)
    mod_contacts = contacts(mod, mod_t, mod_b, mod_off)
    if not ref_contacts:
        return {"error": "no contacts in reference"}
    fnat = len(ref_contacts & mod_contacts) / len(ref_contacts)
    fnonnat = (1 - len(ref_contacts & mod_contacts) / len(mod_contacts)) if mod_contacts else 1.0

    # LRMSD: superpose on target, measure binder
    Pt, Qt = common_atoms(ref[ref_t], mod[mod_t], ref_off, mod_off)
    R, t, _ = kabsch(Pt, Qt)
    Pb, Qb = common_atoms(ref[ref_b], mod[mod_b], 0, 0)
    lrms = float(np.sqrt((((R @ Qb.T).T + t - Pb) ** 2).sum(1).mean())) if len(Pb) else float("nan")

    # iRMSD: superpose on interface residues of both chains
    itgt = {c[0] for c in ref_contacts}
    ibnd = {c[1] for c in ref_contacts}
    Pi1, Qi1 = common_atoms(ref[ref_t], mod[mod_t], ref_off, mod_off, itgt)
    Pi2, Qi2 = common_atoms(ref[ref_b], mod[mod_b], 0, 0, ibnd)
    Pi, Qi = np.vstack([Pi1, Pi2]), np.vstack([Qi1, Qi2])
    _, _, irms = kabsch(Pi, Qi)

    score = (fnat + 1 / (1 + (lrms / 8.5) ** 2) + 1 / (1 + (irms / 1.5) ** 2)) / 3
    quality = ("high" if score >= 0.80 else "medium" if score >= 0.49
               else "acceptable" if score >= 0.23 else "incorrect")
    return {"dockq": round(score, 3), "fnat": round(fnat, 3), "fnonnat": round(fnonnat, 3),
            "lrmsd_A": round(lrms, 2), "irmsd_A": round(irms, 2),
            "n_ref_contacts": len(ref_contacts), "n_model_contacts": len(mod_contacts),
            "quality": quality}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reference", type=Path, help="design-time model")
    ap.add_argument("model", type=Path, help="independent prediction")
    ap.add_argument("--ref-target", default="B")
    ap.add_argument("--ref-binder", default="X")
    ap.add_argument("--mod-target", default="A")
    ap.add_argument("--mod-binder", default="X")
    ap.add_argument("--ref-offset", type=int, default=0)
    ap.add_argument("--mod-offset", type=int, default=0)
    args = ap.parse_args()

    res = dockq(load(args.reference), load(args.model), args.ref_target, args.ref_binder,
                args.mod_target, args.mod_binder, args.ref_offset, args.mod_offset)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
