"""Interface analysis for a designed binder-target complex.

Computes, for a two-chain (or multi-chain) model: interface residues on both sides,
buried surface area, polar contacts (H-bond / salt-bridge geometry by distance),
shape complementarity proxy, and overlap of the binder footprint with a reference
epitope definition.

Usage:
    python analyze_interface.py complex.cif --target-chain A --binder-chain X \
        --epitope-json il6_target.json
"""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley

SR = ShrakeRupley()

ACIDIC = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
BASIC = {"ARG": ["NH1", "NH2", "NE"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"]}
DONOR_ACCEPTOR = set("NO")

KD = {  # Kyte-Doolittle
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure("m", str(path))[0]


def chain_atoms(chain):
    return [a for r in chain for a in r if r.id[0] == " "]


def sasa_of(entities):
    """Run Shrake-Rupley over a list of chains treated as one body; returns dict."""
    import copy
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model

    s = Structure("t")
    m = Model(0)
    s.add(m)
    for ch in entities:
        m.add(copy.deepcopy(ch))
    SR.compute(m, level="R")
    return {(ch.id, r.id[1]): r.sasa for ch in m for r in ch if r.id[0] == " "}


def analyze(model, target_chain, binder_chain, contact_cut=5.0):
    tgt, bnd = model[target_chain], model[binder_chain]
    ta, ba = chain_atoms(tgt), chain_atoms(bnd)
    tc = np.array([a.coord for a in ta])
    bc = np.array([a.coord for a in ba])

    d = np.linalg.norm(tc[:, None, :] - bc[None, :, :], axis=-1)
    close = d <= contact_cut

    tgt_iface, bnd_iface, pairs = {}, {}, []
    for i, j in zip(*np.where(close)):
        at, ab = ta[i], ba[j]
        rt, rb = at.get_parent(), ab.get_parent()
        tgt_iface.setdefault(rt.id[1], rt.get_resname())
        bnd_iface.setdefault(rb.id[1], rb.get_resname())
        pairs.append((rt, at, rb, ab, float(d[i, j])))

    # polar contacts: N/O to N/O within 3.5 A
    hbonds, salt = [], []
    for rt, at, rb, ab, dist in pairs:
        if dist <= 3.5 and at.element in DONOR_ACCEPTOR and ab.element in DONOR_ACCEPTOR:
            hbonds.append((f"{rt.get_resname()}{rt.id[1]}:{at.get_id()}",
                           f"{rb.get_resname()}{rb.id[1]}:{ab.get_id()}", round(dist, 2)))
        if dist <= 4.0:
            for a_res, a_at, b_res, b_at in ((rt, at, rb, ab), (rb, ab, rt, at)):
                if a_res.get_resname() in ACIDIC and a_at.get_id() in ACIDIC[a_res.get_resname()] \
                   and b_res.get_resname() in BASIC and b_at.get_id() in BASIC[b_res.get_resname()]:
                    salt.append((f"{a_res.get_resname()}{a_res.id[1]}",
                                 f"{b_res.get_resname()}{b_res.id[1]}", round(dist, 2)))

    # buried surface area
    bound = sasa_of([tgt, bnd])
    free_t = sasa_of([tgt])
    free_b = sasa_of([bnd])
    bsa_t = sum(max(0.0, free_t[k] - bound[k]) for k in free_t)
    bsa_b = sum(max(0.0, free_b[k] - bound[k]) for k in free_b)

    seq_b = "".join(THREE2ONE.get(r.get_resname(), "X") for r in bnd if r.id[0] == " ")
    hydroph = np.mean([KD.get(c, 0) for c in seq_b])

    return {
        "target_interface_residues": {int(k): v for k, v in sorted(tgt_iface.items())},
        "binder_interface_residues": {int(k): v for k, v in sorted(bnd_iface.items())},
        "n_target_interface": len(tgt_iface),
        "n_binder_interface": len(bnd_iface),
        "n_atom_pairs_within_5A": int(close.sum()),
        "n_atom_pairs_within_4A": int((d <= 4.0).sum()),
        "min_heavy_atom_distance": round(float(d.min()), 2),
        "bsa_target_A2": round(bsa_t, 1),
        "bsa_binder_A2": round(bsa_b, 1),
        "bsa_total_A2": round(bsa_t + bsa_b, 1),
        "hbond_like_contacts": sorted(set(hbonds), key=lambda x: x[2]),
        "salt_bridges": sorted({(a, b, c) for a, b, c in salt}, key=lambda x: x[2]),
        "binder_sequence": seq_b,
        "binder_length": len(seq_b),
        "binder_mean_kd_hydropathy": round(float(hydroph), 2),
    }


def epitope_overlap(result, epitope_resnums, label):
    hit = sorted(set(result["target_interface_residues"]) & set(epitope_resnums))
    frac = len(hit) / max(1, len(result["target_interface_residues"]))
    return {
        f"{label}_residues_contacted": hit,
        f"{label}_recall": round(len(hit) / len(epitope_resnums), 2),
        f"{label}_footprint_fraction": round(frac, 2),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("complex", type=Path)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain", default="X")
    ap.add_argument("--offset", type=int, default=0,
                    help="add to target residue numbers to reach 1P9M numbering")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    model = load(args.complex)
    res = analyze(model, args.target_chain, args.binder_chain)

    # IL-6 reference epitopes, 1P9M chain-B numbering (mature IL-6)
    site1 = [30, 33, 54, 61, 66, 69, 73, 74, 75, 78, 172, 175, 178, 179, 180, 182, 183]
    site2 = [19, 24, 27, 28, 30, 31, 34, 110, 111, 113, 114, 117, 118, 121, 124, 125, 128]
    shifted = {k + args.offset: v for k, v in res["target_interface_residues"].items()}
    tmp = dict(res, target_interface_residues=shifted)
    res.update(epitope_overlap(tmp, site1, "site_I_IL6Ra"))
    res.update(epitope_overlap(tmp, site2, "site_II_gp130"))

    print(json.dumps(res, indent=2))
    if args.out:
        args.out.write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
