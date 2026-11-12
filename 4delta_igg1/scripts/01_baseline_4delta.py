#!/usr/bin/env python3
"""Task 1b: baseline conformational assessment of the 4-Delta reference peptide.

Reference structure: PDB 6R2X, solution-NMR ensemble (15 models) of human
chromogranin A F39-D63  = FETLRGDERILSILRHQNLLKELQD.

4-Delta = FETLRGDLRILSILRHQNLLKEL  = CgA(39-61) carrying E46L.
So 6R2X is the closest experimental template: identical except E46->L and the
C-terminal Q62-D63 truncation.

Outputs per-residue helix fraction, RMSF, RGD exposure and the graft-critical
end-to-end span that any host loop must bridge.
"""
import json
import sys

import mdtraj as md
import numpy as np

PDB = sys.argv[1] if len(sys.argv) > 1 else "pdb/6R2X.pdb"
OUT = sys.argv[2] if len(sys.argv) > 2 else "results/01_baseline_4delta.json"

NATIVE = "FETLRGDERILSILRHQNLLKELQD"   # 6R2X construct, CgA 39-63
FOURD = "FETLRGDLRILSILRHQNLLKEL"      # 4-Delta, CgA 39-61 + E46L
CGA_FIRST = 39                          # residue 1 of the construct is CgA F39

t = md.load(PDB)
top = t.topology
seq = "".join(r.code for r in top.residues if r.is_protein)
assert seq == NATIVE, f"unexpected construct sequence: {seq}"
n_res = len(seq)
print(f"loaded {t.n_frames} NMR models, {n_res} residues: {seq}")

# ---------------------------------------------------------------- secondary structure
dssp = md.compute_dssp(t, simplified=False)           # (frames, residues)
helix_codes = {"H", "G", "I"}
helix_frac = np.array([
    np.mean([code in helix_codes for code in dssp[:, i]]) for i in range(n_res)
])
alpha_frac = np.array([np.mean(dssp[:, i] == "H") for i in range(n_res)])

# ---------------------------------------------------------------- flexibility (RMSF on CA)
ca = top.select("name CA")
t_ca = t.atom_slice(ca)
t_ca.superpose(t_ca, 0)
mean_xyz = t_ca.xyz.mean(axis=0)
rmsf = np.sqrt(((t_ca.xyz - mean_xyz) ** 2).sum(axis=2).mean(axis=0)) * 10.0  # A

# superposition restricted to the structured helix core (CgA E46-N56 per the
# published ensemble alignment) gives a fairer picture of RGD mobility
core_idx = [i for i in range(n_res) if 46 <= CGA_FIRST + i <= 56]
t_core = t.atom_slice(ca)
t_core.superpose(t_core, 0, atom_indices=core_idx)
mean_core = t_core.xyz.mean(axis=0)
rmsf_core = np.sqrt(((t_core.xyz - mean_core) ** 2).sum(axis=2).mean(axis=0)) * 10.0

# ---------------------------------------------------------------- SASA
sasa_res = md.shrake_rupley(t, mode="residue") * 100.0   # nm^2 -> A^2
sasa_mean = sasa_res.mean(axis=0)
sasa_sd = sasa_res.std(axis=0)

# Tripeptide reference values for relative exposure of the free peptide
MAX_ACC = {"A": 106, "R": 248, "N": 157, "D": 163, "C": 135, "Q": 198, "E": 194,
           "G": 84, "H": 184, "I": 169, "L": 164, "K": 205, "M": 188, "F": 197,
           "P": 136, "S": 130, "T": 142, "W": 227, "Y": 222, "V": 142}
rel_sasa = np.array([sasa_mean[i] / MAX_ACC[seq[i]] for i in range(n_res)])

# ---------------------------------------------------------------- graft geometry
xyz_ca = t.atom_slice(ca).xyz * 10.0            # (frames, res, 3) in A
# span the host loop must bridge if 4-Delta (res 1..23 of the construct) is
# inserted as a contiguous cassette
i_n, i_c = 0, 22                                 # F39 and L61 -> 4-Delta termini
end_to_end = np.linalg.norm(xyz_ca[:, i_c] - xyz_ca[:, i_n], axis=1)

# RGD = CgA R43-G44-D45 -> construct indices 4,5,6
rgd_idx = [43 - CGA_FIRST, 44 - CGA_FIRST, 45 - CGA_FIRST]
# pharmacophore distance: Arg guanidinium CZ to Asp carboxyl CG
def atom_xyz(resid, name):
    sel = top.select(f"resid {resid} and name {name}")
    if len(sel) == 0:
        return None
    return t.xyz[:, sel[0], :] * 10.0

arg_cz = atom_xyz(rgd_idx[0], "CZ")
asp_cg = atom_xyz(rgd_idx[2], "CG")
rgd_pharmacophore = np.linalg.norm(arg_cz - asp_cg, axis=1)

# how far the RGD tip projects away from the helix core centroid
helix_centroid = xyz_ca[:, core_idx, :].mean(axis=1)
rgd_tip = xyz_ca[:, rgd_idx, :].mean(axis=1)
rgd_projection = np.linalg.norm(rgd_tip - helix_centroid, axis=1)

# ---------------------------------------------------------------- amphipathicity
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

def hydrophobic_moment(segment):
    """Eisenberg hydrophobic moment per residue, 100 deg/residue."""
    vx = vy = 0.0
    for k, aa in enumerate(segment):
        ang = np.deg2rad(100.0 * k)
        vx += KD[aa] * np.cos(ang)
        vy += KD[aa] * np.sin(ang)
    return float(np.hypot(vx, vy) / len(segment))

helix_seg_4d = FOURD[7:21]      # L46..K59 equivalent in 4-Delta
mu_h = hydrophobic_moment(helix_seg_4d)
mu_h_native = hydrophobic_moment(NATIVE[7:21])

report = {
    "reference_structure": {
        "pdb": "6R2X", "method": "solution NMR", "models": int(t.n_frames),
        "construct": NATIVE, "construct_range": "CgA F39-D63",
    },
    "four_delta": {
        "sequence": FOURD, "length": len(FOURD),
        "relationship_to_template": "CgA(39-61) with E46L; template differs by E46L + Q62/D63 truncation",
        "rgd_position_in_peptide": "R5-G6-D7",
    },
    "per_residue": [
        {
            "i": i + 1, "cga": CGA_FIRST + i, "aa": seq[i],
            "in_4delta": (i < len(FOURD)),
            "helix_frac": round(float(helix_frac[i]), 3),
            "alpha_frac": round(float(alpha_frac[i]), 3),
            "rmsf_A": round(float(rmsf[i]), 2),
            "rmsf_corealigned_A": round(float(rmsf_core[i]), 2),
            "sasa_A2": round(float(sasa_mean[i]), 1),
            "sasa_sd_A2": round(float(sasa_sd[i]), 1),
            "rel_sasa": round(float(rel_sasa[i]), 3),
        }
        for i in range(n_res)
    ],
    "graft_geometry": {
        "end_to_end_CA_A_mean": round(float(end_to_end.mean()), 2),
        "end_to_end_CA_A_sd": round(float(end_to_end.std()), 2),
        "end_to_end_CA_A_min": round(float(end_to_end.min()), 2),
        "end_to_end_CA_A_max": round(float(end_to_end.max()), 2),
        "rgd_pharmacophore_ArgCZ_AspCG_A_mean": round(float(rgd_pharmacophore.mean()), 2),
        "rgd_pharmacophore_A_sd": round(float(rgd_pharmacophore.std()), 2),
        "rgd_tip_projection_from_helix_core_A_mean": round(float(rgd_projection.mean()), 2),
        "rgd_tip_projection_A_sd": round(float(rgd_projection.std()), 2),
    },
    "amphipathicity": {
        "helix_segment_4delta": helix_seg_4d,
        "hydrophobic_moment_per_residue_4delta": round(mu_h, 3),
        "hydrophobic_moment_per_residue_native": round(mu_h_native, 3),
    },
    "derived": {
        "helical_core_cga": "E46-K59" ,
        "helix_frac_mean_core": round(float(helix_frac[core_idx].mean()), 3),
        "rgd_mean_rel_sasa": round(float(rel_sasa[rgd_idx].mean()), 3),
        "rgd_mean_rmsf_corealigned_A": round(float(rmsf_core[rgd_idx].mean()), 2),
    },
}

with open(OUT, "w") as fh:
    json.dump(report, fh, indent=2)

print("\n== per-residue (helix fraction / core-aligned RMSF / rel SASA) ==")
for r in report["per_residue"]:
    bar = "#" * int(round(r["helix_frac"] * 20))
    mark = "  <- RGD" if r["cga"] in (43, 44, 45) else ("" if r["in_4delta"] else "  (not in 4D)")
    print(f"{r['i']:>3} CgA{r['cga']:<4}{r['aa']}  H={r['helix_frac']:.2f} {bar:<20} "
          f"RMSF={r['rmsf_corealigned_A']:>5.2f}  relSASA={r['rel_sasa']:.2f}{mark}")

g = report["graft_geometry"]
print(f"\nend-to-end CA span (4D termini) = {g['end_to_end_CA_A_mean']} +/- {g['end_to_end_CA_A_sd']} A "
      f"(range {g['end_to_end_CA_A_min']}-{g['end_to_end_CA_A_max']})")
print(f"RGD pharmacophore Arg-CZ..Asp-CG = {g['rgd_pharmacophore_ArgCZ_AspCG_A_mean']} +/- {g['rgd_pharmacophore_A_sd']} A")
print(f"RGD tip projection from helix core = {g['rgd_tip_projection_from_helix_core_A_mean']} +/- {g['rgd_tip_projection_A_sd']} A")
print(f"hydrophobic moment/res (4D helix) = {report['amphipathicity']['hydrophobic_moment_per_residue_4delta']}")
print(f"wrote {OUT}")
