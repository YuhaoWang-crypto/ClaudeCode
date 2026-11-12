#!/usr/bin/env python3
"""Task 3c: did the graft damage the antibody scaffold?

Three things decide whether a grafted construct is still an antibody:
  1. the host domain must keep its fold          -> scaffold RMSD vs the WT model
  2. CH3 must still homodimerise                 -> native interface contacts kept
     (an Fc that cannot dimerise is not an Fc), and for Fv the VH/VL pairing
  3. the model must actually be confident where  -> per-residue pLDDT read from
     it matters                                     the CIF B_iso column

Scaffold RMSD uses only the segments that are byte-identical to wild type on
either side of the graft, so no alignment ambiguity is introduced by the insert.
"""
import glob
import json
import os
import sys

import mdtraj as md
import numpy as np

MODELS = sys.argv[1] if len(sys.argv) > 1 else "models"
OUT = "results/09_scaffold_integrity.json"

cand = json.load(open("results/03_candidates.json"))
by_id = {c["id"]: c for c in cand["candidates"] + cand["controls"]}
CH3 = cand["scaffold"]["CH3_341_447"]
VH = cand["scaffold"]["VH"]
VL = cand["scaffold"]["VL"]
CDRH3 = "WGGDGFYAMDY"
CDR3_START = VH.index(CDRH3)
SITES = {"CH3-AB": (355, 362), "CH3-CD": (383, 391), "CH3-EF": (413, 422)}
CH3_LO = 341
sites_json = json.load(open("results/02_graft_sites.json"))
IFACE = set(sites_json["CH3"]["dimer_interface_residues_eu"])


def flanks(cid):
    """WT-identical sequence segments on either side of the graft."""
    c = by_id[cid]
    site = c.get("site")
    if site == "CDRH3":
        return VH[:CDR3_START], VH[CDR3_START + len(CDRH3):], "Fv"
    if site in SITES:
        lo, hi = SITES[site]
        return CH3[: lo - CH3_LO], CH3[hi - CH3_LO + 1 :], "CH3"
    if cid == "WT_Fv":
        return VH[:CDR3_START], VH[CDR3_START + len(CDRH3):], "Fv"
    if cid in ("WT_CH3", "BM_CH3_Cterm_fusion"):
        lo, hi = SITES["CH3-EF"]
        return CH3[: lo - CH3_LO], CH3[hi - CH3_LO + 1 :], "CH3"
    return None, None, None


def plddt_by_residue(path):
    """(auth_asym_id, auth_seq_id) -> pLDDT, from the CIF B_iso column."""
    cols, vals, in_loop = [], {}, False
    for line in open(path):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1])
            in_loop = True
            continue
        if in_loop and line.startswith(("ATOM", "HETATM")):
            f = line.split()
            if len(f) != len(cols):
                continue
            rec = dict(zip(cols, f))
            try:
                key = (rec["auth_asym_id"], int(rec["auth_seq_id"]))
                vals.setdefault(key, []).append(float(rec["B_iso_or_equiv"]))
            except (KeyError, ValueError):
                pass
        elif in_loop and line.startswith("#"):
            in_loop = False
    return {k: float(np.mean(v)) for k, v in vals.items()}


def chain_residues(t):
    return {ch.index: [r for r in ch.residues if r.is_protein] for ch in t.topology.chains}


def load(cid):
    out = []
    for f in sorted(glob.glob(os.path.join(MODELS, cid, "*.cif"))):
        t = md.load(f)
        t = t.atom_slice(t.topology.select("protein"))
        if t.n_atoms:
            out.append((f, t))
    return out


def seg_ca(rs, seq, segment, at_start):
    """CA indices for a WT-identical flank inside a predicted chain."""
    s = "".join(r.code or "X" for r in rs)
    pos = s.find(segment) if at_start else s.rfind(segment)
    if pos < 0:
        return None
    out = []
    for r in rs[pos:pos + len(segment)]:
        a = next((x for x in r.atoms if x.name == "CA"), None)
        if a is None:
            return None
        out.append(a.index)
    return out


def kabsch_rmsd(m, t):
    mc, tc = m.mean(0), t.mean(0)
    m, t = m - mc, t - tc
    u, s, vt = np.linalg.svd(m.T @ t)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1, 1, d]) @ vt
    return float(np.sqrt(((m @ r - t) ** 2).sum(1).mean()))


# ------------------------------------------------------------------ WT references
# The flanks differ per graft site, so the WT reference coordinates have to be
# extracted with THAT candidate's own flanks - one global reference would only
# match the site it was built from.
wt_models = {}
for wt in ("WT_CH3", "WT_Fv"):
    ms = load(wt)
    if not ms:
        print(f"! {wt} missing")
        continue
    kind = "CH3" if wt == "WT_CH3" else "Fv"
    wt_models[kind] = ms[0][1]                    # sample_0 as the reference frame
    print(f"{wt}: reference frame loaded ({kind})")


def wt_flank_xyz(kind, pre, post):
    """CA coordinates of these two flanks inside the WT model of that module."""
    t = wt_models.get(kind)
    if t is None:
        return None
    rs = chain_residues(t)[min(chain_residues(t))]
    a = seg_ca(rs, None, pre, True)
    b = seg_ca(rs, None, post, False)
    if not a or not b:
        return None
    return t.xyz[0][a + b] * 10.0


def dimer_contacts(t, native_only=True):
    """Heavy-atom contacts (<4.5 A) between the two host chains."""
    top = t.topology
    chs = chain_residues(t)
    if len(chs) < 2:
        return None
    ids = sorted(chs)[:2]
    xyz = t.xyz[0] * 10.0
    sel = []
    for ci in ids:
        sel.append([a.index for a in top.chain(ci).atoms if a.element.symbol != "H"])
    d = np.linalg.norm(xyz[sel[0]][:, None, :] - xyz[sel[1]][None, :, :], axis=2)
    return int((d < 4.5).sum())


if __name__ == "__main__":
    out = {}
    for cid, c in by_id.items():
        ms = load(cid)
        if not ms:
            continue
        pre, post, kind = flanks(cid)
        payload = c.get("payload") or (
            "FETLRGDLRILSILRHQNLLKEL" if cid in ("REF_4delta_free", "BM_CH3_Cterm_fusion") else None)
        rows = []
        for f, t in ms:
            chs = chain_residues(t)
            host = min(chs)
            rs = chs[host]
            pl = plddt_by_residue(f)
            # pLDDT over the graft and over the scaffold
            s = "".join(r.code or "X" for r in rs)
            g_plddt = s_plddt = rgd_plddt = None
            if payload:
                p0 = s.find(payload)
                if p0 >= 0:
                    keys = []
                    for r in rs[p0:p0 + len(payload)]:
                        a = next(iter(r.atoms))
                        keys.append(pl.get((a.residue.chain.chain_id
                                            if hasattr(a.residue.chain, "chain_id") else "A",
                                            r.resSeq)))
                    got = [x for x in keys if x is not None]
                    if got:
                        g_plddt = float(np.mean(got))
                        rgd_plddt = float(np.mean(got[4:7])) if len(got) > 6 else None
            allv = list(pl.values())
            s_plddt = float(np.mean(allv)) if allv else None

            rmsd = None
            wt_xyz = wt_flank_xyz(kind, pre, post) if (pre and post) else None
            if wt_xyz is not None:
                a = seg_ca(rs, None, pre, True)
                b = seg_ca(rs, None, post, False)
                if a and b and len(a) + len(b) == len(wt_xyz):
                    rmsd = kabsch_rmsd(t.xyz[0][a + b] * 10.0, wt_xyz)

            rows.append({
                "file": os.path.basename(f),
                "graft_plddt": None if g_plddt is None else round(g_plddt, 1),
                "rgd_plddt": None if rgd_plddt is None else round(rgd_plddt, 1),
                "model_plddt": None if s_plddt is None else round(s_plddt, 1),
                "scaffold_rmsd_vs_wt_A": None if rmsd is None else round(rmsd, 2),
                "interchain_contacts": dimer_contacts(t),
            })
        def ag(k):
            v = [r[k] for r in rows if r.get(k) is not None]
            return None if not v else {"mean": round(float(np.mean(v)), 2),
                                       "sd": round(float(np.std(v)), 2),
                                       "min": round(float(np.min(v)), 2),
                                       "max": round(float(np.max(v)), 2)}
        out[cid] = {"n_samples": len(rows), "module": kind,
                    "metrics": {k: ag(k) for k in
                                ("graft_plddt", "rgd_plddt", "model_plddt",
                                 "scaffold_rmsd_vs_wt_A", "interchain_contacts")},
                    "per_sample": rows}
        m = out[cid]["metrics"]
        print(f"{cid:<22} graft-pLDDT {str(m['graft_plddt'] and m['graft_plddt']['mean']):>6}  "
              f"scaf-RMSD {str(m['scaffold_rmsd_vs_wt_A'] and m['scaffold_rmsd_vs_wt_A']['mean']):>6}  "
              f"interchain {str(m['interchain_contacts'] and m['interchain_contacts']['mean']):>8}")

    json.dump(out, open(OUT, "w"), indent=2)
    print(f"\nwrote {OUT}")
