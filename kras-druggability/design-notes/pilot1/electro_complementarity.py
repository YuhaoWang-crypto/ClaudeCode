#!/usr/bin/env python3
"""
界面静电互补分析器（Pilot 1：复现"为什么选 CypA"）。
量化一个蛋白-蛋白界面的电荷互补性，用于比较不同 presenter（CypA vs FKBP12）
对 KRAS 效应面的互补程度。纯 numpy，无第三方依赖。

度量：
  - 两侧界面残基 patch 的净形式电荷（Asp/Glu=-1, Lys/Arg=+1, His=+0.5）
  - 跨界面盐桥数（酸性O ↔ 碱性N < 4.0 Å）
  - 跨界面同电荷冲突数（O↔O 或 N↔N < 4.0 Å）
  - 静电互补分 EC = Σ_cross (−q_i·q_j / r_ij)，跨界面带电原子对 < 8 Å；正=互补
用法: python3 electro_complementarity.py <cif> --presenter <chain> --kras <chain> [--label L]
"""
import sys, argparse, math
import numpy as np

# 带电侧链原子的部分电荷（概念性，用于相对比较，非精确力场）
POS = {"LYS": {"NZ": 1.0},
       "ARG": {"NH1": 1/3, "NH2": 1/3, "NE": 1/3},
       "HIS": {"ND1": 0.05, "NE2": 0.05}}
NEG = {"ASP": {"OD1": -0.5, "OD2": -0.5},
       "GLU": {"OE1": -0.5, "OE2": -0.5}}
RES_FORMAL = {"ASP": -1, "GLU": -1, "LYS": +1, "ARG": +1, "HIS": +0.5}

def parse_cif(path):
    atoms = []; lines = open(path).read().splitlines(); i = 0
    while i < len(lines):
        if lines[i].startswith("_atom_site."):
            cols = []
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip().split(".", 1)[1]); i += 1
            idx = {c: k for k, c in enumerate(cols)}
            def col(*names):
                for n in names:
                    if n in idx: return idx[n]
                return None
            ci = {"grp": col("group_PDB"), "atom": col("label_atom_id"),
                  "comp": col("label_comp_id"), "ch": col("auth_asym_id", "label_asym_id"),
                  "seq": col("auth_seq_id", "label_seq_id"),
                  "x": col("Cartn_x"), "y": col("Cartn_y"), "z": col("Cartn_z"),
                  "model": col("pdbx_PDB_model_num")}
            while i < len(lines):
                s = lines[i].strip()
                if s == "" or s[0] == "#" or s.startswith(("loop_", "_", "data_")): break
                p = s.split()
                if len(p) < len(cols): i += 1; continue
                if ci["model"] is not None and p[ci["model"]] not in (".", "1"): i += 1; continue
                try:
                    atoms.append({"grp": p[ci["grp"]], "atom": p[ci["atom"]].strip('"'),
                                  "comp": p[ci["comp"]], "ch": p[ci["ch"]], "seq": p[ci["seq"]],
                                  "xyz": np.array([float(p[ci["x"]]), float(p[ci["y"]]), float(p[ci["z"]])])})
                except (ValueError, IndexError): pass
                i += 1
            continue
        i += 1
    return atoms

def prot_atoms(atoms, ch):
    return [a for a in atoms if a["ch"] == ch and a["grp"] == "ATOM"]

def interface_residues(A, B, cut=4.5):
    XA = np.array([a["xyz"] for a in A]); XB = np.array([a["xyz"] for a in B])
    d2 = ((XA[:, None, :] - XB[None, :, :]) ** 2).sum(-1)
    ia, ib = np.where(d2 <= cut * cut)
    resA = {(A[k]["comp"], A[k]["seq"]) for k in ia}
    resB = {(B[k]["comp"], B[k]["seq"]) for k in ib}
    return resA, resB

def patch_net_charge(resset):
    return sum(RES_FORMAL.get(c, 0) for c, s in resset)

def charged_atoms(prot, resset):
    """返回带电原子 (xyz, q)，仅限界面残基。"""
    want = {(c, s) for c, s in resset}
    out = []
    for a in prot:
        key = (a["comp"], a["seq"])
        if key not in want: continue
        for tbl in (POS, NEG):
            if a["comp"] in tbl and a["atom"] in tbl[a["comp"]]:
                out.append((a["xyz"], tbl[a["comp"]][a["atom"]]))
    return out

def analyze(cif, pres_ch, kras_ch, label):
    atoms = parse_cif(cif)
    P = prot_atoms(atoms, pres_ch); K = prot_atoms(atoms, kras_ch)
    print(f"\n########## {label} ##########")
    print(f"presenter链 {pres_ch}: {len(P)} 原子 | KRAS链 {kras_ch}: {len(K)} 原子")
    if not P or not K:
        print("  !! 缺链，跳过"); return None
    rP, rK = interface_residues(P, K)
    qP = patch_net_charge(rP); qK = patch_net_charge(rK)
    print(f"界面残基: presenter {len(rP)} (净电荷 {qP:+.0f}) | KRAS {len(rK)} (净电荷 {qK:+.0f})")

    cP = charged_atoms(P, rP); cK = charged_atoms(K, rK)
    # 盐桥 / 冲突（跨界面，<4.0 Å，用带电侧链极性原子）
    saltbridge = clash = 0
    EC = 0.0
    for xp, qp in cP:
        for xk, qk in cK:
            r = math.sqrt(((xp - xk) ** 2).sum())
            if r < 8.0:
                EC += (-qp * qk) / max(r, 1.0)
            if r < 4.0:
                if qp * qk < 0: saltbridge += 1
                elif qp * qk > 0: clash += 1
    npair = max(len(rP) + len(rK), 1)
    print(f"跨界面盐桥(<4Å): {saltbridge} | 同电荷冲突(<4Å): {clash}")
    print(f"静电互补分 EC = {EC:+.2f}  (每界面残基 {EC/npair:+.3f}; 正=互补)")
    print(f"电荷互补判读: KRAS面 {'酸性' if qK<0 else '碱性' if qK>0 else '中性'} ↔ "
          f"presenter面 {'碱性' if qP>0 else '酸性' if qP<0 else '中性'} "
          f"→ {'互补 ✓' if qP*qK<0 else '不互补'}")
    return {"label": label, "qP": qP, "qK": qK, "iface_P": len(rP), "iface_K": len(rK),
            "saltbridge": saltbridge, "clash": clash, "EC": round(EC, 2),
            "EC_per_res": round(EC/npair, 3)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cif"); ap.add_argument("--presenter", required=True)
    ap.add_argument("--kras", required=True); ap.add_argument("--label", default="")
    a = ap.parse_args()
    analyze(a.cif, a.presenter, a.kras, a.label or a.cif)
