#!/usr/bin/env python3
"""
三元复合物界面分析：CYPA : 分子胶 : KRAS。
- 跨链接触（KRAS <-> CYPA 残基对，最近重原子 <= CUT）
- 分子胶是否同时桥接两个蛋白（neo-PPI 的充分条件）
- 蛋白-蛋白埋藏表面积 BSA（自实现 Shrake-Rupley SASA，纯 numpy）
依赖仅 numpy。距离/面积单位 Å / Å²。

用法:
  python3 tri_analyze.py <cif> --kras A --cypa B [--glue G]
"""
import sys, argparse, math
import numpy as np

VDW = {"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80,"F":1.47,
       "CL":1.75,"MG":1.73,"NA":2.27,"K":2.75,"ZN":1.39,"H":1.10}
PROBE = 1.40
CUT = 4.5

def elem_of(atom_name, type_symbol):
    if type_symbol and type_symbol.strip() not in (".","?"):
        return type_symbol.strip().upper()
    s = atom_name.strip('"').lstrip("0123456789")
    return (s[:2].upper() if s[:2].upper() in VDW else s[:1].upper())

def parse_cif(path):
    atoms=[]; lines=open(path).read().splitlines()
    i=0
    while i<len(lines):
        if lines[i].startswith("_atom_site."):
            cols=[]
            while i<len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].strip().split(".",1)[1]); i+=1
            idx={c:k for k,c in enumerate(cols)}
            def col(name,*alts):
                for n in (name,)+alts:
                    if n in idx: return idx[n]
                return None
            ci={"grp":col("group_PDB"),"atom":col("label_atom_id"),
                "comp":col("label_comp_id"),
                "ch":col("auth_asym_id","label_asym_id"),
                "seq":col("auth_seq_id","label_seq_id"),
                "sym":col("type_symbol"),
                "x":col("Cartn_x"),"y":col("Cartn_y"),"z":col("Cartn_z"),
                "model":col("pdbx_PDB_model_num")}
            while i<len(lines):
                s=lines[i].strip()
                if s=="" or s[0]=="#" or s.startswith("loop_") or s.startswith("_") or s.startswith("data_"):
                    break
                p=s.split()
                if len(p)<len(cols): i+=1; continue
                if ci["model"] is not None and p[ci["model"]] not in (".","1"):
                    i+=1; continue
                try:
                    atoms.append({
                        "grp":p[ci["grp"]],
                        "atom":p[ci["atom"]].strip('"'),
                        "comp":p[ci["comp"]],
                        "ch":p[ci["ch"]],
                        "seq":p[ci["seq"]],
                        "elem":elem_of(p[ci["atom"]], p[ci["sym"]] if ci["sym"] is not None else ""),
                        "xyz":np.array([float(p[ci["x"]]),float(p[ci["y"]]),float(p[ci["z"]])]),
                    })
                except (ValueError,IndexError): pass
                i+=1
            continue
        i+=1
    return atoms

def sasa(atoms, n_sphere=192):
    """Shrake-Rupley。atoms: list of dict with xyz, elem. 返回每原子 SASA 之和。"""
    if not atoms: return 0.0
    xyz=np.array([a["xyz"] for a in atoms])
    r=np.array([VDW.get(a["elem"],1.70)+PROBE for a in atoms])
    # golden-spiral 球面采样点
    k=np.arange(n_sphere); phi=math.pi*(3-math.sqrt(5))
    z=1-2*(k+0.5)/n_sphere; rad=np.sqrt(np.maximum(0,1-z*z)); th=phi*k
    sph=np.stack([rad*np.cos(th), rad*np.sin(th), z],axis=1)  # (S,3)
    # 邻居预筛
    from math import inf
    rmax=r.max()
    total=0.0
    # 构建简单网格邻居（原子数几千，直接 O(N^2) 对距离平方可接受）
    d2=((xyz[:,None,:]-xyz[None,:,:])**2).sum(-1)
    for i in range(len(atoms)):
        pts=xyz[i]+r[i]*sph             # (S,3)
        neigh=np.where((d2[i]<(r[i]+r)**2)&(np.arange(len(atoms))!=i))[0]
        if len(neigh)==0:
            total+=4*math.pi*r[i]**2; continue
        acc=np.ones(n_sphere,dtype=bool)
        for j in neigh:
            dd=((pts-xyz[j])**2).sum(-1)
            acc&=dd>r[j]**2
        total+=(acc.sum()/n_sphere)*4*math.pi*r[i]**2
    return total

def protein_atoms(atoms, ch):
    return [a for a in atoms if a["ch"]==ch and a["grp"]=="ATOM" and a["elem"]!="H"]

def chain_atoms_any(atoms, ch):
    return [a for a in atoms if a["ch"]==ch and a["elem"]!="H"]

def residue_min_contacts(aA, aB, cut=CUT):
    """返回 A 侧接触残基集合、B 侧接触残基集合、接触原子对数。"""
    if not aA or not aB: return set(),set(),0
    XA=np.array([a["xyz"] for a in aA]); XB=np.array([a["xyz"] for a in aB])
    d2=((XA[:,None,:]-XB[None,:,:])**2).sum(-1)
    mask=d2<=cut*cut
    resA=set(); resB=set(); pairs=int(mask.sum())
    ia,ib=np.where(mask)
    for k in range(len(ia)):
        resA.add((aA[ia[k]]["comp"],aA[ia[k]]["seq"]))
        resB.add((aB[ib[k]]["comp"],aB[ib[k]]["seq"]))
    return resA,resB,pairs

def fmt_res(resset):
    def key(t):
        try: return int(t[1])
        except: return 0
    return ", ".join(f"{c}{s}" for c,s in sorted(resset,key=key))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--kras",required=True); ap.add_argument("--cypa",required=True)
    ap.add_argument("--glue",default=None, help="胶所在链(取该链 HETATM)")
    ap.add_argument("--glue-comp",dest="glue_comp",default=None, help="按组分名选胶(如 A1AOD)，优先于 --glue")
    ap.add_argument("--label",default="")
    a=ap.parse_args()
    atoms=parse_cif(a.cif)
    K=protein_atoms(atoms,a.kras); C=protein_atoms(atoms,a.cypa)
    print(f"\n########## {a.label or a.cif} ##########")
    print(f"KRAS链 {a.kras}: {len(K)} 重原子 | CYPA链 {a.cypa}: {len(C)} 重原子")

    # 1) KRAS<->CYPA 直接界面
    rK,rC,pairs=residue_min_contacts(K,C)
    print(f"\n[CYPA–KRAS 蛋白-蛋白界面]  接触原子对={pairs}  界面残基: KRAS {len(rK)} / CYPA {len(rC)}")
    if rK: print(f"  KRAS 界面残基: {fmt_res(rK)}")
    if rC: print(f"  CYPA 界面残基: {fmt_res(rC)}")

    # 2) BSA（仅蛋白-蛋白）
    s_k=sasa(K); s_c=sasa(C); s_kc=sasa(K+C)
    bsa=(s_k+s_c-s_kc)/2.0
    print(f"\n[BSA 埋藏表面积 (CYPA–KRAS, 仅两蛋白)]")
    print(f"  SASA(KRAS)={s_k:.0f}  SASA(CYPA)={s_c:.0f}  SASA(复合)={s_kc:.0f}")
    print(f"  => 界面 BSA ≈ {bsa:.0f} Å²  (每侧埋藏 ≈ {bsa:.0f} Å²)")

    # 3) 分子胶桥接
    if a.glue_comp or a.glue:
        if a.glue_comp:
            G=[x for x in atoms if x["comp"]==a.glue_comp and x["elem"]!="H"]
        else:
            G=[x for x in atoms if x["ch"]==a.glue and x["grp"]=="HETATM" and x["elem"]!="H"]
        gK_K,_,pK=residue_min_contacts(K,G)
        gC_C,_,pC=residue_min_contacts(C,G)
        print(f"\n[分子胶桥接检验]  胶链 {a.glue}: {len(G)} 重原子")
        print(f"  胶 ↔ KRAS: 接触 {len(gK_K)} 残基 ({pK} 原子对): {fmt_res(gK_K)}")
        print(f"  胶 ↔ CYPA: 接触 {len(gC_C)} 残基 ({pC} 原子对): {fmt_res(gC_C)}")
        bridges = len(gK_K)>0 and len(gC_C)>0
        print(f"  => 胶同时接触两个蛋白: {'是 ✓ (molecular glue 几何成立)' if bridges else '否'}")
        # KRAS 第12位到胶（共价靶点/热点）
        pos12=[a2 for a2 in K if a2["seq"] in ("12",)]
        if pos12 and G:
            XK=np.array([x["xyz"] for x in pos12]); XG=np.array([x["xyz"] for x in G])
            d=math.sqrt(((XK[:,None,:]-XG[None,:,:])**2).sum(-1).min())
            print(f"  KRAS 第12位({pos12[0]['comp']}12) ↔ 胶 最近 = {d:.2f} Å")

if __name__=="__main__":
    main()
