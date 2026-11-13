#!/usr/bin/env python3
"""
降解分支的"诱导邻近"分析器（与遮盖分支的静电互补判据不同）。
对 CRBN : PROTAC : STAT3 三元，评估：
  1) PROTAC 两端是否分别锚定两蛋白（双功能桥是否成立）
  2) CRBN–STAT3 诱导界面（接触残基数、最近重原子）
  3) STAT3 赖氨酸递呈：每个 Lys Nζ 到 CRBN 的距离，统计"朝 CRBN 一侧"的 Lys
     —— 泛素化需要把靶蛋白的 Lys 递到 E2~Ub 可及区；这里用"朝 CRBN 递呈的可及 Lys"做方向性代理。
纯 numpy。用法: python3 proximity_analyze.py <cif> --crbn C --target S --protac L
"""
import sys, argparse, math
import numpy as np

def parse_cif(path):
    atoms=[]; L=open(path).read().splitlines(); i=0
    while i<len(L):
        if L[i].startswith("_atom_site."):
            cols=[]
            while i<len(L) and L[i].startswith("_atom_site."): cols.append(L[i].strip().split(".",1)[1]); i+=1
            idx={c:k for k,c in enumerate(cols)}
            def c(*n):
                for x in n:
                    if x in idx: return idx[x]
            ci={"g":c("group_PDB"),"a":c("label_atom_id"),"comp":c("label_comp_id"),
                "ch":c("auth_asym_id","label_asym_id"),"seq":c("auth_seq_id","label_seq_id"),
                "el":c("type_symbol"),"x":c("Cartn_x"),"y":c("Cartn_y"),"z":c("Cartn_z"),"m":c("pdbx_PDB_model_num")}
            while i<len(L):
                s=L[i].strip()
                if s=="" or s[0]=="#" or s.startswith(("loop_","_","data_")): break
                p=s.split()
                if len(p)>=len(cols) and (ci["m"] is None or p[ci["m"]] in (".","1")):
                    try: atoms.append({"g":p[ci["g"]],"a":p[ci["a"]].strip('"'),"comp":p[ci["comp"]],
                                       "ch":p[ci["ch"]],"seq":p[ci["seq"]],
                                       "el":(p[ci["el"]].strip() if ci["el"] is not None else ""),
                                       "xyz":np.array([float(p[ci["x"]]),float(p[ci["y"]]),float(p[ci["z"]])])})
                    except: pass
                i+=1
            continue
        i+=1
    return atoms

def heavy(atoms,ch,grp=None):
    return [a for a in atoms if a["ch"]==ch and a["el"]!="H" and (grp is None or a["g"]==grp)]

def mindist(A,B):
    XA=np.array([a["xyz"] for a in A]); XB=np.array([a["xyz"] for a in B])
    return math.sqrt(((XA[:,None,:]-XB[None,:,:])**2).sum(-1).min())

def contacts(A,B,cut=4.5):
    XA=np.array([a["xyz"] for a in A]); XB=np.array([a["xyz"] for a in B])
    d2=((XA[:,None,:]-XB[None,:,:])**2).sum(-1)
    ia,ib=np.where(d2<=cut*cut)
    rA={(A[k]["comp"],A[k]["seq"]) for k in ia}; rB={(B[k]["comp"],B[k]["seq"]) for k in ib}
    return rA,rB,int((d2<=cut*cut).sum())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("cif"); ap.add_argument("--crbn",required=True)
    ap.add_argument("--target",required=True); ap.add_argument("--protac",required=True)
    ap.add_argument("--label",default="")
    a=ap.parse_args()
    atoms=parse_cif(a.cif)
    C=heavy(atoms,a.crbn,"ATOM"); T=heavy(atoms,a.target,"ATOM")
    Lig=[x for x in atoms if x["ch"]==a.protac and x["el"]!="H"]
    print(f"\n########## {a.label or a.cif} ##########")
    print(f"CRBN链{a.crbn}:{len(C)} | 靶(STAT3)链{a.target}:{len(T)} | PROTAC链{a.protac}:{len(Lig)}原子")
    if not (C and T and Lig): print("  !!缺链"); return

    # 1) PROTAC 双端锚定
    cC,_,pC=contacts(C,Lig); cT,_,pT=contacts(T,Lig)
    print(f"\n[1] PROTAC 双功能桥：")
    print(f"  PROTAC↔CRBN: 接触 {len(cC)} 残基 ({pC} 原子对)  最近 {mindist(C,Lig):.2f} Å")
    print(f"  PROTAC↔STAT3: 接触 {len(cT)} 残基 ({pT} 原子对)  最近 {mindist(T,Lig):.2f} Å")
    bridges=len(cC)>0 and len(cT)>0
    print(f"  => 双端同时锚定: {'是 ✓ (诱导邻近成立)' if bridges else '否 ✗'}")
    # CRBN 芳香袋(Trp)是否接触谷氨酰亚胺端
    trp=[x for x in C if x["comp"]=="TRP"]
    if trp: print(f"  CRBN Trp 口袋 ↔ PROTAC 最近 {mindist(trp,Lig):.2f} Å (IMiD 识别袋)")
    # STAT3 Arg609 磷酸袋是否接触膦酸端
    r609=[x for x in T if x["seq"]=="609"]
    if r609: print(f"  STAT3 Arg609 ↔ PROTAC 最近 {mindist(r609,Lig):.2f} Å (pTyr袋/膦酸端)")

    # 2) CRBN–STAT3 诱导界面
    rC,rT,pp=contacts(C,T)
    print(f"\n[2] CRBN–STAT3 诱导界面：接触 {len(rC)}/{len(rT)} 残基, {pp} 原子对, 最近 {mindist(C,T):.2f} Å")

    # 3) STAT3 赖氨酸递呈
    nz=[x for x in T if x["comp"]=="LYS" and x["a"]=="NZ"]
    Ccen=np.mean([x["xyz"] for x in C],axis=0); Tcen=np.mean([x["xyz"] for x in T],axis=0)
    print(f"\n[3] STAT3 赖氨酸递呈 (共 {len(nz)} 个 Lys)：")
    rows=[]
    for k in nz:
        dC=mindist([k],C)                     # Nζ 到 CRBN 最近重原子
        side = np.linalg.norm(k["xyz"]-Ccen) < np.linalg.norm(k["xyz"]-Tcen)  # 更靠CRBN?
        rows.append((int(k["seq"]),dC,side))
    rows.sort(key=lambda r:r[1])
    zone=[r for r in rows if r[1]<=25.0]      # 朝CRBN一侧、25Å内 = 泛素化可及区代理
    for s,dC,side in rows:
        tag=" ★朝CRBN侧" if side else ""
        z=" [可及区<25Å]" if dC<=25 else ""
        print(f"  Lys{s}: Nζ→CRBN {dC:.1f} Å{z}{tag}")
    print(f"  => 朝 CRBN 递呈且在泛素化可及区(<25Å)的 Lys: {len([r for r in zone if r[2]])} 个 "
          f"(最近 Lys{rows[0][0]} @ {rows[0][1]:.1f} Å)")
    print(f"  => STAT3 表面 Lys 总数 {len(nz)} (泛素化'把手'丰度)")

if __name__=="__main__": main()
