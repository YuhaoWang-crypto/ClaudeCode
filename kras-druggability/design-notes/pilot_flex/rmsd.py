import sys, numpy as np
def parse(path):
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
                "x":c("Cartn_x"),"y":c("Cartn_y"),"z":c("Cartn_z"),"m":c("pdbx_PDB_model_num")}
            while i<len(L):
                s=L[i].strip()
                if s=="" or s[0]=="#" or s.startswith(("loop_","_","data_")): break
                p=s.split()
                if len(p)>=len(cols):
                    if ci["m"] is None or p[ci["m"]] in (".","1"):
                        try: atoms.append((p[ci["g"]],p[ci["a"]].strip('"'),p[ci["comp"]],p[ci["ch"]],p[ci["seq"]],
                                          np.array([float(p[ci["x"]]),float(p[ci["y"]]),float(p[ci["z"]])])))
                        except: pass
                i+=1
            continue
        i+=1
    return atoms
def ca(atoms,ch):
    d={}
    for g,a,comp,c,seq,xyz in atoms:
        if g=="ATOM" and a=="CA" and c==ch:
            try: d[int(seq)]=xyz
            except: pass
    return d
def kabsch(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(Pc.T@Qc)
    d=np.sign(np.linalg.det(Wt.T@V.T))
    D=np.diag([1,1,d]); R=Wt.T@D@V.T
    Pr=Pc@R.T
    return np.sqrt(((Pr-Qc)**2).sum(1).mean()), R, P.mean(0), Q.mean(0)
apo=parse("apo_cypa_2cpl.cif"); bound=parse("../pilot1/../../demo/tricomplex/structures/9BFX.cif")
# pick chains
apo_ch=sorted({c for g,a,comp,c,s,x in apo if g=="ATOM"})
print("apo CypA chains:",apo_ch)
A=ca(apo,apo_ch[0]); B=ca(bound,"D")  # 9BFX CypA = chain D
common=sorted(set(A)&set(B))
print(f"共有 Cα 残基数: {len(common)} (范围 {common[0]}–{common[-1]})")
P=np.array([A[k] for k in common]); Q=np.array([B[k] for k in common])
rmsd_all,R,pm,qm=kabsch(P,Q)
print(f"全局 Cα RMSD (apo CypA vs 结合态CypA) = {rmsd_all:.2f} Å")
# interface rim residues
rim=[55,57,69,71,73,121,122,125,148,149,151]
rim=[r for r in rim if r in A and r in B]
# apply global alignment then measure per-rim deviation
Pc=P-pm; Pr=Pc@R.T+qm
dev={k:np.linalg.norm((A[k]-pm)@R.T+qm - B[k]) for k in rim}
print(f"界面 rim 残基 Cα 偏差 (叠合后): 均值 {np.mean(list(dev.values())):.2f} Å, 最大 {max(dev.values()):.2f} Å")
for k in rim: print(f"   残基{k}: {dev[k]:.2f} Å")
