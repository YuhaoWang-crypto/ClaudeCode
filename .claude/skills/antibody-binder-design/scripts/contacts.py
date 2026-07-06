#!/usr/bin/env python3
"""Epitope/paratope contact analysis + PTM-occlusion check from a co-fold CIF.

Given a predicted complex (Boltz/Chai mmCIF), compute for each antigen residue the
minimum heavy-atom distance to the antibody, list the epitope (<= cutoff), and
optionally check whether a modification anchor atom (e.g. a lipidated Lys Nz) is
buried by the paratope -- the failure mode where a naive unmodified-peptide model
looks fine but the real modified drug would clash.

Usage:
  python contacts.py complex.cif --ab H,L --antigen P [--cutoff 4.5] \
        [--occlusion P:20:NZ]   # chain:auth_seq:atom_name of the modification anchor
"""
import argparse, math, json

def read_atoms(fn):
    # returns list of (auth_asym, auth_seq, comp, atom, x, y, z)
    cols=None; out=[]
    with open(fn) as f:
        order=[]
        for line in f:
            s=line.strip()
            if s.startswith("_atom_site."):
                order.append(s.split(".")[1]); continue
            if line.startswith(("ATOM","HETATM")):
                p=line.split()
                if not order or len(p)<len(order): 
                    # fallback fixed Boltz layout
                    out.append((p[15],p[7],p[5],p[3],float(p[10]),float(p[11]),float(p[12]))); continue
                d=dict(zip(order,p))
                out.append((d.get("auth_asym_id"),d.get("auth_seq_id"),d.get("auth_comp_id",d.get("label_comp_id")),
                            d.get("label_atom_id"),float(d["Cartn_x"]),float(d["Cartn_y"]),float(d["Cartn_z"])))
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("cif"); ap.add_argument("--ab",required=True,help="antibody chains, comma sep e.g. H,L")
    ap.add_argument("--antigen",required=True,help="antigen chain id e.g. P")
    ap.add_argument("--cutoff",type=float,default=4.5)
    ap.add_argument("--occlusion",help="chain:auth_seq:atom_name of modification anchor, e.g. P:20:NZ")
    ap.add_argument("--json")
    a=ap.parse_args()
    ab=set(a.ab.split(",")); at=read_atoms(a.cif)
    abo=[x for x in at if x[0] in ab]; ag=[x for x in at if x[0]==a.antigen]
    if not abo or not ag: raise SystemExit("no atoms for given chains; check --ab/--antigen ids")
    import collections
    g=collections.defaultdict(list)
    for x in abo: g[(int(x[4]//5),int(x[5]//5),int(x[6]//5))].append(x)
    def nearest(px,py,pz):
        kx,ky,kz=int(px//5),int(py//5),int(pz//5); m=1e9
        for dx in(-1,0,1):
         for dy in(-1,0,1):
          for dz in(-1,0,1):
            for b in g.get((kx+dx,ky+dy,kz+dz),()):
                d=math.dist((px,py,pz),b[4:7]); m=min(m,d)
        return m
    perres={}
    for x in ag:
        d=nearest(x[4],x[5],x[6]); r=x[1]
        if r not in perres or d<perres[r][0]: perres[r]=(d,x[2])
    epi=sorted([(int(r),perres[r][1],round(perres[r][0],2)) for r in perres if perres[r][0]<=a.cutoff],key=lambda t:t[0])
    res={"epitope_cutoff":a.cutoff,"n_epitope":len(epi),
         "epitope":[f"{c}{r}({d}A)" for r,c,d in epi]}
    if a.occlusion:
        ch,seq,atom=a.occlusion.split(":")
        anc=[x for x in ag if x[0]==ch and x[1]==seq and x[3]==atom]
        if anc:
            dmin=min(nearest(x[4],x[5],x[6]) for x in anc)
            res["occlusion"]={"anchor":a.occlusion,"min_dist_A":round(dmin,2),
                "buried":dmin<5.0,
                "verdict":"BURIED -> real modification would clash" if dmin<5.0 else "solvent-exposed -> modification tolerated"}
        else:
            res["occlusion"]={"anchor":a.occlusion,"error":"anchor atom not found"}
    print(json.dumps(res,indent=1))
    if a.json: json.dump(res,open(a.json,"w"),indent=1)

if __name__=="__main__": main()
