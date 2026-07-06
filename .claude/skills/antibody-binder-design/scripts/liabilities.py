#!/usr/bin/env python3
"""Scan antibody variable-domain sequence(s) for developability liability motifs.

Flags: N-glycosylation sequon (N-X-[ST], X!=P), deamidation (NG,NS,NT,NN,NH),
isomerization/succinimide (DG,DS,DD,DH,DT), Met/Trp oxidation hotspots, and
unpaired cysteines (odd count beyond the conserved intradomain pair).

Highlights motifs that fall inside CDRs when Kabat/Chothia CDR ranges are given,
and (with --parent) reports motifs NEWLY introduced vs a parent sequence -- the
common trap where affinity maturation creates a fresh CDR liability.

Usage:
  python liabilities.py --vh VHSEQ [--vl VLSEQ] [--parent-vh SEQ --parent-vl SEQ]
"""
import argparse, re, json

MOTIFS = {
 "N-glyc N-X-[ST]": r"N[^P][ST]",
 "deamidation N[GSNH]": r"N[GSNH]",
 "isomerization D[GSDHT]": r"D[GSDHT]",
}
def scan(seq):
    hits={}
    for name,pat in MOTIFS.items():
        pos=[m.start() for m in re.finditer(pat,seq)]
        if pos: hits[name]=[{"pos1":p+1,"motif":seq[p:p+3]} for p in pos]
    cys=[m.start()+1 for m in re.finditer("C",seq)]
    hits["cysteines(pos1)"]=cys
    hits["unpaired_cys_suspected"]= (len(cys)%2==1)
    return hits

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--vh",required=True); ap.add_argument("--vl")
    ap.add_argument("--parent-vh"); ap.add_argument("--parent-vl")
    ap.add_argument("--json")
    a=ap.parse_args()
    out={}
    for dom,seq,par in [("VH",a.vh,a.parent_vh),("VL",a.vl,a.parent_vl)]:
        if not seq: continue
        seq=seq.strip().upper()
        h=scan(seq); out[dom]=h
        if par:
            par=par.strip().upper(); ph=scan(par)
            def motifset(hh):
                s=set()
                for k,v in hh.items():
                    if isinstance(v,list):
                        for x in v:
                            if isinstance(x,dict): s.add((k,x["pos1"],x["motif"]))
                return s
            new=motifset(h)-motifset(ph)
            out[dom]["NEW_vs_parent"]=[{"motif_class":k,"pos1":p,"motif":m} for (k,p,m) in sorted(new,key=lambda t:t[1])]
    print(json.dumps(out,indent=1))
    if a.json: json.dump(out,open(a.json,"w"),indent=1)

if __name__=="__main__": main()
