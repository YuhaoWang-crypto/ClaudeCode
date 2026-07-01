#!/usr/bin/env python3
"""Sequence-based developability / liability flags for peptide candidates.
Heuristics only (no external DB): net charge, GRAVY hydrophobicity, aromatic
content, Cys parity (for disulfide staples), and common sequence liabilities:
  - NG / NS  -> Asn deamidation
  - DG / DP  -> Asp isomerization / cleavage
  - Met, Trp -> oxidation-prone
  - N-x-S/T  -> potential N-glycosylation sequon (recombinant expression)
  - runs of >4 hydrophobic -> aggregation risk
"""
import json, re

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,
'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,
'W':-0.9,'Y':-1.3,'V':4.2}
POS='KR'; NEG='DE'; AROM='FWY'; HPHOB='AILMFVWY'

def flags(seq):
    f=[]
    for m in re.finditer('N[GS]', seq): f.append(f"deamid@{m.start()+1}")
    for m in re.finditer('D[GP]', seq): f.append(f"isomer@{m.start()+1}")
    for m in re.finditer('N[^P][ST]', seq): f.append(f"glycosite@{m.start()+1}")
    if seq.count('M'): f.append(f"Met_ox({seq.count('M')})")
    for m in re.finditer('['+HPHOB+']{5,}', seq): f.append(f"hydrophobic_run@{m.start()+1}")
    return f

def analyze(seq):
    n=len(seq)
    pos=sum(seq.count(a) for a in POS); neg=sum(seq.count(a) for a in NEG)
    cys=seq.count('C')
    return {
      "len":n,
      "net_charge": pos-neg,
      "frac_pos": round(pos/n,3), "frac_neg": round(neg/n,3),
      "gravy": round(sum(KD[a] for a in seq)/n,2),
      "aromatic_frac": round(sum(seq.count(a) for a in AROM)/n,3),
      "cys": cys, "cys_paired": cys%2==0,
      "liabilities": flags(seq),
      "n_liabilities": len(flags(seq)),
    }

if __name__=="__main__":
    lib=json.load(open("design/library_v1.json"))
    out=[]
    for c in lib:
        d=analyze(c['seq']); d.update({k:c[k] for k in ('lib_id','id','seq','family','note')})
        out.append(d)
    json.dump(out, open("design/library_developability.json","w"), indent=1)
    # quick summary
    import statistics as st
    print("candidates:", len(out))
    print("mean len:", round(st.mean(x['len'] for x in out),1))
    print("mean aromatic frac:", round(st.mean(x['aromatic_frac'] for x in out),3))
    print("with >=1 liability:", sum(1 for x in out if x['n_liabilities']>0))
    print("odd-Cys (unpaired):", sum(1 for x in out if not x['cys_paired']))
