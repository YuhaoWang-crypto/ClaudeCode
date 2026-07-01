#!/usr/bin/env python3
"""
Combine Boltz-2 screen metrics with local developability to produce the final
wet-lab candidate panel.

Screen metrics per candidate:
  iptm                 interface pTM (0-1, higher better; >0.8 confident geom)
  binding_confidence   Boltz-2 binder head (0-1, higher = predicted binder)
  structure_confidence monomer/complex confidence
  min_interaction_pae  interface PAE, Angstrom (lower better)

Composite (z-scored within the screened set), weighted toward the two most
trustworthy signals for a peptide target (ipTM + interface PAE), with
binding_confidence as a strong tie-breaker and a developability penalty.
"""
import json, sys, statistics as st

def load_screen(path):
    d=json.load(open(path)); out={}
    for r in d['results']:
        m=r['metrics']
        out[r['external_id']]={
          "iptm":m.get('iptm'),"bind_conf":m.get('binding_confidence'),
          "struct_conf":m.get('structure_confidence'),
          "min_ipae":m.get('min_interaction_pae'),
          "helix":m.get('helix_fraction'),"sheet":m.get('sheet_fraction'),
          "loop":m.get('loop_fraction'),
          "cif":r['artifacts']['structure']['url'] if r.get('artifacts') else None,
        }
    return out

def zscores(rows, key, direction=1):
    vals=[r[key] for r in rows if isinstance(r.get(key),(int,float))]
    m=st.mean(vals); s=st.pstdev(vals) or 1.0
    for r in rows:
        r[f"z_{key}"]=direction*((r[key]-m)/s) if isinstance(r.get(key),(int,float)) else 0.0

def composite(rows):
    zscores(rows,"iptm",+1); zscores(rows,"bind_conf",+1)
    zscores(rows,"min_ipae",-1); zscores(rows,"struct_conf",+1)
    W={"z_iptm":2.0,"z_min_ipae":1.5,"z_bind_conf":1.5,"z_struct_conf":0.5}
    for r in rows:
        base=sum(W[k]*r[k] for k in W)/sum(W.values())
        pen=0.10*r.get("n_liabilities",0) + (0.5 if not r.get("cys_paired",True) else 0)
        r["score"]=round(base-pen,4)

if __name__=="__main__":
    new_screen=sys.argv[1]      # results json for 100 new
    dev=json.load(open("design/library_developability.json"))
    devmap={d['lib_id']:d for d in dev}
    sc=load_screen(new_screen)
    rows=[]
    for lib_id,metrics in sc.items():
        d=devmap.get(lib_id,{})
        row={**metrics,**{k:d.get(k) for k in
             ('seq','family','note','n_liabilities','liabilities','cys_paired',
              'net_charge','gravy','aromatic_frac','len')}}
        row['lib_id']=lib_id
        rows.append(row)
    composite(rows)
    rows.sort(key=lambda r:-r['score'])
    json.dump(rows, open("results/new_library_ranked.json","w"), indent=1, default=str)
    print(f"{'rank':>4} {'id':7} {'fam':16} {'iptm':>5} {'bind':>5} {'mPAE':>5} {'struct':>6} {'liab':>4} {'score':>6}")
    for i,r in enumerate(rows[:30],1):
        print(f"{i:>4} {r['lib_id']:7} {str(r['family'])[:16]:16} {r['iptm']:5.2f} {r['bind_conf']:5.2f} {r['min_ipae']:5.2f} {r['struct_conf']:6.2f} {r['n_liabilities']:>4} {r['score']:6.2f}")
