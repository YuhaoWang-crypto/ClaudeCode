import sys, json, re, csv, time
import numpy as np, torch, esm
from scipy.stats import spearmanr

D="/home/user/ClaudeCode/protein_design_benchmark/data"
RES="/home/user/ClaudeCode/protein_design_benchmark/results"
MODEL=sys.argv[1] if len(sys.argv)>1 else "esm2_t12_35M_UR50D"
MAXMUT=int(sys.argv[2]) if len(sys.argv)>2 else 4000   # cap single-mutants scored per assay (CPU)

seqs=json.load(open(f"{D}/seqs.json"))
assays={
 "SPG1_STRSG_Olson_2014":"GB1 (IgG-Fc binding)",
 "BLAT_ECOLX_Stiffler_2015":"beta-lactamase (Amp resistance)",
 "GFP_AEQVI_Sarkisyan_2016":"GFP (fluorescence)",
}
print(f"loading {MODEL} ...",flush=True)
model,alphabet=esm.pretrained.__dict__[MODEL]()
model.eval(); bc=alphabet.get_batch_converter()
tok=lambda s:[alphabet.get_idx(a) for a in s]

def logits_full(seq):
    # one forward pass, return log-softmax logits [L, vocab] over the sequence (no mask)
    _,_,toks=bc([("wt",seq)])
    with torch.no_grad():
        lg=model(toks)["logits"][0]          # [L+2, vocab] incl BOS/EOS
    return torch.log_softmax(lg,dim=-1)

def masked_logprobs(seq, positions):
    # masked-marginal: for each needed 0-based residue position, mask it and read log-softmax
    out={}
    midx=alphabet.mask_idx
    uniq=sorted(set(positions))
    for i in range(0,len(uniq),32):
        batch=uniq[i:i+32]
        rows=[]
        _,_,base=bc([("wt",seq)])
        for p in batch:
            t=base.clone(); t[0,p+1]=midx; rows.append(t)
        toks=torch.cat(rows,0)
        with torch.no_grad():
            lg=model(toks)["logits"]
        lg=torch.log_softmax(lg,dim=-1)
        for j,p in enumerate(batch):
            out[p]=lg[j,p+1]
    return out

results={}
for aid,label in assays.items():
    seq=seqs[aid]["seq"]
    rows=[r for r in csv.DictReader(open(f"{D}/{aid}.csv"))]
    # single mutants only, well-formed, WT matches
    singles=[]
    for r in rows:
        m=r["mutant"]
        if ":" in m or ";" in m: continue
        mm=re.match(r"^([A-Z])(\d+)([A-Z])$",m)
        if not mm: continue
        wt,pos,mut=mm.group(1),int(mm.group(2)),mm.group(3)
        if pos<1 or pos>len(seq) or seq[pos-1]!=wt: continue
        try: y=float(r["DMS_score"])
        except: continue
        singles.append((wt,pos-1,mut,y))
    if len(singles)>MAXMUT:
        import random; random.seed(0); singles=random.sample(singles,MAXMUT)
    positions=[s[1] for s in singles]
    t0=time.time()
    lw=logits_full(seq)                      # wildtype-marginal
    lm=masked_logprobs(seq,positions)        # masked-marginal
    wt_m,mk_m,y=[],[],[]
    for wt,pos,mut,yy in singles:
        wi,mi=alphabet.get_idx(wt),alphabet.get_idx(mut)
        wt_m.append((lw[pos+1,mi]-lw[pos+1,wi]).item())
        mk_m.append((lm[pos][mi]-lm[pos][wi]).item())
        y.append(yy)
    rw=spearmanr(wt_m,y).correlation
    rm=spearmanr(mk_m,y).correlation
    results[aid]={"label":label,"n":len(singles),"len":len(seq),
                  "spearman_wt_marginal":round(rw,3),"spearman_masked_marginal":round(rm,3)}
    print(f"{aid:28s} n={len(singles):5d} L={len(seq):4d} "
          f"rho_wt={rw:.3f} rho_masked={rm:.3f}  ({time.time()-t0:.0f}s)",flush=True)

out=f"{RES}/esm_{MODEL}.json"
json.dump(results,open(out,"w"),indent=1)
print("written",out)
