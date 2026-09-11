#!/usr/bin/env python
"""Does the published loss weighting starve the proteome-prediction task?

Why this exists
---------------
`scripts/optimize_model.py` runs at the published setting, fixed lambda = 0.8,
which the Supplementary Information confirms is what the reported results used.
At that setting every variant scores a delta correlation near zero, and an
earlier version of this repository concluded from that "nothing learns the
direction of the response".

That conclusion was wrong, and this script is what shows it. lambda = 0.8 puts
only 0.2 weight on Loss1, the proteome reconstruction, while the ridge control it
was being compared against optimises the trajectory alone. The comparison was
therefore not like-for-like. Sweeping lambda separates two explanations:

  (a) the architecture cannot represent the response, or
  (b) the published objective barely trains it.

The answer is (b), with a residue of (a): see docs/OPTIMIZATION.md.
"""
import sys, os, copy, numpy as np, torch, torch.nn as nn
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteintalks import SyntheticPerturbationProteome, make_splits, trajectory_metrics
from proteintalks.improved import ProteinTalksR
torch.set_num_threads(os.cpu_count() or 4)
sim=SyntheticPerturbationProteome(n_proteins=200,n_cell_lines=18,n_drugs=62,seed=0)
ds=sim.generate(); tr,va,te=make_splits(ds,setting=1,seed=0)
def pack(i): return (torch.tensor(ds.p0[i][...,None]),torch.tensor(ds.pert[i][...,None]),
                     torch.tensor(ds.drug_feats[i][:,0][...,None]),torch.tensor(ds.drug_feats[i][:,1][...,None]),
                     torch.tensor(ds.p_future[i]),torch.tensor(ds.label[i]))
Xtr,Xva=pack(tr),pack(va)
print(f"{'lambda':>7s} {'flags':22s} {'skill':>9s} {'delta_r':>9s}")
print('-'*52)
for lam in [0.8, 0.5, 0.0]:
    for tag,cfg in [('released',dict(delta_decode=False,real_time=False,time_conditioned=False,pert_conditioned=False,coupling_rank=0,head_rank=0,substeps=1)),
                    ('all five changes',dict(delta_decode=True,real_time=True,time_conditioned=True,pert_conditioned=True,coupling_rank=8,head_rank=64,substeps=2))]:
        sk,dr=[],[]
        for seed in [0,1]:
            torch.manual_seed(seed); np.random.seed(seed)
            m=ProteinTalksR(pro_feats=ds.n_proteins,hidden_feats=32,drug_feature_feats=935,dropout=0.0,**cfg)
            opt=torch.optim.AdamW(m.parameters(),lr=5e-4,weight_decay=1e-4)
            mse,bce=nn.MSELoss(),nn.BCEWithLogitsLoss()
            best,bs,bad=np.inf,None,0
            for ep in range(60):
                m.train(); perm=torch.randperm(len(tr))
                for i in range(0,len(perm),64):
                    b=perm[i:i+64]; x,p,a,c,y,lab=(t[b] for t in Xtr)
                    opt.zero_grad(); pf,lg=m(x,p,a,c)
                    ((1-lam)*mse(pf,y)+lam*bce(lg,lab)).backward()
                    torch.nn.utils.clip_grad_norm_(m.parameters(),1.0); opt.step()
                m.eval()
                with torch.no_grad():
                    x,p,a,c,y,lab=Xva; pf,lg=m(x,p,a,c)
                    v=(1-lam)*mse(pf,y).item()+lam*bce(lg,lab).item()
                if v<best-1e-6: best,bad,bs=v,0,copy.deepcopy(m.state_dict())
                else:
                    bad+=1
                    if bad>=20: break
            m.load_state_dict(bs); m.eval()
            with torch.no_grad():
                x,p,a,c,_,_=pack(te); pf,_=m(x,p,a,c)
            t_=trajectory_metrics(ds.p_future[te],pf.numpy(),ds.p0[te])
            sk.append(t_['skill_vs_nochange']); dr.append(t_['delta_pearson_r'])
        print(f"{lam:7.1f} {tag:22s} {np.mean(sk):9.3f} {np.nanmean(dr):9.3f}")
print('\nridge (trajectory only, no ODE):        +0.971    +0.948   [48h]')
