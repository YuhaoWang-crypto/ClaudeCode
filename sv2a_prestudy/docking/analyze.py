import sys, json, numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFreeSASA
from meeko import PDBQTMolecule, RDKitMolCreate
from Bio.PDB import MMCIFParser, Superimposer
out='run_8JLC'; S='/tmp/claude-0/-home-user-ClaudeCode/ddb10227-2e96-5c4a-85ab-f3cb88392cf9/scratchpad/pdb/'
sc=pd.read_csv(out+'/scores.csv'); cal=pd.read_csv('../sv2a_calibration_set_all.csv')
def pose(cid,n=0):
    pm=PDBQTMolecule.from_file(f'{out}/{cid}_out.pdbqt',skip_typing=True)
    m=RDKitMolCreate.from_pdbqt_mol(pm)[0]; m=Chem.RemoveHs(m,sanitize=False)
    m.RemoveAllConformers() if False else None
    return m
def rmsd_to_xtal(cid,smiles,xtal_pdb,coords_xf=None):
    pm=PDBQTMolecule.from_file(f'{out}/{cid}_out.pdbqt',skip_typing=True)
    mol=RDKitMolCreate.from_pdbqt_mol(pm)[0]
    tmpl=Chem.MolFromSmiles(smiles)
    x=Chem.MolFromPDBFile(xtal_pdb,removeHs=True,sanitize=False)
    x=Chem.AssignBondsFromTemplate(tmpl,x)
    if coords_xf is not None:
        c=x.GetConformer(); 
        for i in range(x.GetNumAtoms()):
            p=np.array(c.GetAtomPosition(i)); q=coords_xf(p); c.SetAtomPosition(i,q.tolist())
    mol=Chem.RemoveHs(mol); res=[]
    for cid_ in range(mol.GetNumConformers()):
        mm=Chem.Mol(mol,confId=cid_) if False else Chem.Mol(mol); mm.RemoveAllConformers(); mm.AddConformer(mol.GetConformer(cid_),assignId=True)
        res.append(rdMolAlign.CalcRMS(mm,x))   # symmetry-aware, no alignment
    return res
# LEV redock
lev_smi='CC[C@@H](C(N)=O)N1CCCC1=O'
r_lev=rmsd_to_xtal('CHEMBL1286',lev_smi,'8JLC_lig_xtal.pdb')
print('LEV redock RMSD per mode (A):',np.round(r_lev,2).tolist())
# BRV: transform 8K77 -> 8JLC frame by CA superposition of TM residues
TM=[(170,190),(206,226),(234,254),(263,283),(295,315),(335,355),(448,468),(599,619),(627,647),(652,672),(686,708),(713,731)]
a=MMCIFParser(QUIET=True).get_structure('a',S+'8JLC.cif')[0]['A']; b=MMCIFParser(QUIET=True).get_structure('b',S+'8K77.cif')[0]['A']
ids=[i for lo,hi in TM for i in range(lo,hi+1) if i in a and i in b and 'CA' in a[i] and 'CA' in b[i]]
sup=Superimposer(); sup.set_atoms([a[i]['CA'] for i in ids],[b[i]['CA'] for i in ids])
print('8K77->8JLC TM CA RMSD %.2f A (n=%d)'%(sup.rms,len(ids)))
R,t=sup.rotran
brv_smi='CCC[C@@H]1CC(=O)N([C@@H](CC)C(N)=O)C1'
r_brv=rmsd_to_xtal('CHEMBL607400',brv_smi,'8K77_lig_xtal.pdb',coords_xf=lambda p:p@R+t)
print('BRV cross-dock RMSD vs 8K77 pose, per mode:',np.round(r_brv,2).tolist())
json.dump(dict(lev_redock=r_lev,brv_crossdock=r_brv,tm_ca_rmsd=sup.rms),open('pose_rmsd.json','w'))
# ranking
d=sc.merge(cal.groupby(['assay_chembl_id','standard_type','molecule_chembl_id'],as_index=False).pAct.first(),on='molecule_chembl_id')
rows=[]
for (aid,typ),g in d.groupby(['assay_chembl_id','standard_type']):
    g=g.dropna(subset=['best']); 
    if len(g)<5: continue
    rho,p=spearmanr(-g.best,g.pAct); tau,_=kendalltau(-g.best,g.pAct); rh,ph=spearmanr(g.n_heavy,g.pAct)
    rows.append(dict(assay=aid,type=typ,n=len(g),spearman_vina=rho,p=p,kendall=tau,spearman_heavyatoms=rh,p_ha=ph,pAct_span=g.pAct.max()-g.pAct.min()))
rk=pd.DataFrame(rows); print(rk.round(3).to_string()); rk.to_csv('ranking_stats.csv',index=False)
d.to_csv('docking_vs_activity.csv',index=False)
for c,n in [('CHEMBL1286','LEV'),('CHEMBL607400','BRV')]:
    print(n, sc[sc.molecule_chembl_id==c].best.tolist())
