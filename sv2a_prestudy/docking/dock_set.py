"""Dock calibration set into SV2A (8JLC) with AutoDock Vina 1.2.5.
usage: dock_set.py <receptor_prefix> <vina_bin> <out_dir>"""
import sys, os, json, subprocess, numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
pref, vina, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
box=json.load(open(pref+'_box.json'))['center']
cal=pd.read_csv('../sv2a_calibration_set_all.csv')
ref=pd.DataFrame([dict(molecule_chembl_id='CHEMBL607400',smiles='CCC[C@@H]1CC(=O)N([C@@H](CC)C(N)=O)C1')])
mols=pd.concat([cal.drop_duplicates('molecule_chembl_id')[['molecule_chembl_id','smiles']],ref]).drop_duplicates('molecule_chembl_id')
res=[]
for r in mols.itertuples():
    m=Chem.AddHs(Chem.MolFromSmiles(r.smiles))
    if AllChem.EmbedMolecule(m,randomSeed=7)!=0: print('embed fail',r.molecule_chembl_id); continue
    AllChem.MMFFOptimizeMolecule(m)
    pre=MoleculePreparation(); setup=pre.prepare(m)[0]
    pq,ok,err=PDBQTWriterLegacy.write_string(setup)
    f=f'{out}/{r.molecule_chembl_id}'; open(f+'.pdbqt','w').write(pq)
    cmd=[vina,'--receptor',pref+'_rec.pdbqt','--ligand',f+'.pdbqt','--center_x',str(box[0]),'--center_y',str(box[1]),'--center_z',str(box[2]),'--size_x','22','--size_y','22','--size_z','22',
         '--exhaustiveness','32','--num_modes','9','--seed','11','--cpu','4','--out',f+'_out.pdbqt']
    o=subprocess.run(cmd,capture_output=True,text=True).stdout
    sc=[float(l.split()[1]) for l in o.splitlines() if l.strip().startswith(('1 ','2 ','3 '))][:3]
    res.append(dict(molecule_chembl_id=r.molecule_chembl_id,best=sc[0] if sc else np.nan,n_heavy=m.GetNumHeavyAtoms()))
    print(r.molecule_chembl_id,sc[:1],flush=True)
pd.DataFrame(res).to_csv(out+'/scores.csv',index=False)
