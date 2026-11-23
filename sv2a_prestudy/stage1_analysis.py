"""Stage 1 pre-study: SV2A (Q7L0J3) residue classification from 8JLC / 8K77.
Usage: python stage1_analysis.py <cif> <ligand_resname> <out_prefix>
Membrane: no OPM/PPM entry available offline -> slab fitted from UniProt TM annotation
(axis = mean TM-helix axis, centre = TM CA centroid, hydrophobic half-thickness 15 A). Approximation; replace by PPM3 output."""
import sys, json, numpy as np, pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
cif, lig, out = sys.argv[1:4]
s = MMCIFParser(QUIET=True).get_structure('x', cif)[0]
A = s['A']
up = json.load(open('../../../../tmp/claude-0/-home-user-ClaudeCode/ddb10227-2e96-5c4a-85ab-f3cb88392cf9/scratchpad/up.json')) if False else None
TM = [(170,190),(206,226),(234,254),(263,283),(295,315),(335,355),(448,468),(599,619),(627,647),(652,672),(686,708),(713,731)]
LD4 = (469,598)
CYT = [(1,169),(227,233),(284,294),(356,447),(620,626),(673,685),(732,742)]
MAXASA = dict(ALA=129,ARG=274,ASN=195,ASP=193,CYS=167,GLN=225,GLU=223,GLY=104,HIS=224,ILE=197,LEU=201,LYS=236,MET=224,PHE=240,PRO=159,SER=155,THR=172,TRP=285,TYR=263,VAL=174)
res = [r for r in A if r.id[0]==' ']
ca = {r.id[1]: r['CA'].coord for r in res if 'CA' in r}
# membrane normal
axes=[]
for a,b in TM:
    pts=np.array([ca[i] for i in range(a,b+1) if i in ca])
    if len(pts)<15: continue
    u,sv,vt=np.linalg.svd(pts-pts.mean(0)); ax=vt[0]
    if ax@(pts[-1]-pts[0])<0: ax=-ax
    axes.append(ax)
n = np.mean(axes,0); n/=np.linalg.norm(n)
# helices alternate direction (up/down) -> use sign-aligned axes
axes=[ax if ax@axes[0]>0 else -ax for ax in axes]
n=np.mean(axes,0); n/=np.linalg.norm(n)
tmca=np.array([ca[i] for a,b in TM for i in range(a,b+1) if i in ca]); c0=tmca.mean(0)
# SASA: protein chain A only (lipid exposure proxy), residue level
sr=ShrakeRupley(probe_radius=1.4, n_points=200)
from Bio.PDB.Structure import Structure; from Bio.PDB.Model import Model
import copy
mA=copy.deepcopy(A); 
for r in [r for r in mA if r.id[0]!=' ']: mA.detach_child(r.id)
sr.compute(mA, level='R')
sasa={r.id[1]:r.sasa for r in mA}
# ligand + toxin
L=[r for r in s.get_residues() if r.get_resname()==lig]
latoms=[a for r in L for a in r if a.element!='H']
B=[a for a in s['B'].get_atoms() if a.element!='H'] if 'B' in s else []
glyc=[a for ch in s for r in ch if r.get_resname() in('NAG','FUC','BMA','MAN','NDG') for a in r]
allA=[a for r in res for a in r if a.element!='H']
nsA=NeighborSearch(allA)
def dmin(r,atoms):
    ns=NeighborSearch(atoms) if atoms else None
    best=99
    for a in r:
        if a.element=='H': continue
        for b in atoms:
            d=np.linalg.norm(a.coord-b.coord)
            if d<best: best=d
    return best
rows=[]
for r in res:
    i=r.id[1]; name=r.get_resname()
    sc=[a for a in r if a.get_id() not in('N','C','O','CA') and a.element!='H'] or [r['CA']]
    cen=np.mean([a.coord for a in sc],0); z=float((cen-c0)@n)
    rsasa=sasa[i]/MAXASA.get(name,200)
    dl=dmin(r,latoms) if latoms else np.nan
    dt=dmin(r,B) if B else np.nan
    intm=any(a<=i<=b for a,b in TM)
    if intm: seg='TM%d'%([k for k,(a,b) in enumerate(TM) if a<=i<=b][0]+1)
    elif LD4[0]<=i<=LD4[1]: seg='LD4/lumenal'
    elif any(a<=i<=b for a,b in CYT): seg='cytoplasmic'
    else: seg='lumenal loop'
    if dl<=5: cls='pocket'
    elif dt<=4.5 and not intm: cls='toxin_interface'
    elif intm and dl<=8: cls='pocket_shell2'
    elif intm and abs(z)<=15 and rsasa>=0.20: cls='lipid_facing'
    elif intm and abs(z)<=15 and rsasa>=0.08: cls='lipid_borderline'
    elif intm and abs(z)<=15: cls='core_buried'
    elif intm: cls='TM_end_(out_of_slab)'
    elif seg=='cytoplasmic': cls='cytoplasmic_loop'
    else: cls='lumenal_loop/LD4'
    rows.append(dict(resnum=i,aa=name,seg=seg,z_A=round(z,1),rSASA=round(rsasa,2),d_lig=round(dl,1),d_toxin=round(dt,1),cls=cls))
df=pd.DataFrame(rows); df.to_csv(out+'_residues.csv',index=False)
print(cif,'resolved',len(df),'range',df.resnum.min(),df.resnum.max(),'ligand atoms',len(latoms))
print(df.cls.value_counts().to_string())
print('\nPocket residues (<=5 A):',' '.join(f"{r.aa}{r.resnum}" for r in df[df.cls=='pocket'].itertuples()))
print('\nToxin-interface (<=4.5 A, non-TM):',' '.join(f"{r.aa}{r.resnum}" for r in df[df.cls=='toxin_interface'].itertuples()))
# glycan contacts to toxin
for r in s['A']:
    if r.id[0]!=' ': pass
print('\nLipid-facing TM composition:',df[df.cls=='lipid_facing'].aa.value_counts().to_dict())
print('\nmissing resnums inside 1..742:',[i for i in range(df.resnum.min(),743) if i not in set(df.resnum)][:0] or 'n=%d'%(743-df.resnum.min()-len(df)))
