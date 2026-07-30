"""
In silico toxicological qualification of two vancomycin B impurities.
Structure-based + 3D shape read-across analysis using RDKit.
"""
import json, math
from rdkit import Chem
from rdkit.Chem import (Descriptors, rdMolDescriptors, AllChem, Draw, QED,
                        Crippen, rdMolAlign, rdShapeHelpers)
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

OUT = "/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"

MOLS = {
 "Parent": ("Vancomycin B (parent API)",
  "OC1=CC([C@@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1"),
 "RRT087": ("RRT 0.87  (7S,26R)-Vancomycin B",
  "OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1"),
 "RRT075": ("RRT 0.75  (26R)-DMe-DeLI",
  "OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](N)[C@H](CC)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1"),
}

def mk(smi):
    m = Chem.MolFromSmiles(smi)
    return m

mset = {k:(name, mk(smi), smi) for k,(name,smi) in MOLS.items()}

# ---------- 1. Physicochemical descriptors ----------
def physchem(m):
    d = {}
    d["MolFormula"] = rdMolDescriptors.CalcMolFormula(m)
    d["MW"] = round(Descriptors.MolWt(m),2)
    d["ExactMW"] = round(Descriptors.ExactMolWt(m),4)
    d["cLogP_Crippen"] = round(Crippen.MolLogP(m),2)
    d["TPSA"] = round(rdMolDescriptors.CalcTPSA(m),1)
    d["HBD"] = rdMolDescriptors.CalcNumHBD(m)
    d["HBA"] = rdMolDescriptors.CalcNumHBA(m)
    d["RotBonds"] = rdMolDescriptors.CalcNumRotatableBonds(m)
    d["AromaticRings"] = rdMolDescriptors.CalcNumAromaticRings(m)
    d["Rings"] = rdMolDescriptors.CalcNumRings(m)
    d["HeavyAtoms"] = m.GetNumHeavyAtoms()
    d["FractionCsp3"] = round(rdMolDescriptors.CalcFractionCSP3(m),3)
    d["QED"] = round(QED.qed(m),3)
    centers = Chem.FindMolChiralCenters(m, useLegacyImplementation=False, includeUnassigned=True)
    d["StereoCenters"] = len(centers)
    # Lipinski / Veber (informational only for a glycopeptide)
    d["Lipinski_violations"] = sum([d["MW"]>500, d["cLogP_Crippen"]>5, d["HBD"]>5, d["HBA"]>10])
    d["Veber_pass"] = (d["RotBonds"]<=10 and d["TPSA"]<=140)
    return d

# ---------- 2. Structural alerts ----------
def alerts(m):
    res = {}
    for cat_name, cat_enum in [("PAINS",FilterCatalogParams.FilterCatalogs.PAINS),
                               ("BRENK",FilterCatalogParams.FilterCatalogs.BRENK),
                               ("NIH",FilterCatalogParams.FilterCatalogs.NIH)]:
        p = FilterCatalogParams(); p.AddCatalog(cat_enum)
        cat = FilterCatalog(p)
        hits=[]
        for e in cat.GetMatches(m):
            hits.append(e.GetDescription())
        res[cat_name]=sorted(set(hits))
    return res

# ---------- 3. Functional-group / toxicophore inventory ----------
FG = {
 "phenol":"[OX2H][cX3]:[c]",
 "aromatic_chloride(ArCl)":"[cX3][Cl]",
 "primary_aliphatic_amine":"[NX3;H2;!$(NC=O)][CX4]",
 "secondary_aliphatic_amine(NHMe/NHR)":"[NX3;H1;!$(NC=O)]([CX4])[CX4]",
 "primary_amide":"[NX3H2][CX3]=[OX1]",
 "secondary_amide":"[NX3H1]([#6])[CX3]=[OX1]",
 "carboxylic_acid":"[CX3](=O)[OX2H1]",
 "aromatic_ring":"c1ccccc1",
 "catechol/ortho-diphenol":"[OX2H]c1ccccc1[OX2H]",
 "amino_sugar":"[NX3;H2,H1][CX4]([OX2])",
 "benzylic_alcohol":"[OX2H][CX4][c]",
 "aryl_ether":"[cX3][OX2][#6]",
}
def fg_inventory(m):
    out={}
    for name,sma in FG.items():
        q=Chem.MolFromSmarts(sma)
        if q is None:
            out[name]="BAD_SMARTS"; continue
        out[name]=len(m.GetSubstructMatches(q, uniquify=True))
    return out

# ---------- 4. 2D similarity to parent ----------
from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
def morgan(m): return AllChem.GetMorganFingerprintAsBitVect(m,2,2048)
parent_m = mset["Parent"][1]
p_ecfp = morgan(parent_m); p_maccs = MACCSkeys.GenMACCSKeys(parent_m)

results={}
for k,(name,m,smi) in mset.items():
    r={"name":name,"smiles":smi}
    r["physchem"]=physchem(m)
    r["alerts"]=alerts(m)
    r["functional_groups"]=fg_inventory(m)
    r["ECFP4_Tanimoto_vs_parent"]=round(DataStructs.TanimotoSimilarity(p_ecfp,morgan(m)),4)
    r["MACCS_Tanimoto_vs_parent"]=round(DataStructs.TanimotoSimilarity(p_maccs,MACCSkeys.GenMACCSKeys(m)),4)
    results[k]=r

# ---------- 5. 3D embedding + shape/pharmacophore read-across ----------
def embed3d(m, seed=0xC0FFEE, n=1):
    mh=Chem.AddHs(m)
    params=AllChem.ETKDGv3(); params.randomSeed=seed; params.useSmallRingTorsions=True
    cids=AllChem.EmbedMultipleConfs(mh, numConfs=n, params=params)
    ok=[]
    for cid in cids:
        try:
            ff=AllChem.MMFFGetMoleculeForceField(mh, AllChem.MMFFGetMoleculeProperties(mh), confId=cid)
            if ff: ff.Minimize(maxIts=200); ok.append(cid)
        except Exception:
            ok.append(cid)
    return mh, (cids[0] if len(cids) else -1)

def descr3d(mh, cid):
    d={}
    try:
        d["RadiusOfGyration"]=round(rdMolDescriptors.CalcRadiusOfGyration(mh,confId=cid),3)
        d["Asphericity"]=round(rdMolDescriptors.CalcAsphericity(mh,confId=cid),3)
        d["SpherocityIndex"]=round(rdMolDescriptors.CalcSpherocityIndex(mh,confId=cid),3)
        d["Eccentricity"]=round(rdMolDescriptors.CalcEccentricity(mh,confId=cid),3)
        npr1=rdMolDescriptors.CalcNPR1(mh,confId=cid); npr2=rdMolDescriptors.CalcNPR2(mh,confId=cid)
        d["NPR1"]=round(npr1,3); d["NPR2"]=round(npr2,3)
        d["InertialShapeFactor"]=round(rdMolDescriptors.CalcInertialShapeFactor(mh,confId=cid),4)
    except Exception as e:
        d["error"]=str(e)
    return d

print("Embedding 3D conformers (this can take ~1-2 min for a 1450 Da glycopeptide)...")
emb={}
for k,(name,m,smi) in mset.items():
    mh,cid=embed3d(m)
    emb[k]=(mh,cid)
    results[k]["descr3d"]=descr3d(mh,cid)
    print(f"  {k}: conf {cid}  RoG={results[k]['descr3d'].get('RadiusOfGyration')}")

# O3A alignment (Crippen) + shape Tanimoto vs parent, on 3D confs
pmh,pcid=emb["Parent"]
for k,(name,m,smi) in mset.items():
    mh,cid=emb[k]
    o3a_score=None; shape_tani=None
    try:
        o3a=rdMolAlign.GetCrippenO3A(mh, pmh, prbCid=cid, refCid=pcid)
        o3a.Align(); o3a_score=round(o3a.Score(),2)
        # ShapeTanimoto = 1 - ShapeTanimotoDist
        shape_tani=round(1.0-rdShapeHelpers.ShapeTanimotoDist(mh,pmh,confId1=cid,confId2=pcid),3)
    except Exception as e:
        o3a_score=f"err:{e}"
    results[k]["CrippenO3A_score_vs_parent"]=o3a_score
    results[k]["ShapeTanimoto_vs_parent_(aligned)"]=shape_tani
    print(f"  {k}: O3A={o3a_score} ShapeTanimoto={shape_tani}")

with open(f"{OUT}/analysis_results.json","w") as f:
    json.dump(results,f,indent=2)
print("Saved analysis_results.json")
