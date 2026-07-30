"""Render structures, MCS core-based 3D alignment, and summary figures."""
import json
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS, rdMolAlign
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
OUT="/tmp/claude-0/-home-user-ClaudeCode/0e1c2284-3690-5d68-a96a-480ebe908c1c/scratchpad"

SM = {
 "Parent (Vancomycin B)":"OC1=CC([C@@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.87 (7S,26R)-Vanco B":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](NC)CC(C)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
 "RRT 0.75 (26R)-DMe-DeLI":"OC1=CC([C@H](C(O)=O)NC2=O)=C(C3=C(O)C=CC([C@@H](NC([C@H](NC([C@H](CC(N)=O)NC([C@@H]([C@H](O)C(C=C4Cl)=CC=C4O5)NC([C@H](N)[C@H](CC)C)=O)=O)=O)C6=CC5=C(O[C@H]7[C@H](O[C@]8([H])O[C@@H](C)[C@@H](O)[C@](N)(C)C8)[C@@H](O)[C@H](O)[C@@H](CO)O7)C9=C6)=O)C(N[C@H]2[C@H](O)C%10=CC(Cl)=C(O9)C=C%10)=O)=C3)C(O)=C1",
}
mols=[Chem.MolFromSmiles(s) for s in SM.values()]
names=list(SM.keys())
for m in mols: AllChem.Compute2DCoords(m)

# grid of 2D structures
img=Draw.MolsToGridImage(mols, legends=names, molsPerRow=3, subImgSize=(560,520))
img.save(f"{OUT}/fig1_structures.png")
print("saved fig1_structures.png")

# Highlight the N-terminal residue-1 difference in RRT 0.75 vs parent
# residue-1 fragment: N-methyl-leucine (parent) vs isoleucine primary amine (075)
patt_parent = Chem.MolFromSmarts("[NX3H1][CX4]CC(C)C")     # NHMe-CH(...)-CH2CH(CH3)2 (N-Me-Leu)
patt_075    = Chem.MolFromSmarts("[NX3H2][CX4][CX4](CC)C")  # NH2-CH(...)-CH(CH2CH3)CH3 (Ile)
def draw_hl(m, patt, fname, title):
    hit = m.GetSubstructMatch(patt)
    d = rdMolDraw2D.MolDraw2DCairo(900,750)
    d.drawOptions().legendFontSize=22
    rdMolDraw2D.PrepareAndDrawMolecule(d, m, highlightAtoms=list(hit), legend=title)
    d.FinishDrawing()
    open(f"{OUT}/{fname}","wb").write(d.GetDrawingText())
    return hit
h1=draw_hl(mols[0], patt_parent, "fig2_parent_residue1.png", "Parent: N-methyl-D-leucine (residue 1) - secondary amine")
h3=draw_hl(mols[2], patt_075, "fig2_rrt075_residue1.png", "RRT 0.75: des-Me isoleucine (residue 1) - primary amine")
print("residue-1 highlight hits:", h1, h3)

# ---- MCS-based 3D core alignment (fair read-across shape metric) ----
def embed(m, seed=1):
    mh=Chem.AddHs(m)
    p=AllChem.ETKDGv3(); p.randomSeed=seed
    AllChem.EmbedMolecule(mh,p)
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=200)
    except Exception: pass
    return mh
res={}
pm=embed(mols[0])
for i in (1,2):
    dm=embed(mols[i])
    mcs=rdFMCS.FindMCS([Chem.RemoveHs(pm),Chem.RemoveHs(dm)],
                       bondCompare=rdFMCS.BondCompare.CompareOrder,
                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                       ringMatchesRingOnly=True, completeRingsOnly=True, timeout=60)
    q=Chem.MolFromSmarts(mcs.smartsString)
    pa=Chem.RemoveHs(pm).GetSubstructMatch(q)
    da=Chem.RemoveHs(dm).GetSubstructMatch(q)
    amap=list(zip(da,pa))
    rmsd=rdMolAlign.AlignMol(dm, pm, atomMap=amap)
    res[names[i]]={"MCS_numAtoms":mcs.numAtoms,"MCS_numBonds":mcs.numBonds,
                   "core_align_RMSD_Angstrom":round(rmsd,3),
                   "core_atoms_matched":len(amap)}
    print(names[i], res[names[i]])
json.dump(res, open(f"{OUT}/mcs_core_alignment.json","w"), indent=2)
print("saved mcs_core_alignment.json")
