"""Extract SV2A chain A (protein only) + crystal ligand from a cif; write rec.pdb, lig_xtal.pdb, box.json."""
import sys, json, numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
cif, lig, out = sys.argv[1:4]
m = MMCIFParser(QUIET=True).get_structure('x', cif)[0]
class P(Select):
    def accept_residue(self, r): return r.get_parent().id=='A' and r.id[0]==' '
class L(Select):
    def accept_residue(self, r): return r.get_resname()==lig
io=PDBIO(); io.set_structure(m); io.save(out+'_rec.pdb',P()); io.save(out+'_lig_xtal.pdb',L())
xyz=np.array([a.coord for r in m.get_residues() if r.get_resname()==lig for a in r])
json.dump(dict(center=xyz.mean(0).round(2).tolist(),n_atoms=len(xyz)),open(out+'_box.json','w'))
print(out,json.load(open(out+'_box.json')))
