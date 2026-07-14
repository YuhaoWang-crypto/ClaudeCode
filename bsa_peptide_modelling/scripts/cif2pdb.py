import sys, warnings
warnings.filterwarnings('ignore')
from Bio.PDB import MMCIFParser, PDBIO, Select
p=MMCIFParser(QUIET=True)
inp,out=sys.argv[1],sys.argv[2]
s=p.get_structure('x',inp)
io=PDBIO()
io.set_structure(s)
io.save(out)
# report chains/residues
m=list(s.get_models())[0]
for ch in m:
    resns=[r.get_resname() for r in ch]
    print('chain',ch.id,'nres',len(resns),'first3',resns[:3],'last3',resns[-3:])
