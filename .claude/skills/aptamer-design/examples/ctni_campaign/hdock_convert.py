#!/usr/bin/env python3
import warnings; warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, PDBIO, Select
p = MMCIFParser(QUIET=True)
jobs = [("tni1","tni1_rna.pdb"),("ctni","ctni_prot.pdb"),
        ("tnni2","tnni2_prot.pdb"),("tnni1","tnni1_prot.pdb")]
for name, out in jobs:
    st = p.get_structure(name, f"{name}.cif")
    io = PDBIO(); io.set_structure(st); io.save(out)
    n = sum(1 for l in open(out) if l.startswith("ATOM"))
    print(f"{name}: {n} atoms -> {out}")
