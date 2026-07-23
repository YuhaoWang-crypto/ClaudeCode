#!/usr/bin/env python3
"""Design an RNA aptamer library biased to the ORM1-divergent epitope (~110-118,
'MLAFDVND', mildly acidic/hydrophobic surface loop; away from AGP's 5 glycans).
Diverse scaffolds; ViennaRNA fold-filtered. Deterministic."""
import random, RNA
random.seed(30903)  # 3KQ0-themed seed

# recognition-loop motifs: mixed to read a hydrophobic/H-bonding loop (not pure charge)
MOTIFS = ["GGAUC","GAUGG","UGGAC","ACGGA","GGACU","CUGGA","AUGGC","GGUAC",
          "GACUG","UGACG","CAGGU","GGCAU","AUCGG","GGAAC","UCGGA","GAGGU"]
STEMS = [("GGCAG","CUGCC"),("GCGGC","GCCGC"),("GGGAC","GUCCC"),("CGGAG","CUCCG")]

def make(idx):
    s5,s3 = random.choice(STEMS)
    nloop = random.choice([2,3])
    loop = "".join(random.choice(MOTIFS) for _ in range(nloop))
    pad = "".join(random.choice("ACGU") for _ in range(max(0,23-len(loop))))
    seq = (s5+loop+pad+s3)[:33].ljust(33,"A")
    return seq

lib, seen = [], set()
tries=0
while len(lib)<10 and tries<6000:
    tries+=1
    s=make(len(lib))
    if s in seen: continue
    st,mfe=RNA.fold(s)
    npair=st.count("(")
    if npair>=6 and mfe<=-6.5:
        seen.add(s); lib.append((f"agp{len(lib)+1}",s,st,mfe,npair))

import json
out={}
print(f"AGP/ORM1 aptamer library ({len(lib)} candidates, 33 nt):\n")
for name,s,st,mfe,npair in lib:
    print(f"{name}: {s}")
    print(f"      {st}  MFE={mfe:.1f}  {npair}bp")
    out[name]=s
# scramble decoy (shuffle agp1) for the calibrated gate baseline
base=lib[0][1]; sc=list(base); random.shuffle(sc); scr="".join(sc)
out["scramble"]=scr
print(f"\nscramble (decoy): {scr}")
json.dump(out,open("agp_lib.json","w"),indent=1)
print("-> agp_lib.json")
