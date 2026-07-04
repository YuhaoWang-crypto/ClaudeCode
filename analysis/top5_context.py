import json, numpy as np
from collections import defaultdict

NATIVE = ("MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEP"
          "VAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRH"
          "QALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKT"
          "SAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
          "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA")
TOP5 = [168, 193, 231, 248, 257]
CATALYTIC = {11, 18, 20, 60, 105}

atoms = []
for line in open("data/8WT9.pdb"):
    if line.startswith(("ATOM", "HETATM")):
        atoms.append(dict(ch=line[21], resn=line[17:20].strip(),
                          resi=int(line[22:26]), an=line[12:16].strip(),
                          xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                          rec=line[:6].strip()))

RNA, DNA = set("EF"), set("GHIJ")
Aatoms = defaultdict(list)
for a in atoms:
    if a["ch"] == "A" and a["rec"] == "ATOM":
        Aatoms[a["resi"]].append(a)

def near(resi, pool_pred, cutoff=5.0):
    mine = Aatoms[resi]
    hits = []
    for a in atoms:
        if not pool_pred(a):
            continue
        dmin = min(np.linalg.norm(a["xyz"] - m["xyz"]) for m in mine)
        if dmin < cutoff:
            hits.append((round(dmin, 2), a["ch"], a["resn"], a["resi"], a["an"]))
    hits.sort()
    return hits

for resi in TOP5:
    wt = NATIVE[resi - 1]
    print(f"\n===== residue {wt}{resi} =====")
    # nearest RNA
    rna = near(resi, lambda a: a["ch"] in RNA, 5.0)
    dna = near(resi, lambda a: a["ch"] in DNA and a["resn"] != "MG", 5.0)
    prot_other = near(resi, lambda a: a["ch"] in set("BCD"), 5.0)
    prot_self = near(resi, lambda a: a["ch"] == "A" and a["resi"] not in (resi, resi-1, resi+1), 4.5)
    print(f"  nearest RNA atoms (<5A): {rna[:4]}")
    print(f"  nearest DNA atoms (<5A): {dna[:4]}")
    print(f"  nearest other-protomer atoms (<5A): {prot_other[:3]}")
    print(f"  packing neighbors chainA (<4.5A): {[(h[2]+str(h[3]),h[0]) for h in prot_self[:6]]}")
