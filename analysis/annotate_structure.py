import numpy as np, json
from collections import defaultdict

CATALYTIC = {11, 18, 20, 60, 105}  # RuvC-like DEDD active-site acidic residues (chain A)

atoms = []
for line in open("data/8WT9.pdb"):
    if line.startswith(("ATOM", "HETATM")):
        atoms.append(dict(ch=line[21], resn=line[17:20].strip(), resi=int(line[22:26]),
                          an=line[12:16].strip(),
                          xyz=[float(line[30:38]), float(line[38:46]), float(line[46:54])],
                          rec=line[:6].strip(), el=line[76:78].strip()))

RNA_CH, DNA_CH = set("EF"), set("GHIJ")
def sel(pred): return np.array([a["xyz"] for a in atoms if pred(a)]) if any(pred(a) for a in atoms) else np.empty((0, 3))

rna = sel(lambda a: a["ch"] in RNA_CH)
dna = sel(lambda a: a["ch"] in DNA_CH and a["resn"] != "MG")
other_prot = sel(lambda a: a["ch"] in set("BCD"))
cat_atoms = sel(lambda a: a["ch"] == "A" and a["resi"] in CATALYTIC)

# chain A residues
byres = defaultdict(list)
ca_coords = {}
for a in atoms:
    if a["ch"] == "A" and a["rec"] == "ATOM":
        byres[a["resi"]].append(a["xyz"])
        if a["an"] == "CA":
            ca_coords[a["resi"]] = np.array(a["xyz"])

def mindist(pts, target):
    if len(target) == 0 or len(pts) == 0:
        return None
    P = np.asarray(pts); d = np.linalg.norm(P[:, None, :] - target[None, :, :], axis=-1)
    return float(d.min())

ca_arr = np.array(list(ca_coords.values()))
ann = {}
for resi, pts in byres.items():
    P = np.array(pts)
    # burial proxy: number of CA within 10 A of this residue's CA
    burial = None
    if resi in ca_coords:
        burial = int((np.linalg.norm(ca_arr - ca_coords[resi], axis=1) < 10).sum() - 1)
    ann[resi] = {
        "d_catalytic": None if resi in CATALYTIC else mindist(pts, cat_atoms),
        "d_rna": mindist(pts, rna),
        "d_dna": mindist(pts, dna),
        "d_interface": mindist(pts, other_prot),
        "cn10": burial,  # CA contact number within 10A (higher = more buried/core)
        "is_catalytic": resi in CATALYTIC,
    }

# environment label per residue
for resi, d in ann.items():
    tags = []
    if d["is_catalytic"]: tags.append("CATALYTIC")
    if d["d_rna"] is not None and d["d_rna"] < 4.0: tags.append("RNA-contact")
    if d["d_dna"] is not None and d["d_dna"] < 4.0: tags.append("DNA-contact")
    if d["d_interface"] is not None and d["d_interface"] < 4.0: tags.append("oligomer-interface")
    if d["d_catalytic"] is not None and d["d_catalytic"] < 6.0: tags.append("active-site-shell")
    if d["cn10"] is not None:
        tags.append("core" if d["cn10"] >= 17 else ("surface" if d["cn10"] <= 11 else "intermediate"))
    d["tags"] = tags

json.dump(ann, open("data/structure_annotation.json", "w"))
# quick summary
from collections import Counter
c = Counter(t for d in ann.values() for t in d["tags"])
print("residues annotated:", len(ann))
print("tag counts:", dict(c))
print("catalytic:", sorted(CATALYTIC))
print("RNA-contact residues:", sorted(r for r, d in ann.items() if "RNA-contact" in d["tags"]))
print("DNA-contact residues:", sorted(r for r, d in ann.items() if "DNA-contact" in d["tags"]))
