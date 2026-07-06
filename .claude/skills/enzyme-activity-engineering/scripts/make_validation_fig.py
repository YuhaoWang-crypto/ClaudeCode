import numpy as np, warnings, json
from Bio.PDB import MMCIFParser, PDBParser
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.simplefilter("ignore")

def ca_cif(f):
    for m in MMCIFParser(QUIET=True).get_structure("x", f):
        for ch in m:
            c = {r.id[1]: r['CA'].coord for r in ch if 'CA' in r}
            if len(c) > 200: return c
        break
def ca_pdb(f, cid):
    d = {}
    for m in PDBParser(QUIET=True).get_structure("r", f):
        for ch in m:
            if ch.id != cid: continue
            for r in ch:
                if 'CA' in r and r.id[0] == ' ': d[r.id[1]] = r['CA'].coord
        break
    return d
def krmsd(A, B, idx):
    idx = [i for i in idx if i in A and i in B]
    P = np.array([A[i] for i in idx]); Q = np.array([B[i] for i in idx])
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc); d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1, 1, d]) @ Wt
    return np.sqrt(((Pc @ U - Qc) ** 2).sum(1).mean())

wt = ca_cif("mdqm/boltz_wt.cif"); mut = ca_cif("mdqm/boltz_5mut.cif"); ref = ca_pdb("data/8WT9.pdb", "A")
Ndom, Cdom = range(12, 126), range(150, 323)
groups = {
    "WT vs 8WT9\n(catalytic N-dom)": krmsd(wt, ref, Ndom),
    "5-mut vs 8WT9\n(catalytic N-dom)": krmsd(mut, ref, Ndom),
    "5-mut vs WT\n(catalytic N-dom)": krmsd(mut, wt, Ndom),
    "5-mut vs WT\n(C-dom, apo)": krmsd(mut, wt, Cdom),
}
local = {p: krmsd(mut, wt, range(p - 6, p + 7)) for p in [168, 193, 231, 248, 257]}
LAB = {168: "L168K", 193: "H193R", 231: "F231Y", 248: "A248K", 257: "V257A"}

INK = "#12141c"; FG = "#e8eaf0"
plt.rcParams.update({"figure.facecolor": INK, "axes.facecolor": INK, "savefig.facecolor": INK,
                     "text.color": FG, "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 5.4), gridspec_kw={"width_ratios": [1.15, 1]})

names = list(groups); vals = [groups[n] for n in names]
cols = ["#38bdf8", "#38bdf8", "#a3e635", "#6b7280"]
a1.bar(range(len(names)), vals, color=cols, edgecolor="#0a0c12")
a1.axhline(2.0, ls="--", lw=1, color="#f59e0b", alpha=.7)
a1.text(len(names)-0.5, 2.05, "2 Å", color="#f59e0b", fontsize=9, ha="right")
for i, v in enumerate(vals): a1.text(i, v+0.05, f"{v:.2f}", ha="center", fontsize=10, color=FG)
a1.set_xticks(range(len(names))); a1.set_xticklabels(names, fontsize=9)
a1.set_ylabel("Cα RMSD (Å)", fontsize=10)
a1.set_title("Fold self-consistency (Boltz-2)\ncatalytic domain preserved; mutant ≈ WT", fontsize=11)

ys = np.arange(len(local))[::-1]
lv = [local[p] for p in LAB]
lc = ["#a3e635" if v < 1 else "#f59e0b" for v in lv]
a2.barh(ys, lv, color=lc, edgecolor="#0a0c12")
for y, p in zip(ys, LAB): a2.text(local[p]+0.05, y, f"{local[p]:.2f} Å", va="center", fontsize=9, color=FG)
a2.set_yticks(ys); a2.set_yticklabels([LAB[p] for p in LAB], fontsize=10, fontfamily="monospace")
a2.set_xlabel("local backbone RMSD, 5-mut vs WT (±6 res)", fontsize=10)
a2.set_title("Local perturbation per mutation\nonly the DNA/RNA-interface loops move", fontsize=11)
for ax in (a1, a2):
    for sp in ax.spines.values(): sp.set_color("#333")
fig.suptitle("IS621 — physics-based validation: fold self-consistency (Boltz-2.1, apo monomer)", fontsize=13, y=1.0)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("figures/fig5_validation.png", dpi=150, bbox_inches="tight")
print("wrote figures/fig5_validation.png")
json.dump({"domain_rmsd": groups, "local_rmsd": local}, open("mdqm/selfconsist_summary.json", "w"))
print(json.dumps({k: round(v,2) for k,v in groups.items()}, indent=0))
