import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

NATIVE = ("MEHELHYIGIDTAKEKLDVDVLRPDGRHRTKKFANTTKGHDELVSWLKGHKIDHAHICIEATGTYMEP"
          "VAECLYDAGYIVSVINPALGKAFAQSEGLRNKTDTVDARMLAEFCRQKRPAAWEAPHPLERALRALVVRH"
          "QALTDMHTQELNRTETAREVQRPSIDAHLLWLEAELKRLEKQIKDLTDDDPDMKHRRKLLESIPGIGEKT"
          "SAVLLAYIGLKDRFAHARQFAAFAGLTPRRYESGSSVRGASRMSKAGHVSLRRALYMPAMVATSKTEWGR"
          "AFRDRLAANGKKGKVILGAMMRKLAQVAYGVLKSGVPFDASRHNPVAA")
CATALYTIC = {11, 18, 20, 60, 105}
TOP5 = {168: "L168K", 193: "H193R", 231: "F231Y", 248: "A248K", 257: "V257A"}
MECH = {168: "core packing → fold stability / yield",
        193: "surface basic patch → stability / yield",
        231: "new H-bond to bridge-RNA U82 → guide affinity",
        248: "new contact to target/donor DNA backbone → substrate affinity",
        257: "relieve steric crowding at RNA:DNA junction → turnover"}

atoms = []
for line in open("data/8WT9.pdb"):
    if line.startswith(("ATOM", "HETATM")):
        atoms.append(dict(ch=line[21], resn=line[17:20].strip(), resi=int(line[22:26]),
                          an=line[12:16].strip(),
                          xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                          rec=line[:6].strip()))
RNA, DNA = set("EF"), set("GHIJ")
ca = {a["resi"]: a["xyz"] for a in atoms if a["ch"] == "A" and a["an"] == "CA"}
def res_atoms(resi): return [a for a in atoms if a["ch"] == "A" and a["resi"] == resi and a["rec"] == "ATOM"]

INK = "#12141c"; FG = "#e8eaf0"
plt.rcParams.update({"figure.facecolor": INK, "axes.facecolor": INK,
                     "savefig.facecolor": INK, "text.color": FG,
                     "axes.labelcolor": FG, "xtick.color": FG, "ytick.color": FG,
                     "font.family": "DejaVu Sans"})

# ---------- Figure 1: whole-complex overview ----------
fig = plt.figure(figsize=(11, 9)); ax = fig.add_subplot(111, projection="3d")
cx = np.array([ca[r] for r in sorted(ca)])
ax.plot(cx[:, 0], cx[:, 1], cx[:, 2], color="#6b7280", lw=1.6, alpha=0.8, label="IS621 protomer (Cα)")
rna = np.array([a["xyz"] for a in atoms if a["ch"] in RNA and a["an"] in ("P", "C1'")])
dna = np.array([a["xyz"] for a in atoms if a["ch"] in DNA and a["an"] in ("P", "C1'") and a["resn"] != "MG"])
ax.scatter(rna[:, 0], rna[:, 1], rna[:, 2], s=28, color="#22d3ee", alpha=0.85, label="bridge RNA")
ax.scatter(dna[:, 0], dna[:, 1], dna[:, 2], s=28, color="#f59e0b", alpha=0.85, label="target/donor DNA")
catxyz = np.array([ca[r] for r in CATALYTIC if r in ca])
ax.scatter(catxyz[:, 0], catxyz[:, 1], catxyz[:, 2], s=170, color="#ef4444",
           edgecolor="white", lw=1.2, label="catalytic DEDD (D11/D18/D20/E60/D105)", depthshade=False)
for r, lab in TOP5.items():
    p = ca[r]
    ax.scatter(*p, s=210, color="#a3e635", edgecolor="black", lw=1.4, depthshade=False)
    ax.text(p[0], p[1], p[2] + 2, lab, color="#d9f99d", fontsize=11, fontweight="bold")
ax.set_title("IS621 bridge recombinase (PDB 8WT9) — 5 engineering candidates in context",
             color=FG, fontsize=13, pad=8)
ax.legend(loc="upper left", fontsize=8.5, facecolor="#1c1f2b", edgecolor="#333", labelcolor=FG)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
ax.grid(False)
for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
    pane.pane.set_facecolor((0, 0, 0, 0)); pane.pane.set_edgecolor((1, 1, 1, 0.05))
fig.tight_layout(); fig.savefig("figures/fig1_overview.png", dpi=150); plt.close(fig)
print("wrote figures/fig1_overview.png")

# ---------- Figure 2: top-5 local environment panels ----------
fig = plt.figure(figsize=(17, 10.5))
elem_color = {"C": "#cbd5e1", "N": "#60a5fa", "O": "#f87171", "S": "#fbbf24", "P": "#fb923c"}
for k, (resi, lab) in enumerate(TOP5.items()):
    ax = fig.add_subplot(2, 3, k + 1, projection="3d")
    center = ca[resi]
    # local backbone
    loc = [r for r in sorted(ca) if np.linalg.norm(ca[r] - center) < 13]
    lc = np.array([ca[r] for r in loc])
    ax.plot(lc[:, 0], lc[:, 1], lc[:, 2], color="#4b5563", lw=1.4, alpha=0.7)
    # WT sidechain sticks
    for a in res_atoms(resi):
        c = elem_color.get(a["an"][0], "#cbd5e1")
        ax.scatter(*a["xyz"], s=70, color=c, edgecolor="black", lw=0.4, depthshade=False)
    # nearby nucleic acid atoms
    for a in atoms:
        if a["ch"] in RNA and np.linalg.norm(a["xyz"] - center) < 8:
            ax.scatter(*a["xyz"], s=34, color="#22d3ee", alpha=0.9, depthshade=False)
        elif a["ch"] in DNA and a["resn"] != "MG" and np.linalg.norm(a["xyz"] - center) < 8:
            ax.scatter(*a["xyz"], s=34, color="#f59e0b", alpha=0.9, depthshade=False)
    ax.scatter(*center, s=240, facecolor="none", edgecolor="#a3e635", lw=2.4, depthshade=False)
    ax.set_title(f"{lab}\n{MECH[resi]}", color="#d9f99d", fontsize=10.5, pad=2)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([]); ax.grid(False)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor((0, 0, 0, 0)); pane.pane.set_edgecolor((1, 1, 1, 0.05))
# legend panel
ax = fig.add_subplot(2, 3, 6); ax.axis("off")
handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markersize=11, lw=0, label=l)
           for l, c in [("carbon", "#cbd5e1"), ("nitrogen", "#60a5fa"), ("oxygen", "#f87171"),
                        ("sulfur", "#fbbf24"), ("bridge RNA", "#22d3ee"), ("target/donor DNA", "#f59e0b")]]
handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
                          markeredgecolor="#a3e635", markersize=13, lw=0, label="mutation site (Cα)"))
ax.legend(handles=handles, loc="center", fontsize=12, facecolor="#1c1f2b",
          edgecolor="#333", labelcolor=FG, title="Local environment (atoms <8 Å)",
          title_fontsize=12)
fig.suptitle("Top-5 IS621 mutations — local structural environment (real coordinates, PDB 8WT9)",
             color=FG, fontsize=14, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig("figures/fig2_top5_panels.png", dpi=150); plt.close(fig)
print("wrote figures/fig2_top5_panels.png")
