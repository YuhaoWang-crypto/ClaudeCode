#!/usr/bin/env python3
"""Report figures. Validated categorical palette; no dual-axis anywhere."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

OUT = "figures"
os.makedirs(OUT, exist_ok=True)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8982"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
SITE_COLOR = {"CDRH3": S1, "CH3-AB": S2, "CH3-CD": S3, "CH3-EF": S4}
SITE_MARKER = {"CDRH3": "o", "CH3-AB": "s", "CH3-CD": "^", "CH3-EF": "D"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE, "text.color": INK,
    "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": "#d8d7d0", "axes.linewidth": 0.8,
    "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.spines.top": False,
    "axes.spines.right": False, "grid.color": "#ebeae4", "grid.linewidth": 0.7,
    "legend.frameon": False, "figure.dpi": 150,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"  wrote {OUT}/{name}")


# ============================================================ Fig 1 baseline 4-Delta
base = json.load(open("results/01_baseline_4delta.json"))
pr = base["per_residue"]
x = [r["i"] for r in pr]
helix = [r["helix_frac"] for r in pr]
rmsf = [r["rmsf_corealigned_A"] for r in pr]
labels = [f"{r['aa']}{r['cga']}" for r in pr]
in4d = [r["in_4delta"] for r in pr]
rgd = [i for i, r in enumerate(pr) if r["cga"] in (43, 44, 45)]

fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.5, 5.2), sharex=True,
                             gridspec_kw={"hspace": 0.18})
cols = [S1 if ok else MUTED for ok in in4d]
a1.bar(x, helix, color=cols, width=0.72)
a1.axvspan(x[rgd[0]] - 0.5, x[rgd[-1]] + 0.5, color=S2, alpha=0.13, zorder=0)
a1.text(x[rgd[1]], 1.07, "RGD", ha="center", color=S2, fontweight="bold", fontsize=9)
a1.set_ylim(0, 1.18)
a1.set_ylabel("helix fraction\n(15 NMR models)")
a1.set_title("a   4-Delta reference conformation from PDB 6R2X: a mobile RGD on a rigid helix")
a1.grid(axis="y", zorder=0)
a1.set_axisbelow(True)
a1.annotate("fully helical core\nCgA E46–N56", xy=(9.6, 1.0), xytext=(1.0, 0.70),
            color=INK2, fontsize=8.5, va="center",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8,
                            connectionstyle="arc3,rad=-0.15"))

a2.plot(x, rmsf, color=S2, lw=2, marker="o", ms=4.5, mfc=SURFACE, mew=1.4)
a2.axvspan(x[rgd[0]] - 0.5, x[rgd[-1]] + 0.5, color=S2, alpha=0.13, zorder=0)
a2.set_ylabel("backbone RMSF (Å)\nhelix-core aligned")
a2.set_xlabel("residue (4-Delta position / chromogranin A number)")
a2.set_xticks(x)
a2.set_xticklabels(labels, rotation=90, fontsize=7)
a2.grid(axis="y", zorder=0)
a2.set_axisbelow(True)
a2.set_title("b   the RGD end stays mobile while the helix is rigid", fontsize=10)
for i in rgd:
    a2.annotate(f"{rmsf[i]:.1f}", (x[i], rmsf[i]), textcoords="offset points",
                xytext=(0, 7), ha="center", fontsize=7.5, color=S2, fontweight="bold")
fig.text(0.008, -0.03, "grey bars = residues outside 4-Delta (template extends to CgA D63)",
         fontsize=7.5, color=MUTED)
save(fig, "fig1_baseline_4delta.png")

# ============================================================ Fig 2 span reconciliation
sites = json.load(open("results/02_graft_sites.json"))
span = sites["span_reconciliation"]
order = ["CDRH3", "CH3-AB", "CH3-CD", "CH3-EF"]
vals = [span["sites"][s]["anchor_CA_A"] for s in order]
PEP = span["fourdelta_end_to_end_A"]

fig, ax = plt.subplots(figsize=(8.6, 2.9))
y = np.arange(len(order))
ax.barh(y, vals, color=[SITE_COLOR[s] for s in order], height=0.58)
for i, v in enumerate(vals):
    ax.text(v + 0.6, i, f"{v:.1f} Å", va="center", fontsize=9,
            color=INK, fontweight="bold")
    ax.text(PEP - 0.8, i, f"needs +{PEP - v:.1f} Å", va="center", ha="right",
            fontsize=8, color=INK2)
ax.axvline(PEP, color=INK, lw=1.6, ls=(0, (4, 2)))
ax.text(PEP + 0.5, len(order) - 0.35, f"4-Delta bioactive\nend-to-end {PEP:.1f} Å",
        fontsize=8.5, color=INK, va="top", fontweight="bold")
ax.set_yticks(y)
ax.set_yticklabels(order)
ax.set_xlim(0, 44)
ax.set_xlabel("distance (Å)")
ax.set_title("Every requested site is a closed loop far too short to span 4-Delta")
ax.grid(axis="x", zorder=0)
ax.set_axisbelow(True)
ax.invert_yaxis()
save(fig, "fig2_span_reconciliation.png")

# ============================================================ Fig 3 exit geometry
try:
    xv = json.load(open("results/07_exit_vector.json"))
    stj = json.load(open("results/06_integrin_steric.json"))["results"]
    cand = json.load(open("results/03_candidates.json"))
    meta = {c["id"]: c for c in cand["candidates"] + cand["controls"]}
    cone = xv["escape_cone"]["free_direction_fraction"] * 100
    nat = xv.get("natural_exit_vectors", {})
    nat_ang = min([v["angle_to_nearest_open_direction_deg"] for v in nat.values()],
                  default=1.8)

    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    seen = set()
    pts = []
    for cid, d in xv["candidates"].items():
        if cid not in stj or cid not in meta:
            continue
        site = meta[cid].get("site")
        if site not in SITE_COLOR:
            continue
        ang = d["angle_to_open_deg"]["min"]
        fit = stj[cid]["motif_fit_rmsd_A"]["min"]
        ax.scatter(ang, fit, s=120, color=SITE_COLOR[site],
                   marker=SITE_MARKER[site], edgecolor=SURFACE, linewidth=1.6,
                   zorder=3, label=site if site not in seen else None)
        seen.add(site)
        pts.append((ang, fit, cid.split("_")[0]))
    # fan out labels for points that sit close together
    pts.sort(key=lambda p: (p[0], p[1]))
    placed = []
    for ang, fit, lab in pts:
        dy, dx, ha = -14, 0, "center"
        for pa, pf in placed:
            if abs(pa - ang) < 4.0 and abs(pf - fit) < 0.055:
                dy, dx, ha = 6, 9, "left"
                break
        ax.annotate(lab, (ang, fit), textcoords="offset points",
                    xytext=(dx, dy), ha=ha, fontsize=7, color=INK2)
        placed.append((ang, fit))
    ax.margins(x=0.09, y=0.16)
    ax.axhline(2.08, color=INK, lw=1.3, ls=(0, (4, 2)), zorder=1)
    ax.text(ax.get_xlim()[1], 2.095, "free peptide (6R2X) vs bound pose = 2.08 Å  ",
            fontsize=8, color=INK, va="bottom", ha="right", fontweight="bold")
    ax.axvline(nat_ang, color=S3, lw=1.3, ls=(0, (2, 2)), zorder=1)
    ax.text(nat_ang + 1.2, ax.get_ylim()[0], f" natural pro-TGF-β1 exit = {nat_ang}°",
            fontsize=8, color=S3, va="bottom", fontweight="bold")
    ax.set_xlabel("exit-vector angle off the nearest open direction (°, best of 5 samples)")
    ax.set_ylabel("RGDLxxL motif fit to the receptor-bound pose (Å, best of 5)")
    ax.set_title(f"Only {cone:.1f}% of directions out of the αvβ6 cleft fit a domain-sized body\n"
                 f"lower-left = productive geometry", fontsize=10.5)
    ax.grid(zorder=0)
    ax.set_axisbelow(True)
    ax.legend(title="graft site", loc="upper right", fontsize=8.5, title_fontsize=8.5)
    save(fig, "fig3_exit_geometry.png")
except (FileNotFoundError, KeyError) as e:
    print(f"  skip fig3: {e}")

# ============================================================ Fig 4 + 5 ranking
try:
    rk = json.load(open("results/10_ranking.json"))
    designs = rk["designs_ranked"]

    # ---- Fig 4: metric matrix, sequential single hue
    keys = [("exit_angle_best_deg", "exit angle (°)", True),
            ("motif_fit_best_A", "motif fit (Å)", True),
            ("helix_frac_core", "helix core", False),
            ("rmsd_to_6R2X_A", "RMSD vs 6R2X (Å)", True),
            ("rgd_rel_sasa", "RGD rel. SASA", False),
            ("scaffold_rmsd_vs_wt_A", "scaffold RMSD (Å)", True),
            ("interface_ratio_vs_wt", "interface vs WT", None),
            ("graft_plddt", "graft pLDDT", False)]
    mat, raw = [], []
    for d in designs:
        row_n, row_r = [], []
        for k, _, lower_better in keys:
            v = d["metrics"].get(k)
            row_r.append(v)
            row_n.append(np.nan if v is None else float(v))
        mat.append(row_n)
        raw.append(row_r)
    mat = np.array(mat, dtype=float)
    norm = np.zeros_like(mat)
    for j, (k, _, lb) in enumerate(keys):
        col = mat[:, j]
        ok = ~np.isnan(col)
        if not ok.any():
            continue
        if lb is None:                      # closeness to 1.0 is best
            d_ = np.abs(col - 1.0)
            rng = np.nanmax(d_) - np.nanmin(d_) or 1
            norm[:, j] = 1 - (d_ - np.nanmin(d_)) / rng
        else:
            rng = np.nanmax(col) - np.nanmin(col) or 1
            z = (col - np.nanmin(col)) / rng
            norm[:, j] = 1 - z if lb else z

    fig, ax = plt.subplots(figsize=(9.6, 7.0))
    im = ax.imshow(norm, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([k[1] for k in keys], rotation=35, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(designs)))
    ax.set_yticklabels([f"{d['rank']:>2}. {d['id']}" for d in designs], fontsize=8)
    for i in range(len(designs)):
        for j in range(len(keys)):
            v = raw[i][j]
            if v is None:
                txt = "–"
            elif abs(v) >= 100:
                txt = f"{v:.0f}"
            elif abs(v) >= 10:
                txt = f"{v:.1f}"
            else:
                txt = f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.8,
                    color="#f7f7f5" if norm[i, j] > 0.58 else INK2)
    ax.set_title("Evaluation matrix, ranked (cell shading = better within column; "
                 "values printed)", fontsize=10.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.42, pad=0.02)
    cb.set_label("relative quality", fontsize=8)
    cb.outline.set_visible(False)
    ax.set_xticks(np.arange(-.5, len(keys), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(designs), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=1.6)
    ax.tick_params(which="minor", length=0)
    save(fig, "fig4_evaluation_matrix.png")

    # ---- Fig 5: ranking bars
    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    y = np.arange(len(designs))
    cols = [SITE_COLOR.get(d["site"], MUTED) for d in designs]
    ax.barh(y, [d["score_total"] for d in designs], color=cols, height=0.62)
    for i, d in enumerate(designs):
        ax.text(d["score_total"] + 0.7, i, f"{d['score_total']:.1f}", va="center",
                fontsize=8, color=INK, fontweight="bold")
        if not d["gates_passed"]:
            ax.text(1.2, i, "gate fail", va="center", fontsize=7,
                    color="#f7f7f5", fontweight="bold")
    if len(designs) > 12:
        ax.axhline(11.5, color=INK, lw=1.4, ls=(0, (4, 2)))
        ax.text(ax.get_xlim()[1], 11.5, "top-12 cut  ", ha="right", va="bottom",
                fontsize=8.5, color=INK, fontweight="bold")
    bm = next((c for c in rk["controls"] if c["id"] == "BM_CH3_Cterm_fusion"), None)
    if bm:
        ax.axvline(bm["score_total"], color=S3, lw=1.6, ls=(0, (2, 2)))
        ax.text(bm["score_total"] + 0.5, len(designs) - 0.4,
                f"C-terminal fusion\nbenchmark {bm['score_total']:.1f}",
                fontsize=8.5, color=S3, va="top", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{d['rank']:>2}. {d['id']}" for d in designs], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("total score (structural 0-100 minus programme deductions)")
    ax.set_title("Candidate ranking and the down-selected 12")
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(facecolor=SITE_COLOR[s], label=s) for s in order],
              title="graft site", loc="lower right", fontsize=8.5, title_fontsize=8.5)
    save(fig, "fig5_ranking.png")
except (FileNotFoundError, KeyError) as e:
    print(f"  skip fig4/5: {e}")

# ============================================================ Fig 6 binding-mode ladder
try:
    bm = json.load(open("results/08_binding_mode.json"))
    cr = bm["crystal_reference"]
    order6 = [("INT_LAP_methodctrl", "LAP  GRGDLGRL\n(method control)", S1),
              ("INT_4delta_free", "4Δ  RGDLRIL\n(free peptide)", S3),
              ("INT_4delta_RGE_neg", "4Δ RGE mutant\n(specificity control)", S2)]
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(10.6, 4.6),
                                 gridspec_kw={"wspace": 0.28})
    for ax, key, lab, ref in (
            (b1, "acid_to_MIDAS_metal_A", "acidic carboxylate → MIDAS metal (Å)",
             cr.get("asp_to_MIDAS_metal_A")),
            (b2, "arg_to_alphaV_Asp218_A", "RGD Arg → αv Asp218 (Å)",
             cr.get("arg_to_alphaV_Asp218_A"))):
        for i, (cid, nm, col) in enumerate(order6):
            r = bm["results"].get(cid)
            if not r:
                continue
            vals = [s[key] for s in r["per_sample"] if s.get(key) is not None]
            jitter = np.linspace(-0.13, 0.13, len(vals))
            ax.scatter([i + j for j in jitter], vals, s=85, color=col,
                       edgecolor=SURFACE, linewidth=1.4, zorder=3)
            # label the BEST (= lowest) value, next to the point it describes
            ax.annotate(f"best {min(vals):.2f}", (i + jitter[int(np.argmin(vals))],
                                                  min(vals)),
                        textcoords="offset points", xytext=(13, -3), ha="left",
                        fontsize=8, color=INK, fontweight="bold")
        if ref:
            ax.axhline(ref, color=INK, lw=1.4, ls=(0, (4, 2)), zorder=1)
            ax.text(2.45, ref, f"crystal {ref:.2f} Å ", ha="right", va="bottom",
                    fontsize=8, color=INK, fontweight="bold")
        ax.set_xticks(range(len(order6)))
        ax.set_xticklabels([n for _, n, _ in order6], fontsize=8.5)
        ax.set_ylabel(lab)
        ax.grid(axis="y", zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.5, 2.5)
    b1.set_title("a   metal coordination", fontsize=10)
    b2.set_title("b   the αv salt bridge", fontsize=10)
    fig.suptitle("Binding-mode validation: the method reproduces the co-crystal, and the "
                 "dead mutant fails\n(each dot = one of 5 predicted samples)",
                 fontsize=11, fontweight="bold", x=0.009, ha="left", y=1.06)
    fig.text(0.009, -0.07,
             "ipTM cannot separate these: 0.943 (4Δ) vs 0.944 (dead RGE mutant). "
             "Only the geometry does.", fontsize=8.5, color=INK2, style="italic")
    save(fig, "fig6_binding_mode.png")
except (FileNotFoundError, KeyError) as e:
    print(f"  skip fig6: {e}")

print("figures done")
