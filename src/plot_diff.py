"""
plot_diff.py
------------
Render the conformer-difference profile for a target:

    python src/plot_diff.py AdK

Top panel : per-residue Cα displacement (open->closed), raw + smoothed, with
            distal candidate regions shaded, active-site residues marked, and the
            top hinge residues flagged.
Bottom    : distance of each residue to the active-center centroid, with the
            distal cutoff line. Together they show *where* the structure responds
            and which responses are far enough from the active site to be
            engineerable allosteric handles.
"""
from __future__ import annotations
import os
import sys
import yaml
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import fetch_structures as fetch          # noqa: E402
import analyze_conformers as ana          # noqa: E402


def build(key):
    with open(os.path.join(ROOT, "targets.yaml")) as fh:
        doc = yaml.safe_load(fh)
    cfg, spec = dict(doc["defaults"]), doc["targets"][key]
    pair, offline = spec["pair"], spec.get("offline_files") or {}
    ref_path = fetch.resolve_structure(pair["ref"]["pdb"], offline.get("ref"))
    alt_path = fetch.resolve_structure(pair["alt"]["pdb"], offline.get("alt"))
    res = ana.analyze(ref_path, pair["ref"]["chain"], alt_path, pair["alt"]["chain"],
                      spec.get("active_site", []), cfg)
    return spec, cfg, res


def plot(key, spec, cfg, res, out):
    rn = res["resnums"]
    disp, disp_s = res["displacement"], res["displacement_smooth"]
    rdist = res["r_to_active"]
    aset = set(spec.get("active_site", []))

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(11, 6.2), sharex=True,
                                   gridspec_kw={"height_ratios": [2.4, 1]})

    # domains (if defined) as colored bands along the top
    domain_colors = {"CORE": "#9e9e9e", "NMP": "#1f77b4", "LID": "#d62728"}
    domains = spec.get("domains") or {}
    ymax = disp.max() * 1.12
    for dname, spans in domains.items():
        for s, e in spans:
            ax0.axvspan(s, e, ymin=0.94, ymax=1.0,
                        color=domain_colors.get(dname, "#cccccc"), alpha=0.9)
            ax0.text((s + e) / 2, ymax * 0.985, dname, ha="center", va="top",
                     fontsize=8, color="white", fontweight="bold")

    # distal candidate regions shaded
    for i, reg in enumerate(res["regions"]):
        ax0.axvspan(reg["start"], reg["end"], color="gold", alpha=0.28,
                    label="distal candidate region" if i == 0 else None)

    ax0.plot(rn, disp, lw=0.8, color="#888", alpha=0.7, label="Cα displacement (raw)")
    ax0.plot(rn, disp_s, lw=2.0, color="#222", label="smoothed")

    # active site markers
    asite_x = [r for r in rn if int(r) in aset]
    asite_y = [disp[k] for k, r in enumerate(rn) if int(r) in aset]
    ax0.scatter(asite_x, asite_y, marker="v", s=42, color="#2ca02c",
                zorder=5, label="active-site residue")

    # top hinge residues
    hs = res["hinge_score"]
    top = [k for k in hs.argsort()[::-1] if hs[k] > 0][:8]
    ax0.scatter(rn[top], disp[top], marker="*", s=150, color="crimson",
                edgecolor="k", linewidth=0.4, zorder=6, label="top hinge (clamp point)")
    for k in top[:5]:
        ax0.annotate(int(rn[k]), (rn[k], disp[k]), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=7.5, color="crimson")

    ax0.axhline(res["disp_threshold"], ls="--", lw=1, color="orange",
                label=f"responsive P{int(res['displacement_pctl'])} = {res['disp_threshold']:.1f} Å")
    ax0.set_ylim(0, ymax)
    ax0.set_ylabel("Cα displacement\nopen→closed (Å)")
    ax0.set_title(f"{key} — {spec['name']}   |   {spec['pair']['ref']['pdb']} ⇄ "
                  f"{spec['pair']['alt']['pdb']}   "
                  f"(core RMSD {res['rmsd_core']:.1f} Å, all-CA {res['rmsd_all']:.1f} Å)")
    ax0.legend(loc="upper left", fontsize=7.5, ncol=2, framealpha=0.9)

    # bottom: distance to active center
    ax1.plot(rn, rdist, lw=1.5, color="#1f77b4")
    ax1.axhline(res["distal_cutoff"], ls="--", color="red", lw=1,
                label=f"distal cutoff {res['distal_cutoff']:.0f} Å")
    ax1.fill_between(rn, res["distal_cutoff"], rdist,
                     where=rdist > res["distal_cutoff"], color="#1f77b4", alpha=0.15)
    ax1.set_ylabel("dist → active\ncenter (Å)")
    ax1.set_xlabel("residue number")
    ax1.legend(loc="upper right", fontsize=8)
    for reg in res["regions"]:
        ax1.axvspan(reg["start"], reg["end"], color="gold", alpha=0.28)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("wrote", out)


if __name__ == "__main__":
    key = sys.argv[1] if len(sys.argv) > 1 else "AdK"
    spec, cfg, res = build(key)
    out = os.path.join(ROOT, "results", f"{key}_displacement.png")
    plot(key, spec, cfg, res, out)
