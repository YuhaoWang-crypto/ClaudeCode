"""
scoring_psmos.py — the six-layer PSMOS model (G/D/N/R/E/X) for the Hippo–YAP
regeneration panel, with computed layers overlaid on curated priors.

    PSMOS(s,P) = wG·G + wD·D + wN·N + wR·R + wE·E + wX·X     (weights per goal)

Layers (0..1):
  G  core-gene orthology + sequence conservation   <- CURATED prior, BLENDED with
                                                       computed Evo2 constraint
                                                       (live: 5 species) and
                                                       informed by computed
                                                       Compara paralogue counts
  D  domain / active-site conservation              <- curated
  N  network topology + feedback fidelity           <- curated
  R  regulatory grammar (promoter/enhancer/TF/splice)<- AlphaGenome (human/mouse
                                                       only); curated elsewhere
  E  spatiotemporal expression match                <- curated
  X  experimental utility (tools/throughput/cost)   <- curated prior

Every computed overlay is labelled. The composite uses the Evo2-blended G (Evo2
constraint IS part of "sequence conservation" in the G definition), so the live
model data actually moves the ranking — but the blend is explicit, not hidden.

Curated baseline = the Hippo comparative-biology priors from the PSMOS
regeneration dashboard.
"""

import json
import os

from .pathways import SPECIES, get_pathway

DATA = os.path.join(os.path.dirname(__file__), "data")

# Curated Hippo baseline: layers (0..1), arch-role fit (0..1), fidelity/utility,
# clade, note. From the PSMOS regeneration dashboard.
CURATED_HIPPO = {
    "zebrafish":      dict(G=.85, D=.88, N=.90, R=.80, E=.90, X=.87, fid=.867, util=.867, clade="aqua",
        arch=dict(mechanistic=.756, imaging=.910, genetic=.858, translational=.256, comparative=.093),
        note="心脏/鳍再生由YAP驱动;透明胚胎→活体成像无敌;CRISPR、通量高。硬骨鱼WGD使旁系同源偏多。"),
    "mouse":          dict(G=.97, D=.96, N=.97, R=.90, E=.60, X=.57, fid=.889, util=.570, clade="mammal",
        arch=dict(mechanistic=.767, imaging=.360, genetic=.635, translational=.856, comparative=.0),
        note="哺乳类近人架构,肝再生等有限再生;转化外推首选,但通量/成像弱。"),
    "human":          dict(G=1., D=1., N=1., R=1., E=.50, X=.34, fid=.910, util=.341, clade="ref",
        arch=dict(mechanistic=.0, imaging=.0, genetic=.0, translational=.0, comparative=.0),
        note="参考/靶标(非模型)。人有纤维化(再生失败)背景,但活体研究再生几乎不可行。"),
    "axolotl":        dict(G=.82, D=.85, N=.88, R=.75, E=.92, X=.50, fid=.845, util=.496, clade="aqua",
        arch=dict(mechanistic=.685, imaging=.544, genetic=.487, translational=.220, comparative=.090),
        note="脊椎动物肢体/多组织再生冠军;工具改善中(大基因组),通量与成本受限。"),
    "naked_mole_rat": dict(G=.95, D=.95, N=.88, R=.80, E=.45, X=.31, fid=.818, util=.314, clade="mammal",
        arch=dict(mechanistic=.630, imaging=.300, genetic=.395, translational=.721, comparative=.079),
        note="Hippo核心保守,但高分子量透明质酸/ECM使上游机械-接触抑制输入天然重连;抗癌背景。"),
    "fly":            dict(G=.70, D=.82, N=.92, R=.65, E=.35, X=.87, fid=.700, util=.869, clade="invert",
        arch=dict(mechanistic=.900, imaging=.670, genetic=.873, translational=.192, comparative=.168),
        note="Hippo通路的发现物种;单拷贝核心、冗余最低、遗传学无敌;但无器官再生对应。"),
    "planaria":       dict(G=.60, D=.70, N=.70, R=.55, E=.95, X=.76, fid=.696, util=.757, clade="invert",
        arch=dict(mechanistic=.730, imaging=.670, genetic=.772, translational=.162, comparative=.264),
        note="全身再生(neoblast);Hippo控干细胞/体型,核心保守但布线以neoblast为中心、远离人;RNAi高通量。"),
}

# Curated Wnt/β-catenin baseline (regeneration ↔ fibrosis). Same layer scheme.
CURATED_WNT = {
    "human":          dict(G=1., D=1., N=1., R=1., E=.50, X=.34, fid=.910, util=.341, clade="ref",
        arch=dict(mechanistic=.0, imaging=.0, genetic=.0, translational=.0, comparative=.0),
        note="参考/靶标。人 IPF/器官纤维化中 Wnt/β-catenin 持续异常激活;活体再生研究几乎不可行。"),
    "mouse":          dict(G=.97, D=.95, N=.95, R=.90, E=.60, X=.60, fid=.870, util=.600, clade="mammal",
        arch=dict(mechanistic=.750, imaging=.360, genetic=.700, translational=.880, comparative=.0),
        note="肝再生 + 肺/肾/肝纤维化模型齐全;哺乳类近人架构,纤维化↔修复转化首选。"),
    "zebrafish":      dict(G=.85, D=.86, N=.88, R=.80, E=.88, X=.88, fid=.850, util=.880, clade="aqua",
        arch=dict(mechanistic=.750, imaging=.920, genetic=.860, translational=.280, comparative=.100),
        note="鳍/心脏无瘢痕再生由 Wnt 驱动;透明胚胎活体成像;CRISPR 高通量。硬骨鱼WGD冗余偏高。"),
    "axolotl":        dict(G=.80, D=.83, N=.86, R=.72, E=.90, X=.48, fid=.820, util=.480, clade="aqua",
        arch=dict(mechanistic=.680, imaging=.550, genetic=.480, translational=.220, comparative=.120),
        note="附肢/多组织再生冠军,Wnt 为再生必需;工具改善中(大基因组),通量与成本受限。"),
    "naked_mole_rat": dict(G=.94, D=.93, N=.86, R=.80, E=.45, X=.30, fid=.800, util=.300, clade="mammal",
        arch=dict(mechanistic=.600, imaging=.300, genetic=.380, translational=.700, comparative=.300),
        note="抗纤维化 / 抗癌;高分子量透明质酸-ECM 使损伤反应天然重连;哺乳类比较/对照。"),
    "fly":            dict(G=.68, D=.80, N=.90, R=.62, E=.32, X=.88, fid=.700, util=.880, clade="invert",
        arch=dict(mechanistic=.900, imaging=.660, genetic=.880, translational=.180, comparative=.200),
        note="arm/pan 是 Wnt 遗传学的发现体系;冗余低、遗传学无敌;但无器官纤维化对应。"),
    "planaria":       dict(G=.58, D=.68, N=.68, R=.52, E=.95, X=.78, fid=.680, util=.780, clade="invert",
        arch=dict(mechanistic=.750, imaging=.660, genetic=.800, translational=.150, comparative=.420),
        note="全身再生(neoblast);Wnt 定头尾极性(β-catenin 梯度),布线以极性为中心、远离人;RNAi 高通量。最佳阴性/对照。"),
    "yeast":          dict(G=.0, D=.0, N=.0, R=.0, E=.0, X=.90, fid=.05, util=.90, clade="fungus",
        arch=dict(mechanistic=.0, imaging=.0, genetic=.0, translational=.0, comparative=.0),
        note="阴性对照:真菌无 Wnt/β-catenin–TCF 通路(门控检索零命中)。"),
}

CURATED = {"Hippo": CURATED_HIPPO, "Wnt": CURATED_WNT}

WEIGHTS = {
    "default (PSMOS)":            dict(G=.18, D=.17, N=.18, R=.17, E=.15, X=.15),
    "regeneration":              dict(G=.12, D=.12, N=.20, R=.14, E=.27, X=.15),
    "fibrosis":                  dict(G=.10, D=.12, N=.25, R=.18, E=.20, X=.15),
    "drug-target":               dict(G=.15, D=.28, N=.17, R=.10, E=.10, X=.20),
    "enhancer/cis-reg":          dict(G=.10, D=.10, N=.15, R=.35, E=.15, X=.15),
    "development (spatiotemporal)": dict(G=.12, D=.12, N=.16, R=.15, E=.30, X=.15),
    "signal-dynamics":           dict(G=.10, D=.12, N=.33, R=.15, E=.15, X=.15),
}

ROLES = ["mechanistic", "imaging", "genetic", "translational", "comparative"]


def _load(path):
    return json.load(open(path)) if os.path.exists(path) else None


def _gate_ok(pw, curated):
    """Hard gate from the ortholog search, honest about annotation gaps
    (mirrors scoring.compute). Returns {species_key: (gate_ok, provenance)}."""
    orths = _load(os.path.join(DATA, "orthologs", f"{pw.key}_gates.json")) or []
    present = {}
    for o in orths:
        if o["found"]:
            present.setdefault(o["species_key"], set()).add(o["family"])
    gate_fams = {f.key for f in pw.gate_families}
    searched = {o["species_key"] for o in orths}
    out = {}
    for sk in curated:
        found = gate_fams & present.get(sk, set())
        cur_absent = sum(curated[sk][l] for l in ("G", "D", "N")) == 0
        if sk not in searched:
            out[sk] = (not cur_absent, "curated")
        elif len(found) == len(gate_fams):
            out[sk] = (True, "computed:uniprot")
        elif len(found) == 0 and cur_absent:
            out[sk] = (False, "computed:uniprot (confirmed absent)")
        elif len(found) == 0:
            out[sk] = (not cur_absent, "curated (no hit — likely query miss)")
        else:
            out[sk] = (True, "computed-partial (annotation gap → curated)")
    return out


def compute_hippo():
    return compute_psmos("Hippo")


def compute_psmos(pathway_key="Hippo"):
    pw = get_pathway(pathway_key)
    curated = CURATED[pathway_key]
    gate = _gate_ok(pw, curated)

    # computed overlay: Evo2 constraint (normalized across live species)
    evo2 = _load(os.path.join(DATA, "evo2", f"{pathway_key}_scores.json")) or {}
    meanLL = {sk: sum(v.values()) / len(v) for sk, v in evo2.items() if v}
    norm = {}
    if meanLL:
        lo, hi = min(meanLL.values()), max(meanLL.values())
        norm = {sk: (0.5 if hi == lo else (ll - lo) / (hi - lo)) for sk, ll in meanLL.items()}

    # computed overlay: Compara paralogue-based redundancy
    try:
        from .paralogs import species_simplicity
        _, compara_simp = species_simplicity(pw, use_cache=True)
    except Exception:
        compara_simp = {}

    # computed overlay: AlphaGenome R (human↔mouse regulatory concordance).
    # New format: {"R": {family: score}, "detail": {...}}. Mouse's R layer
    # becomes the mean per-family concordance; human is the reference (R=1).
    ag_raw = _load(os.path.join(DATA, "alphagenome", f"{pathway_key}_R.json")) or {}
    ag_R = ag_raw.get("R", ag_raw) if isinstance(ag_raw, dict) else {}
    ag_mouse_R = None
    if ag_R:
        vals = [v for v in ag_R.values() if isinstance(v, (int, float))]
        if vals:
            ag_mouse_R = round(sum(vals) / len(vals), 3)

    rows = []
    for sk, cur in curated.items():
        gate_ok, gate_prov = gate.get(sk, (True, "curated"))
        layers = dict(G=cur["G"], D=cur["D"], N=cur["N"], R=cur["R"], E=cur["E"], X=cur["X"])
        prov = {k: "curated" for k in layers}

        evo2_c = norm.get(sk)
        if evo2_c is not None:
            # G := blend of curated ortholog-conservation prior and computed
            # Evo2 constraint (both are "sequence conservation" in G's definition).
            layers["G"] = round(0.5 * cur["G"] + 0.5 * evo2_c, 3)
            prov["G"] = "computed:evo2+curated"

        # R: AlphaGenome gives a human↔mouse concordance. Human is the reference
        # (R stays 1.0); mouse's R layer becomes the computed concordance.
        if sk == "mouse" and ag_mouse_R is not None:
            layers["R"] = ag_mouse_R
            prov["R"] = "computed:alphagenome"

        rows.append(dict(
            key=sk, name=SPECIES[sk].name, clade=cur["clade"], note=cur["note"],
            layers=layers, prov=prov, arch=cur["arch"],
            fidelity=cur["fid"], utility=cur["util"],
            gate_ok=gate_ok, gate_prov=gate_prov,
            evo2_constraint=(round(evo2_c, 3) if evo2_c is not None else None),
            evo2_meanLL=(round(meanLL[sk], 4) if sk in meanLL else None),
            compara_simplicity=(round(compara_simp[sk], 3) if sk in compara_simp else None),
        ))

    # PSMOS score per weight profile (gate zeroes a disqualified species).
    profiles = {}
    for prof, w in WEIGHTS.items():
        scored = sorted(rows, key=lambda r: -(r["gate_ok"] * sum(w[k] * r["layers"][k] for k in w)))
        profiles[prof] = [dict(key=r["key"], name=r["name"], clade=r["clade"], gate=r["gate_ok"],
                               psmos=round(100 * r["gate_ok"] * sum(w[k] * r["layers"][k] for k in w), 1))
                          for r in scored]

    # role winners (from curated arch fit)
    winners = {}
    for role in ROLES:
        best = max(rows, key=lambda r: r["arch"][role])
        winners[role] = dict(key=best["key"], name=best["name"],
                             score=round(100 * best["arch"][role], 1))

    ag_status = ("live" if ag_R else
                 "needs-key (adapter ready; human/mouse intervals resolved)")
    return dict(pathway=pw.label, rows=rows, profiles=profiles, winners=winners,
                weights=WEIGHTS, evo2_live=sorted(norm.keys()),
                alphagenome_status=ag_status, ag_R=ag_R, ag_mouse_R=ag_mouse_R)


if __name__ == "__main__":
    D = compute_hippo()
    print(f"{D['pathway']}\nEvo2 live: {D['evo2_live']}\nAlphaGenome R: {D['alphagenome_status']}\n")
    print("regeneration ranking (PSMOS):")
    for r in D["profiles"]["regeneration"]:
        print(f"  {r['psmos']:5.1f}  {r['name']}")
    print("\nrole winners:")
    for role, w in D["winners"].items():
        print(f"  {role:14s} -> {w['name']} ({w['score']})")
    print("\ncomputed overlays (G blends Evo2; redundancy from Compara):")
    for r in D["rows"]:
        e = f"Evo2 {r['evo2_constraint']}" if r["evo2_constraint"] is not None else "Evo2 —"
        s = f"Compara-simp {r['compara_simplicity']}" if r["compara_simplicity"] is not None else "Compara —"
        print(f"  {r['key']:15s} G={r['layers']['G']:.2f}[{r['prov']['G']:22s}] {e:11s} {s}")
