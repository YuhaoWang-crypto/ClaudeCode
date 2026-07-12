"""
scoring.py — PSMOS scoring with explicit provenance on every axis.

The whole point of wiring Modal/Evo2 in is to move specific axes from *curated
proxy* to *computed*. This module keeps those two classes strictly separated and
labels every value:

  provenance = "computed:uniprot"   presence / hard gate (ortholog search)
             = "computed:ensembl"   CDS availability
             = "computed:evo2"      sequence constraint (Evo2 log-likelihood)  <-- the live layer
             = "curated"            architecture / redundancy / tractability priors

Design decision (honesty over flash): Evo2's log-likelihood measures sequence
*constraint / naturalness across the tree of life*, NOT "similarity to human".
So we do NOT silently overwrite the curated fidelity axis with Evo2 and pretend
the composite is now "model-measured". Instead:

  1. The hard GATE (is the pathway even present?) becomes EMPIRICAL — driven by
     the UniProt/Ensembl ortholog search, not by assertion.
  2. A NEW computed axis, `evo2_constraint`, is reported alongside the curated
     axes. It is the headline "now live" quantity.
  3. Composite goal scores still use the transparent 6-axis model; we report how
     evo2_constraint *corroborates* the curated fidelity axis (correlation),
     rather than laundering one into the other.

Curated baseline axes are the Notch comparative-genomics priors from the pilot
dashboard (component presence, paralogue counts, divergence-time proxy,
tractability) — unchanged, just now clearly fenced off as `curated`.
"""

import json
import os
from dataclasses import dataclass, field

from .pathways import SPECIES, get_pathway

DATA = os.path.join(os.path.dirname(__file__), "data")

# --------------------------------------------------------------------------- #
# Curated baseline (Notch): the comparative-genomics priors. CURATED.
# axes: completeness, fidelity(divergence proxy), simplicity(low-redundancy),
#       arch_match, tract, thru, mean_par
# --------------------------------------------------------------------------- #
CURATED_NOTCH = {
    "human":       dict(completeness=1.00, fidelity=1.000, simplicity=0.31, arch_match=1.00, mean_par=4.14),
    "mouse":       dict(completeness=1.00, fidelity=0.805, simplicity=0.31, arch_match=1.00, mean_par=4.14),
    "rat":         dict(completeness=1.00, fidelity=0.805, simplicity=0.31, arch_match=1.00, mean_par=4.14),
    "chicken":     dict(completeness=1.00, fidelity=0.450, simplicity=0.33, arch_match=0.93, mean_par=3.86),
    "frog":        dict(completeness=1.00, fidelity=0.415, simplicity=0.33, arch_match=0.93, mean_par=3.86),
    "zebrafish":   dict(completeness=1.00, fidelity=0.342, simplicity=0.26, arch_match=0.79, mean_par=5.00),
    "tunicate":    dict(completeness=1.00, fidelity=0.185, simplicity=1.00, arch_match=0.31, mean_par=1.29),
    "fly":         dict(completeness=1.00, fidelity=0.136, simplicity=0.64, arch_match=0.48, mean_par=2.00),
    "worm":        dict(completeness=0.875, fidelity=0.136, simplicity=0.64, arch_match=0.48, mean_par=2.00),
    "sea_anemone": dict(completeness=1.00, fidelity=0.127, simplicity=1.00, arch_match=0.31, mean_par=1.29),
    "yeast":       dict(completeness=0.00, fidelity=0.063, simplicity=0.00, arch_match=0.00, mean_par=0.0),
    "plant":       dict(completeness=0.00, fidelity=0.024, simplicity=0.00, arch_match=0.00, mean_par=0.0),
}

# Transparent goal weight vectors over the 6 axes (sum to 1 each).
GOALS = {
    "G1 机制解析 (mechanism)": dict(
        note="拆解机制 → 通路完整 + 遗传冗余低(单敲除即出表型)+ 好操作。",
        w=dict(completeness=0.18, fidelity=0.10, simplicity=0.27, arch_match=0.10, tract=0.20, thru=0.15)),
    "G2 疾病建模 (disease model)": dict(
        note="模拟人类疾病 → 对人保真度高 + 架构接近人 + 哺乳在体相关。",
        w=dict(completeness=0.15, fidelity=0.34, simplicity=0.05, arch_match=0.26, tract=0.12, thru=0.08)),
    "G3 药物/表型筛选 (screening)": dict(
        note="药物/表型筛选 → 高通量 + 好操作 + 靶点保守。",
        w=dict(completeness=0.15, fidelity=0.10, simplicity=0.18, arch_match=0.07, tract=0.20, thru=0.30)),
}


@dataclass
class SpeciesScore:
    key: str
    name: str
    clade: str
    gate_ok: bool
    gate_provenance: str
    axes: dict            # {axis: value}
    provenance: dict      # {axis: provenance-string}
    evo2: dict = field(default_factory=dict)   # per-gene LL + mean + normalized
    goal_scores: dict = field(default_factory=dict)


def _load_json(path):
    return json.load(open(path)) if os.path.exists(path) else None


def compute(pathway_key="Notch"):
    pw = get_pathway(pathway_key)

    # ---- computed layer 1: presence / hard gate (UniProt) ----------------- #
    orths = _load_json(os.path.join(DATA, "orthologs", f"{pathway_key}_gates.json")) or []
    present = {}   # species -> set(families with ortholog found)
    for o in orths:
        if o["found"]:
            present.setdefault(o["species_key"], set()).add(o["family"])
    gate_families = {f.key for f in pw.gate_families}

    # ---- computed layer 2: Evo2 constraint (log-likelihood) --------------- #
    evo2_raw = _load_json(os.path.join(DATA, "evo2", f"{pathway_key}_scores.json")) or {}
    # evo2_raw: {species_key: {family: mean_LL, ...}}  (natural-log per-token LL)
    per_species_meanLL = {}
    for sk, fams in evo2_raw.items():
        vals = [v for v in fams.values() if v is not None]
        if vals:
            per_species_meanLL[sk] = sum(vals) / len(vals)
    # normalize mean-LL across species with live data -> 0..1 (higher LL = 1)
    norm = {}
    if per_species_meanLL:
        lo, hi = min(per_species_meanLL.values()), max(per_species_meanLL.values())
        for sk, ll in per_species_meanLL.items():
            norm[sk] = 0.5 if hi == lo else (ll - lo) / (hi - lo)

    searched = {o["species_key"] for o in orths}

    scores = []
    for sk, sp in SPECIES.items():
        cur = CURATED_NOTCH[sk]
        # Hard gate, done honestly. Absence of a UniProt hit is not the same as
        # a gene being absent — it can be an annotation gap. So:
        #   * all gate families found            -> PASS  (computed)
        #   * none found AND curated prior agrees -> FAIL  (confirmed absent:
        #                                            the true negative controls)
        #   * partial hits (some found, some not) -> annotation gap; trust the
        #                                            curated prior, flag it
        #   * not searched                        -> curated fallback
        n_gate = len(gate_families)
        found_gate = gate_families & present.get(sk, set())
        if sk not in searched:
            gate_ok, gate_prov = cur["completeness"] > 0, "curated"
        elif len(found_gate) == n_gate:
            gate_ok, gate_prov = True, "computed:uniprot"
        elif len(found_gate) == 0:
            if cur["completeness"] == 0:
                gate_ok, gate_prov = False, "computed:uniprot (confirmed absent)"
            else:
                gate_ok, gate_prov = cur["completeness"] > 0, "curated (no hit — likely query miss)"
        else:
            gate_ok = cur["completeness"] > 0
            gate_prov = "computed-partial (annotation gap → curated)"

        axes = dict(
            completeness=cur["completeness"],
            fidelity=cur["fidelity"],
            simplicity=cur["simplicity"],
            arch_match=cur["arch_match"],
            tract=sp.tract,
            thru=sp.thru,
        )
        prov = dict(
            completeness="curated", fidelity="curated", simplicity="curated",
            arch_match="curated", tract="curated", thru="curated",
        )
        # If we empirically confirmed both gate families, completeness of the
        # gate is computed (2/2). Keep the full 7-family completeness curated.
        if gate_prov == "computed:uniprot":
            prov["gate"] = "computed:uniprot"

        evo2 = {}
        if sk in evo2_raw:
            evo2 = dict(
                per_gene=evo2_raw[sk],
                mean_LL=round(per_species_meanLL.get(sk, float("nan")), 4),
                constraint=round(norm.get(sk, float("nan")), 3),
                provenance="computed:evo2",
            )

        # goal scores (gate zeroes everything)
        gs = {}
        for gname, gdef in GOALS.items():
            s = sum(gdef["w"][a] * axes[a] for a in gdef["w"])
            gs[gname] = round(100 * s * (1.0 if gate_ok else 0.0), 1)

        scores.append(SpeciesScore(
            key=sk, name=sp.name, clade=sp.clade,
            gate_ok=gate_ok, gate_provenance=gate_prov,
            axes={k: round(v, 3) for k, v in axes.items()},
            provenance=prov, evo2=evo2, goal_scores=gs,
        ))
    return pw, scores


def evo2_vs_curated_correlation(scores):
    """Pearson r between computed Evo2 constraint and curated fidelity, over the
    species that have live Evo2 data — reported, not laundered."""
    xs, ys = [], []
    for s in scores:
        if s.evo2 and s.gate_ok:
            xs.append(s.evo2["constraint"])
            ys.append(s.axes["fidelity"])
    if len(xs) < 3:
        return None, len(xs)
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) ** 2 for x in xs) ** 0.5
    vy = sum((y - my) ** 2 for y in ys) ** 0.5
    if vx == 0 or vy == 0:
        return None, n
    return cov / (vx * vy), n


if __name__ == "__main__":
    pw, scores = compute("Notch")
    print(f"pathway: {pw.label}  (human mean paralogs/family = {pw.human_mean_paralogs:.2f})\n")
    live = [s for s in scores if s.evo2]
    print(f"species with LIVE Evo2 constraint: {len(live)}")
    for s in sorted(scores, key=lambda x: -(x.evo2.get('mean_LL', -1e9) if x.evo2 else -1e9)):
        tag = "GATE-FAIL" if not s.gate_ok else "ok"
        e = f"Evo2 meanLL={s.evo2['mean_LL']:.4f} constraint={s.evo2['constraint']:.2f}" if s.evo2 else "— (no live Evo2)"
        print(f"  {s.key:12s} gate={tag:9s}[{s.gate_provenance:16s}] {e}")
    r, n = evo2_vs_curated_correlation(scores)
    if r is not None:
        print(f"\nEvo2 constraint vs curated fidelity: Pearson r = {r:+.2f} (n={n} live species)")
