"""
campaign.py -- success-oriented design orchestrator for a NEW analyte.

Encodes the optimization thinking of the biosensor literature (see
reference/optimization-playbook.md) as an explicit decision workflow: given an
analyte, it (1) triages the analyte, (2) sources/accepts a receptor and checks
the affinity match, (3) selects an architecture + reporter, (4) designs a
pocket-safe circular-permutation library, (5) plans in-silico validation,
(6) plans the site×linker optimization, and (7) emits a transparent SUCCESS
SCORECARD with specific fixes.

Deterministic steps run for real (library construction, scorecard). Boltz/MD
steps are *planned* (commands emitted) so nothing costs money without consent.

    python3 -m biosensor_pipeline.campaign            # cortisol demo (new analyte)
"""

from __future__ import annotations
import json
import math
import os
from dataclasses import dataclass, field

from .design import choose_cp_linker
from .systems import Receptor, Reporter, System, REPORTERS, TEM1, GDH, NANOLUC

OUT = os.path.join(os.getcwd(), "biosensor_out")


# ---------------------------------------------------------------------------
# specs
# ---------------------------------------------------------------------------
@dataclass
class AnalyteSpec:
    name: str
    smiles: str = None
    analyte_class: str = "small_molecule"        # small_molecule | peptide | protein
    physiological_range_nM: tuple = (10.0, 1000.0)  # operating concentration window
    desired_readout: str = "any"                 # colorimetric | electrochemical | luminescent | any
    in_cell: bool = False


@dataclass
class ReceptorSpec:
    name: str
    seq: str
    pdb: str = None
    known_Kd_nM: float = None                     # if known (from PDB paper / ITC)
    loop_sites: list = field(default_factory=list)
    pocket_indices: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# (1) analyte triage
# ---------------------------------------------------------------------------
def triage_analyte(a: AnalyteSpec) -> dict:
    lo, hi = a.physiological_range_nM
    target_kd = round(math.sqrt(lo * hi), 1)      # geometric center of the range
    # candidate architectures by class + readout
    arch = []
    if a.analyte_class == "small_molecule":
        arch = ["insertion"]
    elif a.analyte_class in ("peptide", "protein"):
        arch = ["insertion", "split_complementation"]   # PPI/sandwich possible
    return {"target_Kd_nM": target_kd, "operating_window_nM": [lo, hi],
            "candidate_architectures": arch}


# ---------------------------------------------------------------------------
# (2) affinity-match check  (Binkowski lesson)
# ---------------------------------------------------------------------------
def affinity_match(receptor_Kd_nM, target_Kd_nM, window_nM) -> dict:
    lo, hi = window_nM
    if receptor_Kd_nM is None:
        return {"status": "unknown", "note": "measure/predict Kd; aim near the "
                f"operating center ~{target_Kd_nM} nM"}
    if receptor_Kd_nM < lo / 3:
        return {"status": "too_tight", "note": "Kd well below the operating window "
                "→ SATURATION risk; use a weaker binder variant or accept a narrow range"}
    if receptor_Kd_nM > hi * 3:
        return {"status": "too_weak", "note": "Kd above the operating window → low "
                "signal; add an AUXILIARY binding domain (offset CP penalty) or use avidity"}
    return {"status": "matched", "note": "receptor affinity brackets the operating range"}


# ---------------------------------------------------------------------------
# (3) architecture + reporter selection
# ---------------------------------------------------------------------------
def select_reporter(readout: str, in_cell: bool) -> Reporter:
    if readout == "electrochemical":
        return GDH
    if readout == "luminescent" or in_cell:
        return NANOLUC
    return TEM1                                    # colorimetric default


# ---------------------------------------------------------------------------
# (4) library design (pocket-safe CP × linker grid)
# ---------------------------------------------------------------------------
def design_plan(receptor: ReceptorSpec, reporter: Reporter) -> dict:
    from .screen import pocket_safe_sites
    safe, dropped = pocket_safe_sites(receptor.loop_sites, receptor.pocket_indices)
    # auto-scaled base linker + a small grid to scan (the DR knob)
    base = choose_cp_linker()
    linker_grid = ["G", "GGG", "GGGGG", "GSGSG"]   # flank linker series
    return {
        "reporter": reporter.name, "readout": reporter.readout, "topology": reporter.mode,
        "candidate_sites": receptor.loop_sites, "pocket_safe_sites": safe,
        "pocket_dropped_sites": dropped, "cp_gs_linker": base,
        "flank_linker_grid": linker_grid,
        "core_library_size": min(len(safe) if safe else len(receptor.loop_sites), 9),
        "grid_size_if_optimizing": (len(safe) or len(receptor.loop_sites)) * len(linker_grid),
    }


# ---------------------------------------------------------------------------
# (7) success scorecard  (transparent, ⚠️ heuristic)
# ---------------------------------------------------------------------------
def scorecard(a: AnalyteSpec, r: ReceptorSpec, tri: dict, aff: dict,
              rep: Reporter, plan: dict, boltz=None) -> dict:
    axes = {}
    recs = []

    # receptor size / rigidity
    L = len(r.seq)
    axes["receptor_size"] = (0.20, 1.0 if L <= 150 else 0.6 if L <= 250 else 0.3)
    if L > 250:
        recs.append(f"receptor is large ({L} aa) — switching is harder; prefer a smaller/de-novo binder")

    # affinity match
    axes["affinity_match"] = (0.20, {"matched": 1.0, "unknown": 0.6,
                                     "too_weak": 0.5, "too_tight": 0.4}[aff["status"]])
    if aff["status"] not in ("matched", "unknown"):
        recs.append(aff["note"])

    # pocket-safe sites
    ns = len(plan["pocket_safe_sites"]) if plan["pocket_safe_sites"] else len(plan["candidate_sites"])
    axes["pocket_safe_sites"] = (0.15, 1.0 if ns >= 3 else 0.5 if ns >= 1 else 0.0)
    if ns < 3:
        recs.append("few pocket-safe CP sites — widen the loop search or redesign the binder")

    # architecture fit
    fit = (a.desired_readout in ("any", rep.readout))
    axes["architecture_fit"] = (0.15, 1.0 if fit else 0.5)
    if not fit:
        recs.append(f"requested readout {a.desired_readout} ≠ reporter {rep.readout}; switch reporter")

    # DR-tunability (always available here)
    axes["DR_tunability"] = (0.10, 1.0)

    # Boltz-derived axes (only if a validation was run)
    if boltz and "binding_retention" in boltz:
        br = boltz["binding_retention"]
        axes["binding_retention"] = (0.12, min(1.0, br))
        if br < 0.85:
            recs.append(f"binding retention {br:.2f} < 0.85 — scan sites/linker or graft elsewhere")
        axes["active_site_intact"] = (0.08, 1.0 if boltz.get("active_site_intact") else 0.3)
        if not boltz.get("active_site_intact", True):
            recs.append("active site perturbed in silico — try another insertion loop")
    else:
        recs.append("run Boltz holo/apo/control to score binding retention + active-site integrity")

    wsum = sum(w for w, _ in axes.values())
    readiness = round(sum(w * v for w, v in axes.values()) / wsum, 3)
    return {
        "readiness": readiness,
        "axes": {k: round(v, 2) for k, (_, v) in axes.items()},
        "recommendations": recs,
        "label": "⚠️ readiness is a transparent prioritization, NOT a success probability "
                 "or a dynamic range; ground truth is a wet-lab titration.",
    }


# ---------------------------------------------------------------------------
# orchestrator
# ---------------------------------------------------------------------------
def plan_campaign(analyte: AnalyteSpec, receptor: ReceptorSpec, boltz=None) -> dict:
    tri = triage_analyte(analyte)
    aff = affinity_match(receptor.known_Kd_nM, tri["target_Kd_nM"], analyte.physiological_range_nM)
    rep = select_reporter(analyte.desired_readout, analyte.in_cell)
    plan = design_plan(receptor, rep)
    card = scorecard(analyte, receptor, tri, aff, rep, plan, boltz)
    return {"analyte": analyte.name, "triage": tri, "receptor": receptor.name,
            "affinity_match": aff, "design_plan": plan, "scorecard": card,
            "next_steps": [
                f"build library: register a System(receptor={receptor.name}, reporter={rep.name}) "
                f"and run `run_repro`",
                "submit Boltz holo/apo/control (reference/boltz-validation.md) → analyze_boltz",
                "optimize the marginal variants over the site×linker grid (tune_linker)",
                "shortlist top 2-3 for synthesis + wet-lab titration",
            ]}


def _demo_cortisol() -> dict:
    """New-analyte demo: cortisol (paper's flagship), receptor mined = hcy129 (8UQF)."""
    analyte = AnalyteSpec(
        name="cortisol", smiles="CC12CCC(=O)C=C1CCC3C2C(CC4(C3CCC4(C(=O)CO)O)C)O",
        analyte_class="small_molecule",
        physiological_range_nM=(10.0, 600.0),      # free/total serum cortisol span
        desired_readout="electrochemical",          # e.g. a wearable/point-of-care sensor
        in_cell=False)
    # receptor sourced by Tier-A mining (discover.py finds 8UQF for cortisol)
    hcy = ("MSGTSSEEAKEAIIAMLKEWYDAMNEGDMEKLRSLVDPDASFVDARTNQVYDKDQFLQMIKEALEQDLKV"
           "EVKSIDIEQQPDGDVVIVKVKVRATMVRNGQEHVFEVVDTYEFRRKGDSWKIVKLVSEITQLGSG")
    receptor = ReceptorSpec(
        name="hcy129", seq=hcy, pdb="8UQF", known_Kd_nM=None,
        # loop sites: try structure-derived; fall back to a spread of interior residues
        loop_sites=_loop_sites_or_fallback("8UQF", hcy),
        pocket_indices=[])
    return plan_campaign(analyte, receptor)


def _loop_sites_or_fallback(pdb, seq):
    try:
        from .structure import annotate
        _, _, _, centers = annotate(pdb, "A", seq)
        if centers:
            return centers
    except Exception:
        pass
    # heuristic fallback: evenly spread interior residues
    n = len(seq)
    return [int(n * f) for f in (0.2, 0.35, 0.5, 0.65, 0.8)]


def main():
    res = _demo_cortisol()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "campaign_cortisol.json"), "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"=== DESIGN CAMPAIGN — new analyte: {res['analyte']} ===")
    t = res["triage"]
    print(f"[1] triage: class=small_molecule  operating window={t['operating_window_nM']} nM  "
          f"→ target Kd ~{t['target_Kd_nM']} nM  archs={t['candidate_architectures']}")
    print(f"[2] receptor: {res['receptor']} (mined 8UQF)  affinity_match={res['affinity_match']['status']}"
          f"  — {res['affinity_match']['note']}")
    p = res["design_plan"]
    print(f"[3] reporter: {p['reporter']} ({p['readout']}, {p['topology']})")
    print(f"[4] library: {p['core_library_size']} pocket-safe CP variants; "
          f"site×linker grid = {p['grid_size_if_optimizing']}  (linkers {p['flank_linker_grid']})")
    c = res["scorecard"]
    print(f"[7] SUCCESS SCORECARD  readiness={c['readiness']}")
    for k, v in c["axes"].items():
        print(f"      {k:20s} {v}")
    print("    recommendations (what to fix / do next):")
    for r in c["recommendations"]:
        print(f"      • {r}")
    print("    " + c["label"])
    print(f"\nwrote {OUT}/campaign_cortisol.json")
    return res


if __name__ == "__main__":
    main()
