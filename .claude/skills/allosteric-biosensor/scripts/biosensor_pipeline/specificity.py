"""
specificity.py -- design receptors that DISCRIMINATE closely-related metabolites.

Problem: vitamin-D metabolites differ by one/two hydroxyls (D3 -> 25(OH)D3 ->
1,25(OH)2D3 -> 24,25(OH)2D3). A useful sensor must bind ITS target and REJECT
the others. That is a negative-design / counter-selection problem, not a
single-affinity problem.

Strategy (two routes + one shared scoring core):
  A. PDB route -- mine proteins already selective for a metabolite (discover.py);
     e.g. VDR nuclear receptor is naturally selective for 1,25(OH)2D3; DBP for
     25(OH)D3; the CDL designed binders for D3.
  B. Generative route -- if no selective binder exists, design one against the
     target with RFdiffusion(AA)/Boltz protein design, then COUNTER-SELECT it in
     silico against the off-target metabolites (Boltz protein-design + this
     specificity filter). See generative_design_plan().

  Shared core -- CROSS-PANEL CO-FOLDING: co-fold every candidate receptor against
  EVERY metabolite (Boltz), build a receptor x metabolite matrix of interface
  confidence / predicted affinity, and rank receptors by a SPECIFICITY score
  (on-target minus the best off-target). This is the number that matters.

Honesty: the specificity MATRIX is built from Boltz metrics (✅ numbers) but
discrimination of a single -OH is at the edge of what co-folding resolves;
treat the ranking as a ⚠️ prioritization and confirm selectivity with a
competition assay (the wet-lab ground truth).
"""

from __future__ import annotations
from dataclasses import dataclass, field


# ---------------------------------------------------------------------------
# the metabolite panel (SMILES from PubChem; discriminating feature annotated)
# ---------------------------------------------------------------------------
METABOLITE_PANEL = {
    "D3":        {"name": "cholecalciferol",  "smiles": "CC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C",
                  "feature": "3-OH only (A-ring)"},
    "25OHD3":    {"name": "calcifediol (25(OH)D3)", "smiles": "CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C",
                  "feature": "+25-OH (side-chain terminus)"},
    "1,25OH2D3": {"name": "calcitriol (1,25(OH)2D3)", "smiles": "CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CC(C3=C)O)O)C",
                  "feature": "+1a-OH (A-ring) AND +25-OH — the active hormone"},
    "24,25OH2D3":{"name": "24,25(OH)2D3", "smiles": "CC(CCC(C(C)(C)O)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C",
                  "feature": "+24-OH and +25-OH"},
    "D2":        {"name": "ergocalciferol (D2)", "smiles": "CC(C)C(C)C=CC(C)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C",
                  "feature": "different side chain (24-methyl, 22,23-ene)"},
}

# natural/selective receptor templates found by discover.py (PDB route)
SELECTIVE_TEMPLATES = {
    "1,25OH2D3": [("VDR nuclear receptor LBD", "1DB1", 259),
                  ("VDR", "7QPI", 287)],
    "25OHD3":    [("vitamin-D binding protein DBP/GC", "1J78", 458),  # >max-len; the natural 25OH carrier
                  ("designed CDL binder", "5IEN", 137)],
    "D3":        [("CDL2.2 designed binder", "5IEN", 137),
                  ("beta-lactoglobulin", "2GJ5", 162)],
}


# ---------------------------------------------------------------------------
# specificity scoring core
# ---------------------------------------------------------------------------
@dataclass
class SpecificityRow:
    receptor: str
    target: str
    affinity: dict                      # metabolite_key -> score (ligand_iptm or binding_prob)
    on_target: float = 0.0
    best_off_target: float = 0.0
    best_off_name: str = ""
    specificity: float = 0.0            # on_target - best_off_target


def specificity_score(receptor: str, target: str, affinity_by_metabolite: dict) -> SpecificityRow:
    """on-target minus the strongest off-target = the discrimination margin."""
    on = affinity_by_metabolite.get(target, 0.0)
    offs = {k: v for k, v in affinity_by_metabolite.items() if k != target}
    if offs:
        bo_name = max(offs, key=offs.get)
        bo = offs[bo_name]
    else:
        bo_name, bo = "", 0.0
    return SpecificityRow(receptor=receptor, target=target, affinity=affinity_by_metabolite,
                          on_target=round(on, 3), best_off_target=round(bo, 3),
                          best_off_name=bo_name, specificity=round(on - bo, 3))


def specificity_matrix(rows: list) -> dict:
    """Rank candidate receptors by discrimination margin; flag the winners."""
    ranked = sorted(rows, key=lambda r: r.specificity, reverse=True)
    return {
        "n_receptors": len(rows),
        "ranking": [{"receptor": r.receptor, "target": r.target,
                     "on_target": r.on_target, "best_off": f"{r.best_off_name}={r.best_off_target}",
                     "specificity_margin": r.specificity,
                     "verdict": "selective" if r.specificity >= 0.15 else
                                "weakly selective" if r.specificity >= 0.05 else "non-selective"}
                    for r in ranked],
        "label": "⚠️ margins from co-folding interface confidence; resolving a single "
                 "-OH is at the edge of the method — confirm with a competition assay.",
    }


def boltz_specificity_plan(receptors_smiles_pairs) -> dict:
    """Emit the Boltz co-folding jobs for a full receptor x metabolite matrix.

    receptors_smiles_pairs: [(receptor_name, receptor_seq), ...].
    One holo prediction per (receptor, metabolite); read ligand_iptm /
    binding_probability into specificity_score.
    """
    jobs = []
    for rname, _ in receptors_smiles_pairs:
        for mkey, m in METABOLITE_PANEL.items():
            jobs.append({"receptor": rname, "metabolite": mkey, "smiles": m["smiles"],
                         "idempotency_key": f"spec-{rname}-{mkey}"})
    return {"n_jobs": len(jobs), "cost_usd_estimate": round(0.05 * len(jobs), 2),
            "jobs": jobs,
            "note": "submit each as a Boltz holo (protein + metabolite SMILES, "
                    "ligand_protein_binding); fill specificity_score() with ligand_iptm."}


# ---------------------------------------------------------------------------
# generative (de-novo) route with counter-selection
# ---------------------------------------------------------------------------
def generative_design_plan(target: str, off_targets: list) -> dict:
    """Plan a de-novo selective binder when no PDB template is selective enough."""
    t = METABOLITE_PANEL[target]
    return {
        "target": target, "target_feature": t["feature"],
        "off_targets": {o: METABOLITE_PANEL[o]["feature"] for o in off_targets},
        "steps": [
            "1. SCAFFOLD: RFdiffusion-AllAtom (RFdiffusionAA) or Boltz protein-design "
            "against the target metabolite as the ligand; generate ~10^2-10^3 backbones "
            "with a pocket complementary to the DISCRIMINATING group "
            f"({t['feature']}).",
            "2. SEQUENCE: LigandMPNN (ligand-aware) to design pocket sequences that "
            "H-bond the distinguishing -OH; avoid making room for the off-target groups.",
            "3. FILTER (positive): co-fold designs vs the TARGET (Boltz/AF3) — keep high "
            "interface confidence + correct pocket contact to the discriminating -OH.",
            "4. COUNTER-SELECT (negative): co-fold survivors vs EVERY off-target; KEEP only "
            "those with a large specificity margin (specificity_score) — this is the step "
            "that buys discrimination, and the one single-target design skips.",
            "5. Feed the selective binder into the biosensor pipeline (chimera -> reporter).",
        ],
        "tools": ["RFdiffusionAA / Boltz protein-design (MCP: boltz_start_protein_design)",
                  "LigandMPNN", "Boltz-2 co-fold (this package's boltz_io + specificity core)"],
        "label": "⚠️ de-novo selective small-molecule binders are hard; expect to design "
                 "and screen many. The counter-selection filter is what maximizes the odds.",
    }


def main():
    import json, os
    OUT = os.path.join(os.getcwd(), "biosensor_out")
    print("=== Vitamin-D metabolite discrimination — design framework ===\n")
    print("Panel (what must be told apart):")
    for k, m in METABOLITE_PANEL.items():
        print(f"  {k:11s} {m['name']:26s} {m['feature']}")
    print("\nPDB route — naturally-selective templates (discover.py):")
    for tgt, tmpls in SELECTIVE_TEMPLATES.items():
        print(f"  {tgt:11s} <- " + "; ".join(f"{n}({p})" for n, p, _ in tmpls))
    # illustrative matrix (⚠️) — real values come from the Boltz co-fold run
    illus = {
        "VDR-LBD(1DB1)": ("1,25OH2D3", {"D3": .45, "25OHD3": .55, "1,25OH2D3": .90, "24,25OH2D3": .50, "D2": .40}),
        "CDL2.2(5IEN)": ("D3", {"D3": .85, "25OHD3": .80, "1,25OH2D3": .60, "24,25OH2D3": .62, "D2": .78}),
        "anti-25OH(design)": ("25OHD3", {"D3": .50, "25OHD3": .82, "1,25OH2D3": .55, "24,25OH2D3": .58, "D2": .48}),
    }
    rows = [specificity_score(r, t, a) for r, (t, a) in illus.items()]
    mat = specificity_matrix(rows)
    print("\nSpecificity ranking (ILLUSTRATIVE ⚠️ — run the Boltz matrix for real numbers):")
    for r in mat["ranking"]:
        print(f"  {r['receptor']:20s} target={r['target']:10s} margin={r['specificity_margin']:+.2f} "
              f"[{r['verdict']}]  (best off {r['best_off']})")
    plan = boltz_specificity_plan([(r, "") for r in illus])
    print(f"\nBoltz co-fold plan for the REAL matrix: {plan['n_jobs']} jobs (~${plan['cost_usd_estimate']})")
    gp = generative_design_plan("25OHD3", ["D3", "1,25OH2D3"])
    print(f"Generative fallback for 25OHD3: counter-select vs {list(gp['off_targets'])} "
          f"(steps in specificity.generative_design_plan)")
    os.makedirs(OUT, exist_ok=True)
    json.dump({"panel": METABOLITE_PANEL, "templates": SELECTIVE_TEMPLATES,
               "illustrative_ranking": mat, "boltz_plan": plan, "generative_plan": gp},
              open(os.path.join(OUT, "specificity_vitd.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/specificity_vitd.json")


if __name__ == "__main__":
    main()
