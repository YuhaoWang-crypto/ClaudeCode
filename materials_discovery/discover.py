"""
discover.py -- the GNoME active-learning discovery loop.

    generate candidates (structural + compositional)
      -> cheap prescreen: score every candidate, rank by predicted E_hull
      -> verify: (re)score the top-k with the expensive back-end (DFT)   [hook]
      -> decide: keep candidates on/near the convex hull (E_hull <= tol)
      -> (retrain the cheap model on the new labels; repeat)             [hook]

Offline, prescreen and verify both use the surrogate; a real run passes
scorer="uma" for prescreen and verifies with qe_modal_energy. The convex-hull
decision is exact either way.
"""
from __future__ import annotations
from dataclasses import dataclass
from .composition import Composition
from .hull import PhaseDiagram, Entry
from .generate import structural, compositional
from .score import SCORERS


@dataclass
class Discovery:
    formula: str
    prototype: str
    source: str
    formation_energy_per_atom: float
    e_above_hull: float
    stable: bool
    decomposes_to: list


def build_reference_pd(reference_energies: dict) -> PhaseDiagram:
    """reference_energies: {formula: formation_energy_per_atom}."""
    entries = [Entry(Composition.from_formula(f), e, label=Composition.from_formula(f).reduced_formula())
               for f, e in reference_energies.items()]
    return PhaseDiagram(entries)


def screen(candidates: list, pd: PhaseDiagram, scorer="surrogate",
           tol: float = 0.02) -> list:
    """Score candidates, compute E_hull vs the reference PD, flag stable ones."""
    fn = SCORERS[scorer]
    results = []
    for c in candidates:
        comp = c.composition
        e = fn(comp)
        eh = pd.e_above_hull(comp, e, exclude_self=True)
        results.append(Discovery(
            formula=comp.reduced_formula(), prototype=c.prototype, source=c.source,
            formation_energy_per_atom=round(e, 4), e_above_hull=round(eh, 4),
            stable=bool(eh <= tol),
            decomposes_to=pd.decomposition(comp, exclude_label=comp.reduced_formula())))
    return sorted(results, key=lambda d: d.e_above_hull)


def discover(role_candidates: dict, prototype: str, reference_energies: dict,
             scorer="surrogate", tol: float = 0.02, top_k: int = 20) -> dict:
    """One structural-pipeline discovery pass over a prototype."""
    pd = build_reference_pd(reference_energies)
    cands = structural(prototype, role_candidates)
    ranked = screen(cands, pd, scorer=scorer, tol=tol)
    stable = [d for d in ranked if d.stable]
    return {
        "prototype": prototype, "scorer": scorer, "tol": tol,
        "n_candidates": len(cands), "n_stable": len(stable),
        "hit_rate": round(len(stable) / max(1, len(cands)), 3),
        "top": ranked[:top_k],
        "stable": stable,
    }


def main():
    import json, os
    from .prototypes import DEMO_REFERENCE_ENERGIES
    OUT = os.path.join(os.getcwd(), "materials_out")
    os.makedirs(OUT, exist_ok=True)
    print("=== GNoME-style discovery demo: alkaline-earth ABO3 perovskite substitutions ===")
    print("⚠️ energies are the OFFLINE SURROGATE (ionicity heuristic), NOT DFT.\n")
    # substitute into ABO3 perovskite over the alkaline-earth titanate/zirconate space
    role_candidates = {"A": ["Ca", "Sr", "Ba"],
                       "B": ["Ti", "Zr"],
                       "X": ["O"]}
    res = discover(role_candidates, "perovskite", DEMO_REFERENCE_ENERGIES,
                   scorer="surrogate", tol=0.05, top_k=12)
    print(f"prototype=perovskite  candidates={res['n_candidates']}  "
          f"charge-neutral+tolerance-formable; predicted-stable={res['n_stable']} "
          f"(hit rate {res['hit_rate']})\n")
    print(f"{'formula':10s} {'Eform/atom':>11s} {'E_hull':>8s}  verdict   decomposes_to")
    for d in res["top"]:
        v = "STABLE" if d.stable else "-"
        print(f"{d.formula:10s} {d.formation_energy_per_atom:11.3f} {d.e_above_hull:8.3f}  "
              f"{v:8s}  {d.decomposes_to}")
    out = {k: (v if not isinstance(v, list) else [d.__dict__ for d in v])
           for k, v in res.items()}
    json.dump(out, open(os.path.join(OUT, "discovery_perovskite_demo.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/discovery_perovskite_demo.json")
    print("\nReal run: pass scorer='uma' for the prescreen and verify the top hits with "
          "qe_modal_energy (DFT on Modal) -- swap the surrogate for real energies.")


if __name__ == "__main__":
    main()
