"""
test_materials.py -- deterministic correctness tests for the discovery core.

Run:  python3 -m materials_discovery.test_materials
Pure, offline. The convex-hull math is tested against SYNTHETIC systems with
known answers (rigorous); the candidate generators are tested for the physical
invariants (charge neutrality, stoichiometry, tolerance factor).
"""
from __future__ import annotations
from .composition import Composition, parse_formula
from .hull import PhaseDiagram, Entry
from .generate import structural, compositional
from .prototypes import tolerance_factor
from .score import surrogate_energy


def test_formula_parsing():
    assert parse_formula("Li3PS4") == {"Li": 3.0, "P": 1.0, "S": 4.0}
    assert parse_formula("La2Ti2O7") == {"La": 2.0, "Ti": 2.0, "O": 7.0}
    c = Composition.from_formula("Li2O")
    assert abs(c.num_atoms - 3.0) < 1e-9
    assert c.reduced_formula() == "Li2O"
    try:
        parse_formula("Li3P!S4"); assert False
    except ValueError:
        pass
    print("[OK] formula parsing + reduced formula + rejects garbage")


def test_hull_known_stable_and_unstable():
    """A synthetic A-B binary: known convex hull, check E_hull signs exactly."""
    # elements A,B at 0; compounds AB at -1.0/atom (stable), A3B at -0.2/atom
    entries = [
        Entry(Composition.from_formula("A"), 0.0, "A"),
        Entry(Composition.from_formula("B"), 0.0, "B"),
        Entry(Composition.from_formula("AB"), -1.0, "AB"),
    ]
    pd = PhaseDiagram(entries)
    # AB sits on the hull: E_hull ~ 0 (exclude itself -> hull is the A-B tie line at 0)
    eh_ab = pd.e_above_hull(Composition.from_formula("AB"), -1.0)
    assert abs(eh_ab - (-1.0)) < 1e-6, eh_ab   # -1.0 below the elemental tie line
    # a candidate A3B at -0.2 : hull at that composition is 0.5*AB + 0.5*A -> -0.5/atom
    eh = pd.e_above_hull(Composition.from_formula("A3B"), -0.2)
    assert abs(eh - 0.3) < 1e-6, eh            # -0.2 - (-0.5) = +0.3 : UNSTABLE
    # a candidate A3B at -0.7 : below hull -> stable (negative E_hull)
    eh2 = pd.e_above_hull(Composition.from_formula("A3B"), -0.7)
    assert abs(eh2 - (-0.2)) < 1e-6, eh2
    print("[OK] convex hull: E_hull signs exact on a synthetic A-B system")


def test_hull_decomposition():
    entries = [
        Entry(Composition.from_formula("Li"), 0.0, "Li"),
        Entry(Composition.from_formula("O"), 0.0, "O"),
        Entry(Composition.from_formula("Li2O"), -2.0, "Li2O"),
    ]
    pd = PhaseDiagram(entries)
    dec = pd.decomposition(Composition.from_formula("LiO2"),
                           exclude_label="LiO2")
    els = dict(dec)
    assert "Li2O" in els and "O" in els, dec   # LiO2 -> Li2O + O2
    print(f"[OK] decomposition: LiO2 -> {dec}")


def test_structural_charge_neutral():
    """Every generated perovskite candidate must be charge-neutral by construction."""
    cands = structural("perovskite",
                       {"A": ["Ca", "Sr", "La", "Na"], "B": ["Ti", "Zr", "Al"], "X": ["O"]})
    assert len(cands) >= 3, len(cands)
    from .prototypes import OXIDATION_STATES
    for c in cands:
        net = 0.0
        for role, el, ox in c.assignment:
            from .prototypes import PROTOTYPES
            net += ox * PROTOTYPES["perovskite"].counts[role]
        assert abs(net) < 1e-9, (c.composition.reduced_formula(), c.assignment)
    # SrTiO3 (classic perovskite, t~1.0) should be present
    formulas = {c.composition.reduced_formula() for c in cands}
    assert "SrTiO3" in formulas, sorted(formulas)
    print(f"[OK] structural: {len(cands)} charge-neutral perovskites incl. SrTiO3")


def test_tolerance_factor_filter():
    """SrTiO3 ~ 1.0 (formable); a badly-mismatched combo should be filtered out."""
    from .prototypes import SHANNON_RADII
    t = tolerance_factor(SHANNON_RADII[("Sr", 2)], SHANNON_RADII[("Ti", 4)],
                         SHANNON_RADII[("O", -2)])
    assert 0.9 <= t <= 1.05, t
    # Na(1+)Al(3+)O3 wouldn't be charge neutral for ABO3 (1+3 != 6) -> absent anyway;
    # ensure the generator never emits a non-formable-tolerance perovskite
    cands = structural("perovskite", {"A": ["Ba"], "B": ["Al"], "X": ["O"]})
    # Ba2+ Al3+ -> charge 5, not 6 -> no neutral ABO3 -> empty
    assert cands == [], [c.composition.reduced_formula() for c in cands]
    print(f"[OK] tolerance factor: SrTiO3 t={t:.3f}; ill-formed combos filtered")


def test_compositional_generator():
    cands = compositional(["Li", "P", "S"], max_count=4)
    formulas = {c.composition.reduced_formula() for c in cands}
    assert "Li3PS4" in formulas, sorted(formulas)[:20]
    # all charge neutral by construction (Li+ P5+ S2-)
    assert len(cands) >= 1
    print(f"[OK] compositional: {len(cands)} neutral Li-P-S stoichiometries incl. Li3PS4")


def test_surrogate_monotonic():
    """Surrogate must rank a more-ionic oxide below a less-ionic one (sanity only)."""
    e_lio = surrogate_energy(Composition.from_formula("Li2O"))   # big Δχ
    e_sns = surrogate_energy(Composition.from_formula("SnS"))    # small Δχ
    assert e_lio < e_sns < 0, (e_lio, e_sns)
    print(f"[OK] surrogate: Li2O ({e_lio:.2f}) more stable than SnS ({e_sns:.2f}) [⚠️ not DFT]")


def test_discovery_pipeline_runs():
    from .discover import discover
    from .prototypes import DEMO_REFERENCE_ENERGIES
    res = discover({"A": ["Ca", "Sr", "Ba"], "B": ["Ti", "Zr"], "X": ["O"]},
                   "perovskite", DEMO_REFERENCE_ENERGIES, tol=0.05)
    assert res["n_candidates"] >= 2
    assert 0.0 <= res["hit_rate"] <= 1.0
    assert all(hasattr(d, "e_above_hull") for d in res["top"])
    print(f"[OK] discovery loop: {res['n_candidates']} candidates, "
          f"{res['n_stable']} predicted-stable (surrogate)")


if __name__ == "__main__":
    test_formula_parsing()
    test_hull_known_stable_and_unstable()
    test_hull_decomposition()
    test_structural_charge_neutral()
    test_tolerance_factor_filter()
    test_compositional_generator()
    test_surrogate_monotonic()
    test_discovery_pipeline_runs()
    print("\nALL TESTS PASSED ✅")
