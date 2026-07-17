"""
test_design.py -- reproducibility + correctness tests for the design engine.

Run:  python3 -m biosensor_pipeline.test_design
Pure, deterministic, offline.  No network, no Boltz.
"""

from __future__ import annotations
from .design import circular_permute, build_chimera, verify_chimera
from .screen import build_library
from .systems import SYSTEMS, TEM1


def test_circular_permutation_conserves_residues():
    seq = "ACDEFGHIKL"
    site = 4  # remove 'F'
    cp = circular_permute(seq, site, "GS")
    assert cp == "GHIKL" + "GS" + "ACDE", cp
    body = cp.replace("GS", "", 1)
    assert sorted(body) == sorted(seq[:site] + seq[site + 1:])
    print("[OK] circular permutation removes exactly the site residue, conserves the rest")


def test_insertion_preserves_reporter():
    reporter = "MREPORTER" * 3
    receptor = "ABCDEFGHIK"
    ch = build_chimera(
        name="t", reporter_name="R", reporter_seq=reporter,
        receptor_name="X", receptor_seq=receptor,
        insertion_index=5, permutation_site=4, gs_linker="GGS", flank_linker="G",
    )
    assert ch.sequence.startswith(reporter[:5])
    assert ch.sequence.endswith(reporter[5:])
    chk = verify_chimera(ch, reporter, receptor)
    assert chk["all_ok"], chk
    print("[OK] insertion preserves reporter head/tail and conserves receptor residues")


def test_all_library_constructions_valid():
    for key, sysm in SYSTEMS.items():
        lib = build_library(sysm)
        # insertion: <10 per paper; terminal fusion: <10 sites × 2 orientations
        cap = 10 if sysm.reporter.mode == "insertion" else 20
        assert 1 <= len(lib) < cap, f"{key}: library size {len(lib)} out of range"
        for c in lib:
            chk = verify_chimera(c, sysm.reporter.seq, sysm.receptor.seq)
            assert chk["all_ok"], (key, c.name, chk)
        print(f"[OK] {key}: {len(lib)} valid variants ({sysm.reporter.readout})")


def test_determinism():
    """Same inputs -> identical sequence, always."""
    sysm = SYSTEMS["dig"]
    a = build_library(sysm)
    b = build_library(sysm)
    assert [c.sequence for c in a] == [c.sequence for c in b]
    print("[OK] construction is deterministic (identical output across runs)")


def test_catalytic_residues_present():
    """The TEM-1 catalytic serines are where the motif says they are."""
    assert TEM1.seq[TEM1.catalytic["S70"]] == "S"
    assert TEM1.seq[TEM1.catalytic["S130"]] == "S"
    assert TEM1.seq[TEM1.catalytic["E166"]] == "E"
    print("[OK] TEM-1 catalytic residues verified by position")


def test_reporter_functional_residues():
    """Integrated reporters: functional residues match their labels."""
    from .systems import REPORTERS
    gdh = REPORTERS["PQQ-GDH"]
    for lab, idx in gdh.catalytic.items():
        assert gdh.seq[idx] == lab[0], (lab, gdh.seq[idx])
    print("[OK] PQQ-GDH functional residues (W346/T348/R406/R408) verified by position")


def test_terminal_fusion_topology():
    """NanoLuc CP-reporter + terminal-fusion library builds and verifies."""
    from .systems import get_system
    lib = build_library(get_system("dig-nluc"))
    assert len(lib) >= 2
    assert all(c.mode == "cp_reporter_terminal" for c in lib)
    assert {c.orientation for c in lib} == {"N", "C"}
    for c in lib:
        assert verify_chimera(c, get_system("dig-nluc").reporter.seq,
                              get_system("dig-nluc").receptor.seq)["all_ok"]
    print(f"[OK] NanoLuc terminal-fusion: {len(lib)} variants (N+C), all valid")


def test_split_complementation():
    """NanoBiT: both constructs carry their binder + fragment."""
    from .systems import NANOBIT, PROXIMITY_PAIRS
    from .architectures import build_split_complementation, verify_split
    p = PROXIMITY_PAIRS["rapamycin"]
    (an, aseq), (bn, bseq) = p["binder_a"], p["binder_b"]
    sp = build_split_complementation(
        analyte="rapamycin", binder_a_name=an, binder_a_seq=aseq,
        binder_b_name=bn, binder_b_seq=bseq,
        lgbit_seq=NANOBIT["LgBiT"], smbit_seq=NANOBIT["SmBiT"])
    chk = verify_split(sp, aseq, bseq, NANOBIT["LgBiT"], NANOBIT["SmBiT"])
    assert chk["all_ok"], chk
    print("[OK] NanoBiT split pair: LgBiT+SmBiT constructs carry both binders")


def test_and_gate_logic():
    """Two orthogonal ligand-gated modules -> AND gate with correct truth table."""
    from .systems import TEM1, LIGAND_GATED_MODULES
    from .architectures import build_and_gate, truth_table
    A, B = LIGAND_GATED_MODULES["TrpR"], LIGAND_GATED_MODULES["MetJ"]
    gate = build_and_gate(
        name="t", reporter_name=TEM1.name, reporter_seq=TEM1.seq,
        modules=[("TrpR", A["ligand"], A["seq"], len(A["seq"]) // 2, TEM1.insertion_sites["41"]),
                 ("MetJ", B["ligand"], B["seq"], len(B["seq"]) // 2, TEM1.insertion_sites["197"])])
    # reporter residues remain a subsequence (inserts are interspersed, not deleted)
    it = iter(gate.sequence)
    assert all(c in it for c in TEM1.seq), "reporter residues not preserved in order"
    tt = truth_table(gate)
    assert tt[0]["reporter_ON"] == 0 and tt[-1]["reporter_ON"] == 1
    assert sum(r["reporter_ON"] for r in tt) == 1, "AND gate should be ON in exactly 1 of 4 states"
    print("[OK] AND gate: reporter preserved; truth table ON only when both ligands present")


def test_campaign_orchestrator():
    """The design campaign returns a well-formed plan + scorecard, and encodes
    the optimization logic (readout→reporter, physiological range→target Kd)."""
    from .campaign import AnalyteSpec, ReceptorSpec, plan_campaign
    a = AnalyteSpec(name="X", analyte_class="small_molecule",
                    physiological_range_nM=(10.0, 1000.0), desired_readout="electrochemical")
    r = ReceptorSpec(name="rX", seq="M" + "AEKLQ" * 20, loop_sites=[20, 40, 60, 80], pocket_indices=[])
    res = plan_campaign(a, r)
    assert res["design_plan"]["reporter"] == "PQQ-GDH"          # electrochemical → GDH
    assert res["triage"]["target_Kd_nM"] == round((10.0 * 1000.0) ** 0.5, 1)  # geometric center
    assert 0.0 <= res["scorecard"]["readiness"] <= 1.0
    assert isinstance(res["scorecard"]["recommendations"], list)
    print("[OK] campaign: readout→reporter + range→target-Kd logic; scorecard well-formed")


if __name__ == "__main__":
    test_circular_permutation_conserves_residues()
    test_insertion_preserves_reporter()
    test_all_library_constructions_valid()
    test_determinism()
    test_catalytic_residues_present()
    test_reporter_functional_residues()
    test_terminal_fusion_topology()
    test_split_complementation()
    test_and_gate_logic()
    test_campaign_orchestrator()
    print("\nALL TESTS PASSED ✅")
