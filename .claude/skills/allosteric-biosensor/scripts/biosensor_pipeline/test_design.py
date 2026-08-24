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



def test_specificity_scoring():
    """on-target minus best off-target margin; ranking flags selective receptors."""
    from .specificity import specificity_score, specificity_matrix, METABOLITE_PANEL
    assert "1,25OH2D3" in METABOLITE_PANEL and "25OHD3" in METABOLITE_PANEL
    sel = specificity_score("VDR", "1,25OH2D3",
                            {"D3": 0.4, "25OHD3": 0.5, "1,25OH2D3": 0.9})
    non = specificity_score("promiscuous", "D3",
                            {"D3": 0.7, "25OHD3": 0.72, "1,25OH2D3": 0.71})
    assert sel.specificity > 0.15 and non.specificity < 0.05
    mat = specificity_matrix([sel, non])
    assert mat["ranking"][0]["receptor"] == "VDR"          # most selective ranks first
    print("[OK] specificity: discrimination margin + ranking (VDR selective, promiscuous not)")


def test_pocket_redesign_plan():
    """LigandMPNN-style pocket redesign: value string keeps length + designs only the
    requested windows; campaign emits a valid design job + counter-select plan."""
    from .specificity import build_design_value, pocket_redesign_plan
    seq = "ACDEFGHIKLMNPQRST"
    val = build_design_value(seq, [(2, 2), (10, 3)])   # design idx 2-3 and 10-12
    assert val == "AC" + "2" + "FGHIKL" + "3" + "QRST", val
    assert len(val.replace("2", "").replace("3", "")) == len(seq) - 5
    plan = pocket_redesign_plan("VDR", "M" + "AEKLQ" * 30, [(10, 3), (40, 3)],
                                "1,25OH2D3", ["25OHD3", "D3"], num_proteins=12)
    assert plan["n_designed_residues"] == 6
    dj = plan["design_job"]
    assert dj["binder_specification"]["type"] == "no_template"
    assert dj["target"]["entities"][0]["type"] == "ligand_smiles"
    assert "25OHD3" in plan["counter_select"]["metabolites"] and \
           "1,25OH2D3" in plan["counter_select"]["metabolites"]
    print("[OK] pocket redesign: value keeps length, designs windows; campaign well-formed")


def test_transducer_split_bundle():
    """Split-bundle + cpGFP sensor with the paper's background-suppression knobs."""
    from .transducer import (build_split_bundle_sensor, verify_split_bundle,
                             fragment_duplication, interface_weakening, truncate_terminus)
    binder = "".join("HELIX%d____" % i for i in range(5))  # 5 pseudo-helices, 50 aa
    boundaries = [0, 10, 20, 30, 40]                        # each helix starts every 10
    s = build_split_bundle_sensor("serotonin", "SRO_b", binder, boundaries,
                                  split_after_helix=2,
                                  duplicate_segment=(0, 5), weaken_positions=[1, 3])
    chk = verify_split_bundle(s)
    assert chk["all_ok"], chk
    # duplication lengthens the N-fragment; both fragments still present in constructs
    assert len(s.frag_n) > (boundaries[2])            # grew by the duplicated 5 residues
    assert "fragment_duplication" in s.knobs and "interface_weakening" in s.knobs
    # pure ops
    assert fragment_duplication("ABCDEF", (1, 3)) == "ABC" + "BC" + "DEF"
    assert interface_weakening("LLLL", [0, 2]) == "SLSL"
    assert truncate_terminus("ABCDEF", 2, "C") == "ABCD"
    print("[OK] transducer: split-bundle+cpGFP sensor + duplication/weakening/truncation knobs")


def test_transducer_induced_folding_and_metal():
    from .transducer import (build_induced_folding_sensor, metal_coordination_plan,
                             recommend_transducer)
    fr = build_induced_folding_sensor("DTG", "DTG_b", "MKACDEFGH", reporter="FRET")
    assert "mTurquoise2" in fr.sensor_seq and fr.reporter == "FRET"
    zn = metal_coordination_plan("Zn")
    assert zn["geometry"] == "tetrahedral" and zn["n_coordinating"] == 4
    assert zn["residue_types"].count("H") == 2
    rec = recommend_transducer("metal_ion", "helical_bundle")
    assert "coordination" in rec["recommended"]
    rec2 = recommend_transducer("small_molecule", "helical_bundle", "fluorescent")
    assert "split-bundle" in rec2["recommended"]
    print("[OK] transducer: induced-folding + metal-coordination plans + recommender")


def test_end_to_end_sensor_demo():
    """The full small-molecule -> binder-spec -> split+cpGFP sensor -> checklist demo
    assembles and self-verifies."""
    from .transducer import helical_bundle_binder_plan
    from .demo_sensor import run, _placeholder_bundle
    plan = helical_bundle_binder_plan("serotonin", "C1=CC2=C(C=C1O)C(=CN2)CCN",
                                      privileged_contacts=["5-OH", "amine"])
    assert plan["topology"].endswith("bundle") and len(plan["steps"]) == 5
    bundle, boundaries = _placeholder_bundle()
    assert len(boundaries) == 5 and bundle[:5] == "EIAAL"
    run()  # writes biosensor_out/demo_sensor_serotonin.json; raises if anything breaks
    print("[OK] end-to-end demo: binder spec -> split+cpGFP sensor + knobs -> checklist")


def test_25ohd3_blueprint_demo():
    """The 25(OH)D3 blueprint: 24-construct split library + gates assemble correctly."""
    from .transducer import enumerate_split_sensor_library
    from .demo_sensor_25ohd3 import run, _placeholder_6helix
    bundle, boundaries = _placeholder_6helix()
    assert len(boundaries) == 6
    lib = enumerate_split_sensor_library(bundle, boundaries, splits=[1, 3, 5],
                                         copied_repeat_aa=[0, 7], truncation_aa=[0, 3],
                                         interfaces=["WT", "weak1"])
    assert len(lib) == 24 and lib[0]["id"] == "D3S-001" and lib[-1]["id"] == "D3S-024"
    # repeat=7 lengthens the N-fragment vs repeat=0 at the same split
    s15_r0 = next(c for c in lib if c["split"] == "1+5" and c["copied_repeat_aa"] == 0
                  and c["truncation_aa"] == 0 and c["interface"] == "WT")
    s15_r7 = next(c for c in lib if c["split"] == "1+5" and c["copied_repeat_aa"] == 7
                  and c["truncation_aa"] == 0 and c["interface"] == "WT")
    assert s15_r7["frag_n_len"] == s15_r0["frag_n_len"] + 7
    run()  # writes biosensor_out/demo_sensor_25ohd3.json
    print("[OK] 25(OH)D3 blueprint: 24-construct library + A-E gates + neighbour panel")


def test_calibrated_triage():
    """Boltz-score calibration: expected hit-rate curve + within-panel ranking."""
    from .calibration import expected_hit_rate, calibrated_rank, BASELINE_HIT_RATE
    assert expected_hit_rate(1.0) == BASELINE_HIT_RATE          # no filter = baseline
    assert expected_hit_rate(0.05) >= expected_hit_rate(0.5)    # tighter keep = higher hit rate
    assert expected_hit_rate(0.3) == 38.9                        # reproduced point
    cands = [{"construct": "a", "ligand_iptm": 0.6},
             {"construct": "b", "ligand_iptm": 0.9},
             {"construct": "c", "ligand_iptm": 0.75},
             {"construct": "d", "ligand_iptm": None}]
    r = calibrated_rank(cands, "ligand_iptm")
    assert [x["construct"] for x in r[:3]] == ["b", "c", "a"]   # sorted desc by score
    assert r[0]["rank"] == 1 and r[0]["within_panel_percentile"] == 1.0
    assert r[0]["expected_hit_rate_pct"] >= r[2]["expected_hit_rate_pct"]  # top ranked >= lower
    assert r[-1]["construct"] == "d" and r[-1]["rank"] is None  # None score sorts last
    assert "affinity" in r[0]["triage_note"].lower()           # feasibility-not-affinity flag
    print("[OK] calibrated triage: enrichment curve + within-panel ranking + affinity caveat")


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
    test_specificity_scoring()
    test_pocket_redesign_plan()
    test_transducer_split_bundle()
    test_transducer_induced_folding_and_metal()
    test_end_to_end_sensor_demo()
    test_25ohd3_blueprint_demo()
    test_calibrated_triage()
    print("\nALL TESTS PASSED ✅")
