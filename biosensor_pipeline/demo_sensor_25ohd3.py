"""
demo_sensor_25ohd3.py -- END-TO-END sensor design blueprint for 25-hydroxyvitamin D3
(25(OH)D3 / calcifediol), implementing the uploaded design blueprint
(25OHD3_split_cpGFP_end_to_end_demo_CN).

Discipline (from the blueprint): every un-validated sequence/structure here is a
CANDIDATE DESIGN, never a "binder" or "sensor" that has been obtained. Nothing here
is wet-lab validated. The pipeline, the 24-construct matrix, the knobs and the A-E
gates are real; a purified-protein fluorescence titration is the ground truth.

Run:  python3 -m biosensor_pipeline.demo_sensor_25ohd3
"""
from __future__ import annotations
import json, os
from .transducer import (helical_bundle_binder_plan, enumerate_split_sensor_library,
                         recommend_transducer, CPGFP)

ANALYTE = "25(OH)D3"
# isomeric SMILES (blueprint 1.1); stereochemistry MUST be preserved in real design
SMILES = r"C[C@H](CCCC(C)(C)O)[C@H]1CC[C@@H]\2[C@@]1(CCC/C2=C\C=C/3\C[C@H](CCC3=C)O)C"
POLAR_ANCHORS = ["3b-OH (A-ring)", "25-OH (side-chain terminus)"]

# 1.3 neighbour panel that must be designed against (negative design + assay)
NEIGHBORS = {
    "VitaminD3":     "no 25-OH -> tests the side-chain polar anchor",
    "25OHD2":        "D2 side chain (24-Me, 22,23-ene) -> D2/D3 discrimination",
    "1,25OH2D3":     "extra 1a-OH (A-ring) -> the active hormone",
    "24R,25OH2D3":   "extra 24-OH on the side chain",
    "3-epi-25OHD3":  "3-OH stereochemistry flipped -> hardest stereoselectivity test",
}

MEASUREMENT = {
    "measurand": "sensor-accessible 25(OH)D3 in buffer (define free vs total)",
    "conc_window": "10 nM - 10 uM; target mid-point 0.1-1 uM",
    "readout": "cpGFP fluorescence-ON, ex ~488 nm",
    "dynamic_range_gate": "round-1 dF/F0 >= 2; optimize toward >= 5",
    "background_gate": "apo baseline low, but absolute brightness >= 10x blank",
    "response_time_goal": "<=5 min to 90% (design goal, must be measured)",
    "reversibility": "signal returns on dilution / competitor displacement",
    "logP_risk": "XLogP ~6.2 -> aqueous delivery + nonspecific adsorption are first-order risks",
}

CDL2_ROLE = {
    "template": "CDL2.2 / PDB 5IEN (designed 25(OH)D3 binder, ~300 nM; 2.09 A)",
    "use_for": ["pocket geometry / hydrophobic burial reference", "binding-assay positive control",
                "local negative-design parent (6-10 pocket sites) as a fast selectivity baseline"],
    "do_not_claim": ["that it is a 6-helix bundle", "that it splits directly",
                     "that it has strong apo/holo switching", "that it already is selective"],
    "read": "5IEN holo pocket is reproducible but under-discriminates D3 / 25OHD2 / 3-epi; "
            "prove binding+selectivity FIRST, then add fluorescent transduction.",
}


def _placeholder_6helix(helix_len=20, loop="GSGG"):
    """CANDIDATE-DESIGN placeholder 6-helix bundle standing in for a ProBuilder/
    RFdiffusionAA output, so the 24-construct matrix is concrete. NOT a real binder."""
    helix = "EIAALKQEIAALKKENAALK"[:helix_len]
    seq, boundaries, cur = "", [], 0
    for i in range(6):
        boundaries.append(cur)
        seq += helix
        cur += helix_len
        if i < 5:
            seq += loop
            cur += len(loop)
    return seq, boundaries


def run():
    OUT = os.path.join(os.getcwd(), "biosensor_out")
    os.makedirs(OUT, exist_ok=True)
    print("=== 25(OH)D3 -> 6-helix bundle -> split+cpGFP : end-to-end blueprint ===")
    print("    (all sequences are CANDIDATE DESIGNS; nothing is wet-lab validated)\n")

    # 1. binder spec (6-helix bundle; anchor 3b-OH + 25-OH; negative design vs neighbours)
    plan = helical_bundle_binder_plan(ANALYTE, SMILES, n_helices=6, length=130,
                                      privileged_contacts=POLAR_ANCHORS)
    plan["affinity_targets"] = {"Kd_25OHD3": "50-500 nM (do not chase max affinity)",
                                "compute_gate": "ddG(decoy-target) >= +0.5 kcal/mol to advance",
                                "assay_selectivity": "each neighbour < 10% of target response",
                                "split_ready": ">=3 inter-helix split sites that don't cut the pocket"}
    plan["negative_design"] = {
        "VitaminD3": "focus on the 25-OH anchor: no water/side-chain rescue when it's absent",
        "3-epi-25OHD3": "focus on 3-OH direction + local steric crowding, not bulk hydrophobicity"}
    print(f"[1] BINDER SPEC (candidate) - 6-helix bundle ~130 aa; anchors {POLAR_ANCHORS}")
    print(f"     negative design vs {list(NEIGHBORS)}")
    print(f"     Kd target {plan['affinity_targets']['Kd_25OHD3']}; "
          f"advance gate {plan['affinity_targets']['compute_gate']}")

    # 2. transducer + 3. the 24-construct first-round library
    rec = recommend_transducer("small_molecule", "helical_bundle", "fluorescent")
    bundle, boundaries = _placeholder_6helix()
    lib = enumerate_split_sensor_library(
        bundle, boundaries,
        splits=[1, 3, 5],                 # 1+5, 3+3, 5+1
        copied_repeat_aa=[0, 7],
        truncation_aa=[0, 3],
        interfaces=["WT", "weak1"],
        linker="GG", id_prefix="D3S")
    assert len(lib) == 24
    print(f"\n[2] TRANSDUCER - {rec['recommended']}; single-chain "
          f"frag_N-GG-cpGFP-GG-frag_C")
    print(f"[3] FIRST-ROUND LIBRARY - {len(lib)} constructs "
          f"(3 splits x repeat[0,7] x trunc[0,3] x iface[WT,weak1]):")
    for c in lib[:4]:
        print(f"     {c['id']}: split {c['split']}, repeat {c['copied_repeat_aa']}aa, "
              f"trunc {c['truncation_aa']}aa, {c['interface']} -> construct {c['construct_len']} aa")
    print(f"     ... ({len(lib)} total; D3S-001..024)")
    plate_controls = ["cpGFP-only", "full binder-cpGFP", "pocket-dead", "no-ligand",
                      "saturating target", "DMSO vehicle", *[f"neighbour:{k}" for k in NEIGHBORS]]
    print(f"     plate controls: {plate_controls}")

    # 4. validation gates A-E (blueprint section 6)
    gates = {
        "A_chemistry_sample": ["identity/stereochem/purity/lot logged", "stock conc confirmed orthogonally",
                               "no precipitate; plastic-adsorption + stability checked",
                               "define added vs accessible/free conc"],
        "B_binder_itself": ["soluble, single SEC peak (SEC-MALS monomer)",
                            "measured 25(OH)D3 Kd (not predicted)", "pocket-dead loses binding",
                            "all 5 neighbours full dose curves; orthogonal method",
                            "2-4 binders pass before cpGFP fusion"],
        "C_sensor_basic": ["cpGFP-only + pocket-dead give no target response",
                          "absolute F >= 10x blank", "8-12 pt dose-response ~1 decade around EC50",
                          "independent protein preps; triplicate", "round-1 dF/F0 >= 2 with CI/Hill"],
        "D_background_mechanism": ["test >=3 sensor concs (rule out conc-driven self-assembly)",
                                  "compare repeat/weak1/truncation main+interaction effects",
                                  "reversibility by dilution/competition", "response+decay time constants",
                                  "pH/temp/salt/DMSO effect on cpGFP itself"],
        "E_selectivity_matrix": ["each neighbour dose + competition vs target",
                                "any neighbour < 10% target response at use conc",
                                "buffer -> delipidated -> real matrix (no skipping)",
                                "spike-recovery 80-120%, intra-CV <=15%",
                                "storage/freeze-thaw/photobleach/lot logged"],
    }
    print("\n[4] VALIDATION GATES A-E")
    for g, items in gates.items():
        print(f"     Gate {g}: {len(items)} checks")

    stop_go = {
        "no measured binder selectivity": "pause sensor library; return to pocket negative design",
        "selective binder but all 24 high background": "increase repeat/interface weakening; add 2+4,4+2 splits",
        "all low background but no signal": "reduce background suppression; scan junction/cpGFP geometry",
        "responds but neighbours also respond": "is it binder selectivity or cpGFP nonspecific env response?",
        "buffer works but matrix fails": "fix free-conc/adsorption/matrix-binding/optics before more sequence work",
    }

    out = {"analyte": ANALYTE, "smiles": SMILES, "polar_anchors": POLAR_ANCHORS,
           "measurement": MEASUREMENT, "neighbour_panel": NEIGHBORS,
           "cdl2_role": CDL2_ROLE, "binder_spec": plan, "transducer": rec,
           "first_round_library": lib, "plate_controls": plate_controls,
           "validation_gates": gates, "stop_go": stop_go,
           "milestone": "purified-protein prototype: soluble, reversible, fluorescence-ON to "
                        "25(OH)D3, discriminates key neighbours in defined buffer (NOT serum yet)",
           "label": "⚠️ CANDIDATE DESIGNS only; backbone needs ProBuilder/RFdiffusionAA; "
                    "dynamic range/selectivity are hypotheses until the wet-lab titration."}
    json.dump(out, open(os.path.join(OUT, "demo_sensor_25ohd3.json"), "w"), indent=2)
    print(f"\n[5] MILESTONE: {out['milestone']}")
    print(f"\nwrote {OUT}/demo_sensor_25ohd3.json")


if __name__ == "__main__":
    run()
