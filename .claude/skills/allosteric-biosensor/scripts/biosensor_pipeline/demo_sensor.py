"""
demo_sensor.py -- END-TO-END sensor design demo, tying together the two papers we
learned:

  small molecule (SMILES)
    -> helical-bundle BINDER SPEC          (Cao/ProBuilder: grow bundle around ligand
                                            from RIF anchors; fold+bind scored jointly)
    -> SPLIT-BUNDLE + cpGFP SENSOR         (transducer coupling)
       + background-suppression KNOBS      (fragment duplication / interface weakening)
    -> VALIDATION CHECKLIST                (what must be measured, in-silico + wet-lab)

Runnable, deterministic, offline:
    python3 -m biosensor_pipeline.demo_sensor

Honesty: the binder SPEC is a recipe (backbone generation needs ProBuilder/RFdiffusionAA,
not run here) so the split machinery is demonstrated on a ⚠️ PLACEHOLDER idealized
helical-bundle sequence. The pipeline structure, the knobs, and the checklist are real.
"""
from __future__ import annotations
import json, os
from .transducer import (helical_bundle_binder_plan, build_split_bundle_sensor,
                         verify_split_bundle, recommend_transducer, CPGFP)


# example analyte: serotonin (5-HT), a very small neurotransmitter the paper sensed
ANALYTE = "serotonin"
SMILES = "C1=CC2=C(C=C1O)C(=CN2)CCN"     # serotonin
PRIVILEGED = ["5-OH (indole)", "primary amine (ethylamine)"]  # key polar groups to H-bond


def _placeholder_bundle(n_helices=5, helix_len=20, loop="GSG"):
    """⚠️ Idealized amphipathic helix bundle standing in for a ProBuilder/RFdiffusionAA
    output, so the split + knob machinery runs end-to-end. NOT a real binder."""
    helix = "EIAALKQEIAALKKENAALK"[:helix_len]     # heptad-repeat amphipathic helix
    seq, boundaries, cur = "", [], 0
    for i in range(n_helices):
        boundaries.append(cur)
        seq += helix
        cur += helix_len
        if i < n_helices - 1:
            seq += loop
            cur += len(loop)
    return seq, boundaries


def run():
    OUT = os.path.join(os.getcwd(), "biosensor_out")
    os.makedirs(OUT, exist_ok=True)
    print(f"=== END-TO-END sensor design demo: {ANALYTE} ===\n")

    # 1. binder spec (the tractable de-novo recipe)
    plan = helical_bundle_binder_plan(ANALYTE, SMILES, n_helices=5, length=120,
                                      privileged_contacts=PRIVILEGED)
    print(f"[1] BINDER SPEC — {plan['topology']}, {plan['length_aa']} aa; "
          f"anchor privileged H-bonds to {PRIVILEGED}")
    for s in plan["steps"]:
        print("     " + s)

    # 2. transducer choice
    rec = recommend_transducer("small_molecule", "helical_bundle", "fluorescent")
    print(f"\n[2] TRANSDUCER — {rec['recommended']}; knobs: {rec['background_knobs']}")

    # 3. build the split-bundle + cpGFP sensor WITH background-suppression knobs
    bundle, boundaries = _placeholder_bundle()
    sensor = build_split_bundle_sensor(
        ANALYTE, f"{ANALYTE}_binder", bundle, boundaries, split_after_helix=2,
        duplicate_segment=(0, 8),           # duplicate 8 res of helix-1 -> raise self-assoc barrier
        weaken_positions=[3, 7])            # weaken 2 interface positions
    chk = verify_split_bundle(sensor)
    print(f"\n[3] SPLIT+cpGFP SENSOR (⚠️ placeholder bundle) — split after helix 2; "
          f"knobs applied: {list(sensor.knobs)}")
    print(f"     construct_N: frag_N({len(sensor.frag_n)} aa) + linker + cpGFP  "
          f"[len {len(sensor.construct_n)}]")
    print(f"     construct_C: cpGFP + linker + frag_C({len(sensor.frag_c)} aa)  "
          f"[len {len(sensor.construct_c)}]")
    print(f"     verify: {chk['all_ok']}")

    # 4. validation checklist
    checklist = [
        "in-silico: apo bundle folds (Boltz/AF2 pLDDT); holo retains ligand contact "
        "(Boltz ligand_iptm) + the privileged H-bonds present.",
        "in-silico: split fragments do NOT co-fold without ligand (low interface conf), "
        "DO co-fold with ligand (proximity restored) -> the switch.",
        "wet-lab: express both cpGFP fusions; measure fluorescence +/- analyte titration.",
        "wet-lab: BACKGROUND = signal without ligand (fragment-duplication should lower it); "
        "RESPONSE = fold-change with saturating ligand.",
        "wet-lab: SPECIFICITY = titrate structurally-similar off-targets; DYNAMIC RANGE + "
        "apparent Kd from the titration (the ground-truth number).",
        "iterate: tune knobs (more/less duplication, interface weakening, terminal truncation) "
        "to trade background vs response.",
    ]
    print("\n[4] VALIDATION CHECKLIST")
    for i, c in enumerate(checklist, 1):
        print(f"     {i}. {c}")

    out = {"analyte": ANALYTE, "smiles": SMILES, "binder_spec": plan,
           "transducer": rec,
           "sensor": {"split_after_helix": sensor.split_after_helix, "knobs": sensor.knobs,
                      "construct_n_len": len(sensor.construct_n),
                      "construct_c_len": len(sensor.construct_c),
                      "verify": chk, "placeholder_bundle": True},
           "validation_checklist": checklist,
           "label": "⚠️ binder backbone needs ProBuilder/RFdiffusionAA; split shown on a "
                    "placeholder bundle. Pipeline + knobs + checklist are real; dynamic range "
                    "is a hypothesis until the wet-lab titration."}
    json.dump(out, open(os.path.join(OUT, "demo_sensor_serotonin.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/demo_sensor_serotonin.json")


if __name__ == "__main__":
    run()
