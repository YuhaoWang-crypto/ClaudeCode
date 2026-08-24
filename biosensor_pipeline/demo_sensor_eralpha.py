"""
demo_sensor_eralpha.py -- END-TO-END estradiol sensor built on the ERalpha-LBD binder we
VALIDATED in silico (biosensor_out/newproject_ERalpha_specificity.json): real Boltz-2.1
co-folds ranked estradiol (binding_confidence 0.715) far above testosterone (0.334) and
cortisol (0.201). This demo connects that binder to the transducer layer.

Unlike the placeholder-bundle demos, the receptor sequence here is REAL (ESR1 LBD) and its
estradiol preference is grounded on real co-folds. What remains ⚠️ is the transduction:
whether a given fusion switches, and its dynamic range, needs a wet-lab titration.

Two coupling routes (from transducer.py / the Cao paper), primary = the biologically real one:
  PRIMARY  ligand-induced folding: estradiol binding repositions helix-12 (AF-2); read by a
           FRET pair (the classic nuclear-receptor sensor route).
  ALT      split at the H11/H12 (AF-2) junction + cpGFP, with background-suppression knobs.

Run:  python3 -m biosensor_pipeline.demo_sensor_eralpha
"""
from __future__ import annotations
import json, os
from .transducer import (build_induced_folding_sensor, build_split_bundle_sensor,
                         verify_split_bundle, recommend_transducer)

ANALYTE = "estradiol (17b-estradiol)"
SMILES = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O"
# ESR1 ligand-binding domain (UniProt P03372, residues 305-554)
ERA_LBD = ("SLALSLTADQMVSALLDAEPPILYSEYDPTRPFSEASMMGLLTNLADRELVHMINWAKRVPGFVDLTLHDQVHLLE"
           "CAWLEILMIGLVWRSMEHPGKLLFAPNLLLDRNQGKCVEGMVEIFDMLLATSSRFRMMNLQGEEFVCLKSIILLNS"
           "GVYTFLSSTLKSLEEKDHIHRVLDKITDTLIHLMAKAGLTLQQQHQRLAQLLLILSHIRHMSNKGMEHLYSMKCKN"
           "VVPLYDLLLEMLDAHRLHAPTS")

# validated specificity (real Boltz-2.1 co-folds; see newproject_ERalpha_specificity.json)
SPECIFICITY = {
    "estradiol":    {"binding_confidence": 0.715, "role": "target (native agonist)"},
    "testosterone": {"binding_confidence": 0.334, "role": "off-target (androgen)"},
    "cortisol":     {"binding_confidence": 0.201, "role": "off-target (glucocorticoid)"},
}
AF2_MOTIF = "LLEMLDAHR"      # helix-12 / AF-2 (its repositioning is the switch)


def run():
    OUT = os.path.join(os.getcwd(), "biosensor_out")
    os.makedirs(OUT, exist_ok=True)
    print("=== ERalpha estradiol sensor: validated binder -> transducer ===\n")

    margin = (SPECIFICITY["estradiol"]["binding_confidence"]
              - max(v["binding_confidence"] for k, v in SPECIFICITY.items() if k != "estradiol"))
    print(f"[0] BINDER (real, validated): ESR1 LBD 250 aa; estradiol binding_confidence "
          f"{SPECIFICITY['estradiol']['binding_confidence']} >> off-targets; specificity "
          f"margin +{margin:.3f}  ✅ (real co-folds)")

    # PRIMARY: ligand-induced folding -> FRET (helix-12 repositioning)
    fret = build_induced_folding_sensor(ANALYTE, "ERa_LBD", ERA_LBD, reporter="FRET",
                                        fret_pair=("mTurquoise2", "mVenus"))
    print(f"\n[1] PRIMARY sensor - ligand-induced folding (helix-12/AF-2 repositioning) -> FRET")
    print(f"    topology: [mTurquoise2]-ERa_LBD-[mVenus]; {fret.__dict__.get('note','')}")
    print(f"    mechanism note: estradiol binding orders AF-2 ({AF2_MOTIF}); the donor-acceptor")
    print(f"    distance change reports it. (also viable: AF-2 -> coactivator-peptide recruitment)")

    # ALT: split at the H11/H12 (AF-2) junction + cpGFP, with background knobs
    cut = ERA_LBD.find(AF2_MOTIF)                       # split just before AF-2 helix
    boundaries = [0, cut]                               # frag_N = H1..H11 ; frag_C = H12/AF-2
    split = build_split_bundle_sensor(
        ANALYTE, "ERa_LBD", ERA_LBD, boundaries, split_after_helix=1,
        duplicate_segment=(max(0, cut - 8), cut),       # duplicate 8 aa at the junction -> lower background
        weaken_positions=[1, 5])
    chk = verify_split_bundle(split)
    print(f"\n[2] ALT sensor - split at H11/H12 (idx {cut}) + cpGFP, background knobs "
          f"{list(split.knobs)}")
    print(f"    construct_N: frag_N({len(split.frag_n)} aa)+linker+cpGFP [{len(split.construct_n)}]")
    print(f"    construct_C: cpGFP+linker+frag_C({len(split.frag_c)} aa) [{len(split.construct_c)}]  "
          f"verify={chk['all_ok']}")

    # validation gates (specificity gate PARTLY satisfied by the real co-folds)
    gates = {
        "A_binder(real)": "✅ estradiol co-fold >> testosterone/cortisol (margin +%.3f). "
                          "TODO: measured Kd + broader steroid panel (estrone, DES, tamoxifen)." % margin,
        "B_folding": "apo LBD folds (Boltz/AF2 pLDDT); holo orders AF-2 helix-12 "
                     "(compare apo vs holo helix-12 position).",
        "C_transduction(FRET)": "express FRET fusion; measure ratiometric FRET +/- estradiol "
                                "titration; ΔR/R0 and apparent EC50.",
        "C_transduction(split)": "express both cpGFP fusions; fluorescence +/- estradiol; "
                                 "fragment-duplication should suppress apo background.",
        "D_selectivity": "titrate testosterone/cortisol/estrone/tamoxifen; each < 10% of "
                         "estradiol response at use conc; competition assay.",
        "E_reversibility_matrix": "reversibility by dilution/competition; then buffer -> "
                                  "serum-mimic (estradiol is ~pM-nM in serum, largely SHBG-bound "
                                  "-> define free vs total first).",
    }
    print("\n[3] VALIDATION GATES")
    for g, v in gates.items():
        print(f"    {g}: {v}")

    out = {"analyte": ANALYTE, "smiles": SMILES, "receptor": "ESR1 LBD (P03372, 305-554, 250aa)",
           "binder_validation": SPECIFICITY, "specificity_margin_binding_confidence": round(margin, 3),
           "primary_sensor": {"mechanism": "ligand-induced folding (helix-12/AF-2) -> FRET",
                              "topology": fret.sensor_seq, "reporter": "mTurquoise2/mVenus FRET"},
           "alt_sensor": {"mechanism": "split at H11/H12 + cpGFP", "split_index": cut,
                          "knobs": split.knobs, "construct_n_len": len(split.construct_n),
                          "construct_c_len": len(split.construct_c), "verify": chk},
           "validation_gates": gates,
           "label": "✅ binder + its estradiol specificity are real Boltz-2.1 co-folds; "
                    "⚠️ sensor constructs are CANDIDATE DESIGNS - switch behaviour & dynamic "
                    "range need a wet-lab FRET/fluorescence titration."}
    json.dump(out, open(os.path.join(OUT, "demo_sensor_eralpha.json"), "w"), indent=2)
    print(f"\nwrote {OUT}/demo_sensor_eralpha.json")


if __name__ == "__main__":
    run()
