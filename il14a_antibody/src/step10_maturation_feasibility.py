"""Step 10 - Is affinity maturation DEFINED for each candidate line?

Affinity maturation is not a thing you can always attempt. It needs three
things, and if any is missing the exercise produces confident nonsense rather
than an improved binder:

  1. A binding mode that reproduces. Mutations are scored as perturbations of
     a pose; if the pose moves between runs, the perturbation is meaningless.
  2. A scoring function that passes BOTH a positive and a negative control.
     Without a negative control you cannot tell binding from placement.
  3. A signal that exceeds the metric's own noise. If a round's improvement is
     smaller than the spread you get re-measuring one sequence, the round
     taught you nothing.

This applies those three tests to every line now on the table, using each
line's own published numbers. It deliberately does not propose new chemistry -
the output is a go/no-go per line, with the reason.
"""
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS

# --- peptide line: the 5-round groove-registration trajectory (G034x) -------
# round -> best design score
PEP_TRAJECTORY = {"native_36mer_start": 0.407, "R1": 0.755, "R2": 0.723,
                  "R3": 0.820, "R4": 0.792, "R5": 0.792}
# G034x re-measured ONE sequence across rounds and got these three values.
# That is the metric's repeatability, and the yardstick for every delta above.
PEP_REPEATABILITY = [0.820, 0.792, 0.750]
PEP_CONTROLS = {"scrambled_early_rounds": 0.000, "scrambled_late_rounds": 0.643,
                "register_shift": None,   # reported only as ~0.08 below designs
                "native_H3_positive": 1.000}
# leave-one-out effects G034x measured in round 4
PEP_LOO = {"E18R": -0.296, "E23R": -0.116}


def assess(name, reproducible, controls_pass, deltas, noise_sd, notes):
    """deltas: {label: observed improvement}. Returns a verdict dict."""
    resolved = {k: (abs(v) > 2 * noise_sd) for k, v in deltas.items()}
    defined = reproducible and controls_pass and any(resolved.values())
    return {"line": name, "pose_reproducible": reproducible,
            "scorer_passes_controls": controls_pass,
            "noise_sd": round(noise_sd, 4),
            "deltas": {k: round(v, 4) for k, v in deltas.items()},
            "delta_exceeds_2sd": resolved,
            "maturation_defined": defined, "notes": notes}


def main():
    verdicts = []

    # --- peptide line -----------------------------------------------------
    noise = statistics.stdev(PEP_REPEATABILITY)
    traj = PEP_TRAJECTORY
    deltas = {
        "native_start -> R1": traj["R1"] - traj["native_36mer_start"],
        "R1 -> best (R3)": traj["R3"] - traj["R1"],
        "R3 -> final (R5)": traj["R5"] - traj["R3"],
    }
    verdicts.append(assess(
        "peptide (groove registration)", True, True, deltas, noise,
        "Scorer passes its own positive (native H3 = 1.00) and negative "
        "(scrambled = 0.00) controls, unlike ipTM. But the reference interface "
        "is DEFINED by the native STX4 co-fold, so native H3 scoring 1.00 is "
        "near-tautological and the scale is 'overlap with where STX4 binds', "
        "not 'fraction of STX4's affinity'."))

    # --- antibody line, RFantibody (G034x) --------------------------------
    verdicts.append(assess(
        "antibody, RFdiffusion-Ab/RF2 on 502-522", False, False,
        {"best_interaction_pae_vs_threshold": 12.22 - 10.0}, 0.0,
        "0/64 passed interaction_pae<10; binding register undetermined. "
        "G034x correctly declined to mature on an undetermined pose."))

    # --- antibody line, Boltz (this repo) ---------------------------------
    verdicts.append(assess(
        "antibody, Boltz co-fold on 367-382 / 493-523", False, False,
        {"stage3_signal_to_artifact": 1.27 - 1.0}, 0.0,
        "Identical-paratope control failed (framework-only change moved ipTM "
        "0.438); signal-to-artifact 1.27; no design reached "
        "binding_confidence 0.5 across 350."))

    print("=" * 96)
    print("IS AFFINITY MATURATION DEFINED FOR EACH LINE?")
    print("=" * 96)
    for v in verdicts:
        print(f"\n### {v['line']}")
        print(f"  pose reproducible        : {v['pose_reproducible']}")
        print(f"  scorer passes controls   : {v['scorer_passes_controls']}")
        if v["noise_sd"]:
            print(f"  metric noise (1 SD)      : {v['noise_sd']:.4f}"
                  f"   -> resolution limit 2 SD = {2*v['noise_sd']:.3f}")
        for k, d in v["deltas"].items():
            mark = "REAL" if v["delta_exceeds_2sd"][k] else "within noise"
            print(f"    {k:28} {d:+.3f}   {mark}")
        print(f"  => maturation defined    : "
              f"{'YES' if v['maturation_defined'] else 'NO'}")

    print("\n" + "=" * 96)
    print("PEPTIDE LINE - WHERE THE GAIN ACTUALLY CAME FROM")
    print("=" * 96)
    print(f"  metric repeatability from re-measuring one sequence: "
          f"{PEP_REPEATABILITY} -> 1 SD = {noise:.3f}")
    print(f"  The whole resolvable gain is native 0.407 -> R1 0.755 (+0.348,"
          f" {0.348/noise:.1f} SD).")
    print(f"  R1 -> R3 (+0.065) and R3 -> R5 (-0.028) are inside the noise band.")
    print(f"  => Rounds 2-5 produced no improvement this metric can resolve.")
    print(f"     More rounds of the SAME scorer cannot help; its resolution,")
    print(f"     not the chemistry, is now the limit.")
    print(f"\n  Two results from those rounds ARE above noise and do stand:")
    for k, v in PEP_LOO.items():
        print(f"    leave-one-out {k}: {v:+.3f}  ({abs(v)/noise:.1f} SD) - load-bearing")
    print(f"    de-immunisation K35E: 4 strong MHC-II binders -> 0"
          f" (separate assay, not this metric)")

    (RESULTS / "step10_maturation_feasibility.json").write_text(
        json.dumps({"verdicts": verdicts,
                    "peptide_trajectory": PEP_TRAJECTORY,
                    "peptide_repeatability": PEP_REPEATABILITY,
                    "peptide_noise_sd": round(noise, 4),
                    "peptide_leave_one_out": PEP_LOO,
                    "peptide_controls": PEP_CONTROLS}, indent=2))
    print(f"\nwrote results/step10_maturation_feasibility.json")


if __name__ == "__main__":
    main()
