"""Step 7 - QC of the Boltz-2 stage1 screens. RESULT: both screens FAIL.

THE CONTROL
-----------
Two of the nine screened candidates are a matched pair:

    parent_G034-Ab-Ax2-C-05   QVQLVESGGGLVQAGGSLRLSCAAS...CDR3...WGQGTQVTVSS
    humanised_VHH_C05         EVQLVESGGGLVQPGGSLRLSCAAS...CDR3...WGQGTQVTVSS

They carry the IDENTICAL CDR1 (GFTFSSYA), CDR2 (AISGSGGSTYYADSVKG) and CDR3
(VKEFEHVLHMGAIVK). They differ ONLY at framework positions introduced by step4
humanisation. The paratope is the same molecule surface.

A metric that reports the antigen-paratope interface must therefore give these
two nearly the same score. If it does not, the metric is being driven by
something else and cannot be used to rank candidates.

THE RESULT
----------
It does not. On the coiled-coil epitope the same paratope moves ipTM
0.242 -> 0.680 and interaction PAE 15.3 A -> 6.3 A; on the C-terminal epitope
it moves ipTM 0.960 -> 0.771. Two further checks agree:

  * ipTM is ANTI-correlated with binding_confidence on the C-terminal screen.
    The highest-ipTM candidates carry binding_confidence ~3e-5.
  * The DISORDERED epitope (mean pLDDT 41.6) scores far better than the
    well-ordered one (mean pLDDT 97.5). That is backwards for a real interface
    and is the signature of a flexible chain draping over the binder.
  * No candidate in either screen reaches a meaningful binding_confidence
    (max 0.059 over 18 co-folds).

CONSEQUENCE
-----------
The stage2 and stage3 design runs rank their outputs with the same scoring
machinery that just failed this control, so their rankings are not usable
either. Their generated CDR3 sequences remain valid as a DIVERSITY POOL for
experimental screening; their scores do not order that pool.

This reproduces, with a cleaner diagnosis, the earlier campaign's failure
(24 Boltz-2 models, replicate target-contact Jaccard 0.012). Adding epitope
and non-binding constraints did not rescue it. The blocker is the target -
a long, largely low-confidence coiled-coil protein with no experimental
complex - not the constraint setup.
"""
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS

RES = RESULTS / "boltz_results"
# the matched pair: same CDR1/CDR2/CDR3, framework-only difference
CONTROL = ("parent_G034-Ab-Ax2-C-05", "humanised_VHH_C05")
# a framework-only change should not move a genuine interface more than this
IPTM_TOLERANCE = 0.10
PAE_TOLERANCE = 2.0


def pearson(xs, ys):
    n = len(xs)
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    den = (sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys)) ** 0.5
    return num / den if den else float("nan")


def main():
    screens = {}
    for path in sorted(RES.glob("stage1_*.json")):
        d = json.loads(path.read_text())
        screens[d["epitope"]] = {x["id"]: x for x in d["results"]}

    verdicts, report = {}, {}
    print("=" * 92)
    print("STEP 7 - STAGE-1 SCREEN QC")
    print("=" * 92)

    for epi, rows in screens.items():
        a, b = rows[CONTROL[0]], rows[CONTROL[1]]
        d_iptm = abs(b["iptm"] - a["iptm"])
        d_pae = abs(b["min_interaction_pae"] - a["min_interaction_pae"])
        control_pass = d_iptm <= IPTM_TOLERANCE and d_pae <= PAE_TOLERANCE
        ipt = [x["iptm"] for x in rows.values()]
        bc = [x["binding_confidence"] for x in rows.values()]
        r = pearson(ipt, bc)
        verdicts[epi] = control_pass
        report[epi] = {
            "control_pair": list(CONTROL),
            "delta_iptm": round(d_iptm, 4),
            "delta_min_interaction_pae": round(d_pae, 4),
            "control_passes": bool(control_pass),
            "pearson_iptm_vs_binding_confidence": round(r, 4),
            "max_binding_confidence": max(bc),
            "best_iptm": max(ipt),
        }
        print(f"\n### {epi}")
        print(f"  identical-paratope control  delta ipTM {d_iptm:.3f} "
              f"(tol {IPTM_TOLERANCE})  delta PAE {d_pae:.2f} A (tol {PAE_TOLERANCE})")
        print(f"  control: {'PASS' if control_pass else '*** FAIL ***'}")
        print(f"  r(ipTM, binding_confidence) = {r:+.3f}")
        print(f"  max binding_confidence      = {max(bc):.4f}  "
              f"(no candidate predicted to bind)")
        print(f"  best ipTM                   = {max(ipt):.3f}")

    print("\n" + "=" * 92)
    print("VERDICT")
    print("=" * 92)
    for epi, ok in verdicts.items():
        print(f"  {epi:18} {'usable for ranking' if ok else 'NOT usable for ranking'}")
    if not any(verdicts.values()):
        print("\n  Both screens fail. Do NOT rank candidates on these numbers, and do")
        print("  NOT rank the stage2/stage3 design outputs either - they are scored by")
        print("  the same machinery. Treat the designed CDR3s as an unordered")
        print("  diversity pool for experimental screening, or not at all.")
        print("\n  The blocker is the target, not the constraint setup: adding")
        print("  epitope_residues and non_binding_residues did not rescue convergence.")

    (RES / "step7_screen_qc.json").write_text(json.dumps({
        "control_definition": (
            "parent_G034-Ab-Ax2-C-05 and humanised_VHH_C05 share CDR1 GFTFSSYA, "
            "CDR2 AISGSGGSTYYADSVKG and CDR3 VKEFEHVLHMGAIVK exactly; they differ "
            "only at framework positions from humanisation."),
        "tolerances": {"iptm": IPTM_TOLERANCE, "min_interaction_pae": PAE_TOLERANCE},
        "per_epitope": report,
        "overall_usable_for_ranking": any(verdicts.values()),
    }, indent=2))
    print(f"\nwrote {RES / 'step7_screen_qc.json'}")


if __name__ == "__main__":
    main()
