"""Run the basic package (+ optional complete-package stages) and make figures.

    python run_all.py           # basic package only  (~90 s)  -> fig1..fig4
    python run_all.py --full    # + complete package   (~4 min) -> fig5..fig8 + GIF
"""
from __future__ import annotations
import sys
import bands
import topology
import edge_states
import param_scan


def main(full=False):
    print("=" * 70)
    print("Wu-Hu photonic topological insulator — BASIC computation package")
    print("=" * 70)

    print("\n[1/4] Band structure: Dirac cone + gap opening + band inversion")
    diag = bands.run()
    for k, v in diag.items():
        print(f"      {k.splitlines()[0]:22s} -> {v['phase']} "
              f"(lower|l|={v['l_lower']}, upper|l|={v['l_upper']})")

    print("\n[2/4] Berry curvature + spin-Chern number")
    t = topology.run()
    print(f"      expanded : C_spin={t['topo']['C_spin']:+.0f}, "
          f"C_total={t['topo']['C_total']:+.0f}")
    print(f"      shrunken : C_spin={t['triv']['C_spin']:+.0f}, "
          f"C_total={t['triv']['C_total']:+.0f}")

    print("\n[3/4] Helical edge states + disorder robustness")
    edge_states.run()

    print("\n[4/4] Geometry parameter scan / phase diagram")
    param_scan.run()

    figs = ["fig1_bands_dirac_gap.png", "fig2_berry_spin_chern.png",
            "fig3_edge_states_robustness.png", "fig4_param_scan_phase.png"]

    if full:
        print("\n" + "=" * 70)
        print("COMPLETE package stages")
        print("=" * 70)
        import fdfd_edge
        import fdtd_transport
        import phononic_compare
        print("\n[5] Real-rod FDFD edge states (removes the effective-model caveat)")
        fdfd_edge.run()
        print("\n[6] FDTD wave propagation: straight / bend / defect + GIF")
        fdtd_transport.snapshots_figure()
        d = fdtd_transport.energy_delivery_figure()
        print("    delivered/input:",
              {k.replace(chr(10), ' '): round(v, 3) for k, v in d.items()})
        print("\n[7] Photonic vs phononic two-system comparison")
        phononic_compare.run()
        import bott_disorder
        import loss_budget
        print("\n[8] Real-space spin Bott index under disorder (+ TAI)")
        bott_disorder.run()
        print("\n[9] Loss budget: absorption floor vs topology")
        loss_budget.run()
        figs += ["fig5_fdfd_real_edge.png", "fig6_fdtd_snapshots.png",
                 "fdtd_bend.gif", "fig7_energy_delivery.png",
                 "fig8_photonic_vs_phononic.png", "fig9_spin_bott_disorder.png",
                 "fig10_loss_budget.png"]

    print("\nDone. Figures written to ../figures/:")
    for f in figs:
        print("   ", f)


if __name__ == "__main__":
    main(full="--full" in sys.argv)
