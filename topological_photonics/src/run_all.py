"""Run the whole basic package and produce all four figures + a text summary."""
from __future__ import annotations
import bands
import topology
import edge_states
import param_scan


def main():
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

    print("\nDone. Figures written to ../figures/:")
    for f in ("fig1_bands_dirac_gap.png", "fig2_berry_spin_chern.png",
              "fig3_edge_states_robustness.png", "fig4_param_scan_phase.png"):
        print("   ", f)


if __name__ == "__main__":
    main()
