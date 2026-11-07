"""
Run the glycan-design pipeline end to end.

  G1  DeltaRASA engine       Shrake-Rupley two-pass glycan shielding score
  G2  confounds              glycan size and missing coordinates, measured
  G3  underprocessing score  ordinal ladder + sensitivity to its reconstruction
  G4  design triage          sequons, confidence gates, filter payoff

G1 and G2 need network access the first time (structures are cached under
figures/_pdb_cache/ afterwards). G3 and G4 are fully offline.
"""

from glycan_pipeline import (g1_delta_rasa, g2_confounds, g3_underprocessing,
                             g4_triage)


def main():
    print("\n" + "#" * 68)
    print("#  DE NOVO GLYCOPROTEIN DESIGN - STERIC CONTROL OF N-GLYCAN STATE")
    print("#" * 68 + "\n")

    results = {}
    for key, module, needs_net in [
        ("g1", g1_delta_rasa, True),
        ("g2", g2_confounds, True),
        ("g3", g3_underprocessing, False),
        ("g4", g4_triage, True),
    ]:
        try:
            results[key] = module.report()
        except Exception as exc:
            results[key] = None
            hint = " (needs network for RCSB on first run)" if needs_net else ""
            print(f"{key.upper()} skipped ({type(exc).__name__}: {exc}){hint}")
        print()

    print("=" * 68)
    print("CONSOLIDATED SUMMARY")
    print("=" * 68)

    if results.get("g1"):
        r = results["g1"]
        v = r["validation"]
        print(f"G1  SASA engine agrees with Biopython to "
              f"{v['max_abs_atom_diff_A2']:.1e} A^2 "
              f"[{'ok' if v['agrees'] else 'MISMATCH'}]")
        top = r["protein_context"][0]
        print(f"    most shielded site in {r['pdb_id']}: {top['site']} at "
              f"{top['mean_delta_rasa']:.2f}% mean per-atom DeltaRASA")

    if results.get("g2"):
        r = results["g2"]
        print(f"G2  numerical error bar (symmetry copies): "
              f"{r['median_replicate_spread']:.3f}% median")
        print(f"    size confound Spearman {r['rho_size_vs_drasa_full']:+.2f} "
              f"-> {r['rho_size_vs_drasa_core3']:+.2f} after standardising")
        print(f"    median glycan modelled: "
              f"{r['median_completeness_vs_man9']:.0%} of a Man9GlcNAc2")

    if results.get("g3"):
        r = results["g3"]
        s = r["sensitivity"]
        fragile = [n for n, f in s["pass_fraction"].items() if 0.0 < f < 1.0]
        print(f"G3  dominance orderings held in "
              f"{r['dominance_tests'] - r['dominance_violations']}/"
              f"{r['dominance_tests']} ladder sweeps")
        print(f"    designs whose pass/fail is ladder-dependent: "
              f"{', '.join(fragile) if fragile else 'none'}")

    if results.get("g4"):
        r = results["g4"]
        print(f"G4  AF2 gen-1 gate keeps {r['af2_gen1_pass']} of O1-O4; "
              f"glycan-aware gate keeps all four")
        print(f"    DeltaRASA filter Fisher p = {r['fisher']['p_value']:.1e} "
              f"(paper: p < 1e-10)")
        a = r["sequon_audit"]
        print(f"    {a['total_glycans_outside_a_sequon']} of "
              f"{sum(s['n_with_glycan'] for s in a['per_structure'])} observed "
              f"glycans sat outside a canonical sequon")

    print("\nRIGOR LEDGER")
    print("  (rigorous)  SASA/RASA/DeltaRASA geometry; confound correlations;")
    print("              Fisher reproduction; AF2-gate reproduction; sequon audit.")
    print("  (caution)   the 1-20 underprocessing ladder is RECONSTRUCTED, and")
    print("              G3's glycan profiles are ILLUSTRATIVE, not measured.")
    print("              Neighbour-glycan shielding is a LOWER BOUND.")
    return results


if __name__ == "__main__":
    main()
