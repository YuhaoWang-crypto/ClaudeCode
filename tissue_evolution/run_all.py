"""
Run the whole tissue-evolution pipeline in dependency order.

    python3 -m tissue_evolution.run_all              # everything
    python3 -m tissue_evolution.run_all --skip t7 t8 # skip the slow ML/structure steps

T2 (BioMart + Ensembl FTP) and T3 (dN/dS for ~60k orthologue pairs) dominate the
runtime; both cache, so a re-run is cheap. T7 needs ~2.5 GB of ESM-2 weights and
T8 needs the foldseek binary - if either is unavailable the step reports that and
the rest of the pipeline still completes.
"""
import argparse
import traceback

STEPS = [
    ("t1", "human tissue expression, tau, enriched sets", "t1_expression"),
    ("t2", "1:1 orthologues + CDS across the ladder", "t2_orthologs"),
    ("t3", "dN/dS from scratch (NG86)", "t3_dnds"),
    ("t4", "per-organ rates, three estimators", "t4_tissue_rates"),
    ("t5", "composition vs intrinsic tissue effect", "t5_confounders"),
    ("t6", "rank stability across divergence depth", "t6_depth"),
    ("t7", "ESM-2 constraint incl. orthologue-less genes", "t7_esm"),
    ("t8", "structural vs sequence divergence", "t8_structure"),
    ("t9", "lineage-specific rate asymmetry", "t9_lineage"),
    ("t10", "organ-profile conservation vs sequence rate", "t10_expression_divergence"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip", nargs="*", default=[])
    ap.add_argument("--only", nargs="*", default=[])
    args = ap.parse_args()

    import importlib
    failed = []
    for key, label, module in STEPS:
        if key in args.skip or (args.only and key not in args.only):
            print(f"\n### skipping {key.upper()} - {label}")
            continue
        print(f"\n\n### {key.upper()} - {label}")
        try:
            mod = importlib.import_module(f"tissue_evolution.{module}")
            (getattr(mod, "report", None) or mod.main)()
        except Exception as exc:                          # noqa: BLE001
            failed.append((key, exc))
            print(f"!!! {key.upper()} failed: {type(exc).__name__}: {exc}")
            traceback.print_exc()

    print("\n" + "=" * 72)
    if failed:
        print("completed with failures:")
        for key, exc in failed:
            print(f"  {key.upper()}: {type(exc).__name__}: {exc}")
    else:
        print("all steps completed")


if __name__ == "__main__":
    main()
