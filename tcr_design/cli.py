"""Command-line entry point.

    python -m tcr_design run 1AO7 --outdir results/1ao7

Stages can be run individually (``interface``, ``scan``, ``offtarget``) when you
only want part of the pipeline.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .benchmark import (
    fit_excluding,
    format_calibration_report,
    held_out_predictions,
    load_combined,
    run_calibration,
    score_entries,
)
from .calibrated import load_model, save_model
from .epistasis import format_report as format_additivity_report
from .epistasis import run_additivity_analysis
from .handoff import write_handoff
from .liability import apply_gates
from .offtarget import download_proteome, peptide_contact_profile, scan_proteome
from .report import (
    plot_additivity,
    plot_calibration,
    plot_gate_funnel,
    plot_interface,
    plot_mutation_heatmap,
    write_report,
    write_summary_json,
    write_tables,
)
from .scan import ranking_key, scan_complex
from .structure import load_complex
from .validate import format_report, retrospective_check

DEFAULT_CACHE = Path.home() / ".cache" / "tcr_design"


def _parse_chain_override(value: str | None) -> dict[str, str]:
    if not value:
        return {}
    override: dict[str, str] = {}
    for item in value.split(","):
        if "=" not in item:
            raise argparse.ArgumentTypeError(
                f"chain override '{item}' should look like role=chain"
            )
        role, chain = item.split("=", 1)
        override[role.strip()] = chain.strip()
    return override


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="tcr_design",
        description=(
            "Structure-based triage of TCR affinity-maturation candidates. "
            "Produces a shortlist and the inputs to score it properly; it does "
            "not predict affinity."
        ),
    )
    parser.add_argument(
        "stage",
        choices=[
            "run", "interface", "scan", "offtarget", "validate", "benchmark",
            "epistasis",
        ],
        help="which part of the pipeline to run",
    )
    parser.add_argument(
        "pdb_id",
        nargs="?",
        default=None,
        help="PDB id of a TCR-pMHC complex, e.g. 1AO7 (not needed for 'benchmark')",
    )
    parser.add_argument(
        "variant_pdb_id",
        nargs="?",
        default=None,
        help=(
            "for 'validate': PDB id of an engineered variant of the same "
            "receptor, e.g. 'validate 3QDG 4L3E'"
        ),
    )
    parser.add_argument("--outdir", type=Path, default=Path("results"))
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--chains",
        default=None,
        help=(
            "override chain typing, e.g. "
            "'mhc=A,b2m=B,peptide=C,tcr_alpha=D,tcr_beta=E'"
        ),
    )
    parser.add_argument(
        "--contact-cutoff", type=float, default=5.0, help="heavy-atom contact cutoff (A)"
    )
    parser.add_argument(
        "--min-contacts",
        type=int,
        default=3,
        help="minimum partner contacts for a position to be scanned",
    )
    parser.add_argument(
        "--loops",
        default=None,
        help="restrict the scan to these loops, e.g. 'A-CDR3,B-CDR3'",
    )
    parser.add_argument(
        "--z-threshold",
        type=float,
        default=1.0,
        help="z-score above which a clean candidate advances",
    )
    parser.add_argument(
        "--shortlist", type=int, default=12, help="how many candidates to hand off"
    )
    parser.add_argument(
        "--offtarget-top", type=int, default=50, help="off-target hits to report"
    )
    parser.add_argument(
        "--offtarget-min-similarity",
        type=float,
        default=0.72,
        help="minimum contact-weighted similarity for an off-target hit",
    )
    parser.add_argument(
        "--offtarget-identity-positions",
        type=int,
        default=1,
        help=(
            "require exact identity at this many of the most-contacted peptide "
            "positions"
        ),
    )
    parser.add_argument(
        "--no-calibration",
        action="store_true",
        help="ignore the calibrated model and use the raw heuristic",
    )
    parser.add_argument(
        "--calibrated",
        action="store_true",
        help=(
            "for 'validate': rank with the calibrated ddG model, refitted "
            "without the structures being tested"
        ),
    )
    parser.add_argument(
        "--save-model",
        action="store_true",
        help="write the fitted model into the package so the pipeline can use it",
    )
    parser.add_argument(
        "--no-atlas", action="store_true", help="benchmark without ATLAS"
    )
    parser.add_argument(
        "--no-skempi", action="store_true", help="benchmark without SKEMPI"
    )
    parser.add_argument(
        "--skip-offtarget",
        action="store_true",
        help="skip the proteome screen (it downloads ~7 MB and takes about a minute)",
    )
    return parser


def _run_epistasis(args, outdir: Path) -> int:
    """Test whether stacked mutations add up, using measured ATLAS data."""
    import pandas as pd

    from .benchmark import download_atlas

    print("[1/2] loading ATLAS and finding measured decompositions")
    try:
        atlas_path = download_atlas(args.cache / "benchmark")
    except Exception as error:  # noqa: BLE001
        print(f"could not obtain ATLAS: {error}", file=sys.stderr)
        return 1

    results = run_additivity_analysis(atlas_path)
    pooled = results["pooled"]["summary"]
    within = results["within_study"]["summary"]
    print(
        f"      {pooled['n']} decompositions pooled, "
        f"{within['n']} within a single study"
    )

    print("\n[2/2] additivity by number of mutations stacked (within one study)")
    for size, stats in within["by_stack_size"].items():
        print(
            f"      {size:2d} mutations: n={stats['n']:2d}  "
            f"mean error {stats['mean_error']:+.2f} kcal/mol  "
            f"within 0.5: {stats['within_half_kcal']}/{stats['n']}"
        )
    print(f"\n      {results['verdict']}")

    outdir.mkdir(parents=True, exist_ok=True)
    rows = results["pooled"]["decompositions"]
    table = outdir / "additivity.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    figure = plot_additivity(results, outdir / "additivity.png")
    report = outdir / "ADDITIVITY.md"
    report.write_text(format_additivity_report(results))
    (outdir / "additivity.json").write_text(
        json.dumps(results, indent=2, default=str) + "\n"
    )
    print(f"\n      report: {report}")
    print(f"      data:   {table}")
    print(f"      figure: {figure}")
    return 0


def _run_benchmark(args, outdir: Path) -> int:
    """Calibrate the triage score against measured binding free energies."""
    import pandas as pd

    print("[1/3] downloading ATLAS and SKEMPI, then scoring every mutation")
    results, scored = run_calibration(
        args.cache, use_atlas=not args.no_atlas, use_skempi=not args.no_skempi
    )
    for note in results["notes"]:
        print(f"      note: {note}")
    print(
        f"      {results['scored']} mutations scored, "
        f"{results['unmapped']} could not be mapped onto a structure"
    )

    unfitted = results["unfitted"]
    print("\n[2/3] the hand-weighted score as shipped")
    print(
        f"      Spearman {unfitted['spearman']:+.3f} (a working score is "
        f"negative), AUC for finding stabilising mutations "
        f"{unfitted['auc_find_stabilising']:.3f}"
    )

    print("\n[3/3] refitting, held out by structure")
    for key, label in (
        ("fitted_components", "six score terms"),
        ("fitted_extended", "raw structural measurements"),
    ):
        held = results[key]["held_out"]
        print(
            f"      {label:28} held-out Spearman {_fmt_opt(held['spearman'])}, "
            f"RMSE {held['rmse']:.3f} vs {held['rmse_predicting_the_mean']:.3f} "
            f"for predicting the mean"
        )

    outdir.mkdir(parents=True, exist_ok=True)
    best_set = "extended"
    predictions = held_out_predictions(
        scored, results[f"fitted_{best_set}"]["ridge_lambda"], best_set
    )
    rows = []
    for entry, prediction in zip(scored, predictions):
        row = entry.as_row()
        row["ddg_predicted_held_out"] = (
            None if prediction != prediction else round(float(prediction), 3)
        )
        rows.append(row)
    table = outdir / "calibration_data.csv"
    pd.DataFrame(rows).to_csv(table, index=False)

    if args.save_model:
        payload = {
            "feature_set": best_set,
            "coefficients": results[f"fitted_{best_set}"]["coefficients"],
            "intercept": results[f"fitted_{best_set}"]["intercept"],
            "ridge_lambda": results[f"fitted_{best_set}"]["ridge_lambda"],
            "units": "kcal/mol, positive means weaker binding",
            "performance": {
                "held_out": results[f"fitted_{best_set}"]["held_out"],
                "in_sample": results[f"fitted_{best_set}"]["in_sample"],
                "validation": (
                    "leave-one-structure-out: every reported number comes from "
                    "a model that never saw any mutation from the structure it "
                    "was scored on"
                ),
            },
            "provenance": {
                **results["provenance"],
                "scored": results["scored"],
                "sources": ["ATLAS", "SKEMPI 2.0 (TCR/pMHC subset)"],
            },
        }
        saved = save_model(payload)
        print(f"      model saved to {saved}")

    figure = plot_calibration(scored, predictions, results, outdir / "calibration.png")
    report = outdir / "CALIBRATION.md"
    report.write_text(format_calibration_report(results))
    (outdir / "calibration_model.json").write_text(
        json.dumps(results, indent=2, default=str) + "\n"
    )

    print(f"\n      report: {report}")
    print(f"      data:   {table}")
    print(f"      figure: {figure}")
    return 0


def _fmt_opt(value) -> str:
    return "n/a" if value is None else f"{value:+.3f}"


def _run_validate(args, outdir: Path) -> int:
    """Benchmark the triage score against a known engineered variant."""
    if not args.variant_pdb_id:
        print(
            "validate needs two structures: a wild type and an engineered "
            "variant of the same receptor, e.g.\n"
            "  python -m tcr_design validate 3QDG 4L3E",
            file=sys.stderr,
        )
        return 2

    print(
        f"benchmarking: scanning {args.pdb_id} blind, then "
        f"looking up the mutations that distinguish {args.variant_pdb_id}"
    )

    ddg_model = None
    if args.calibrated:
        # Both structures are in the benchmark data, so a model trained on all
        # of it would be marking its own homework. Refit without them.
        entries, _ = load_combined(
            args.cache / "benchmark" / "skempi_v2.csv",
            args.cache / "benchmark" / "atlas_Mutants.csv",
        )
        scored, _ = score_entries(entries, args.cache / "pdb", verbose=False)
        ddg_model = fit_excluding(
            scored, (args.pdb_id, args.variant_pdb_id)
        )
        print(
            f"      using a calibrated model refitted without "
            f"{args.pdb_id}/{args.variant_pdb_id}: "
            f"{ddg_model.provenance['training_mutations']} mutations from "
            f"{ddg_model.provenance['training_structures']} other structures"
        )

    result = retrospective_check(
        args.pdb_id,
        args.variant_pdb_id,
        cache_dir=args.cache / "pdb",
        min_contacts=args.min_contacts,
        z_threshold=args.z_threshold,
        ddg_model=ddg_model,
    )
    if ddg_model is not None:
        result["model"] = ddg_model.provenance
    report = format_report(result)
    path = outdir / f"validation_{args.pdb_id.lower()}_{args.variant_pdb_id.lower()}.md"
    path.write_text(report)
    (outdir / f"validation_{args.pdb_id.lower()}_{args.variant_pdb_id.lower()}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print()
    print(report)
    print(f"written to {path}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.stage == "epistasis":
        return _run_epistasis(args, outdir)
    if args.stage == "benchmark":
        return _run_benchmark(args, outdir)
    if not args.pdb_id:
        print("a PDB id is required for this stage", file=sys.stderr)
        return 2
    if args.stage == "validate":
        return _run_validate(args, outdir)

    print(f"[1/5] loading {args.pdb_id} and measuring the interface")
    complex_obj = load_complex(
        args.pdb_id,
        cache_dir=args.cache / "pdb",
        contact_cutoff=args.contact_cutoff,
        chain_override=_parse_chain_override(args.chains),
    )
    for note in complex_obj.typing_notes:
        print(f"      {note}")
    paratope = complex_obj.paratope()
    buried = sum(r.sasa_free - r.sasa_complex for r in paratope)
    print(
        f"      paratope {len(paratope)} residues, "
        f"{buried:.0f} A^2 of TCR surface buried"
    )
    contact_profile = peptide_contact_profile(complex_obj, args.contact_cutoff)
    print(
        "      peptide contact profile: "
        + ", ".join(
            f"P{i + 1}{a}={c}"
            for i, (a, c) in enumerate(zip(complex_obj.peptide, contact_profile))
        )
    )

    if args.stage == "interface":
        tables = write_tables(complex_obj, [], [], outdir)
        print(f"      wrote {tables['interface']}")
        return 0

    print("[2/5] scanning single mutants over the contacting CDR positions")
    loops = tuple(x.strip() for x in args.loops.split(",")) if args.loops else None
    ddg_model = None if args.no_calibration else load_model()
    if ddg_model is not None:
        print(f"      {ddg_model.summary()}")
    else:
        print(
            "      no calibrated model found: falling back to the uncalibrated "
            "heuristic. Run 'python -m tcr_design benchmark --save-model' to fit one."
        )
    mutations = scan_complex(
        complex_obj,
        min_contacts=args.min_contacts,
        loops=loops,
        ddg_model=ddg_model,
    )
    mutations = apply_gates(complex_obj, mutations, z_threshold=args.z_threshold)
    print(f"      {len(mutations)} mutants scored")

    advance = sorted(
        [m for m in mutations if m.verdict == "advance"], key=ranking_key
    )
    watch = sorted(
        [m for m in mutations if m.verdict == "watch-specificity"], key=ranking_key
    )
    shortlist = (advance + watch)[: args.shortlist]
    print(
        f"      {len(advance)} advance, {len(watch)} watch-specificity; "
        f"shortlisting {len(shortlist)}"
    )
    if not shortlist:
        # Nothing cleared the bar. Hand the best surviving candidates to tier 2
        # anyway: the calibrated model is documented to under-call real
        # improvements, so an empty shortlist is a reason to spend the compute,
        # not to stop. The DMF5 case in `validate` is exactly this failure.
        fallback = sorted(
            [m for m in mutations if m.verdict == "below-threshold"], key=ranking_key
        )[: args.shortlist]
        shortlist = fallback
        print(
            f"      nothing cleared the threshold, so the {len(shortlist)} best "
            "surviving candidates are handed on anyway: the model under-calls "
            "real improvements (see 'validate 3QDG 4L3E --calibrated')"
        )

    if args.stage == "scan":
        tables = write_tables(complex_obj, mutations, [], outdir)
        print(f"      wrote {tables['scan']}")
        return 0

    hits: list = []
    motif: dict = {}
    if not args.skip_offtarget:
        print("[3/5] screening the reviewed human proteome for cross-reactivity risk")
        proteome = download_proteome(args.cache / "proteome")
        hits, motif = scan_proteome(
            complex_obj,
            proteome,
            top_n=args.offtarget_top,
            min_similarity=args.offtarget_min_similarity,
            require_identity_top_n=args.offtarget_identity_positions,
        )
        print(
            f"      scanned {motif['windows_scanned']:,} windows, "
            f"{motif['windows_passing_all_filters']:,} passed the filters, "
            f"{len(hits)} reported"
        )
        if hits:
            top = hits[0]
            print(
                f"      closest human mimic: {top.peptide} "
                f"({top.gene}, similarity {top.similarity:.2f})"
            )
    else:
        print("[3/5] skipping the proteome screen (--skip-offtarget)")
        motif = {
            "target_peptide": complex_obj.peptide,
            "tcr_contacts_per_position": contact_profile,
            "weighted_positions": [],
            "identity_required_at": [],
            "anchor_filter": {},
            "windows_scanned": 0,
            "windows_passing_all_filters": 0,
            "hits_above_threshold": 0,
            "min_similarity": args.offtarget_min_similarity,
            "weights_used": [],
        }

    print("[4/5] writing tier-2 inputs")
    handoff_dir = outdir / "tier2"
    handoff = write_handoff(complex_obj, shortlist, handoff_dir)
    print(f"      {len(handoff.written)} files under {handoff_dir}")

    print("[5/5] writing tables, figures and the report")
    tables = write_tables(complex_obj, mutations, hits, outdir)
    figures = {
        "interface": plot_interface(
            complex_obj, contact_profile, outdir / "interface_profile.png"
        ),
        "heatmap": plot_mutation_heatmap(mutations, outdir / "mutation_heatmap.png"),
        "funnel": plot_gate_funnel(mutations, outdir / "gate_outcomes.png"),
    }
    report = write_report(
        complex_obj,
        mutations,
        hits,
        motif,
        contact_profile,
        figures,
        tables,
        handoff_dir,
        outdir / "REPORT.md",
    )
    write_summary_json(complex_obj, mutations, hits, motif, outdir / "summary.json")
    print(f"      report: {report}")
    if ddg_model is None:
        print(
            "\nReminder: without the calibrated model the ranking is an "
            "uncalibrated heuristic, and the off-target list is a similarity "
            "shortlist. Neither is a prediction. See section 6 of the report."
        )
    else:
        print(
            "\nReminder: predicted ddG comes from a weak model (held-out "
            "Spearman +0.43, RMSE 0.88 kcal/mol) that under-calls real "
            "improvements, and the off-target list is a similarity shortlist. "
            "See section 6 of the report."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
