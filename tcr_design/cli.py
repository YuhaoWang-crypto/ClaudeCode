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

from .handoff import write_handoff
from .liability import apply_gates
from .offtarget import download_proteome, peptide_contact_profile, scan_proteome
from .report import (
    plot_gate_funnel,
    plot_interface,
    plot_mutation_heatmap,
    write_report,
    write_summary_json,
    write_tables,
)
from .scan import scan_complex
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
        choices=["run", "interface", "scan", "offtarget", "validate"],
        help="which part of the pipeline to run",
    )
    parser.add_argument("pdb_id", help="PDB id of a TCR-pMHC complex, e.g. 1AO7")
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
        "--skip-offtarget",
        action="store_true",
        help="skip the proteome screen (it downloads ~7 MB and takes about a minute)",
    )
    return parser


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
        f"benchmarking the triage score: scanning {args.pdb_id} blind, then "
        f"looking up the mutations that distinguish {args.variant_pdb_id}"
    )
    result = retrospective_check(
        args.pdb_id,
        args.variant_pdb_id,
        cache_dir=args.cache / "pdb",
        min_contacts=args.min_contacts,
        z_threshold=args.z_threshold,
    )
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
    mutations = scan_complex(
        complex_obj, min_contacts=args.min_contacts, loops=loops
    )
    mutations = apply_gates(complex_obj, mutations, z_threshold=args.z_threshold)
    print(f"      {len(mutations)} mutants scored")

    advance = sorted(
        [m for m in mutations if m.verdict == "advance"],
        key=lambda m: -(m.peptide_z if m.peptide_z is not None else m.z_score),
    )
    watch = sorted(
        [m for m in mutations if m.verdict == "watch-specificity"],
        key=lambda m: -(m.peptide_z if m.peptide_z is not None else m.z_score),
    )
    shortlist = (advance + watch)[: args.shortlist]
    print(
        f"      {len(advance)} advance, {len(watch)} watch-specificity; "
        f"shortlisting {len(shortlist)}"
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
    print(
        "\nReminder: the triage score is uncalibrated and the off-target list is a "
        "similarity shortlist. Neither is a prediction. See section 6 of the report."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
