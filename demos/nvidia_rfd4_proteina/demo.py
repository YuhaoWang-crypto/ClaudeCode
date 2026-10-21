#!/usr/bin/env python3
"""Demo driver for the hosted NVIDIA RFD4-Proteina NIM endpoint.

Three subcommands:

    check       probe auth + routing without spending GPU time
    generate    unconditional backbone/sequence generation from a template
    scaffold    motif scaffolding: keep a motif, design new flanks

Examples
--------
    python demo.py check
    python demo.py generate --pdb-id 1CRN --replicates 2 --seed 42
    python demo.py scaffold --pdb-id 1CRN --motif A1-15 --n-term 30 --c-term 30

Add --dry-run to print the request without sending it (no key needed).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Any

import proteina_client as pc

HERE = Path(__file__).resolve().parent
INPUT_DIR = HERE / "inputs"
OUTPUT_DIR = HERE / "outputs"
RCSB_TEMPLATE = "https://files.rcsb.org/download/{pdb_id}.cif"

# Metrics worth showing on one screen; the rest land in summary.json.
HEADLINE_METRICS = (
    "clash_metrics.frac_atoms_clashing_1p5A",
    "bond_length_metrics.bond_length.frac_unreasonable",
    "radius_of_gyration_metrics.radius_of_gyration",
    "secondary_structure_metrics.frac_helix",
    "secondary_structure_metrics.frac_strand",
    "secondary_structure_metrics.frac_loop",
)

MOTIF_RE = re.compile(r"^(?P<chain>[A-Za-z0-9]+?)(?P<first>\d+)-(?P<last>\d+)$")


def count_ca_atoms(cif_text: str) -> int | None:
    """Count CA atoms in a returned mmCIF, as a cheap structural self-check.

    The endpoint returns an *annotated* mmCIF whose ``_atom_site`` loop
    carries extra condition columns, so the atom-name column has to be read
    from the loop header rather than assumed to be at a fixed index.
    """
    columns: list[str] = []
    atom_name_index: int | None = None
    count = 0
    for line in cif_text.splitlines():
        if line.startswith("_atom_site."):
            columns.append(line.strip())
            continue
        if line.startswith(("ATOM", "HETATM")):
            if atom_name_index is None:
                if "_atom_site.label_atom_id" not in columns:
                    return None
                atom_name_index = columns.index("_atom_site.label_atom_id")
            fields = line.split()
            if len(fields) > atom_name_index and fields[atom_name_index] == "CA":
                count += 1
    return count if atom_name_index is not None else None


def fetch_template(pdb_id: str) -> Path:
    """Download an mmCIF from RCSB into inputs/, cached by id."""
    pdb_id = pdb_id.upper()
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    target = INPUT_DIR / f"{pdb_id.lower()}.cif"
    if target.is_file() and target.stat().st_size > 0:
        return target

    url = RCSB_TEMPLATE.format(pdb_id=pdb_id)
    print(f"  fetching template {pdb_id} from RCSB ...")
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read()
    target.write_bytes(data)
    return target


def resolve_template(args: argparse.Namespace) -> Path:
    if args.template:
        path = Path(args.template)
        if not path.is_file():
            raise SystemExit(f"template not found: {path}")
        return path
    return fetch_template(args.pdb_id)


def parse_motif(spec: str) -> tuple[str, int, int]:
    """Parse a motif like 'A1-15' into (chain, first, last)."""
    match = MOTIF_RE.match(spec.strip())
    if not match:
        raise SystemExit(f"could not parse --motif {spec!r}; expected e.g. A1-15")
    chain = match.group("chain")
    first, last = int(match.group("first")), int(match.group("last"))
    if last < first:
        raise SystemExit(f"--motif {spec!r}: last residue precedes first")
    return chain, first, last


def report(response: dict[str, Any], out_dir: Path) -> None:
    """Print a short summary and write structures + metrics to out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "response.json").write_text(json.dumps(response, indent=2), encoding="utf-8")

    model = response.get("model", {})
    runtime = response.get("runtime", {})
    schedule = response.get("schedule", {})
    print(
        f"\n  task={response.get('task')}  model={model.get('model_name')} "
        f"v{model.get('package_version')}  gpu={runtime.get('gpu')}"
    )
    print(
        f"  steps={schedule.get('num_steps')} "
        f"recycle={schedule.get('recycle_strategy')} "
        f"trunk_passes={schedule.get('trunk_passes')} "
        f"total={response.get('duration_ms', 0) / 1000:.2f}s"
    )
    for warning in response.get("warnings") or []:
        print(f"  warning: {warning}")

    summary = []
    for design in response.get("designs", []):
        diagnostics = design.get("diagnostics", {})
        metrics = design.get("metrics") or {}
        design_id = design.get("design_id", "design")
        print(f"\n  {design_id}  (replicate {design.get('replicate_index')}, seed {design.get('seed')})")
        print(
            f"    residues={diagnostics.get('num_designed_residues')}"
            f"/{diagnostics.get('num_tokens')} tokens"
            f"  coords_finite={diagnostics.get('all_coordinates_finite')}"
            f"  {diagnostics.get('duration_ms', 0) / 1000:.2f}s"
        )
        sequence = design.get("sequence")
        if sequence:
            print(f"    sequence ({len(sequence)} aa): {sequence}")
        for key in HEADLINE_METRICS:
            if key in metrics and metrics[key] is not None:
                print(f"    {key.split('.')[-1]:<32} {metrics[key]:.4g}")

        structure = design.get("structure") or {}
        content = structure.get("content")
        if content:
            cif_path = out_dir / (structure.get("filename") or f"{design_id}.cif")
            cif_path.write_text(content, encoding="utf-8")
            n_ca = count_ca_atoms(content)
            expected = diagnostics.get("num_tokens")
            if n_ca is not None and expected is not None:
                verdict = "matches token count" if n_ca == expected else f"EXPECTED {expected}"
                print(f"    {'CA atoms in mmCIF':<32} {n_ca} ({verdict})")
            print(f"    structure -> {cif_path}")

        summary.append(
            {
                "design_id": design_id,
                "replicate_index": design.get("replicate_index"),
                "seed": design.get("seed"),
                "sequence": sequence,
                "diagnostics": diagnostics,
                "metrics": metrics,
            }
        )

    (out_dir / "summary.json").write_text(
        json.dumps(
            {
                "task": response.get("task"),
                "model": model,
                "schedule": schedule,
                "duration_ms": response.get("duration_ms"),
                "aggregated_metrics": response.get("aggregated_metrics"),
                "designs": summary,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\n  wrote {out_dir}/response.json and summary.json")


def on_retry(attempt: int, total: int, delay: float, error: pc.ProteinaError) -> None:
    print(f"  attempt {attempt}/{total} failed ({error}); retrying in {delay:.0f}s", flush=True)


def send(payload: dict[str, Any], args: argparse.Namespace, out_name: str) -> int:
    if args.dry_run:
        print(json.dumps(_redact_files(payload), indent=2))
        print("\n  dry run: nothing sent")
        return 0

    api_key = pc.load_api_key(args.env_file)
    print(f"  POST {pc.ENDPOINT}")
    try:
        response = pc.invoke(
            payload,
            api_key=api_key,
            timeout=args.timeout,
            max_attempts=args.max_attempts,
            on_retry=on_retry,
        )
    except pc.ProteinaError as exc:
        print(f"\n  request failed: {exc}", file=sys.stderr)
        return 1

    report(response, OUTPUT_DIR / out_name)
    return 0


def _redact_files(payload: dict[str, Any]) -> dict[str, Any]:
    """Copy of the payload with mmCIF contents replaced by their size."""
    clone = json.loads(json.dumps(payload))
    node = clone.get("input", {})
    for item in ([node.get("file")] if node.get("file") else []) + list(node.get("files") or []):
        if isinstance(item, dict) and "content" in item:
            item["content"] = f"<{len(item['content'])} chars of mmCIF>"
    return clone


def cmd_check(args: argparse.Namespace) -> int:
    api_key = pc.load_api_key(args.env_file)
    try:
        print(f"  {pc.check_credentials(api_key)}")
    except pc.ProteinaError as exc:
        print(f"  check failed: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_generate(args: argparse.Namespace) -> int:
    template = pc.inline_file(resolve_template(args))
    payload = pc.build_generation_request(
        template,
        example_id=args.name,
        replicates=args.replicates,
        seed=args.seed,
        recycle_strategy=args.recycle,
        num_steps=args.steps,
    )
    print(f"  generation from {template['filename']} "
          f"({args.replicates} replicate(s), {args.recycle})")
    return send(payload, args, args.name)


def cmd_scaffold(args: argparse.Namespace) -> int:
    template = pc.inline_file(resolve_template(args))
    chain, first, last = parse_motif(args.motif)
    contig = args.contig or f"{args.n_term}, {chain}{first}-{last}, {args.c_term}"
    payload = pc.build_scaffold_request(
        template,
        contig=contig,
        motif_select=pc.motif_selection(chain, first, last),
        example_id=args.name,
        replicates=args.replicates,
        seed=args.seed,
        recycle_strategy=args.recycle,
    )
    print(f"  motif scaffolding from {template['filename']}  contig=\"{contig}\"")
    return send(payload, args, args.name)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Demo client for NVIDIA's hosted RFD4-Proteina NIM endpoint.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--env-file", help="file holding NVIDIA_API_KEY (default: ./.env)")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common(sp: argparse.ArgumentParser) -> None:
        source = sp.add_mutually_exclusive_group()
        source.add_argument("--template", help="path to a local mmCIF template")
        source.add_argument("--pdb-id", default="1CRN", help="PDB id to fetch from RCSB (default: 1CRN)")
        sp.add_argument("--replicates", type=int, default=1, help="designs to sample (default: 1)")
        sp.add_argument("--seed", type=int, default=42, help="sampling seed (default: 42)")
        sp.add_argument(
            "--recycle",
            default="thinking_low",
            choices=pc.RECYCLE_STRATEGIES,
            help="compute budget preset (default: thinking_low)",
        )
        sp.add_argument("--timeout", type=float, default=900.0, help="per-attempt timeout in seconds")
        sp.add_argument("--max-attempts", type=int, default=6, help="attempts before giving up on 429/5xx")
        sp.add_argument("--dry-run", action="store_true", help="print the request instead of sending it")

    p_check = sub.add_parser("check", help="verify the key and endpoint without running a design")
    p_check.set_defaults(func=cmd_check)

    p_gen = sub.add_parser("generate", help="unconditional generation from an annotated mmCIF")
    add_common(p_gen)
    p_gen.add_argument("--name", default="generate_demo", help="example_id and output folder name")
    p_gen.add_argument("--steps", type=int, help="override the sampler's step count (1-1000)")
    p_gen.set_defaults(func=cmd_generate)

    p_scaf = sub.add_parser("scaffold", help="motif scaffolding via a condition spec")
    add_common(p_scaf)
    p_scaf.add_argument("--name", default="scaffold_demo", help="example_id and output folder name")
    p_scaf.add_argument("--motif", default="A1-15", help="motif to keep, e.g. A1-15 (default: A1-15)")
    p_scaf.add_argument("--n-term", type=int, default=30, help="designed residues before the motif")
    p_scaf.add_argument("--c-term", type=int, default=30, help="designed residues after the motif")
    p_scaf.add_argument("--contig", help="full contig string, overriding --motif/--n-term/--c-term")
    p_scaf.set_defaults(func=cmd_scaffold)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
