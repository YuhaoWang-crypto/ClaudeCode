#!/usr/bin/env python3
"""Regenerate `fftplsr/data/` from a checkout of the paper's own repository.

The vendored CSV/FASTA files are format conversions only; this script is the
record of how. It is not needed to use `fftplsr` -- run it only to re-derive the
data or to check it against the upstream source.

    git clone https://github.com/zjuhaoran/FPFORCOM /tmp/FPFORCOM
    python3 scripts/convert_paper_data.py /tmp/FPFORCOM

Pass `--check` to verify the committed files match what the source produces,
without writing anything.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd

OUT = pathlib.Path(__file__).resolve().parents[1] / "fftplsr" / "data"

SEQUENCES = [("IFRS", "input/seq_IFRS.txt"), ("Com1-IFRS", "input/seq_Com1.txt")]
SHEETS = [
    ("Trainset1", "trainset1_ifrs_singles"),
    ("Trainset2", "trainset2_ifrs_combos"),
    ("Trainset3", "trainset3_com1_singles"),
    ("Trainset4", "trainset4_com1_combos"),
]
PANELS = [
    ("data/IFRS-mutants-data.json", "measured_ifrs_panel"),
    ("data/Com1-IFRS-mutants-data.json", "measured_com1_panel"),
]


def read_sequence(path: pathlib.Path) -> str:
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    )


def as_fasta(name: str, seq: str) -> str:
    body = "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60))
    return f">{name}\n{body}\n"


def panel_table(path: pathlib.Path) -> pd.DataFrame:
    """Pull `(variant name, relative activity)` out of an EnzymeML record.

    Each measured variant is a `proteins[]` entry; its measurement is the matching
    `species_data[]` entry's `initial` field. Entries without a measurement (the
    sfGFP reporter) drop out.
    """
    doc = json.loads(path.read_text())
    names = {protein["id"]: protein["name"] for protein in doc["proteins"]}
    rows = [
        {"Variants": names[entry["species_id"]], "Fitness": float(entry["initial"])}
        for entry in doc["measurements"][0]["species_data"]
        if entry["species_id"] in names and entry.get("initial") is not None
    ]
    return pd.DataFrame(rows).drop_duplicates("Variants").reset_index(drop=True)


def build(source: pathlib.Path) -> dict[str, str]:
    """Return ``{filename: file contents}`` for every vendored artefact."""
    produced: dict[str, str] = {}

    for name, rel in SEQUENCES:
        seq = read_sequence(source / rel)
        if len(seq) != 454:
            raise SystemExit(f"{rel}: expected a 454-residue sequence, got {len(seq)}")
        produced[f"{name}.fasta"] = as_fasta(name, seq)

    book = pd.ExcelFile(source / "input" / "data.xlsx")
    for sheet, stem in SHEETS:
        frame = book.parse(sheet).rename(columns={"Unnamed: 2": "Source"})
        frame = frame[frame["Variants"].notna()]
        produced[f"{stem}.csv"] = frame.to_csv(index=False)

    for rel, stem in PANELS:
        produced[f"{stem}.csv"] = panel_table(source / rel).to_csv(index=False)

    return produced


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=pathlib.Path, help="checkout of zjuhaoran/FPFORCOM")
    parser.add_argument("--check", action="store_true", help="compare instead of writing")
    args = parser.parse_args(argv)

    if not (args.source / "input" / "data.xlsx").exists():
        raise SystemExit(f"{args.source} does not look like a FPFORCOM checkout")

    produced = build(args.source)
    mismatched = []
    for filename, content in sorted(produced.items()):
        target = OUT / filename
        if args.check:
            current = target.read_text() if target.exists() else None
            status = "ok" if current == content else "DIFFERS"
            if status != "ok":
                mismatched.append(filename)
            print(f"  {status:8} {filename}")
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(content)
            print(f"  wrote    {filename} ({len(content.splitlines())} lines)")

    if args.check and mismatched:
        print(f"\n{len(mismatched)} file(s) differ from the source: {', '.join(mismatched)}")
        return 1
    print("\nall files match the source" if args.check else "\ndone")
    return 0


if __name__ == "__main__":
    sys.exit(main())
