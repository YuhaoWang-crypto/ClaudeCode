"""Parse the curated aaRS x ncAA literature spreadsheet into a tidy table.

Input: `dataset_from_paper_v2.0.xlsx`, a hand-curated sheet of published
aaRS-engineering campaigns. It is a reading sheet, not a dataset: section headers
live in data rows ("figure1", "Figure7Phe衍生物"), figure sub-tables are pasted
below the records, several cells hold two concatenated protein sequences, and one
cell can list four ncAAs. This module extracts the campaign-level records.

A record is a row whose `aaRS_sequence` cell actually contains a protein sequence.
Everything else is layout.

What comes out (✅ verified on the v2.0 file): **23 campaign records across 3
scaffolds** --

    MjTyrRS  Methanococcus jannaschii TyrRS   308 aa   17 campaigns  (sheet 2)
    McTyrRS  Methanosaeta concilii TyrRS      319 aa    1 campaign   (sheet 3)
    MbPylRS  Methanosarcina barkeri PylRS     427 aa    5 campaigns  (sheet 1)

Two cells hold two sequences stuck together (737 = 308 + 429, 887 = 427 + 460);
`split_sequences` separates them on the second start-methionine-like boundary and
records both lengths rather than silently truncating.

The `chemotype` column is what routes a new target to a scaffold: an ncAA derived
from Tyr/Phe belongs with a TyrRS, one derived from Lys with PylRS. That mapping
is the single most load-bearing decision in planning a campaign, and it is made
from the ncAA's parent amino acid, not from which paper happened to be read.

Run:  python3 pylrs/literature.py --xlsx dataset_from_paper_v2.0.xlsx --out pylrs/data
"""

from __future__ import annotations

import argparse
import pathlib
import re

import pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")

COLUMNS = ["number", "ncAA", "aaRS spicies", "aaRS_sequence",
           "Library design", "Readout results", "paper"]

#: Known scaffold lengths, used to label records and to split glued sequences.
SCAFFOLDS = {308: "MjTyrRS", 309: "MjTyrRS", 319: "McTyrRS", 427: "MbPylRS", 429: "MbPylRS"}

#: Parent amino acid of an ncAA, inferred from its name. Order matters: the first
#: pattern that matches wins, so the specific cases precede the generic ones.
CHEMOTYPE_PATTERNS = [
    ("Lys", r"lysine|\blys\b|\bk\b-?derivat|pyrrolysin|acetyl-l|azido.*lysine|propargyl-l-lysine|\bprk\b|\back\b"),
    ("Tyr", r"tyrosine|\btyr\b|o-methyl|ome-tyr|phenyllactic"),
    ("Phe", r"phenylalanine|\bphe\b|naphthyl|alanine\b|quinolin|coumarin|bipyridin"),
]


def is_sequence(value) -> bool:
    """True when the cell holds a protein sequence rather than prose or a label."""
    if not isinstance(value, str) or len(value) < 80:
        return False
    letters = re.sub(r"[^A-Za-z]", "", value).upper()
    return len(letters) >= 80 and sum(c in AA for c in letters) / len(letters) > 0.95


#: An aaRS (or its catalytic domain, as curated here) falls in this length range.
#: Used to split glued cells: a candidate split is only accepted when BOTH halves
#: are plausible proteins, which is what distinguishes 427+460 from 308+579 for
#: the 887-residue cell.
PLAUSIBLE_LEN = (250, 520)


def split_sequences(seq: str) -> list[str]:
    """Split a cell that holds two concatenated sequences, else return the one."""
    clean = re.sub(r"[^A-Za-z]", "", seq).upper()
    lo, hi = PLAUSIBLE_LEN
    if len(clean) in SCAFFOLDS or len(clean) <= hi:
        return [clean]

    def plausible(n):
        return lo <= n <= hi

    # Prefer a prefix that is a known scaffold length and leaves a plausible rest.
    for n in sorted(SCAFFOLDS, reverse=True):
        if plausible(n) and plausible(len(clean) - n):
            return [clean[:n], clean[n:]]
    # Otherwise fall back to any split into two plausible halves.
    for n in range(lo, hi + 1):
        if plausible(len(clean) - n):
            return [clean[:n], clean[n:]]
    return [clean]


def chemotype(name: str) -> str:
    """Parent amino acid of an ncAA name; '?' when the name does not say."""
    low = name.lower()
    for label, pattern in CHEMOTYPE_PATTERNS:
        if re.search(pattern, low):
            return label
    return "?"


def split_ncaa(cell: str) -> list[str]:
    """One cell may list several ncAAs; split on newlines and separators."""
    parts = re.split(r"[\n;]+|(?<=\))\s*[,，]\s*", str(cell))
    return [p.strip() for p in parts if p and p.strip() and p.strip().lower() != "nan"] or [str(cell)]


def parse(xlsx: pathlib.Path) -> pd.DataFrame:
    book = pd.ExcelFile(xlsx)
    rows = []
    for sheet in book.sheet_names:
        frame = book.parse(sheet)
        for col in COLUMNS:
            if col not in frame.columns:
                frame[col] = pd.NA
        for idx, row in frame.iterrows():
            if not is_sequence(row["aaRS_sequence"]):
                continue
            seqs = split_sequences(row["aaRS_sequence"])
            lengths = [len(s) for s in seqs]
            scaffold = "/".join(dict.fromkeys(SCAFFOLDS.get(n, f"len{n}") for n in lengths))
            for ncaa in split_ncaa(row["ncAA"]):
                rows.append({
                    "sheet": sheet,
                    "record": str(row["number"]).strip(),
                    "ncAA": ncaa[:120],
                    "chemotype": chemotype(ncaa),
                    "scaffold": scaffold,
                    "seq_lengths": "+".join(str(n) for n in lengths),
                    "species_note": str(row["aaRS spicies"])[:80],
                    "library_design": str(row["Library design"])[:400],
                    "readout": str(row["Readout results"])[:400],
                    "paper": str(row["paper"])[:200],
                    "sequence": seqs[0],
                })
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args(argv)

    table = parse(args.xlsx)
    args.out.mkdir(parents=True, exist_ok=True)
    target = args.out / "literature_campaigns.csv"
    table.to_csv(target, index=False)

    print(f"{len(table)} ncAA x campaign rows from "
          f"{table['record'].nunique()} records -> {target}")
    print("\nby scaffold:")
    print(table.groupby(["scaffold", "chemotype"]).size().to_string())
    print("\nunique scaffold sequences:")
    for scaffold, grp in table.groupby("scaffold"):
        print(f"  {scaffold:<18} {grp['sequence'].nunique()} distinct sequence(s), "
              f"{grp['record'].nunique()} records")
    unknown = table[table["chemotype"] == "?"]
    if len(unknown):
        print(f"\n⚠️  {len(unknown)} rows whose parent amino acid the name does not "
              "state -- classify these by hand before routing:")
        for name in unknown["ncAA"].unique()[:12]:
            print(f"     {name[:90]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
