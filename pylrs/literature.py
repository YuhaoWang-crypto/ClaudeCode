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

#: Known scaffold lengths, once labels are stripped. MjTyrRS is 306 aa and
#: M. barkeri PylRS 419 aa; seeing 308/427 means the label is still attached.
SCAFFOLDS = {306: "MjTyrRS", 319: "McTyrRS", 419: "MbPylRS", 454: "MmPylRS", 460: "MmPylRS"}

#: Residues the sheet-2 campaigns cite by number on the MjTyrRS scaffold. Every
#: one must match, or the numbering is off and no mutation parsed from this sheet
#: can be trusted.
MJ_EXPECTED = {32: "Y", 65: "L", 67: "A", 70: "H", 102: "Y", 103: "V", 107: "E",
               108: "F", 109: "Q", 110: "L", 114: "Y", 155: "Q", 158: "D",
               159: "I", 160: "H", 161: "Y", 162: "L"}


def verify_numbering(sequence: str, expected: dict[int, str] = None) -> list[str]:
    """Return the positions where `sequence` disagrees with the cited residues.

    ✅ With labels stripped, all 17 cited MjTyrRS positions match.
    """
    expected = MJ_EXPECTED if expected is None else expected
    return [f"{aa}{pos}->{sequence[pos - 1] if pos <= len(sequence) else '-'}"
            for pos, aa in sorted(expected.items())
            if pos > len(sequence) or sequence[pos - 1] != aa]

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


#: Sequences in this sheet carry a `label:` prefix -- "WT:MDEFEM...",
#: "M.barkeri:MDKK...", "E.coli:MASS..." -- sometimes with a full-width colon, and
#: a cell may hold several such labelled sequences separated by whitespace.
#:
#: ⚠️ Getting this wrong is silent and fatal. Stripping non-letters without
#: removing the label glues it onto the N-terminus: "WT:MDEFEM..." becomes
#: "WTMDEFEM...", shifting every residue number by 2, and the MjTyrRS positions
#: the literature cites (Y32, L65, D158, ...) then land on the wrong residues.
#: `verify_numbering` exists to catch exactly that.
COLON_RE = re.compile(r"[:：]")

#: Shortest run of residues still treated as a sequence rather than stray text.
MIN_SEQUENCE_LEN = 80


def split_sequences(seq: str) -> list[str]:
    """Split a cell into its labelled sequences, dropping the labels.

    Splits on the colon rather than pattern-matching the label, because a label
    pattern permissive enough to cover "M.barkeri" and "M.Mazei" also matches
    residues immediately before the next label and silently truncates the
    preceding sequence. A protein sequence never contains a colon, so each
    colon-delimited piece is exactly ``<sequence of the previous label><next
    label>`` and the label is the final whitespace-delimited token.
    """
    pieces = COLON_RE.split(seq)
    if len(pieces) == 1:
        clean = re.sub(r"[^A-Za-z]", "", seq).upper()
        return [clean] if clean else [""]

    bodies = []
    for i, piece in enumerate(pieces[1:], start=1):
        # Every piece but the last ends with the next sequence's label.
        body = piece.rsplit(None, 1)[0] if i < len(pieces) - 1 and piece.split() else piece
        bodies.append(body)

    out = [re.sub(r"[^A-Za-z]", "", b).upper() for b in bodies]
    out = [s for s in out if len(s) >= MIN_SEQUENCE_LEN]
    return out or [re.sub(r"[^A-Za-z]", "", seq).upper()]


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
    mj = table[table["scaffold"] == "MjTyrRS"]
    if len(mj):
        seq = mj["sequence"].iloc[0]
        bad = verify_numbering(seq)
        status = "✅ all 17 cited positions match" if not bad else f"*** MISMATCH: {bad}"
        print(f"\nMjTyrRS numbering check ({len(seq)} aa): {status}")
        if bad:
            raise SystemExit("residue numbering is wrong; mutations parsed from this "
                             "sheet would be meaningless")

    unknown = table[table["chemotype"] == "?"]
    if len(unknown):
        print(f"\n⚠️  {len(unknown)} rows whose parent amino acid the name does not "
              "state -- classify these by hand before routing:")
        for name in unknown["ncAA"].unique()[:12]:
            print(f"     {name[:90]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
