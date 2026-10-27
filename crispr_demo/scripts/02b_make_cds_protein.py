#!/usr/bin/env python3
"""
02b_make_cds_protein.py - emit cds.fa / protein.fa so that `snpEff build`
can verify the annotation database instead of skipping the check.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "ref"
OUT = ROOT / "snpeff" / "data" / "DEMO"
OUT.mkdir(parents=True, exist_ok=True)

COMP = str.maketrans("ACGTN", "TGCAN")

# NCBI translation table 11 (bacterial), same table snpEff is told to use
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON_TABLE[b1 + b2 + b3] = AA[i]
            i += 1


def rc(s):
    return s.translate(COMP)[::-1]


def translate(cds):
    prot = []
    for i in range(0, len(cds) - 2, 3):
        prot.append(CODON_TABLE.get(cds[i:i + 3], "X"))
    # bacterial alternative start codons all read as Met
    if prot and cds[:3] in ("TTG", "CTG", "ATT", "ATC", "ATA", "GTG"):
        prot[0] = "M"
    if prot and prot[-1] == "*":
        prot.pop()
    return "".join(prot)


def read_fasta(path):
    chunks = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                chunks.append(line.strip())
    return "".join(chunks)


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def main():
    genome = read_fasta(REF / "DEMO_chr1.fasta")
    records = []
    with open(REF / "DEMO_chr1.gff3") as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            start, end, strand = int(f[3]), int(f[4]), f[6]
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            tx_id = attrs["Parent"]
            seq = genome[start - 1:end]
            if strand == "-":
                seq = rc(seq)
            records.append((tx_id, seq))

    with open(OUT / "cds.fa", "w") as fh:
        for tx, seq in records:
            fh.write(">%s\n%s\n" % (tx, wrap(seq)))
    with open(OUT / "protein.fa", "w") as fh:
        for tx, seq in records:
            fh.write(">%s\n%s\n" % (tx, wrap(translate(seq))))

    internal_stops = sum(1 for _, s in records if "*" in translate(s))
    print("wrote %d CDS / protein records -> %s" % (len(records), OUT))
    print("records with an internal stop codon: %d" % internal_stops)


if __name__ == "__main__":
    main()
