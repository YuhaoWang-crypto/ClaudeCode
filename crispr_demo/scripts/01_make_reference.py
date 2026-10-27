#!/usr/bin/env python3
"""
01_make_reference.py - build the mini bacterial reference genome for the demo.

Outputs:
  ref/DEMO_chr1.fasta     reference sequence (60 kb)
  ref/DEMO_chr1.gff3      gene annotation (gene / mRNA / exon / CDS) for snpEff
  ref/sgRNA.tsv           guide RNAs used in the editing experiment
  ref/design_notes.json   design metadata

Notes:
  * every gene is a clean ORF: ATG ... stop, no internal stop codon
  * near-cognate sgRNA sites (2-4 mismatches + NGG PAM) are planted on purpose
    so the demo has realistic off-target hotspots
  * fixed random seed, fully reproducible
"""

import json
import random
from pathlib import Path

SEED = 20251027
GENOME_LEN = 60000
CHROM = "DEMO_chr1"
ORGANISM = "Escherichia coli str. DEMO-WT"

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "ref"
REF.mkdir(parents=True, exist_ok=True)

rng = random.Random(SEED)

STOPS = ("TAA", "TAG", "TGA")
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def random_seq(n, gc=0.507):
    out = []
    for _ in range(n):
        out.append(rng.choice("GC") if rng.random() < gc else rng.choice("AT"))
    return "".join(out)


def random_codon():
    while True:
        c = random_seq(3, 0.52)
        if c not in STOPS and c != "ATG":
            return c


def make_orf(length_bp):
    n_internal = length_bp // 3 - 2
    return "ATG" + "".join(random_codon() for _ in range(n_internal)) + rng.choice(STOPS)


GENE_NAMES = [
    "thrA", "carB", "dapB", "leuD", "ribF", "fhuA", "dnaK", "aceE", "sucB", "gltA",
    "pheA", "mdlA", "ompF", "nuoF", "ftsZ", "murE", "glnS", "tyrR", "aroG", "pyrD",
    "lacZ", "cysK", "metL", "purF", "glyA", "trpE", "serC", "hisD", "argB", "ilvC",
    "rpoB", "recA", "lexA", "uvrB", "yjgF", "sodB", "katG", "fabI", "acpP", "tolC",
    "ptsG", "crr", "cyaA", "flhD", "motA", "cheY", "fliC", "ompA", "lpp", "slp",
]

# genes that must sit at a known coordinate: name -> (start_1based, length_bp, strand)
ANCHORED = {
    "pheA": (12001, 1167, "+"),
    "lacZ": (26001, 1443, "+"),
    "yjgF": (38001, 519, "+"),
    "rpoB": (49001, 1326, "+"),
}


def build_layout():
    """Place anchored genes first, then fill the gaps with filler genes."""
    placed = []
    for name, (start, length, strand) in ANCHORED.items():
        placed.append({"name": name, "start": start, "end": start + length - 1,
                       "strand": strand, "anchored": True})
    placed.sort(key=lambda g: g["start"])

    filler_pool = [n for n in GENE_NAMES if n not in ANCHORED]
    rng.shuffle(filler_pool)
    fi = 0

    # gaps available for filler genes: before / between / after anchored genes
    gaps = []
    cursor = 500
    for g in placed:
        gaps.append((cursor, g["start"] - 250))
        cursor = g["end"] + 250
    gaps.append((cursor, GENOME_LEN - 500))

    filled = []
    for lo, hi in gaps:
        pos = lo
        while fi < len(filler_pool):
            length = rng.randrange(60, 160) * 9  # 540..1431 bp, multiple of 3
            gap = rng.randrange(80, 320)
            if pos + length - 1 > hi:
                break
            filled.append({"name": filler_pool[fi], "start": pos,
                           "end": pos + length - 1,
                           "strand": rng.choice(["+", "-"]), "anchored": False})
            fi += 1
            pos += length + gap

    genes = sorted(placed + filled, key=lambda g: g["start"])
    for i, g in enumerate(genes, start=1):
        g["locus_tag"] = "DEMO_%04d" % (i * 10)
    return genes


def build_sequence(genes):
    """Concatenate intergenic filler + ORFs into one chromosome string."""
    seq = list(random_seq(GENOME_LEN))
    for g in genes:
        length = g["end"] - g["start"] + 1
        orf = make_orf(length)
        if g["strand"] == "-":
            orf = rc(orf)
        seq[g["start"] - 1:g["end"]] = list(orf)
    return "".join(seq)


def cds_of(seq, gene):
    s = seq[gene["start"] - 1:gene["end"]]
    return s if gene["strand"] == "+" else rc(s)


def has_internal_stop(cds):
    for i in range(0, len(cds) - 3, 3):
        if cds[i:i + 3] in STOPS:
            return True
    return False


SWAP = {"A": "C", "C": "A", "G": "C", "T": "C"}


def repair_orf(seq_list, gene, protect=()):
    """After planting a site inside a CDS, silence any new internal stop codon.

    `protect` holds 0-based genome positions that must not be touched
    (sgRNA protospacer / PAM bases).
    """
    protect = set(protect)
    for _ in range(200):
        cds = cds_of("".join(seq_list), gene)
        bad = None
        for i in range(0, len(cds) - 3, 3):
            if cds[i:i + 3] in STOPS:
                bad = i
                break
        if bad is None:
            return True
        # try the wobble base first, then the other two positions of the codon
        for off in (2, 1, 0):
            cds_pos = bad + off
            if gene["strand"] == "+":
                gpos = gene["start"] - 1 + cds_pos
            else:
                gpos = gene["end"] - 1 - cds_pos
            if gpos in protect:
                continue
            cur = seq_list[gpos]
            seq_list[gpos] = SWAP[cur] if gene["strand"] == "+" else rc(SWAP[rc(cur)])
            break
        else:
            return False
    return False


def main():
    genes = build_layout()
    seq = build_sequence(genes)
    by_name = {g["name"]: g for g in genes}
    seq_list = list(seq)

    # -----------------------------------------------------------------
    # sgRNA-1: protospacer inside pheA, PAM = NGG immediately 3'
    # -----------------------------------------------------------------
    pheA = by_name["pheA"]
    # cut site roughly in the middle of the gene, keep it in frame-agnostic place
    ps1_start = pheA["start"] + 600          # 1-based start of 20 nt protospacer
    # force an NGG PAM right after the protospacer, repair the ORF afterwards
    # 0-based: protospacer = [ps1_start-1 .. ps1_start+18], PAM = [ps1_start+19 .. +21]
    seq_list[ps1_start + 19] = "A"
    seq_list[ps1_start + 20] = "G"
    seq_list[ps1_start + 21] = "G"
    prot1 = set(range(ps1_start - 1, ps1_start + 22))
    repair_orf(seq_list, pheA, protect=prot1)
    seq = "".join(seq_list)
    sg1 = seq[ps1_start - 1:ps1_start + 19]
    pam1 = seq[ps1_start + 19:ps1_start + 22]

    # -----------------------------------------------------------------
    # sgRNA-2: protospacer inside lacZ
    # -----------------------------------------------------------------
    lacZ = by_name["lacZ"]
    ps2_start = lacZ["start"] + 300
    seq_list = list(seq)
    seq_list[ps2_start + 19] = "T"
    seq_list[ps2_start + 20] = "G"
    seq_list[ps2_start + 21] = "G"
    prot2 = set(range(ps2_start - 1, ps2_start + 22))
    repair_orf(seq_list, lacZ, protect=prot2)
    seq = "".join(seq_list)
    sg2 = seq[ps2_start - 1:ps2_start + 19]
    pam2 = seq[ps2_start + 19:ps2_start + 22]

    # -----------------------------------------------------------------
    # plant near-cognate (off-target) sites
    #   OT-1: 3 mismatches vs sgRNA-1, inside yjgF
    #   OT-2: 4 mismatches vs sgRNA-1, intergenic
    #   OT-3: 2 mismatches vs sgRNA-2, inside nuoF  (decoy: no variant there)
    # -----------------------------------------------------------------
    def plant(site_start, guide, mismatch_positions, pam, gene=None, strand="+"):
        """Write guide-with-mismatches + PAM at site_start (1-based)."""
        proto = list(guide)
        for p in mismatch_positions:
            cur = proto[p]
            proto[p] = rng.choice([b for b in "ACGT" if b != cur])
        block = "".join(proto) + pam
        if strand == "-":
            block = rc(block)
        sl = list(seq)
        sl[site_start - 1:site_start - 1 + len(block)] = list(block)
        if gene is not None:
            repair_orf(sl, gene,
                       protect=set(range(site_start - 1, site_start - 1 + len(block))))
        return "".join(sl), "".join(proto)

    yjgF = by_name["yjgF"]
    ot1_start = yjgF["start"] + 210
    seq, ot1_proto = plant(ot1_start, sg1, [1, 4, 9], "AGG", gene=yjgF)

    ot2_start = 33200  # intergenic
    seq, ot2_proto = plant(ot2_start, sg1, [0, 2, 6, 11], "CGG")

    nuoF = by_name["nuoF"]
    ot3_start = nuoF["start"] + 150
    seq, ot3_proto = plant(ot3_start, sg2, [3, 8], "TGG", gene=nuoF,
                           strand=nuoF["strand"])

    # -----------------------------------------------------------------
    # self-check: ORFs intact, PAMs are NGG, planted sites have the
    # intended number of mismatches
    # -----------------------------------------------------------------
    broken = [g["name"] for g in genes if has_internal_stop(cds_of(seq, g))]
    assert not broken, "ORFs with internal stop: %s" % broken

    sg1 = seq[ps1_start - 1:ps1_start + 19]
    pam1 = seq[ps1_start + 19:ps1_start + 22]
    sg2 = seq[ps2_start - 1:ps2_start + 19]
    pam2 = seq[ps2_start + 19:ps2_start + 22]
    for gid, pam in (("sgRNA-1", pam1), ("sgRNA-2", pam2)):
        assert pam[1:] == "GG", "%s PAM is %s, expected NGG" % (gid, pam)

    def mismatches(a, b):
        return sum(1 for x, y in zip(a, b) if x != y)

    for label, start, guide, want_mm, strand in (
            ("OT-1", ot1_start, sg1, 3, "+"),
            ("OT-2", ot2_start, sg1, 4, "+"),
            ("OT-3", ot3_start, sg2, 2, nuoF["strand"])):
        block = seq[start - 1:start + 22]
        if strand == "-":
            block = rc(block)
        proto, pam = block[:20], block[20:23]
        got = mismatches(proto, guide)
        assert pam[1:] == "GG", "%s PAM is %s, expected NGG" % (label, pam)
        assert got == want_mm, "%s has %d mismatches, expected %d" % (label, got, want_mm)

    # -----------------------------------------------------------------
    # write FASTA
    # -----------------------------------------------------------------
    with open(REF / "DEMO_chr1.fasta", "w") as fh:
        fh.write(">%s  %s  mini-chromosome, %d bp\n" % (CHROM, ORGANISM, len(seq)))
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")

    # -----------------------------------------------------------------
    # write GFF3 (with embedded FASTA so snpEff can build a database)
    # -----------------------------------------------------------------
    with open(REF / "DEMO_chr1.gff3", "w") as fh:
        fh.write("##gff-version 3\n")
        fh.write("##sequence-region %s 1 %d\n" % (CHROM, len(seq)))
        for g in genes:
            tag, name = g["locus_tag"], g["name"]
            attrs_gene = "ID=gene-%s;Name=%s;locus_tag=%s;gene_biotype=protein_coding" % (tag, name, tag)
            attrs_mrna = "ID=rna-%s;Parent=gene-%s;Name=%s;locus_tag=%s" % (tag, tag, name, tag)
            attrs_exon = "ID=exon-%s;Parent=rna-%s;Name=%s;locus_tag=%s" % (tag, tag, name, tag)
            attrs_cds = "ID=cds-%s;Parent=rna-%s;Name=%s;locus_tag=%s;protein_id=%s.1" % (tag, tag, name, tag, tag)
            for feat, attrs, phase in (("gene", attrs_gene, "."),
                                       ("mRNA", attrs_mrna, "."),
                                       ("exon", attrs_exon, "."),
                                       ("CDS", attrs_cds, "0")):
                fh.write("\t".join([CHROM, "demo", feat, str(g["start"]), str(g["end"]),
                                    ".", g["strand"], phase, attrs]) + "\n")
        fh.write("##FASTA\n>%s\n" % CHROM)
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")

    # -----------------------------------------------------------------
    # write sgRNA table
    #   blunt Cas9 cut = between protospacer position 17 and 18 (3 bp from PAM)
    # -----------------------------------------------------------------
    guides = [
        {"id": "sgRNA-1", "target_gene": "pheA", "protospacer": sg1, "pam": pam1,
         "strand": "+", "proto_start": ps1_start, "cut_site": ps1_start + 16,
         "purpose": "precise 24 bp knock-in (FLAG-like tag + stop)"},
        {"id": "sgRNA-2", "target_gene": "lacZ", "protospacer": sg2, "pam": pam2,
         "strand": "+", "proto_start": ps2_start, "cut_site": ps2_start + 16,
         "purpose": "knockout by NHEJ/MMEJ indel"},
    ]
    with open(REF / "sgRNA.tsv", "w") as fh:
        fh.write("sgrna_id\ttarget_gene\tprotospacer_20nt\tpam\tstrand\tproto_start\tcut_site\tpurpose\n")
        for g in guides:
            fh.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%s\n" % (
                g["id"], g["target_gene"], g["protospacer"], g["pam"], g["strand"],
                g["proto_start"], g["cut_site"], g["purpose"]))

    notes = {
        "seed": SEED, "chrom": CHROM, "length": len(seq), "organism": ORGANISM,
        "n_genes": len(genes),
        "gc_percent": round(100.0 * sum(seq.count(b) for b in "GC") / len(seq), 2),
        "guides": guides,
        "planted_offtarget_sites": [
            {"name": "OT-1", "start": ot1_start, "guide": "sgRNA-1", "mismatches": 3,
             "locus": "yjgF (CDS)", "sequence": ot1_proto},
            {"name": "OT-2", "start": ot2_start, "guide": "sgRNA-1", "mismatches": 4,
             "locus": "intergenic", "sequence": ot2_proto},
            {"name": "OT-3", "start": ot3_start, "guide": "sgRNA-2", "mismatches": 2,
             "locus": "nuoF (CDS)", "sequence": ot3_proto,
             "comment": "decoy: predicted hotspot but NOT edited in this demo"},
        ],
        "genes": genes,
        "orfs_with_internal_stop": broken,
    }
    with open(REF / "design_notes.json", "w") as fh:
        json.dump(notes, fh, indent=2)

    print("reference : %s (%d bp, GC %.2f%%)" % (CHROM, len(seq), notes["gc_percent"]))
    print("genes     : %d" % len(genes))
    print("sgRNA-1   : %s %s  (pheA, cut @ %d)" % (sg1, pam1, guides[0]["cut_site"]))
    print("sgRNA-2   : %s %s  (lacZ, cut @ %d)" % (sg2, pam2, guides[1]["cut_site"]))
    print("planted OT: OT-1 @%d (3mm), OT-2 @%d (4mm), OT-3 @%d (2mm, decoy)"
          % (ot1_start, ot2_start, ot3_start))
    print("broken ORFs: %s" % (broken or "none"))


if __name__ == "__main__":
    main()
