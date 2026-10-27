#!/usr/bin/env python3
"""
02_simulate_reads.py - apply the editing outcomes and simulate Illumina PE150 reads.

Two samples:
  WT_parent  - the un-edited parental clone that went into the experiment
  EDIT_pool  - the edited population recovered after CRISPR/Cas9 + HDR donor

Both samples carry the same 4 "parental" differences relative to the reference
genome (the lab strain is not identical to the reference).  Only the EDIT_pool
carries the editing outcomes.  That is the whole point of the comparison:
a variant being present in EDIT_pool is NOT enough - it has to be absent from
WT_parent before it can be called an edit.

Outputs:
  rawdata/WT_parent_R1.fastq.gz  / _R2.fastq.gz
  rawdata/EDIT_pool_R1.fastq.gz  / _R2.fastq.gz
  ref/truth_variants.tsv         ground truth (used only to score the demo)
"""

import gzip
import json
import random
from pathlib import Path

SEED = 42
READ_LEN = 150
ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "ref"
RAW = ROOT / "rawdata"
RAW.mkdir(parents=True, exist_ok=True)

COMP = str.maketrans("ACGTN", "TGCAN")
# Illumina TruSeq adapter, read-through into the adapter happens on short inserts
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

rng = random.Random(SEED)


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                chunks.append(line.strip())
    return name, "".join(chunks)


# ---------------------------------------------------------------------------
# The editing outcomes, in REFERENCE coordinates.
#   kind: snv | ins | del
#   pos : 1-based reference position
#     snv -> base at pos is replaced
#     ins -> inserted sequence goes AFTER pos
#     del -> `length` bases starting AT pos are removed
# ---------------------------------------------------------------------------
HDR_CASSETTE = "G" + "TAA" + "GCTAGCGGATCCAAGCTTGA"   # 24 bp: stop + cloning scar

PARENTAL = [
    dict(id="P1", kind="snv", pos=5500, alt="C",
         note="parental SNV, intergenic"),
    dict(id="P2", kind="snv", pos=21330, alt="G",
         note="parental SNV, inside a CDS"),
    dict(id="P3", kind="snv", pos=44120, alt="C",
         note="parental SNV, inside a CDS"),
    dict(id="P4", kind="del", pos=56500, length=3,
         note="parental 3 bp deletion, intergenic"),
]

EDITS = [
    dict(id="A", kind="ins", pos=12617, seq=HDR_CASSETTE, vaf=1.00,
         klass="on-target", guide="sgRNA-1", gene="pheA",
         note="precise HDR knock-in of a 24 bp stop cassette at the sgRNA-1 cut site"),
    dict(id="B", kind="del", pos=26314, length=7, vaf=0.55,
         klass="on-target", guide="sgRNA-2", gene="lacZ",
         note="7 bp deletion at the sgRNA-2 cut site (NHEJ/MMEJ scar) -> frameshift"),
    dict(id="E", kind="del", pos=25900, length=1150, vaf=0.45,
         klass="on-target", guide="sgRNA-2", gene="lacZ",
         note="1150 bp deletion spanning the sgRNA-2 cut site (large on-target deletion)"),
    dict(id="C", kind="ins", pos=38227, seq="T", vaf=1.00,
         klass="off-target", guide="sgRNA-1", gene="yjgF",
         note="1 bp insertion at off-target OT-1 (3 mismatches) -> frameshift"),
    dict(id="D", kind="del", pos=33214, length=5, vaf=0.30,
         klass="off-target", guide="sgRNA-1", gene="intergenic",
         note="5 bp deletion at off-target OT-2 (4 mismatches), sub-clonal"),
    dict(id="F", kind="snv", pos=49656, alt="A", vaf=1.00,
         klass="background", guide=None, gene="rpoB",
         note="spontaneous point mutation, no sgRNA homology nearby"),
]

# B and E are mutually exclusive (same locus, different sub-populations);
# D is independent.  Everything else is fixed in the whole population.
HAPLOTYPES = [
    dict(name="H1", variants=["A", "C", "F", "B", "D"], freq=0.55 * 0.30),
    dict(name="H2", variants=["A", "C", "F", "B"],      freq=0.55 * 0.70),
    dict(name="H3", variants=["A", "C", "F", "E", "D"], freq=0.45 * 0.30),
    dict(name="H4", variants=["A", "C", "F", "E"],      freq=0.45 * 0.70),
]


def apply_variants(ref, variants):
    """Apply a list of variants (reference coordinates) to build a haplotype."""
    seq = ref
    for v in sorted(variants, key=lambda x: x["pos"], reverse=True):
        p = v["pos"]
        if v["kind"] == "snv":
            assert seq[p - 1] != v["alt"], "%s is not a change" % v["id"]
            seq = seq[:p - 1] + v["alt"] + seq[p:]
        elif v["kind"] == "ins":
            seq = seq[:p] + v["seq"] + seq[p:]
        elif v["kind"] == "del":
            seq = seq[:p - 1] + seq[p - 1 + v["length"]:]
        else:
            raise ValueError(v["kind"])
    return seq


def quality_profile(read_len, is_r2):
    """Illumina-like declining quality: high at the start, degrading at the 3' end."""
    base = 37.0 if not is_r2 else 35.5
    qs = []
    for i in range(read_len):
        q = base - 0.045 * i - 0.00013 * i * i
        q += rng.gauss(0, 1.6)
        qs.append(max(2, min(41, int(round(q)))))
    return qs


def sequence_read(template, is_r2, defect=None):
    """Add base-calling errors consistent with the per-base quality string.

    `defect` injects the two artefacts every real Illumina run contains, so the
    QC step has something to actually remove:
      "lowqual" - a cluster that never resolved; the whole read is near-noise
      "polyg"   - the two-colour chemistry reads "no signal" as G, so a dead
                  cluster produces a poly-G tail at high nominal quality
    """
    qs = quality_profile(len(template), is_r2)
    bases = list(template)
    if defect == "lowqual":
        qs = [max(2, min(14, int(round(rng.gauss(9, 3))))) for _ in qs]
    elif defect == "polyg":
        tail = rng.randrange(25, 70)
        start = len(bases) - tail
        for i in range(start, len(bases)):
            bases[i] = "G"
            qs[i] = rng.randrange(11, 20)
    for i, q in enumerate(qs):
        p_err = 10 ** (-q / 10.0)
        if rng.random() < p_err:
            cur = bases[i]
            if rng.random() < 0.02:
                bases[i] = "N"
                qs[i] = 2
            elif cur in "ACGT":
                bases[i] = rng.choice([b for b in "ACGT" if b != cur])
    qual = "".join(chr(33 + q) for q in qs)
    return "".join(bases), qual


def simulate(sample, haplotypes, depth, flowcell, lane, barcode, out_prefix):
    """Write a pair of gzipped FASTQ files for one sample."""
    names = [h["name"] for h in haplotypes]
    weights = [h["freq"] for h in haplotypes]
    seqs = [h["seq"] for h in haplotypes]

    mean_len = sum(len(s) * w for s, w in zip(seqs, weights))
    n_pairs = int(mean_len * depth / (2 * READ_LEN))

    r1_path = RAW / ("%s_R1.fastq.gz" % out_prefix)
    r2_path = RAW / ("%s_R2.fastq.gz" % out_prefix)

    stats = dict(pairs=0, duplicates=0, adapter_readthrough=0,
                 lowqual_reads=0, polyg_reads=0,
                 haplotype_counts={n: 0 for n in names})
    prev_fragment = None

    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        i = 0
        while i < n_pairs:
            # ~4% PCR/optical duplicates: re-sequence the very same fragment
            is_dup = prev_fragment is not None and rng.random() < 0.04
            if is_dup:
                frag, hname = prev_fragment
                stats["duplicates"] += 1
            else:
                hname = rng.choices(names, weights=weights, k=1)[0]
                hap = seqs[names.index(hname)]
                flen = int(rng.gauss(380, 65))
                flen = max(170, min(700, flen))
                # 1.5% of inserts are shorter than the read -> adapter read-through
                if rng.random() < 0.015:
                    flen = rng.randrange(90, 145)
                start = rng.randrange(0, len(hap) - flen)
                frag = hap[start:start + flen]
                prev_fragment = (frag, hname)

            stats["haplotype_counts"][hname] += 1

            t1 = frag[:READ_LEN]
            t2 = rc(frag)[:READ_LEN]
            if len(t1) < READ_LEN:
                stats["adapter_readthrough"] += 1
                t1 = (t1 + ADAPTER_R1)[:READ_LEN]
                t2 = (t2 + ADAPTER_R2)[:READ_LEN]

            # inject run artefacts: ~0.9% dead clusters, ~1.3% poly-G tails
            r = rng.random()
            if r < 0.009:
                d1 = d2 = "lowqual"
                stats["lowqual_reads"] += 2
            elif r < 0.022:
                d1, d2 = None, "polyg"
                stats["polyg_reads"] += 1
            else:
                d1 = d2 = None

            s1, q1 = sequence_read(t1, False, d1)
            s2, q2 = sequence_read(t2, True, d2)

            tile = 1101 + (i // 4000)
            x, y = rng.randrange(1000, 25000), rng.randrange(1000, 25000)
            rid = "%s:1:%s:%d:%d:%d:%d" % (sample, flowcell, lane, tile, x, y)
            f1.write("@%s 1:N:0:%s\n%s\n+\n%s\n" % (rid, barcode, s1, q1))
            f2.write("@%s 2:N:0:%s\n%s\n+\n%s\n" % (rid, barcode, s2, q2))
            stats["pairs"] += 1
            i += 1

    return stats, r1_path, r2_path


def main():
    chrom, ref = read_fasta(REF / "DEMO_chr1.fasta")
    by_id = {v["id"]: v for v in EDITS}

    # ---- WT parent: reference + parental differences only -----------------
    wt_hap = apply_variants(ref, PARENTAL)
    wt_haps = [dict(name="WT", freq=1.0, seq=wt_hap)]

    # ---- Edited pool: parental differences + editing outcomes -------------
    edit_haps = []
    for h in HAPLOTYPES:
        vs = PARENTAL + [by_id[i] for i in h["variants"]]
        edit_haps.append(dict(name=h["name"], freq=h["freq"],
                              seq=apply_variants(ref, vs)))
    total = sum(h["freq"] for h in edit_haps)
    assert abs(total - 1.0) < 1e-9, "haplotype frequencies sum to %f" % total

    print("reference      : %s  %d bp" % (chrom, len(ref)))
    print("WT_parent      : 1 haplotype, %d bp" % len(wt_hap))
    for h in edit_haps:
        print("EDIT_pool %-3s: freq %.3f, %d bp" % (h["name"], h["freq"], len(h["seq"])))

    # Depth matters here: a sub-clonal off-target event at ~30% VAF is only
    # estimable to about +/-0.07 even at 150x, so typical 50-80x bacterial WGS
    # cannot pin down sub-clonal editing outcomes.  See the report's "limits"
    # section - amplicon deep sequencing of the target site is the usual
    # complement to WGS for exactly this reason.
    st_wt, r1w, r2w = simulate("WT", wt_haps, 120, "HXXKMDRXY", 1, "ATCACG", "WT_parent")
    st_ed, r1e, r2e = simulate("ED", edit_haps, 150, "HXXKMDRXY", 2, "CGATGT", "EDIT_pool")

    for label, st, r1 in (("WT_parent", st_wt, r1w), ("EDIT_pool", st_ed, r1e)):
        print("%s: %d pairs, %d dup, %d adapter read-through, %d low-qual reads, %d poly-G reads"
              % (label, st["pairs"], st["duplicates"], st["adapter_readthrough"],
                 st["lowqual_reads"], st["polyg_reads"]))

    # ---- ground truth table ----------------------------------------------
    rows = []
    for v in PARENTAL:
        rows.append(dict(id=v["id"], klass="parental (shared)", pos=v["pos"],
                         kind=v["kind"], detail=v.get("alt") or v.get("seq") or
                         ("%d bp" % v.get("length", 0)),
                         expected_vaf_wt=1.0, expected_vaf_edit=1.0,
                         guide="-", gene="-", note=v["note"]))
    for v in EDITS:
        rows.append(dict(id=v["id"], klass=v["klass"], pos=v["pos"], kind=v["kind"],
                         detail=v.get("alt") or v.get("seq") or ("%d bp" % v.get("length", 0)),
                         expected_vaf_wt=0.0, expected_vaf_edit=v["vaf"],
                         guide=v["guide"] or "-", gene=v["gene"], note=v["note"]))

    cols = ["id", "klass", "pos", "kind", "detail", "expected_vaf_wt",
            "expected_vaf_edit", "guide", "gene", "note"]
    with open(REF / "truth_variants.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: x["pos"]):
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(REF / "simulation_stats.json", "w") as fh:
        json.dump({"WT_parent": st_wt, "EDIT_pool": st_ed,
                   "read_length": READ_LEN, "seed": SEED}, fh, indent=2)

    print("truth table -> %s" % (REF / "truth_variants.tsv"))


if __name__ == "__main__":
    main()
