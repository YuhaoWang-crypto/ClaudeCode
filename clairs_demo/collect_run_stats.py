#!/usr/bin/env python3
"""Harvest the numbers that only exist inside a completed ClairS run directory.

`analyze.py` works off the committed VCFs alone. A few observations need the
run's intermediates instead — the Clair3 germline VCFs and the haplotagged
tumour BAM, both under <run>/output/tmp/ — so they are collected once here and
frozen into data/run_derived.json for the report to cite.

    python3 clairs_demo/collect_run_stats.py ~/demo        # <base>/<platform>/…

Needs pysam (present in the ClairS conda env).
"""
import gzip
import json
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REGION = ("chr17", 80_000_000, 80_100_000)
PLATFORMS = ["ont", "ilmn", "pacbio_hifi"]


def germline_counts(vcf):
    """PASS germline calls in the region, split by genotype."""
    if not os.path.exists(vcf):
        return None
    het = hom = 0
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[0] != REGION[0] or not REGION[1] <= int(f[1]) <= REGION[2]:
                continue
            if f[6] != "PASS":
                continue
            gt = f[9].split(":")[0]
            if gt in ("0/1", "1/2"):
                het += 1
            elif gt == "1/1":
                hom += 1
    return {"het": het, "hom": hom}


def pass_calls(platform):
    out = []
    with gzip.open(os.path.join(HERE, "results", platform, "output.vcf.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if f[6] == "PASS":
                out.append((int(f[1]), f[4], "H" in f[7].split(";")))
    return out


def bam_stats(platform, run):
    """Read lengths, haplotagged fraction, and alt-read haplotype partitioning."""
    import pysam
    raw = os.path.join(run, "input", "HCC1395_tumor_chr17_demo.bam")
    tagged = os.path.join(run, "output", "tmp", "clair3_output",
                          "phased_output", "tumor_chr17.bam")
    res = {}
    if os.path.exists(raw):
        lens = [r.query_length for r in
                pysam.AlignmentFile(raw).fetch(*REGION) if r.query_length]
        res["read_length_median"] = int(st.median(lens))
        res["read_length_p90"] = int(sorted(lens)[int(len(lens) * .9)])
        res["reads_in_region"] = len(lens)
    if not os.path.exists(tagged):
        return res                      # Illumina: no phasing stage at all

    bam = pysam.AlignmentFile(tagged)
    reads = list(bam.fetch(*REGION))
    res["haplotagged_fraction"] = round(
        sum(r.has_tag("HP") for r in reads) / len(reads), 3)

    one_hap = both_hap = 0
    hp_totals = {0: 0, 1: 0, 2: 0}
    for pos, alt, _ in pass_calls(platform):
        alt_by_hp = {0: 0, 1: 0, 2: 0}
        for col in bam.pileup(REGION[0], pos - 1, pos, truncate=True,
                              min_base_quality=0, stepper="samtools"):
            for p in col.pileups:
                if p.is_del or p.is_refskip:
                    continue
                hp = p.alignment.get_tag("HP") if p.alignment.has_tag("HP") else 0
                hp_totals[hp] += 1
                if p.alignment.query_sequence[p.query_position] == alt:
                    alt_by_hp[hp] += 1
        if alt_by_hp[1] and alt_by_hp[2]:
            both_hap += 1
        else:
            one_hap += 1
    res.update({
        "alt_reads_confined_to_one_haplotype": one_hap,
        "alt_reads_split_across_haplotypes": both_hap,
        # ClairS only sets the H flag when BOTH haplotypes carry reads at the
        # site (all_hp1 * all_hp2 > 0 in src/haplotype_filtering.py).
        "reads_by_hp_untagged_hp1_hp2": [hp_totals[0], hp_totals[1], hp_totals[2]],
    })
    return res


def main(base):
    out = {}
    for p in PLATFORMS:
        run = os.path.join(base, p)
        if not os.path.isdir(run):
            print(f"[skip] {run} not found")
            continue
        c3 = os.path.join(run, "output", "tmp", "clair3_output")
        entry = {
            "germline_normal": germline_counts(
                os.path.join(c3, "clair3_normal_output", "merge_output.vcf.gz")),
            "germline_tumour": germline_counts(
                os.path.join(c3, "clair3_tumor_output", "merge_output.vcf.gz")),
            "h_flagged_calls": sum(1 for _, _, h in pass_calls(p) if h),
            "n_pass_calls": len(pass_calls(p)),
        }
        try:
            entry.update(bam_stats(p, run))
        except ImportError:
            print("[warn] pysam missing — BAM stats skipped")
        out[p] = entry

    path = os.path.join(HERE, "data", "run_derived.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(json.dumps(out, indent=2))
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/demo"))
