#!/usr/bin/env python3
"""
05_make_report.py - assemble the demo report from the pipeline's real output.

Reads results/analysis.json plus the intermediate-data snapshots the pipeline
saved, and injects them into report/template.html.  Nothing in the report is
typed by hand: every number and every data excerpt comes from the run.
"""

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REF = ROOT / "ref"
WORK = ROOT / "work"
RES = ROOT / "results"
REPORT = ROOT / "report"
SNAP = WORK / "snapshots"
REPORT.mkdir(parents=True, exist_ok=True)

PLACEHOLDER = "/*__DATA__*/null/*__END__*/"


def read_text(path, limit=None):
    try:
        t = Path(path).read_text().rstrip("\n")
    except OSError:
        return ""
    if limit:
        lines = t.split("\n")
        if len(lines) > limit:
            t = "\n".join(lines[:limit]) + "\n..."
    return t


def tool_version(cmd, pattern=None):
    try:
        out = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                             timeout=30)
        text = (out.stdout + out.stderr)
        for line in text.splitlines():
            if pattern is None or pattern.lower() in line.lower():
                return line.strip()
    except Exception:
        pass
    return "n/a"


def main():
    analysis = json.load(open(RES / "analysis.json"))
    design = json.load(open(REF / "design_notes.json"))
    simstats = json.load(open(REF / "simulation_stats.json"))

    truth = []
    with open(REF / "truth_variants.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            truth.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))

    snapshots = {
        "raw_fastq": read_text(SNAP / "00_raw_EDIT_pool_R1.txt", 8),
        "raw_fastq_r2": read_text(SNAP / "00_raw_EDIT_pool_R2.txt", 4),
        "bam_header": read_text(SNAP / "02_bam_header.txt", 3),
        "bam_records": read_text(SNAP / "02_bam_records.txt", 3),
        "bam_insertion": read_text(SNAP / "02_bam_cigar_insertion.txt", 3),
        "bam_split": read_text(SNAP / "02_bam_split_reads.txt", 3),
        "pileup_pheA": read_text(SNAP / "03_pileup_on-target_pheA.txt", 12),
        "pileup_lacZ": read_text(SNAP / "03_pileup_on-target_lacZ.txt", 12),
        "vcf_raw": read_text(SNAP / "04_vcf_raw.txt", 14),
        "isec": read_text(SNAP / "06_isec.txt", 20),
        "delly": read_text(SNAP / "07_delly_sv.txt", 4),
        "snpeff": read_text(SNAP / "08_snpeff_annotation.txt", 12),
        "filtered_out": read_text(WORK / "vcf" / "filtered_out.txt", 6),
    }

    # a cleaned read, to show next to the raw one
    clean = subprocess.run(
        "zcat work/qc/EDIT_pool_R1.clean.fastq.gz 2>/dev/null | head -4",
        shell=True, cwd=ROOT, capture_output=True, text=True).stdout.rstrip("\n")
    snapshots["clean_fastq"] = clean

    fastp = {}
    for s in ("WT_parent", "EDIT_pool"):
        d = json.load(open(WORK / "qc" / ("%s.fastp.json" % s)))
        fastp[s] = dict(
            before=d["summary"]["before_filtering"]["total_reads"],
            after=d["summary"]["after_filtering"]["total_reads"],
            low_quality=d["filtering_result"]["low_quality_reads"],
            too_many_n=d["filtering_result"]["too_many_N_reads"],
            too_short=d["filtering_result"]["too_short_reads"],
            adapter=d["adapter_cutting"]["adapter_trimmed_reads"],
            q30_before=round(100 * d["summary"]["before_filtering"]["q30_rate"], 2),
            q30_after=round(100 * d["summary"]["after_filtering"]["q30_rate"], 2),
            insert_peak=d.get("insert_size", {}).get("peak"),
        )

    # how many variants each calling model found (the haploid-vs-pooled point)
    def nvar(p):
        try:
            return int(subprocess.run("bcftools view -H %s | wc -l" % p, shell=True,
                                      cwd=ROOT, capture_output=True,
                                      text=True).stdout.strip())
        except Exception:
            return None

    models = dict(
        haploid_wt=nvar("work/vcf/WT_parent.haploid.vcf.gz"),
        haploid_edit=nvar("work/vcf/EDIT_pool.haploid.vcf.gz"),
        pooled_wt=nvar("work/vcf/WT_parent.vcf.gz"),
        pooled_edit=nvar("work/vcf/EDIT_pool.vcf.gz"),
        joint_raw=nvar("work/vcf/joint.raw.vcf.gz"),
        joint_pass=nvar("work/vcf/joint.pass.vcf.gz"),
    )

    tools = [
        ("fastp", tool_version("fastp --version", "fastp"), "质控 / 接头与低质量碱基切除"),
        ("BWA-MEM", tool_version("bwa 2>&1 | grep Version", "version"), "短读段比对"),
        ("samtools", tool_version("samtools --version | head -1", "samtools"), "排序 / 去重 / 覆盖度"),
        ("bcftools", tool_version("bcftools --version | head -1", "bcftools"), "SNV / InDel 检测与过滤"),
        ("Delly", tool_version("delly --version 2>&1 | head -1", "version"), "结构变异 (大片段缺失)"),
        ("SnpEff", tool_version("snpEff -version 2>/dev/null", "snpeff"), "变异功能注释"),
    ]

    payload = dict(
        meta=dict(
            chrom=analysis["chrom"],
            genome_len=analysis["genome_len"],
            organism=design["organism"],
            n_genes=design["n_genes"],
            gc=design["gc_percent"],
            seed=design["seed"],
            generated=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            read_length=simstats["read_length"],
        ),
        qc=analysis["qc"],
        fastp=fastp,
        models=models,
        guides=analysis["design"],
        genes=[{k: g[k] for k in ("name", "start", "end", "strand")}
               for g in design["genes"]],
        variants=analysis["variants"],
        svs=analysis["svs"],
        offtargets=analysis["offtargets"],
        planted=design["planted_offtarget_sites"],
        coverage=analysis["coverage"],
        validation=analysis["validation"],
        truth=truth,
        snapshots=snapshots,
        tools=[dict(name=n, version=v, role=r) for n, v, r in tools],
    )

    template = (REPORT / "template.html").read_text()
    if PLACEHOLDER not in template:
        raise SystemExit("placeholder %r missing from template.html" % PLACEHOLDER)
    html = template.replace(PLACEHOLDER, json.dumps(payload, ensure_ascii=False))
    out = REPORT / "index.html"
    out.write_text(html)
    print("report -> %s  (%.0f KB)" % (out, out.stat().st_size / 1024))
    print("embedded: %d variants, %d SVs, %d off-target sites, %d coverage bins"
          % (len(payload["variants"]), len(payload["svs"]),
             len(payload["offtargets"]), len(payload["coverage"]["ratio"])))


if __name__ == "__main__":
    main()
