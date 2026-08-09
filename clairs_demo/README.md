# ClairS demo — HCC1395 / HCC1395BL somatic SNVs on chr17:80.0–80.1 Mb

Reproducible run of [ClairS](https://github.com/HKU-BAL/ClairS) (HKU-BAL),
installed **from source without Docker**, on the official quick-demo data, with
the calls benchmarked against the SEQC2 truth set.

Everything here is driven by the scripts in
[`.claude/skills/clairs-install-demo/scripts/`](../.claude/skills/clairs-install-demo/scripts/).

## Reproduce

```bash
bash .claude/skills/clairs-install-demo/scripts/install_clairs.sh ~/clairs   # ~10 min, ~8 GB
bash .claude/skills/clairs-install-demo/scripts/run_demo.sh ont  ~/demo/ont  4
bash .claude/skills/clairs-install-demo/scripts/run_demo.sh ilmn ~/demo/ilmn 4
bash .claude/skills/clairs-install-demo/scripts/run_demo.sh pacbio_hifi ~/demo/pacbio_hifi 4
```

## The data

| | |
|---|---|
| Sample | HCC1395 (breast-cancer cell line, tumour) / HCC1395BL (matched normal lymphoblastoid) |
| Region | `chr17:80,000,000-80,100,000` (100 kb) |
| Truth | SEQC2 high-confidence somatic SNVs v1.2 ([Fang et al., Nat. Biotechnol. 2021](https://www.nature.com/articles/s41587-021-00993-6)), restricted to the SEQC2 high-confidence BED |
| Source | `https://www.bio8.cs.hku.hk/clairs/quick_demo/{ont,ilmn,pacbio_hifi}/` |

Three independent sequencing platforms over the same cell-line pair, so the
same 29 truth SNVs are called from ONT R10.4.1, Illumina NovaSeq and PacBio
Revio reads.

## Results

Somatic SNVs, `PASS` calls only, benchmarked with `clairs.py compare_vcf`
against 29 truth SNVs inside the SEQC2 high-confidence BED:

| Dataset | `--platform` | Precision | Recall | F1 | TP | FP | FN | wall |
|---|---|---|---|---|---|---|---|---|
| ONT R10.4.1, Guppy 6.1.5 (~70×/45×) | `ont_r10_guppy` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m40s |
| Illumina NovaSeq 6000 (~50×/40×) | `ilmn` | 1.0 | 1.0 | 1.0 | 29 | 0 | 0 | 0m25s |
| PacBio Revio HiFi (~60×/40×) | `hifi_revio` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m38s |

**All three reproduce upstream's published expected output exactly** — which is
the point of the exercise: it shows the from-source install behaves identically
to the official Docker image.

The one missed variant is the same site on both long-read platforms —
`chr17:80,094,483 T>C`, SEQC2 TVAF ≈ 0.096, the lowest-VAF truth variant in the
window. Neither long-read run emits it at all (default `--snv_min_af 0.05`);
Illumina calls it at `AF=0.0642` with `NAF=0.0000` in the matched normal.
Per-run allele-fraction floors: ONT 0.139, PacBio 0.167, Illumina 0.064. With
n=29 in one 100 kb window and differing coverage, that is a statement about
this region, not a platform ranking.

## Analysis

[**ANALYSIS.md**](ANALYSIS.md) digs into what the calls say — cross-platform
concordance, why per-site AF agreement is unmeasurable in this window (and why
Pearson *r* misleads there), an APOBEC mutational signature at 5.1× enrichment,
and loss of heterozygosity across the whole window that explains the `H`-flag
difference between ONT and PacBio.

```bash
python3 clairs_demo/analyze.py            # -> analysis.json + figures/
python3 clairs_demo/collect_run_stats.py ~/demo   # numbers from a live run dir
python3 clairs_demo/make_report.py               # -> self-contained report.html
```

The same write-up is published as a standalone HTML report:
[claude.ai/code/artifact/cde4334b](https://claude.ai/code/artifact/cde4334b-43ea-4d1c-bc66-dbbfa08c16ff).
`report_template.html` is its source; `make_report.py` inlines the figures as
data URIs so the page needs no external requests.

## What is in this directory

```
ANALYSIS.md                        the write-up
analyze.py                         analysis from the committed VCFs (stdlib + matplotlib)
collect_run_stats.py               numbers that need a completed run directory
analysis.json                      every computed value
figures/                           four figures used by ANALYSIS.md
report_template.html               source of the published HTML report
make_report.py                     inlines figures -> self-contained report.html
data/truth_chr17_*.vcf             SEQC2 truth, subset to the demo window
data/truth_context.tsv             trinucleotide context per truth SNV
data/run_derived.json              germline/LOH + haplotagging measurements
results/<platform>/output.vcf.gz   ClairS somatic calls (PASS = somatic)
results/<platform>/benchmark.txt   compare_vcf metrics
results/<platform>/fn.vcf          the missed truth variants
results/<platform>/run_clairs.log  full pipeline log
```

## Notes

- ClairS ran CPU-only on 4 threads; no GPU is needed.
- Installed from source (micromamba) because the sandbox has the `docker`
  client but no daemon. The key deviation from upstream's Dockerfile is pinning
  `clair3=1.2.0` — see
  [`reference/installation.md`](../.claude/skills/clairs-install-demo/reference/installation.md).
- ClairS `0.4.5`, commit `a87db4a` (2026-07-17); Clair3 1.2.0; PyTorch 2.1.2 +
  TensorFlow 2.15 (CPU).
