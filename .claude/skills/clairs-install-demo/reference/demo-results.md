# Quick-demo results (measured)

All numbers below were produced by `scripts/run_demo.sh` on a 4-core CPU-only
container, ClairS 0.4.5 (commit `a87db4a`, 2026-07-17) installed from source via
`scripts/install_clairs.sh` — **not** from the `hkubal/clairs` Docker image.

## Setup

| | |
|---|---|
| Sample | HCC1395 tumour / HCC1395BL matched normal |
| Region | `chr17:80,000,000-80,100,000` |
| Truth | SEQC2 high-confidence somatic SNVs v1.2, inside the SEQC2 high-confidence BED |
| Benchmark | `clairs.py compare_vcf --input_filter_tag PASS` |
| Threads | 4 |

## Measured

| Dataset | `--platform` | Precision | Recall | F1 | TP | FP | FN | ClairS wall |
|---|---|---|---|---|---|---|---|---|
| ONT R10.4.1, Guppy 6.1.5 (~70×/45×) | `ont_r10_guppy` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m40s |
| Illumina NovaSeq 6000 (~50×/40×) | `ilmn` | 1.0 | 1.0 | 1.0 | 29 | 0 | 0 | 0m25s |
| PacBio Revio HiFi (~60×/40×) | `hifi_revio` | 1.0 | 0.9655 | 0.9825 | 28 | 0 | 1 | 1m38s |

All three match upstream's documented expected output for these demos exactly.
29 truth SNVs in the window; no records fell outside the high-confidence BED.

## Reading the misses

Both long-read runs miss the *same* single site, and it is the lowest-VAF truth
variant in the window:

```
chr17:80,094,483  T>C   SEQC2 TVAF=0.096  NVAF=0.002  (PacBio-measured TVAF=0.066)
```

Neither the ONT nor the HiFi run emits it at all — not even as `LowQual` — so it
never reached the models: at that depth the support falls under the default
`--snv_min_af 0.05` candidate threshold. Illumina recovers it as a `PASS` call
with `AF=0.0642` (7 alt / 109 reads) and `NAF=0.0000` (0 / 68 in the normal),
QUAL 10.86.

The allele-fraction floor of the PASS calls says the same thing:

| Dataset | min AF | median AF | max AF | n |
|---|---|---|---|---|
| ONT | 0.139 | 0.213 | 0.337 | 28 |
| Illumina | 0.064 | 0.196 | 0.315 | 29 |
| PacBio HiFi | 0.167 | 0.205 | 0.441 | 28 |

Read this as "where each run's calls bottom out in this 100 kb window", not as a
platform ranking — n=29, one region, and the three datasets differ in coverage.
The transferable point is that the misses live in the low-VAF tail, so
`fn.vcf` plus the `AF` column is the first place to look before blaming a model
or an install. If you genuinely need that tail, `--snv_min_af` is the knob, at a
cost in precision and runtime.

## Why this matters as an install check

The from-source route swaps upstream's baked-in Clair3 for a conda `clair3=1.2.0`
build, and compiles `realigner` / `debruijn_graph` against a different toolchain
and Boost version than the Ubuntu 16.04 image uses. Matching upstream's published
expected output exactly, on all three platforms, is the evidence that none of
those substitutions changed the result. If your own install lands on different
TP/FP/FN here, treat it as a broken install, not as run-to-run noise — ClairS is
deterministic on fixed inputs.

## Beyond the metrics

`clairs_demo/ANALYSIS.md` in this repo analyses what the calls themselves say, and
three of its findings generalise beyond the demo:

- **Per-site AF agreement is often unmeasurable.** In this window the truth VAFs
  have an IQR of 0.011 while the binomial sampling error on a single AF at ~100×
  is ~0.039. Pearson *r* against a truth VAF then reports 0.46–0.80 while the
  rank correlation reports ≈0 — the *r* is carried by one leverage point. Check
  the rank correlation, and the VAF spread, before believing an AF-concordance
  number. Aggregate bias (median AF ÷ truth VAF) is the statistic that survives.
- **The `H` flag is unreliable in LOH regions.** ClairS sets it only when both
  haplotypes carry reads at the site. Where the tumour has lost heterozygosity
  there is nothing to phase against, every read lands in one haplotype, and the
  flag silently never fires — here ONT tagged 24/28 calls and PacBio 0/28 with
  identical underlying phasing evidence. Absence of `H` is not evidence against
  a call.
- **`ilmn` is a different pipeline, not just a different model.** It skips the
  Clair3 germline and phasing stages entirely, which is why it runs ~4× faster
  and never emits `H`.

## Scaling expectations

This is a 100 kb toy region. Whole-genome ONT tumour-normal at 50×/25× is hours
on tens of cores and is dominated by the two Clair3 germline passes and
full-alignment calling. Published whole-genome accuracy (ClairS v0.4.0, ONT
R10.4.1, 50×/25×, SEQC2 truth): SNV recall/precision 86.86 % / 93.01 % at
VAF ≥ 0.05, rising to 94.65 % / 96.63 % at VAF ≥ 0.2. Do not extrapolate the
demo's F1 to a genome — a 100 kb high-confidence window is far easier than the
whole genome.
