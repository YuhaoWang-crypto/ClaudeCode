---
name: crispr-guide-scoring
description: Score a fixed CRISPR-Cas9 guide list for on-target activity (Doench 2016 / Rule Set) and off-target CFD specificity, then filter to PASS/FAIL at a CFD threshold. Prepares the tool inputs (30-mer contexts, spacer list) and merges FlashFry or GuideScan2 output back onto the library. Use after designing guides and before ordering, to fill on/off-target QC, or when someone asks to check guides for off-targets/specificity.
---

# CRISPR guide on-/off-target scoring

Fills the two QC fields a designed guide list needs before ordering: on-target
activity and off-target CFD specificity. Two interchangeable tool paths; the
Python glue (input prep, parse, merge, filter) is self-contained and self-tested.

## Inputs expected

A guide table with an `id`, `gene`, `spacer` (20 nt), `pam`, and a 30-mer
`context_30mer` column (4 nt up + 20 protospacer + 3 PAM + 3 nt down). The
**kinome-crispr-library** skill emits guides; build the 30-mer context from the
library's context column (Brunello has one) or from genomic flanks.

## Path A — FlashFry (self-contained: one JAR + genome FASTA)

Computes Doench 2016 on-target **and** CFD off-target in one tool, building its
own index — no prebuilt database download.

```bash
# needs: java, FlashFry.jar, hg38.fa
python3 scripts/flashfry_pipeline.py run \
  --genome hg38.fa --context rs3_context.tsv \
  --out guides_scored.tsv --ff /opt/FlashFry.jar --cfd-min 0.2
```

`prep` (build the target FASTA) and `merge` (parse + filter) run without any
genome/JAR and are unit-tested:

```bash
python3 scripts/flashfry_pipeline.py prep  --context rs3_context.tsv --fasta targets.fa
python3 scripts/flashfry_pipeline.py merge --context rs3_context.tsv --scored ff_out.tsv --out guides_scored.tsv
```

## Path B — rs3 (Rule Set 3) + GuideScan2

Higher-quality on-target (Rule Set 3) with a separate off-target tool.

```bash
python3 scripts/score_guides.py ontarget --context rs3_context.tsv --out ontarget_rs3.tsv   # pip install rs3
# run GuideScan2/CRISPOR on the spacer list -> a spacer-keyed specificity TSV, then:
python3 scripts/score_guides.py merge \
  --context rs3_context.tsv --ontarget ontarget_rs3.tsv \
  --offtarget offtarget_guidescan.tsv --out guides_scored.tsv --cfd-min 0.2
python3 scripts/score_guides.py selftest   # verify plumbing, no genome needed
```

## Output & finishing

`guides_scored.tsv`: per guide — on-target score, CFD specificity, off-target
count, and `pass_filter` (PASS if CFD specificity ≥ `--cfd-min`, default 0.2).

Finish: within each gene keep the top N (e.g. 6) PASS guides by on-target score;
flag any gene left with < N for manual review or a relaxed threshold. Guides from
a pre-filtered library (e.g. Brunello) already passed off-target QC — only
de-novo guides strictly require this step.

## Note

Both paths need external tools/genome that don't ship with the skill. The pure-
Python prep/merge/filter steps run anywhere; the scoring engines run on any
machine (or Latch Pod) with the genome and tool installed.
