# Latch workflow — Kinome Guide Off-Target Scoring (FlashFry)

Packages the guide off-target/on-target scoring step (from `scripts/run_scoring.sh`)
into a registerable Latch workflow, so it runs on Latch instead of a local Pod.

## What it does

For each guide in the design (`data/rs3_context.tsv`, 2,054 guides), it:
1. builds a **FlashFry** database from the genome FASTA (self-contained: one JAR + Java),
2. computes **Doench 2016 on-target** (Rule Set 2) and **CFD off-target specificity**,
3. merges back to our library and flags **PASS/FAIL** at `CFD specificity >= cfd_min`
   (default 0.2), writing `guides_scored.tsv`.

> FlashFry is used because it is fully self-contained (no prebuilt index download).
> Rule Set 3 is an optional refinement — `rs3` is installed in the image and can be
> swapped in for on-target if preferred.

## Inputs

| Parameter | Value |
|-----------|-------|
| `context_tsv` | `data/rs3_context.tsv` (upload to Latch Data) |
| `genome_fasta` | hg38 FASTA (upload to Latch Data, or point at an existing copy) |
| `cfd_min` | 0.2 |
| `output_dir` | e.g. `latch://42942.account/kinome_screen/offtarget` |

## Register (one-time — needs the `latch` CLI + Docker; NOT doable over MCP)

```bash
python3 -m pip install latch      # Latch SDK + CLI
latch login                       # auth into your account (workspace 42942)
cd latch_workflow
latch register . --remote         # builds the container and registers the workflow
```

`latch register` builds and pushes the Docker image, so it must run somewhere with
Docker and Latch credentials. The MCP connection in this session can **launch** and
**monitor** workflows but cannot **register** them (no build/push capability).

## After registration

Tell me it's registered. I will:
1. `list_workflows` (search "off-target") to get the new `workflow_id`;
2. `get_workflow_schema` to confirm parameters;
3. `launch_workflow` with `context_tsv` + `genome_fasta` + `output_dir` into workspace 42942;
4. `get_execution` to poll, `get_task_logs` on failure, then pull `guides_scored.tsv`.

## Files

| File | Purpose |
|------|---------|
| `wf/__init__.py` | Latch `@workflow` + `@large_task` definition + metadata |
| `scripts/flashfry_pipeline.py` | prep / run / merge (pure-Python parts self-tested) |
| `Dockerfile` | Latch base + Java + FlashFry JAR + pipeline code |
| `version` | workflow version tag |
