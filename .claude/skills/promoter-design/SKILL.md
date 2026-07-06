---
name: promoter-design
description: >
  Design stimulus-responsive AND cell-type-selective synthetic promoters by
  running Evo2 (DNA generation) + AlphaGenome (cell-type RNA-seq scoring),
  self-hosted on Modal via proto-tools/evodesign. Use when the user wants to
  design/optimize a promoter or regulatory sequence for a given stimulus (IFN,
  hypoxia, NF-kB, glucocorticoid, estrogen, oxidative/ER stress, cAMP) and a
  target cell type, or to run the Evo2+AlphaGenome design loop on Modal. Captures
  the known-good build recipe and every fix needed to get the stack running.
---

# Stimulus-responsive, cell-type-selective promoter design (Evo2 + AlphaGenome on Modal)

Design synthetic promoters that turn on **only for a given stimulus** (via a fixed
response element) **and only in a target cell type** (optimized by AlphaGenome
contrastive scoring). Evo2 writes the designable region; AlphaGenome scores
predicted RNA-seq expression in the target cell vs off-target cells; MCMC keeps
the best. Everything runs on the user's Modal GPUs — no PROTO_API_KEY, weights
come from HuggingFace.

## When to use
- "Design a promoter that responds to <stimulus> and is specific to <cell>."
- "Run the Evo2 + AlphaGenome design loop / proto-tools pipeline on Modal."
- Extending to stimulus-AND-cell-type gates, or new stimulus/cell pairs.

## Prerequisites (verify first)
- **Modal** auth present: `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET` in env, and the
  client installed with proxy support: `pip install 'modal[api-proxy-support]' python-socks`.
  Confirm with `modal app list`.
- **HuggingFace**: `HF_TOKEN` whose account has **accepted the gated AlphaGenome
  license** at huggingface.co/google/alphagenome-all-folds. Create the Modal
  secret: `modal secret create huggingface HF_TOKEN=$HF_TOKEN`. (Evo2 weights
  `arcinstitute/evo2_*` are Apache-2.0, ungated.)
- The pipeline code lives in `promoter_design/` (this repo): `elements.py`,
  `design_pipeline.py`, `build_constructs.py`, `dual_and_designs.py`,
  `modal_design.py`. The Modal app is `modal_design.py`.

## Quickstart
```bash
cd promoter_design
modal deploy modal_design.py                    # register app (no compute)
modal run modal_design.py::assemble             # CPU sanity check: builds FASTA
# full design (Evo2 generate + AlphaGenome cell-type score):
modal run modal_design.py::full_design --stimulus interferon_typeII --target THP1
# add cell-type specificity (penalize off-target cells):
modal run modal_design.py::full_design --stimulus interferon_typeII --target THP1 --contrastive
# bigger search + larger designable region:
modal run modal_design.py::full_design --contrastive --num-steps 40 --num-results 3 --spacer-len 96
```
First `full_design` builds the Evo2 + AlphaGenome standalone envs into the
`proto-cache` Volume (~30–60 min, one time); later runs reuse them and are fast.

## Architecture
`[insulator] [stimulus RE ×N] [designable region] [minimal promoter] [Kozak] [handle]`
1. **Assemble** — fixed biology (`make_segments` in `design_pipeline.py`).
2. **Generate** — Evo2 fills the designable region (`Evo2Generator`, BF16).
3. **Score** — `alphagenome_interval_track_constraint`: embed cassette in a
   16,384-bp window, predict RNA_SEQ in `ontology_terms` (target) and penalize
   `contrastive_ontology_terms` (off-targets).
4. **Optimize** — `MCMCOptimizer` maximizes target−off-target margin over `num_steps`.

## Known-good config — DO NOT regress these (each was a real failure)
The `modal_design.py` `gpu_image` and `full_design` already bake in all of these.
See `references/fixes.md` for the full table. Critical ones:
- **CUDA 12.4** pinned in Evo2 `setup.sh` (unpinned pulls CUDA 13.x → breaks flash-attn/torch2.6).
- **torch index cu124** forced (driver-derived cu128 has no torch 2.6.0).
- **uv cache on local disk** (`UV_CACHE_DIR=/tmp/uv_cache`, `UV_LINK_MODE=copy`) —
  Modal Volumes can't persist uv temp files (EPERM). Needed for BOTH Evo2 and AlphaGenome.
- **FP8 disabled → BF16** (patch vortex `fp8_autocast(enabled=False)`); FP8 needs
  compute ≥8.9 and its GEMM rejects short sequences. BF16 runs on A100.
- **AlphaGenome context**: pad left/right so left+target+right == `context_length` (16384).
- **Verified ontology terms** (see below); a wrong CL term returns no RNA_SEQ tracks.
- **Sanitize model output to ACGT** before AlphaGenome (see Known issue).

## Verified AlphaGenome RNA_SEQ ontology terms (cell → CL)
Only terms with real RNA_SEQ tracks work. Mapping used (`elements.CELL_CONTEXTS`):
| Cell | ontology | biosample |
|------|----------|-----------|
| THP1 | CL:0001054 | CD14-positive monocyte |
| HepG2 | CL:0000182 | hepatocyte |
| Jurkat | CL:0000084 | T-cell |
| MCF7 | CL:0002327 | mammary epithelial cell |
| HEK293 | CL:0002518 | kidney epithelial cell |
| SHSY5Y | CL:0002319 | neural cell |
Full list of RNA_SEQ-valid terms: `references/rna_seq_ontologies.md`. To discover
more, run `modal run modal_design.py::dump_ag_outputs` (parses AlphaGenome's
example metadata for output_type→ontology→biosample).

## Key parameters (`full_design` / `design`)
- `--stimulus` one of `elements.ELEMENTS` (interferon_typeI/II, inflammation,
  hypoxia, camp, estrogen, glucocorticoid, oxidative_stress, er_stress_*).
- `--target` a cell in `elements.CELL_CONTEXTS`.
- `--contrastive` on = cell-type-specificity objective (target vs 5 off-targets);
  off = maximize target expression only.
- `--spacer-len` size of the designable region (8 = quick; 96+ = real leverage).
- `--num-steps`, `--num-results` = MCMC search size.
- `--lineage <myeloid|hepatocyte|neuronal|...>` = stimulus-AND-cell-type gate.

## Known issue (open)
Long Evo2 generations occasionally emit a non-ACGT token. AlphaGenome's validator
rejects anything but ACGT (it rejects `N` too, despite the message). The constraint
is patched to sanitize composed sequences; **map stray chars to `A`, not `N`**
(mapping to N still fails). Small `--spacer-len` (8) runs are unaffected and verified.

## Troubleshooting
- `setup.sh failed (exit 2)` in a build → read the log; it's almost always one of
  the CUDA / uv / torch-index issues above. The setup scripts are patched in the
  image build; if you change proto-tools versions, re-verify the patches apply.
- `Unable to extract numeric values for 'RNA_SEQ'` → the ontology term has no
  RNA_SEQ track; use a verified CL term.
- `sequence must only contain A/C/G/T/N` → model emitted a stray char; see Known issue.
- Model reloads every MCMC step (`DeviceManager: LRU eviction`) on a single GPU —
  expected; use H100 (both models resident) if you need speed.

## Diagnostics (CPU, cheap) in `modal_design.py`
`inspect_runtime`, `dump_weight_sources`, `dump_standalone_specs`,
`dump_ag_outputs`, `check_hf_access`, `validate_design` (validates the design
graph on CPU with no GPU spend). Use these before spending GPU.
