# Every fix needed to self-host Evo2 + AlphaGenome on Modal

Each row was a real, distinct build/runtime failure. All are already applied in
`promoter_design/modal_design.py`. If proto-tools/Evo2/AlphaGenome versions change,
re-verify each patch still applies (the string-replace patches `assert` their target).

| # | Symptom | Root cause | Fix (where) |
|---|---------|-----------|-------------|
| 1 | `evo2 setup.sh failed (exit 2)`, nvcc build 13.3 | `micromamba create ... cuda-toolkit` unpinned → CUDA 13.3; flash-attn/torch2.6 need 12.x | Patch evo2 `setup.sh` to add `"cuda-version=12.4"` (gpu_image run_commands) |
| 2 | `Could not persist temporary file .../uv_cache ... Operation not permitted` | Modal Volume can't persist uv temp files (EPERM) | Patch setup.sh: `export UV_CACHE_DIR=/tmp/uv_cache; export UV_LINK_MODE=copy` (Evo2 **and** AlphaGenome) |
| 3 | `No solution ... no version of torch==2.6.0` (found on cu128) | proto-tools derives torch index from driver CUDA (12.8→cu128); torch 2.6.0 ships cu124/cu126 | Patch evo2 setup.sh: `export RECOMMENDED_TORCH_INDEX=https://download.pytorch.org/whl/cu124` |
| 4 | `Device compute capability 8.9 or higher required for FP8` | Evo2 transformer-engine FP8 needs SM≥8.9; A100 is 8.0 | (transient) — real fix is #5 |
| 5 | `cuBLAS ... unsupported value or parameter` in FP8 GEMM (even on H100) | FP8 cuBLASLt rejects short sequences | Patch vortex `layers.py`: `fp8_autocast(enabled=False)` → BF16 (in `full_design`, `disable_fp8`). Lets it run on A100. |
| 6 | `Insufficient AlphaGenome context for length N: need 16376 flanking bp` | AlphaGenome scores a fixed 16384-bp window; target needs full flanks | `make_segments` pads left/right with `neutral_flank` so left+target+right == context_length |
| 7 | `Unable to extract numeric values for 'RNA_SEQ'` | ontology term (e.g. CL:0000576 monocyte) has no RNA_SEQ track | Use verified CL terms (`references/rna_seq_ontologies.md`); `CELL_CONTEXTS` updated |
| 8 | `sequence must only contain A/C/G/T/N` (lowercase) | elements.py uses lowercase for spacer regions; Evo2 also emits lowercase | `.upper()` all assembled parts in `make_segments`; patch constraint to uppercase composed seqs |
| 9 | `sequence must only contain A/C/G/T/N` (still, after uppercase) | Evo2 long generations emit a non-ACGT token; AlphaGenome rejects even `N` | **OPEN**: sanitize composed seqs mapping stray chars to **`A`** (not N). Small spacer_len(8) unaffected. |

## Environment facts (verified)
- Evo2 model_checkpoint default `evo2_7b`; use `evo2_1b_base` (smallest) for cheap runs.
- Evo2 weights Apache-2.0 (ungated). AlphaGenome `google/alphagenome-all-folds` is
  hf-gated, non-commercial research only; the run HF account must accept it.
- AlphaGenome loads from `ALPHAGENOME_CHECKPOINT_PATH` (auto-downloaded to the Volume).
- proto-tools builds per-tool micromamba standalone venvs under `PROTO_HOME`
  (mounted to the `proto-cache` Volume); it downloads its own micromamba.
- Single-GPU runs reload Evo2↔AlphaGenome each MCMC step (LRU eviction) — slow but works.
- No `PROTO_API_KEY` needed for local execution; that key is only for evodesign's
  hosted backend (device='cloud').
