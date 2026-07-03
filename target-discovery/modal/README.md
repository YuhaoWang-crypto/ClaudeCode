# Modal GPU runner for Geneformer

Runs the Geneformer hypothesis engine (Step 3) on a Modal GPU. Modal is the
right home for this stage — it needs a custom GPU container with the geneformer
package + checkpoint, which the Boltz MCP (used for structure/binding) doesn't
cover.

## Status

- ✅ **`smoke` — verified working on a Modal T4.** Loads `Geneformer-V2-104M`
  (104.4M params, vocab 20275) and runs a real forward pass
  (`hidden_state_shape (1, 16, 768)`), with `import geneformer` succeeding.
- ⏳ `perturb` (in-silico deletion → goal-shift) — next: tokenize an IPF atlas
  and run `InSilicoPerturber`.

## Hard-won setup notes (all fixed in `geneformer_app.py`)

1. **Proxy**: this environment routes through a proxy; the Modal client needs
   `pip install python-socks` or it can't reach Modal.
2. **Geneformer install**: `pip install git+https://…` fails because pip's
   default `--filter=blob:none` partial clone chokes on HF's git promisor
   remote. Fix: full `git clone` (LFS smudge skipped), then `pip install` the
   local dir.
3. **Checkpoint download**: `git lfs pull` against HF is unreliable here
   (fetches nothing, leaving pointer files → "safetensors header too large").
   Fix: `huggingface_hub.snapshot_download(allow_patterns=[...])` over the HF
   CDN, cached in a Modal Volume (`gf-models`).
4. **transformers version**: Geneformer does `from transformers import
   SpecialTokensMixin`, dropped from the top-level namespace in newer
   transformers. Pin `transformers==4.46.3`.

## Run

```bash
pip install modal python-socks
export MODAL_TOKEN_ID=... MODAL_TOKEN_SECRET=...
modal run modal/geneformer_app.py            # fetch checkpoint + GPU smoke
modal run modal/geneformer_app.py::list_repo # inspect available checkpoints
```

Checkpoints available in the repo: `Geneformer-V1-10M`, `Geneformer-V2-104M`
(default here), `Geneformer-V2-104M_CLcancer`, `Geneformer-V2-316M`.
