# Geneformer — embeddings & in-silico perturbation → Signature

Geneformer is a BERT-style transformer **foundation model** pretrained on tens
of millions of human single-cell transcriptomes. In this skill it serves two
roles: an **orthogonal similarity space** (embeddings) and an **in-silico
perturbation** engine that produces signatures **without wet-lab data** — a
computational complement to the drug (L1000) and CRISPR (screen) modalities.

Sources:
- Weights (Apache-2.0): <https://huggingface.co/ctheodoris/Geneformer>
- BioNeMo model card (architecture/sizes):
  <https://docs.nvidia.com/bionemo-framework/latest/models/geneformer/>

## What it is (checked facts)

- **Architecture:** bidirectional (BERT) encoder, `fill-mask` pretraining
  objective (mask genes, predict them from co-expression context).
- **Input = rank-value encoding.** A cell is tokenised as its genes **ranked by
  expression** (normalised to each gene's corpus median), keeping the top
  ~1024–2048 genes. Tokens are **Ensembl gene IDs** via the packaged
  gene-median dictionary/tokenizer (~25k gene tokens).
- **Sizes (HF checkpoints):** V1-10M, V2-104M, V2-104M_CLcancer, V2-316M. On
  BioNeMo: 10M ≈ 6 layers / 256-dim; 106M ≈ 12 layers / 768-dim; ReLU, bf16.
- **Pretraining corpus** (Genecorpus / CZ CELLxGENE): ~23–95M non-diseased human
  cells. **Biases to remember:** ~9M brain cells, heavy 10x representation,
  donors skewed to <1 year old.

## Two ways to get a Signature

### A. Gene / cell embeddings → similarity (⚠️ interpretive)
Extract hidden states → a dense vector per cell (256-d or 768-d) or per gene.
Use these to cluster cells, or to place a perturbation's before/after cells in
embedding space. Embedding *distance* is a similarity, not a directional gene
signature — good for "are these two states alike?", not directly for WTCS.

### B. In-silico perturbation → a directional Signature (the useful bridge)
This is Geneformer's headline downstream task and how it plugs into
`perturbomics`:

1. Take cells in a starting state; tokenise (rank-value).
2. **Perturb the token sequence** — *delete* a gene token (simulate KO/CRISPRi)
   or *promote/insert* it (simulate activation/OE).
3. Re-embed perturbed vs original cells and measure, **per gene**, the shift its
   embedding/predicted-rank undergoes (Geneformer's toolkit reports a
   `cosine_shift` / rank change per affected gene).
4. Wrap the per-gene shift as a signed vector:

```python
from perturbomics import Signature
# shift: pandas Series indexed by Ensembl gene id, signed (‑ = pushed down)
sig = Signature(shift, name="insilico_KO_TP53", modality="geneformer",
                meta={"model": "V2-104M", "perturbed": "ENSG00000141510"})
```

Now it scores against drug and CRISPR signatures with the **same WTCS**
(`reference/connectivity.md`). Typical uses:
- validate a **CRISPR** screen hit in silico before wet-lab follow-up;
- expand a small experimental library with in-silico knockouts, then run
  `best_combinations` over the mixed set;
- ask whether a **drug** L1000 signature matches an in-silico knockout of its
  putative target (mechanism-of-action check).

## Sketch (real weights; heavy deps `torch`, `transformers`, `geneformer`)

```python
from transformers import AutoModelForMaskedLM
model = AutoModelForMaskedLM.from_pretrained("ctheodoris/Geneformer")
# Use the packaged geneformer.in_silico_perturber.InSilicoPerturber /
# EmbExtractor to tokenise an AnnData, delete/activate a gene token, and
# emit per-gene cosine shifts -> feed the Series into Signature(...) as above.
```

## Rigor labels

- ✅ **rigorous:** the embedding vectors, the cosine shifts, the WTCS computed
  from them — deterministic given fixed weights and input.
- ⚠️ **hypothesis:** *"the in-silico shift predicts the real perturbation"*.
  It's a model extrapolation, corpus-biased (see above). Treat in-silico
  signatures as prioritisation, and confirm in the **relevant tissue** — the
  brain/infant-skewed corpus may not transfer to your disease context.
