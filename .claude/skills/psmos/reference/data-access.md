# Data access — the recipes that actually work behind the agent proxy

Every external service PSMOS uses, with the exact incantation that works in this
environment. All HTTPS egress goes through the agent proxy, which re-terminates
TLS — so tools must trust `/root/.ccr/ca-bundle.crt`.

## Modal (Evo2 compute) — grpclib + python-socks

Modal's client speaks gRPC via **grpclib**, which needs a SOCKS/HTTP-proxy shim
to tunnel through the agent proxy.

```bash
pip install modal python-socks
export SSL_CERT_FILE=/root/.ccr/ca-bundle.crt
# MODAL_TOKEN_ID / MODAL_TOKEN_SECRET are already in the environment
```

Without `python-socks` you get `ConnectionError: Could not connect to the Modal
server` (the underlying cause is `ModuleNotFoundError: No module named
'python_socks'`, swallowed into a generic connection error). Verify with a
trivial `@app.function()` that returns `3+4` before doing anything expensive.

## Evo2 on Modal — the image that loads AND scores

```python
image = (modal.Image
    .from_registry("nvcr.io/nvidia/pytorch:25.04-py3")   # Arc's validated base
    .apt_install("git")
    .pip_install("hf_transfer", "huggingface_hub", "biopython")
    .run_commands("pip install evo2"))                    # pulls vtx; light deps
```

- **Base tag matters.** `25.04-py3` ships flash-attn + transformer-engine
  prebuilt and matches vortex's kernels. `25.01-py3` loads the model fine but
  fails at scoring with `fwd(): incompatible function arguments` (kernel/ABI
  mismatch). Don't guess other tags — this is the one Arc's Dockerfile uses.
- **weights_only.** NGC PyTorch ≥2.6 defaults `torch.load(weights_only=True)`,
  which rejects the checkpoint (`Unsupported global:
  transformer_engine.common.recipe._OverrideLinearPrecision`). Monkeypatch
  `torch.load` to `weights_only=False` around `Evo2(MODEL_NAME)` — the Arc
  checkpoint is trusted.
- **Weights** (`arcinstitute/evo2_7b`, ~15 GB) go in a `modal.Volume` with
  `HF_HOME` pointed at it so they download once.
- The build is **CPU-only** (no GPU-seconds) — nothing compiles; GPU is used
  only by the scoring class. Scoring is `model.score_sequences([cds, ...])` →
  list of mean log-likelihoods. Feed **CDS nucleotides**, not protein.

## UniProt — ortholog protein + the hard gate

```
https://rest.uniprot.org/uniprotkb/search?query=(gene:{sym}) AND (organism_id:{taxon})
        &format=json&size=25&fields=accession,reviewed,protein_name,gene_names,length,sequence
```

Selection: prefer reviewed (Swiss-Prot); then the entry whose length is
**closest to the human ortholog's length** (orthologues are near-equal length —
"take the longest" grabs fusion/mis-annotated isoforms). Zero hits for the
negative controls (yeast/plant on Notch genes) is the *empirical* gate — real,
not asserted.

## Ensembl REST — CDS (Evo2 input) and gene intervals

Vertebrates and metazoa share `https://rest.ensembl.org`.

```
/xrefs/symbol/{species}/{sym}            -> Ensembl gene id
/lookup/id/{gene}?expand=1               -> canonical_transcript (+ Transcript list)
/sequence/id/{transcript}?type=cds       -> ACGT CDS
/lookup/symbol/{species}/{sym}           -> chrom/start/end/strand (interval)
```

- Strip only the `ENS…`/`FB…` version suffix from transcript ids. WormBase ids
  (`R107.8.1`) contain structural dots — don't split on them — and Ensembl
  appends a trailing-dot artifact (`R107.8.1.`), so `rstrip('.')` first.
- Species not on the endpoint (axolotl, planaria, some tunicates) resolve to no
  CDS — an auditable gap, handled as annotation-gap, not an error.

## Ensembl Compara — paralogue counts (redundancy, computed)

```
/homology/symbol/{species}/{sym}?type=paralogues&format=condensed
```

Count distinct within-species paralogue targets. Per-species mean across gate
families → normalise (fewer paralogues = simpler = higher). Recovers the biology
for free: zebrafish 5 paralogues (teleost WGD) = lowest simplicity; fly 0 =
highest.

## AlphaGenome — the R layer (human/mouse regulatory grammar)

API-only model (no open weights → **cannot** run on Modal). Needs a DeepMind key.

```bash
pip install alphagenome --ignore-installed packaging   # don't disturb system packaging
export ALPHAGENOME_API_KEY=...
export GRPC_DEFAULT_SSL_ROOTS_FILE_PATH=/root/.ccr/ca-bundle.crt   # gRPC trusts proxy CA
```

```python
from alphagenome.models import dna_client
from alphagenome.data import genome
model = dna_client.create(API_KEY)
out = model.predict_interval(
    interval=genome.Interval(chromosome="chr11", start=s, end=s+16384),
    organism=dna_client.Organism.HOMO_SAPIENS,            # or MUS_MUSCULUS
    requested_outputs=[dna_client.OutputType.RNA_SEQ,
                       dna_client.OutputType.ATAC,
                       dna_client.OutputType.CHIP_TF],
    ontology_terms=["UBERON:0002107"])                    # e.g. liver
rna = out.rna_seq.values          # (positions, tracks)
```

- **Supported window lengths**: 16 384 / 131 072 / 524 288 / 1 048 576. Centre on
  the TSS.
- **Assemblies**: human = GRCh38 (matches Ensembl). Mouse = **GRCm38/mm10**;
  Ensembl current is GRCm39 → lift over with
  `/map/mouse/GRCm39/{chrom}:{pos}..{pos}/GRCm38`.
- **R metric** used: human↔mouse Pearson concordance of TSS-anchored,
  strand-oriented profiles per track (RNA_SEQ/ATAC/CHIP_TF), mapped `[-1,1]→[0,1]`,
  averaged. Human is the reference (R=1); mouse's R layer = the concordance.
- **protobuf**: installing alphagenome bumps protobuf to 7.x, which breaks
  `modal` (`<7.0`). Fine because Evo2 scores are cached and the two steps run
  separately — but don't import modal and alphagenome in one process.

## Reachability snapshot (verified)

`rest.uniprot.org`, `rest.ensembl.org` (incl. Compara + `/map`), `nvcr.io`
(anonymous pull), `huggingface.co`, Modal API, AlphaGenome API — all reachable.
`biomodels.org`/EBI SBML endpoints are **blocked** (see the network-biomarker
skill for the GitHub-mirror workaround).
