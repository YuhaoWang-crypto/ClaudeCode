# Running AF2BIND: setup, inputs, scaling, cost

## Install

Local side (numpy only, no GPU, no jax):

```bash
pip install numpy modal
```

The GPU side is built by Modal from the image definition in
`af2bind_pipeline/modal_app.py`. Nothing needs installing by hand. The first run
builds the image (roughly 10 minutes, dominated by the CUDA wheels) and downloads
the AlphaFold parameters into a persistent Modal volume named `af2bind-params`.
Both are cached, so later runs start in under a minute.

The pinned stack is ColabDesign v1.1.1 with jax 0.4.28. ColabDesign v1.1.1 is
the commit the paper used; newer jax releases remove APIs it still calls, so do
not float these pins casually.

### Modal behind an egress proxy

If `modal` reports "Could not connect to the Modal server" while `curl` reaches
`api.modal.com` fine, the client has found an `HTTPS_PROXY` but lacks the
optional dependency needed to use it. The real error is only visible with
`MODAL_LOGLEVEL=DEBUG`. Fix:

```bash
pip install 'modal[api-proxy-support]'
```

## Inputs

`--target` accepts three things:

| Input | Example | Source |
|---|---|---|
| Local PDB file | `--target ./my_design.pdb` | used as-is |
| 4-character PDB ID | `--target 6w70` | RCSB |
| UniProt accession | `--target P00520` | AlphaFold DB, model v4 |

Only PDB format is handled. Convert mmCIF upstream if you have one. Pick the
chain with `--chain`; AF2BIND scores one chain at a time, because the 20 baits
must be appended to a single target.

### Apo, holo, and predicted structures

The method's selling point is that it needs neither a ligand nor a homologous
holo structure, so all three input types are legitimate. What differs is what a
high score means:

- **Holo** (ligand present, and stripped before running): the easiest case, and
  the one to use when you are checking that a run is behaving.
- **Apo**: the intended use. Expect somewhat lower scores than on the matched
  holo structure, because the pocket side chains have relaxed.
- **AlphaFold / predicted**: works, but trim disorder first (below).

### Trimming disorder on predicted models

AlphaFold models carry pLDDT in the B-factor column. Long disordered tails and
linkers are not pockets, but they are unpacked and solvent-exposed, which is
exactly the local environment AF2BIND rewards. They generate confident false
positives.

```bash
python -m af2bind_pipeline.run --target P00520 --min-plddt 70 --out results/abl
```

The paper goes further: the numbers behind its web interface were computed on
pLDDT-trimmed structures split into ECOD-defined domains. For a large
multi-domain protein, run each domain separately rather than trimming alone.

## GPU sizing and protein length

Memory is dominated by the (L+20)² × 128 pair representation plus Evoformer
activations, so it grows quadratically. Set the GPU with an environment variable:

```bash
AF2BIND_GPU=A100-40GB python -m af2bind_pipeline.run --target 1abc --chain A
```

| Target length | GPU that comfortably fits |
|---|---|
| up to ~400 | A10G (24 GB, the default) |
| ~400 to ~900 | L40S or A100-40GB |
| above ~900 | A100-80GB, or split into domains |

A single forward pass is seconds of GPU time; the wall clock is dominated by
container start and jax compilation. Splitting a large protein into domains is
usually both cheaper and more accurate than renting a bigger GPU, because a
domain-local pair representation is what the head was trained on.

## Cost

One target is one short GPU call, so cost per target is small change. The two
things that actually cost anything are the one-time image build and the one-time
AlphaFold parameter download into the volume. Batch work by looping over targets
in a single `app.run()` block so the container is reused.

## Re-scoring without a GPU

Every run writes `features.npz`, which holds the (L, 5120) matrix. Changing the
seed, the ensemble, or anything downstream costs no GPU at all:

```bash
python -m af2bind_pipeline.run --features results/6w70/features.npz \
    --target results/6w70/6W70.pdb --seeds 0,1,2,3,4,5,6,7,8,9 \
    --out results/6w70_ensemble
```

The file is a few megabytes per hundred residues. Keep it; re-running the GPU
pass to try a different threshold is pure waste.

## Outputs

| File | Contents |
|---|---|
| `results.csv` | every residue, ranked, with p(bind) and the raw logit |
| `report.json` | summary, top residues, pockets, ligand self-check |
| `p_bind.pdb` | input structure with p(bind) × 100 in the B-factor column |
| `bait_activations.csv` | per-bait logit decomposition for the top residues |
| `features.npz` | cached AF2 features for free re-scoring |

To look at `p_bind.pdb`:

```
load p_bind.pdb
spectrum b, blue_white_red, minimum=0, maximum=100
```

`report.json` also carries ready-made PyMOL selection strings for the top
residues and for each clustered pocket.
