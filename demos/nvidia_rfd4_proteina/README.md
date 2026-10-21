# RFD4-Proteina NIM demo (hosted NVIDIA endpoint)

A small, dependency-free client for NVIDIA's hosted
[RFD4-Proteina](https://build.nvidia.com/ipd/rfd4-proteina) NIM — a generative
protein design model from the Institute for Protein Design. The demo calls the
free hosted endpoint, so it needs no GPU and no model download.

Two call shapes are covered:

| Subcommand | Task | Input dialect | What it does |
|---|---|---|---|
| `generate` | `generation` | `annotated_cif` | Samples a backbone + sequence; length follows the template |
| `scaffold` | `motif_scaffolding` | `condition_spec` | Keeps a motif's backbone, designs new flanks around it |

## Setup

Python 3.9+ and nothing else — the client uses only the standard library.

Get a key from the model page ("Get API Key"), then either export it:

```bash
export NVIDIA_API_KEY=nvapi-...
```

…or copy `.env.example` to `.env` and put the key there. `.env` is gitignored,
and the client never prints the key.

## Usage

```bash
cd demos/nvidia_rfd4_proteina

# 1. Verify the key and the route without spending GPU time.
python3 demo.py check

# 2. Unconditional generation; template is fetched from RCSB and cached.
python3 demo.py generate --pdb-id 1CRN --replicates 2 --seed 42

# 3. Motif scaffolding: keep residues A1-15, design 30 new residues each side.
python3 demo.py scaffold --pdb-id 1CRN --motif A1-15 --n-term 30 --c-term 30

# Inspect the exact request body without sending it (no key needed).
python3 demo.py generate --dry-run
```

Useful flags: `--recycle {thinking_low,thinking_medium,thinking_high}` sets the
compute budget (medium is the server default), `--steps N` overrides the
sampler's step count, `--template path.cif` uses a local file instead of RCSB,
and `--contig` takes a raw contig string for `scaffold`.

Each run writes to `outputs/<name>/`:

- `design_*.cif` — the designed structures, as annotated mmCIF
- `response.json` — the raw endpoint response
- `summary.json` — task, schedule, per-design metrics, aggregated metrics

## Verified run

Checked against the live endpoint on 2026-10-21 (server-side: NVIDIA H200,
`rfd4-proteina` 1.0.0, 60 sampling steps, 110 trunk passes at
`thinking_low`):

| Run | Result |
|---|---|
| `check` | 422 on an empty body → key accepted, route correct |
| `generate --pdb-id 1CRN --replicates 2` | 2 × 46-residue designs, 1.38 s total, 0 clashing atoms, 46 CA atoms per mmCIF |
| `scaffold --motif A1-15 --n-term 30 --c-term 30` | 75-residue design (30+15+30), 0.80 s, 0 clashing atoms, 0 bond-length violations |

Two sequences from the `generate` run, for reference:

```
MIAIILIIIAVILIIVAVKSKNKVVKIILFVVAVILIVVAIVLVIV   (replicate 0)
MSKEEKKKKIIELKKKGLSVKEIAKKLGVSKSYVSKVLKEEEKKKK   (replicate 1)
```

These are smoke tests, not validated designs: `thinking_low` with a single
replicate is the cheapest setting, and the refolding ("oracle") metrics that
would say whether a design actually folds are not served by the hosted
endpoint — they need a self-hosted NIM with a structure-prediction model
mounted.

## Endpoint notes

Things the docs leave implicit, learned by probing:

- **URL** is `POST https://health.api.nvidia.com/v1/biology/ipd/rfd4-proteina/generative`
  — note `/generative`, not `/generate` or `/predict`, and `health.api`, not
  `integrate.api` (the biology NIMs live on a different host from the LLMs).
- **Synchronous.** `NVCF-POLL-SECONDS: 300` asks NVCF to hold the connection
  and return finished designs, so there is no job id to poll.
- **mmCIF only.** Files are passed inline as `{filename, content}`; PDB format
  is rejected, and the file needs a complete `_atom_site` loop.
- **`condition_spec` always needs a structure.** `files` must hold at least one
  entry and `spec.input` must name it, and the contig must reference it with at
  least one preserved segment. A contig that designs everything from scratch
  belongs in the `annotated_cif` dialect instead.
- **Residue ranges in selections must be bracketed**: `sel('A/*/[1-15]/[N,CA,C,O]')`,
  not `sel('A/*/1-15/...')`.
- **Rate limits are tight.** The free endpoint returns `429` on even two
  concurrent requests; the client retries with exponential backoff (20 s, 40 s,
  …) and runs should be issued one at a time.
- **Contig syntax**: comma-separated segments, where a bare number is a
  designed stretch of that length, `A17-131` keeps those residues of chain A,
  and `/0` is a chain break.

Tasks the endpoint accepts beyond the two demoed here: `enzyme_design`,
`small_molecule_binder`, `antibody_design`, `partial_diffusion`,
`inverse_folding`. It also supports reward-guided inference-time search via a
`search` block (`beam_width`, `n_branch`, `step_checkpoints`), which multiplies
cost by roughly `beam_width * n_branch * (len(step_checkpoints) - 1)` rollouts
per replicate and is left out of this demo.

## References

- Model page: <https://build.nvidia.com/ipd/rfd4-proteina>
- API reference: <https://docs.api.nvidia.com/nim/reference/ipd-rfd4-proteina>
- Upstream code: <https://github.com/RosettaCommons/RFD4-Proteina>
- NVIDIA write-up on reward-guided search:
  <https://developer.nvidia.com/blog/design-proteins-with-rfd4-proteina-using-reward-guided-inference-time-search>
