---
name: mlp-modal
description: >-
  Run DFT/AIMD + machine-learning interatomic potential (MLIP) pipelines on the
  user's Modal GPU cloud, driven from a sandboxed Claude Code session. Covers the
  full loop CP2K AIMD sampling (CPU) → DeepMD-kit training set → DeePMD-SE training
  (GPU), plus how to reach Modal through the agent proxy. Use when the user wants to
  train a machine-learning potential (DeepMD/NEP-style), run CP2K/VASP/Quantum-ESPRESSO
  AIMD or single-point DFT on Modal, build a molten-salt / alloy / oxide training set,
  or generally offload materials-science compute (DFT, MD, phonons) to Modal GPUs.
  Encodes the exact image recipes and version pins that took many iterations to find.
---

# ML interatomic potentials on Modal GPU

This skill drives the user's Modal account from a sandbox to run heavy
materials-science compute that the sandbox itself cannot (no DFT engines, no GPU).
It is proven end-to-end for **molten NaCl: CP2K AIMD → DeepMD-kit training → trained
model on a Tesla T4**.

## When to use
- Train an ML interatomic potential (DeepMD-kit `se_e2_a`, etc.) from DFT data.
- Run CP2K AIMD / single-point DFT (or QE/VASP) on Modal GPU/CPU.
- Build a training set for molten salts, alloys, oxides, etc. across temperatures.
- Any "offload this materials calculation to my Modal GPU" request.

## Prerequisites (do these first, once per session)

1. **Install the Modal SDK + proxy support** (the sandbox reaches Modal through a
   CONNECT proxy; Modal's gRPC/blob transport needs the socks helpers):
   ```bash
   pip install modal 'python-socks[asyncio]' aiohttp-socks
   ```
   `python-socks` → gRPC control plane; `aiohttp-socks` → downloading large return
   values (e.g. model bytes). Without these you get "Could not connect to the Modal
   server" or an aiohttp-socks ImportError.

2. **Credentials** come from env (`MODAL_TOKEN_ID`, `MODAL_TOKEN_SECRET`). Verify:
   ```bash
   modal profile current && modal app list
   ```

## Run the pipeline

The working pipeline is `assets/molten_salt_pipeline.py`. Copy it to the working
directory and drive it with `modal run`:

```bash
# Stage 1 only — CP2K AIMD (validate quickly)
modal run molten_salt_pipeline.py --stage aimd --steps 6 --temp 1200

# Full single-temperature loop: AIMD → convert → DeepMD GPU training
modal run molten_salt_pipeline.py --stage all --steps 30 --dp-steps 2000

# Production multi-temperature sweep (sequential CP2K, combined training set)
modal run molten_salt_pipeline.py --stage production \
    --temps "900,1000,1100,1200,1300,1400" --steps 10 --dp-steps 3000
```

Outputs land in `msalt_out/`: `model*.pth` (trained potential), `lcurve*.out`
(learning curve), `arrays*.npz` (DeepMD dataset), CP2K trajectory xyz.

## Critical implementation details

Read `reference/image-recipes.md` **before editing images** — the version pins and
env vars there are load-bearing and non-obvious (torch↔deepmd ABI, MKL threading,
mpich, numpy-in-CP2K-container). Read `reference/scaling.md` to go from demo to
production, and `reference/modal-setup.md` for proxy/connection troubleshooting.

## Honesty
The bundled demo starts from an expanded lattice and heats briefly — it proves the
**pipeline**, not production accuracy. Real density/viscosity numbers need proper
melt equilibration (ps-scale) before sampling. Say so when reporting demo results.
