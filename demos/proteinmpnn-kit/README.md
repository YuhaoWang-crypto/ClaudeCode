# Demo: ProteinMPNN kit from `anthropics/uplifting-biomolecular-modeling`

[anthropics/uplifting-biomolecular-modeling](https://github.com/anthropics/uplifting-biomolecular-modeling) holds 36
drop-in **GPU inference-optimization kits** for open protein / genomics models (AlphaFold 3, Boltz-2, Chai-1,
ProteinMPNN, RFdiffusion, ESM, Evo 2, …). Each kit wraps the unmodified upstream release and adds a `--mode`:
`off` (stock), `exact` (identical output, faster), `fast`, `big`. All kits target an NVIDIA H100.

This demo ran on a **CPU-only** cloud sandbox (4 cores, 15 GB RAM, no GPU), so it uses the smallest kit,
**ProteinMPNN** (inverse folding: backbone structure → amino-acid sequences), whose stock model runs on CPU.

Reproduce: `bash run_demo.sh` (needs `git`, `gcc`, `uv`; ≈1 GB download; a few minutes).

## Results

| step | command | result |
|---|---|---|
| install | `run.sh install --weights ~/ProteinMPNN` | upstream cloned at pin `8907e667`; `WEIGHTS OK: 2/2 pinned files match` |
| dry run | `run.sh check --mode exact` | `gpu=none` → exit 3 (refuses, names the missing CUDA levers) |
| stock | `run.sh design --mode off` | 2 backbones × 8 designs in 9 s on CPU, exit 0 |
| exact | `run.sh design --mode exact` | `ERROR … no CUDA device` → exit 3; nothing designed, no silent fallback |
| partial | `… --mode exact --allow-partial` | 8 of 11 levers on CPU, 4.9 s vs 8.0 s stock wall time, exit 0 |
| tests | `pytest opt/proteinmpnn_opt/tests` | 191 passed |

Inputs: upstream's example monomers `5L33` (106 aa) and `6MRR` (68 aa); `--seed 37 --num_seq_per_target 8`,
temperature 0.1. Sequences are in `results/` (`off_*` = stock, `partial_*` = kit on CPU).

Example (5L33, stock): native sequence, then the first design (score 0.805, 48 % sequence recovery):

```
HMPEEEKAARLFIEALEKGDPELMRKVISPDTRMEDNGREFTGDEVVEYVKEIQKRGEQWHLRRYTKEGNSWRFEVQVDNNGQTEQWEVQIEVRNGRIKRVTITHV
SISPEEKKALDFIKALEKADPELMAKVIEPDTEMEVNGKKYKGEEIVEFVKEIKEKGIKYKLLSYKKEGDKYVFTVEKSYNGKTYIETIEVTVENGKIKKVVIKEE
```

## What the output comparison shows

- **5L33: byte-identical** to stock (`cmp` silent).
- **6MRR: different designs** (similar quality: sequence recovery 0.40–0.54 in the kit vs 0.46–0.54 in stock). The kit's byte-for-byte guarantee
  rests on replaying stock's RNG stream with GPU-only levers (`graph_rng`, `single_graph`, `fused_draw`), which are off
  on CPU. The kit also sorts backbones by length, so 6MRR is drawn first in the kit but second in stock. With the
  GPU levers off, the random draws for the second-in-stock backbone don't line up. That is outside what the kit
  promises: `--allow-partial` is the documented escape hatch, and the run is labeled `partial=…` in its output.

## Caveats

- torch is the CPU build of 2.5.1, not the pinned `2.5.1+cu124`; the kit reports this on its `STACK` line.
- The speed-up figures the kits are built for are H100 numbers. To see them, run on a GPU host (Docker route in
  each kit's README). The larger kits (Boltz-2, AF3, Chai-1, RFdiffusion) need a GPU and weight downloads of several GB.
