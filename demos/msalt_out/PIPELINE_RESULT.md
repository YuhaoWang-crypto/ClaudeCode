# Molten-salt ML potential on Modal GPU — WORKING end-to-end ✅

Full loop runs from the sandbox on the user's Modal account (`wyh-58141`):

```
CP2K AIMD (CPU, 8 cores)  →  DeepMD training set  →  DeepMD training (Tesla T4 GPU)
```

## Demo result (molten NaCl, 64 atoms, 7 frames, 1000 train steps)
- CP2K: PBE/GTH, DZVP-MOLOPT-SR, NVT @ 1200 K, 1 fs, Nosé → `rc=0`, 7 frames
- DeepMD `se_e2_a`, torch backend on **Tesla T4**: `train_rc=0 freeze_rc=0`
- Trained model: `model.pth` (3.1 MB)
- Energy RMSE ≈ **0.9 meV/atom**; force RMSE ≈ 1.1 eV/Å
  (force high only because 7 frames — demo scale; see `lcurve.png`)

## Working image recipes (all the gotchas, solved)

**CP2K (CPU):**
```python
cp2k_image = modal.Image.from_registry("cp2k/cp2k:2024.1", add_python="3.11")
# numpy is NOT imported at module top level — the add_python interpreter
# can't see packages pip-installed into the base image, and CP2K doesn't
# need numpy anyway. Import numpy lazily inside local/GPU functions only.
```

**DeepMD (GPU):**
```python
deepmd_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch==2.10.0", "deepmd-kit", "numpy<2", "mpich")
    .env({"MKL_THREADING_LAYER": "GNU", "MKL_SERVICE_FORCE_INTEL": "1"})
)
```
Three fixes baked in:
1. `MKL_THREADING_LAYER=GNU` — else MKL/libgomp threading conflict at import.
2. `mpich` — deepmd's pt-backend `cxx_op` needs mpich package metadata.
3. `torch==2.10.0` on a **slim** base (not the conda pytorch image) — the
   deepmd-kit 3.1.3 prebuilt op is compiled against torch 2.10.0 with
   CXX11 ABI=1; conda torch is ABI=0 and other torch versions symbol-mismatch.

## Local sandbox deps (for Modal through the agent proxy)
```bash
pip install python-socks[asyncio] aiohttp-socks   # = modal[api-proxy-support]
```
- `python-socks`: grpclib control-plane through the CONNECT proxy.
- `aiohttp-socks`: blob download of large return values (the model bytes).

## Scaling to production (郭硕's spec)
Change only parameters — the code path is identical:
- system size 50–100 atoms → `--nx/--ny/--nz` (already 64; bump to 3×3×3=216 etc.)
- 6 temperature points → loop `run_cp2k` over `temp` in {900,1000,...,1400} K
- ~10 ps equilibration + sampling → `--steps` from 6/30 to ~10000 (1 fs)
- more data → concatenate all T/frames into the DeepMD `set.000` (or multiple systems)
- bigger GPU → `gpu="A100"` / `"H100"` for faster training
- production accuracy needs proper melt equilibration before sampling (this
  demo starts from an expanded lattice and heats — fine for plumbing, not for
  final density/viscosity numbers).
