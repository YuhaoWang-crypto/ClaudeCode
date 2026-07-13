# Modal image recipes — the load-bearing details

These are the exact recipes that work. Each line marked ⚠️ fixes a real failure
that cost an iteration to diagnose.

## CP2K image (CPU AIMD / DFT)

```python
cp2k_image = modal.Image.from_registry("cp2k/cp2k:2024.1", add_python="3.11")
```

- ⚠️ **Do NOT import numpy (or any base-image-installed pkg) at module top level.**
  Modal imports the whole module inside *every* container to hydrate functions.
  The `add_python` standalone interpreter cannot see packages pip-installed into the
  base image's Python, so a top-level `import numpy` crashes the CP2K container even
  though CP2K never needs numpy. → import numpy lazily inside the functions that use it.
- Binary is `cp2k.psmp` at `/opt/cp2k/exe/local/`; data (`BASIS_MOLOPT`,
  `GTH_POTENTIALS`) at `/opt/cp2k/data` (`CP2K_DATA_DIR`). Locate robustly with a
  `find / -name BASIS_MOLOPT` fallback.
- Run single-rank with OMP threads (`OMP_NUM_THREADS=<cpu>`, no mpiexec) — avoids
  MPI launcher flag differences.
- ⚠️ Reserve `memory=` and use a **fresh `tempfile.mkdtemp()` workdir per call** so a
  reused warm container can't collide on output files.
- ⚠️ Concurrency: run temperatures **sequentially** via `.remote()` in a loop by
  default; tune `cpu`/`memory` and account concurrency limits before parallelizing.
- ⚠️ **SCF robustness for AIMD** (the real cause of most `MPI_Abort(1)` aborts — not
  parallelism!): OT/DIIS SCF fails to converge on liquid-like configurations as the
  MD disorders, and CP2K aborts. Fixes, in the `&SCF` block:
  `MINIMIZER CG` + `LINESEARCH 3PNT` (more robust than DIIS), `MAX_SCF 50` with
  `&OUTER_SCF MAX_SCF 30`, and `IGNORE_CONVERGENCE_FAILURE .TRUE.` as a safety net so
  one hard frame can't kill the trajectory. A run that works for 6 steps near the
  starting lattice can still fail at step 10+ — always test with enough steps.

## DeepMD-kit image (GPU training)

```python
deepmd_image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch==2.10.0", "deepmd-kit", "numpy<2", "mpich")
    .env({"MKL_THREADING_LAYER": "GNU", "MKL_SERVICE_FORCE_INTEL": "1"})
)
```

Four fixes, all required:
1. ⚠️ **Slim base, not `pytorch/pytorch`.** The conda PyTorch image is built with
   CXX11 ABI=0; the deepmd-kit prebuilt op needs ABI=1. Getting torch from PyPI on a
   slim base gives ABI=1.
2. ⚠️ **`torch==2.10.0` pinned.** deepmd-kit 3.1.3's `libdeepmd_op_pt.so` is compiled
   against a specific torch; a mismatched torch (e.g. 2.13) raises
   "version of PyTorch used to compile ... is 2.10.0, but runtime is X". Match the
   error's stated version. (Re-check the pin when deepmd-kit updates.)
3. ⚠️ **`mpich`.** The pt-backend `cxx_op` calls `load_mpi_library()` at import and
   needs the `mpich` package metadata, else `PackageNotFoundError: mpich`.
4. ⚠️ **`MKL_THREADING_LAYER=GNU`.** Else "MKL_THREADING_LAYER=INTEL is incompatible
   with libgomp.so.1" at import.

Invoke training with the torch backend: `dp --pt train input.json` then
`dp --pt freeze -o model.pth`. GPU is auto-detected via torch.

## Passing data between stages
Small arrays (coords/forces/energies) ride through the driver as pickled dict args —
no Modal Volume needed for demo scale. For large trajectories, use a `modal.Volume`
mounted in both the CP2K and DeepMD functions.
