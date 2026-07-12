# Running DFT/MD engines on Modal GPU from this sandbox

## Verified working
This sandbox can drive the user's Modal account (`wyh-58141`) and execute on real GPUs.
Proven end-to-end on **Tesla T4** via `modal_gpu_probe.py`
(built a custom image + ran `nvidia-smi` + numpy remotely).

## The one gotcha (solved)
Outbound traffic goes through the agent CONNECT proxy (`HTTPS_PROXY`). Modal's gRPC
(grpclib) supports proxies **natively**, but needs an extra dep:

```bash
pip install 'python-socks[asyncio]'    # a.k.a. modal[api-proxy-support]
```

After that, `modal profile current`, `modal app list`, `modal run`, `modal deploy`
all work. No /etc/hosts hacks or TCP relays needed (an earlier relay attempt was
removed — it was unnecessary).

Credentials come from env: `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET`.

## Feasibility of the materials engines on Modal
| Engine | On Modal GPU | Notes |
|---|---|---|
| CP2K | yes | `Image.from_registry("cp2k/cp2k:2024.1")` or conda-forge; GPU build exists |
| Quantum ESPRESSO | yes | conda-forge `qe` or NGC NVHPC container; OpenACC GPU |
| DeepMD-kit | yes | `deepmodeling/deepmd-kit` image; GPU training is the point |
| phonopy | yes | pip; it's light Python glue — the DFT force calc is the heavy part |
| ORCA | yes (CPU only) | free-but-registered download, cannot redistribute; no GPU benefit |
| VASP | yes, if user provides licensed source | compile in image; cannot redistribute source/license |
```
