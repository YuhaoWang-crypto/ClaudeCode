# Installing TyxonQ — verified recipe

Verified 2026-08 on Linux x86-64, Python 3.11.15, CPU-only container, no GPU.

## The install

```bash
uv venv tq-env --python 3.11
uv pip install --python tq-env/bin/python tyxonq        # -> tyxonq 1.2.0
# plain pip works too:  pip install tyxonq
```

That single command is enough for **everything except real hardware**: circuits,
noise, MPS, VQE/QAOA, and the whole quantum-chemistry stack. No API key, no
compiler toolchain, no CUDA.

Python 3.10+ is claimed by the project; 3.11 is what QCOS requires and what was
tested here — prefer 3.11.

## What comes with it (pulled automatically)

`numpy`, `scipy`, `torch` 2.13 (+`triton`, `torchvision`), `qiskit` 2.5.2,
`pyscf` 2.14, `openfermion`, `renormalizer`, `sympy`/`symengine`, `networkx`,
`rustworkx`, `tqdm`, `requests`, `python-dotenv`, `ruff`.

Consequences worth knowing:

- **PySCF and OpenFermion are already there** — no separate install for the
  chemistry examples, contrary to the README comment `# pip install pyscf`.
- **Torch is a hard dependency**, so the download is large (~2-3 GB with
  CUDA wheels) even if you only want NumPy circuits. Budget disk space; in a
  size-constrained container install `torch` CPU-only first from the PyTorch
  CPU index, then `tyxonq`.
- Qiskit is present, so `HEA.from_qiskit_circuit(...)` and QASM interop work
  out of the box.

## From source (needed for QCOS, or to read/patch internals)

```bash
git clone --depth 1 https://github.com/QureGenAI-Biotech/TyxonQ.git
cd TyxonQ && uv build && uv pip install dist/tyxonq-*.whl
```

The source tree is worth cloning regardless: `examples/` has ~80 runnable
scripts (VQE variants, QAOA, pulse control, readout mitigation, MPS
benchmarks, GPT-QE drug design) that are the real documentation.

## Optional / gated extras

| Feature | Extra requirement | Status here |
|---|---|---|
| TyxonQ QPU (`homebrew_s2`) | API key from https://www.tyxonq.com | **untested** — no key; `api.tyxonq.com` is reachable (HTTP 302) |
| China Mobile QCOS | `wuyue-1.0-py3-none-any.whl` from the ecloud console, Python 3.11, then reinstall tyxonq from source | untested |
| Quafu / IBM providers | provider-specific tokens | untested |
| GPU numerics | `cupynumeric` backend | untested (no GPU) |

Never claim a hardware result without an actual submitted job id.

## Credentials

```python
import tyxonq as tq, os
tq.set_token(os.environ["TYXONQ_API_KEY"], provider="tyxonq")
print(tq.api.list_devices(provider="tyxonq"))
```

Read the key from an environment variable or `getpass`. Do not hard-code a
token into a script, a notebook, or anything that gets committed.

## Sanity check

```bash
tq-env/bin/python .claude/skills/tyxonq-quantum/assets/verify_install.py
```

Expect every line to print `OK`. Anything else means the environment, not your
physics, is what needs fixing.
