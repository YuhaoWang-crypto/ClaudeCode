---
name: tyxonq-quantum
description: >-
  Run quantum circuits, variational algorithms (VQE/QAOA), and quantum
  chemistry (UCCSD / HEA / kUpCCGSD / pUCCD / SQD) with the TyxonQ framework —
  on local simulators (statevector / density-matrix / MPS), with noise and
  readout mitigation, or submitted to real superconducting hardware
  (TyxonQ Homebrew_S2, China Mobile QCOS, Quafu, IBM). Use when asked to
  simulate a quantum circuit, compute a molecular ground-state energy on a
  quantum backend, build a hardware-efficient ansatz, run parameter-shift or
  autograd gradients on a quantum circuit, add shot/gate noise, export QASM,
  or send a job to a QPU. Encodes the verified install recipe, the chain-API
  contract, and the shots/QASM/backend gotchas that silently give wrong
  answers.
---

# TyxonQ — full-stack quantum programming

TyxonQ (太玄量子, QureGen-Biotech) is a modular quantum framework: stable IR →
pluggable compiler → unified device abstraction (simulator **and** real QPU) →
postprocessing. Its distinguishing feature vs. other Python quantum SDKs is that
the *same* circuit object runs on a simulator or on hardware by changing one
`provider=` string, and it ships a PySCF-style quantum-chemistry layer aimed at
drug design.

**Install and version status: verified working.** `pip install tyxonq` →
v1.2.0 on Python 3.11, CPU-only, no GPU or API key needed for everything
except the hardware path. See `reference/install.md`.

## Golden rule — the chain API

Every execution is one chain. Build once, retarget by editing one argument:

```python
import tyxonq as tq
tq.set_backend("numpy")                      # or "pytorch" / "cupynumeric"

c = tq.Circuit(2).h(0).cx(0, 1).measure_z(0).measure_z(1)

res = (c.compile()
        .device(provider="simulator", device="statevector")
        .postprocessing(method=None)
        .run(shots=4096))                    # <-- shots go HERE, see gotcha 1
counts = res[0]["result"]                    # {'00': 2039, '11': 2057}
```

Swap `provider="simulator"` → `provider="tyxonq", device="homebrew_s2"` and the
identical circuit goes to a real QPU. Nothing else changes.

## Four gotchas that silently produce wrong numbers

Each was reproduced on v1.2.0; do not skip them.

1. **`.device(shots=N).run()` is ignored — you get 1024 shots.** Only
   `.run(shots=N)` is honored. The README's own quick-start example has this
   bug: it asks for 4096 and the returned counts sum to 1024. Always pass
   `shots` to `run()`. Shot count drives your statistical error bar, so a
   silent 4× reduction quietly widens every confidence interval.
2. **No `measure_z(...)` ⇒ the compiler auto-adds Z-measurements on all qubits**
   and warns. Fine for a Bell state, wrong the moment you wanted a subset or a
   non-Z basis. Declare measurements explicitly.
3. **`compile(output="qasm")` returns a `Circuit`, not a QASM string.** It sets
   the compilation target; don't `print()` it expecting text.
4. **Gradients need the right backend.** `tq.set_backend("pytorch")` makes
   `c.state()` return a differentiable `torch.Tensor`; with `"numpy"` there is
   no autograd and you must use parameter-shift (`grad="param-shift"` in the
   chem algorithms). Picking the wrong one gives either a crash or a silently
   detached gradient.

## Pick the execution path deliberately

| Goal | Path | Cost |
|---|---|---|
| Exact energy / state, fast optimization | numeric runtime, `set_backend("pytorch")`, autograd on `c.state()` | seconds, ≤ ~20 qubits |
| Realistic NISQ estimate | device runtime, `simulator/statevector`, finite `shots` | shot noise ~1/√N |
| Noise study | `.with_noise("depolarizing", p=0.05)` (routes to `density_matrix`) | ≤ ~12 qubits |
| Large shallow / 1D circuits | `device="matrix_product_state"` | bond-dim limited |
| Real hardware | `provider="tyxonq", device="homebrew_s2"` + API key | queue + real error rates |

Rule: develop and debug on the numeric path, report on the device path, and
state which one produced every number you quote.

## Quantum chemistry in three lines

```python
from tyxonq.applications.chem.algorithms.uccsd import UCCSD
from tyxonq.applications.chem import molecule
ucc = UCCSD(molecule.h2)
e_exact  = ucc.kernel()                                              # -1.137274
e_device = ucc.kernel(shots=4096, provider="simulator", device="statevector")
```

Verified: H₂/STO-3G UCCSD gives **-1.137274 Ha** (numeric) vs. FCI -1.137270,
and **-1.1354 Ha** with 4096 shots — the gap *is* the shot noise, not a bug.
`HEA(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", layers=1, mapping="parity")
.kernel()` gives -1.137284. Full algorithm/mapping/active-space guidance in
`reference/chemistry.md`.

## Before you trust a run

Run `assets/verify_install.py` — it checks the import, the chain API, the shots
contract, noise, MPS, the chemistry stack, and the autograd path in ~30 s, and
prints each result next to its known-correct reference value. If a number
drifts, the framework or a dependency changed; fix that before interpreting
science.

## Honesty labels — apply to every number you report

- ✅ **exact** — numeric runtime, converged optimizer, checked against an
  independent solver (FCI, `scipy.sparse.linalg.eigsh`, analytic value).
- ⚠️ **sampled** — finite shots; quote shots and the ~1/√N error.
- ⚠️ **noisy-model** — depolarizing/damping parameters are *assumed*, not
  device-calibrated, unless taken from a real calibration run.
- ⚠️ **hardware** — real QPU counts; readout mitigation state must be stated.

A VQE energy is an upper bound on the ground state; never present it as "the"
energy without saying which ansatz and which path produced it.

## Reference

- `reference/install.md` — verified install, dependency set, optional extras, hardware credentials.
- `reference/api-map.md` — circuits, gates, devices, noise, postprocessing, cloud submission.
- `reference/chemistry.md` — UCCSD/HEA/kUpCCGSD/pUCCD/SQD, mappings, active spaces, what to trust.
- `reference/variational.md` — the three gradient strategies, ansatz choice, barren plateaus.
- `assets/verify_install.py` — the smoke test described above.
- `assets/vqe_template.py` — working TFIM VQE (torch autograd) validated against exact diagonalization.
