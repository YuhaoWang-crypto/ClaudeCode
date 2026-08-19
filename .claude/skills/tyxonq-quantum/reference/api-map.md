# TyxonQ API map (v1.2.0, all snippets executed)

## Circuit construction

`tq.Circuit(n)` — every gate method returns `self`, so chain freely.

- 1-qubit: `h x y z s sdg t tdg` — `c.h(0)`
- rotations: `rx ry rz` with **keyword** angle — `c.rx(0, theta=0.3)`
- 2-qubit: `cx/cnot cy cz swap iswap` — `c.cx(0, 1)`
- 2-qubit rotations: `rxx ryy rzz` — `c.rzz(0, 1, theta=0.4)`
- arbitrary: `c.unitary(...)`, `c.kraus(...)`, `c.add_calibration(...)`
- structure: `c.measure_z(q)`, `c.reset(q)`, `c.add_barrier()`, `c.compose(other)`,
  `c.inverse()`, `c.remap_qubits(...)`, `c.mid_measurement(...)`

Inspect: `c.draw()` (ASCII, Qiskit-style), `c.gate_count()`, `c.gate_summary()`,
`c.count_flop()`, `c.num_qubits`, `c.ops`, `c.to_json_str()` / `from_json_str`.

## Two execution modes — know which one you are in

**Numeric (exact, no shots):**

```python
psi = c.state()                     # ndarray, or torch.Tensor under set_backend("pytorch")
val = c.expectation((Z, [0]), (Z, [1]))   # (matrix, [qubits]) pairs — NOT ('Z', 0)
```

For a many-term Hamiltonian, build it once instead of calling `expectation()`
per term (each call re-executes the circuit — 10-15× slower):

```python
from tyxonq.libs.quantum_library.kernels.pauli import pauli_string_sum_dense
H = pauli_string_sum_dense([[3, 3]], [1.0])   # 0=I 1=X 2=Y 3=Z ; here Z0Z1
psi = c.state()
e = float(np.real(np.conj(psi) @ (H @ psi)))
```

**Device (sampled, counts):**

```python
res = c.compile().device(provider="simulator", device="statevector").run(shots=4096)
counts = res[0]["result"]                       # {'00': 2039, '11': 2057}
meta   = res[0]["result_meta"]                  # shots actually used, backend, status
from tyxonq.postprocessing import metrics
metrics.expectation(counts, z=[0, 1])           # ⟨Z0Z1⟩ from counts
metrics.expectation(counts, z=[0])              # ⟨Z0⟩
```

`res` is a **list**; take `res[0]`. Check `res[0]["uni_status"] == "completed"`
and `res[0]["error"]` before trusting counts, especially on hardware.

## Devices

| `device=` | provider | Use |
|---|---|---|
| `statevector` | `simulator` | default; exact or sampled |
| `density_matrix` | `simulator` | required for noise channels |
| `matrix_product_state` | `simulator` | shallow / 1D-local circuits, more qubits |
| `homebrew_s2` | `tyxonq` | real superconducting QPU (needs API key) |
| `WuYue-*` | `qcos` | China Mobile QCOS (needs wuyue SDK + access/secret keys) |
| — | `quafu`, `ibm` | third-party providers |

Global defaults: `tq.device(provider=..., device=..., shots=...)` then a bare
`c.run()`; inspect with `tq.get_device_defaults()`. Explicit per-chain
configuration is clearer in scripts you'll re-read.

## Noise

```python
c.with_noise("depolarizing", p=0.05).run(shots=2000)
c.with_noise("amplitude_damping", gamma=0.1).run(shots=1024)
c.with_noise("phase_damping", l=0.1).run(shots=1024)      # 'l', not 'lambda'
c.with_noise("pauli", px=0.01, py=0.01, pz=0.05).run(shots=1024)
```

Verified: a Bell state under `depolarizing p=0.05` over 2000 shots gives
≈ {'00': 898, '11': 921, '01': 103, '10': 78} — the off-diagonal weight is the
noise. Equivalent verbose form:
`c.device(provider="simulator", device="density_matrix", use_noise=True,
noise={"type": "depolarizing", "p": 0.05})`.

Noise parameters you invent are assumptions. Label them ⚠️ noisy-model unless
they come from `examples/noise_t1_t2_calibration.py` against a real device.

## Postprocessing

```python
from tyxonq.postprocessing.readout import ReadoutMit
mit = ReadoutMit()
corrected = mit.apply_readout_mitigation(raw_counts, method="inverse",
                                         qubits=[0, 1], shots=shots)
# or inline in the chain:
c.postprocessing(method="readout_mitigation", cals={0: A0, 1: A1},
                 mitigation="inverse")
```

Also in `tyxonq.postprocessing`: `metrics`, `counts_expval`,
`error_mitigation`, `classical_shadows`, `noise_analysis`, `io`.

## QASM / interop

```python
qasm = c.to_openqasm()          # <-- the actual OPENQASM 2.0 string
c.compile(output="qasm")        # returns a Circuit; sets the compile target only
```

Qiskit is an installed dependency, so Qiskit circuits can be converted in
(`HEA.from_qiskit_circuit`, `real_amplitudes_circuit_template_converter`).

## Hardware submission

```python
import os, tyxonq as tq
tq.set_token(os.environ["TYXONQ_API_KEY"], provider="tyxonq")
print(tq.api.list_devices(provider="tyxonq"))

res = (c.compile(output="qasm")
        .device(provider="tyxonq", device="homebrew_s2")
        .run(shots=4096))
```

Async: `c.submit_task(...)` → `tq.api.get_task_details(task, wait=True)`; also
`c.get_task_details()` / `c.get_result()` / `c.cancel()`. Get a key at
https://www.tyxonq.com. **This path is untested in this environment** (no key);
`api.tyxonq.com` responds, nothing more is verified. Never report a hardware
number without a real job id.

## Pulse level

`c.use_pulse()`, `tyxonq.waveforms`, and `examples/pulse_*.py` (~12 scripts)
cover calibration, TQASM export, three-level qudits, and pulse-level
variational control. Docs: `docs/pulse_support_en.md`.
