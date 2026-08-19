# Variational algorithms — choosing a gradient path

Three ways to get gradients in TyxonQ. Pick on purpose; they differ by orders
of magnitude in cost and in what they model.

## 1. Autograd on the exact statevector (fastest, exact, simulation-only)

```python
import tyxonq as tq, torch
tq.set_backend("pytorch")          # c.state() now returns a differentiable tensor

def energy(params):
    c = tq.Circuit(n)
    ...                            # c.rx(i, theta=params[k]) — params are torch scalars
    return observable(c.state())   # any torch-differentiable function of the state

p = torch.tensor(init, requires_grad=True)
opt = torch.optim.Adam([p], lr=0.08)
for _ in range(200):
    opt.zero_grad(); e = energy(p); e.backward(); opt.step()
```

Verified end to end: 4-qubit TFIM (H = ΣZᵢZᵢ₊₁ − ΣXᵢ), 2-layer
hardware-efficient ansatz, 80 Adam steps → **E = −4.734** vs. exact
diagonalization **−4.7588** (0.5% — an ansatz-expressivity gap, not an
optimizer failure; more layers close it). Runnable copy in
`assets/vqe_template.py`.

Use this for development, ansatz screening, and anything where you want the
answer rather than a hardware forecast. It does not model shot noise, so it is
not a NISQ prediction.

## 2. Parameter-shift on the device runtime (hardware-honest)

The chem algorithms expose `grad="param-shift"` and `energy_and_grad(params)`;
each gradient component costs two extra circuit executions, so cost scales as
2 × n_params × (measurement groups) × shots. This is the only gradient that
also works on a real QPU. Use it when the deliverable is "what would hardware
do", and expect the optimizer to need noise-robust settings.

## 3. Gradient-free

`scipy.optimize.minimize` with COBYLA/Nelder-Mead over a counts-based energy —
see `examples/vqe_noisyopt.py` and `examples/vqe_scipy_optimization.py`.
Sensible when shots are expensive and the parameter count is small; it degrades
badly past ~20 parameters.

## Ansatz choice

| Ansatz | Depth | Parameters | When |
|---|---|---|---|
| Hardware-efficient (RY + CNOT chain) | shallow | (layers+1) × n | real hardware, generic Hamiltonians |
| UCCSD | deep | # screened excitations | chemistry accuracy on a simulator |
| kUpCCGSD | medium | k × O(n²) | chemistry with a depth budget |
| QAOA (cost + mixer layers) | p × (edges + n) | 2p | combinatorial optimization |

QAOA in TyxonQ: build ZZ couplings as `cx(u,v); rz(v, theta=2γw); cx(u,v)`, a
layer of `rx(i, theta=β)`, then `measure_z` everything and score the counts —
`examples/simple_qaoa.py` is the reference implementation.

## Barren plateaus — check before you blame the optimizer

Gradient variance decays exponentially in qubit count for deep random circuits.
Symptoms: gradient norm ~1e-3 from step 1 and flat loss. Mitigations that
actually work: shallower/structured ansatz, layerwise training, chemistry-
informed initialization (`init_method="mp2"`), small random init near zero.
`examples/barren_plateau_benchmark.py` measures it for your ansatz;
`examples/performance_layerwise_optimization.py` shows layerwise training.

## Performance rules

- **Never** loop `circuit.expectation()` over Hamiltonian terms — each call
  re-executes the circuit. Build a dense/sparse Pauli sum once with
  `pauli_string_sum_dense(terms, weights)` and contract against one `state()`
  (10-15× faster, per the framework's own docstring).
- The numeric path is bounded by memory: 2ⁿ complex amplitudes ⇒ ~20 qubits on
  a laptop, ~30 on a big machine. Beyond that use MPS and accept the bond-
  dimension approximation.
- Shot noise on an energy scales as ~1/√N. Quadrupling shots halves the error
  bar; that is often cheaper than more optimizer steps on a noisy estimate.
- `examples/psr_vs_qng_comprehensive_analysis.py`,
  `gradient_benchmark.py`, and `autograd_vs_counts.py` benchmark these paths
  against each other in the upstream repo.

## Validation habit

Every variational result should come with an independent check: exact
diagonalization (`scipy.sparse.linalg.eigsh`) for spin models, `run_fci=True`
for chemistry, or an analytic value for a toy Hamiltonian. A VQE number with
no reference is an upper bound of unknown quality.
