# Quantum chemistry with TyxonQ (all numbers below were computed, not quoted)

The chemistry layer wraps PySCF: you give a molecule, it does RHF, builds the
fermionic Hamiltonian, maps it to qubits, and variationally optimizes an ansatz
on the numeric or device runtime. `pyscf` and `openfermion` ship as
dependencies — nothing extra to install.

## Algorithms

| Class | Import from `tyxonq.applications.chem.algorithms.` | Use it for |
|---|---|---|
| `UCCSD` | `uccsd` | the default: chemically-motivated, MP2/CCSD-initialized, near-FCI in a small active space |
| `HEA` | `hea` | hardware-efficient RY + CNOT-chain; shallow enough for a real QPU; also takes an arbitrary Pauli Hamiltonian |
| `KUPCCGSD` | `kupccgsd` | k-fold paired generalized ansatz; fewer parameters than UCCSD, systematically improvable in `k` |
| `PUCCD` | `puccd` | paired doubles only — cheapest, good for seniority-zero problems |
| `UCC` | `ucc` | base class: supply your own `ex_ops` / `param_ids` |
| SQD, LUCJ | `sqd/`, `lucj/` | sample-based quantum diagonalization, local unitary cluster Jastrow |

## Minimal runs, with verified reference values

```python
import tyxonq as tq
from tyxonq.applications.chem.algorithms.uccsd import UCCSD
from tyxonq.applications.chem import molecule
tq.set_backend("numpy")

ucc = UCCSD(molecule.h2)
ucc.kernel()                                                   # -1.137274 Ha  (FCI -1.137270) ✅
ucc.kernel(shots=4096, provider="simulator", device="statevector")   # -1.1354 Ha  ⚠️ sampled
```

Active space on a real molecule — 8 qubits, 8 parameters, 0.1 s:

```python
u = UCCSD(molecule.water(), active_space=(4, 4), run_fci=True)
u.kernel()      # -74.970570 Ha   vs  u.e_fci = -74.970570 (CASCI in the same space) ✅
                #                 vs  u.e_hf  = -74.963120
```

Hardware-efficient ansatz, two ways:

```python
from tyxonq.applications.chem.algorithms.hea import HEA
HEA(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", layers=1, mapping="parity").kernel()
# -1.137284 Ha ✅

HEA(n_qubits=2, layers=2,
    hamiltonian=[(0.5, [('Z', 0)]), (0.3, [('X', 0), ('X', 1)])]).kernel()
# -0.583095  vs analytic -sqrt(0.5^2+0.3^2) = -0.583095 ✅
```

`KUPCCGSD(molecule.h2, k=2).kernel()` → -1.137274 Ha ✅.

## The `kernel()` shots default is not what you'd guess

`kernel()` with no `shots`:

- simulator/local provider → **shots = 0**, i.e. the analytic path (no sampling
  noise, so the optimizer and the RDMs are clean). This is why the numbers
  above match FCI.
- real-hardware provider → **shots = 2048**.

So `ucc.kernel()` and `ucc.kernel(shots=4096, provider="simulator", ...)` are
answering *different questions*. Say which you ran. Optimization is L-BFGS-B;
options passed to `kernel()` are forwarded to `energy_and_grad`.

## Constructor options that matter

- `active_space=(n_electrons, n_orbitals)` — the main cost lever. Qubit count
  is `2 × n_orbitals` (Jordan-Wigner). Also `active_orbital_indices=[...]` for
  explicit selection instead of energy ordering.
- `init_method` — `"mp2"` (default), `"ccsd"`, `"fe"`, `"zeros"`. MP2/CCSD
  initialization plus amplitude screening is why H₂O converges in 0.1 s;
  `"zeros"` disables screening and is much slower.
- `mapping` (HEA) — `"parity"` (saves 2 qubits for closed-shell molecules),
  `"jordan-wigner"`, `"bravyi-kitaev"`.
- `runtime` — `"device"` (default) or `"numeric"`; `numeric_engine` picks the
  numeric kernel.
- `run_fci=True` — computes the CASCI reference in the same active space. **Do
  this whenever it is affordable**: it is the only cheap way to tell an ansatz
  error from a bug.
- `atom` / `basis` / `unit` / `charge` / `spin` — build the molecule inline
  without touching PySCF.

## Molecule presets

`tyxonq.applications.chem.molecule` provides `h2`, `h_chain(n, d)`,
`h_ring`, `h_square`, `h_cube`, `water()`, `nh3()`, `bh3()`, `hcn()`,
`c2h2()`, `h2co()`, `benzene()`, `indene()`, and more. Anything else: pass a
PySCF `Mole`/`RHF` object or the `atom=` string.

## After optimization

`u.params`, `u.opt_res`, `u.energy()`, `u.energy_and_grad(params)`,
`u.make_rdm1()`, `u.make_rdm2()`, `u.get_circuit(params)`,
`u.print_circuit()`, `u.print_summary()`. RDMs need the chemistry metadata, so
they only exist on molecule-constructed instances.

## Reporting rules

- A VQE/UCCSD energy is a **variational upper bound** for the given ansatz,
  basis, and active space. Always report all four: value, ansatz, basis,
  active space.
- Compare against `run_fci=True` in the same active space before claiming
  accuracy; compare against experiment only after saying the basis set is
  small (STO-3G is not quantitative chemistry).
- "Chemical accuracy" means 1 kcal/mol = 1.6 mHa. Shot noise at 4096 shots was
  ~1.8 mHa for H₂ — i.e. *at* that threshold for the simplest molecule.
  Do not claim chemical accuracy from a sampled run without an error bar.
