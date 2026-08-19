"""VQE template — transverse-field Ising model, validated against exact diagonalization.

Two execution paths on the same ansatz, so you can see exactly what the
hardware-realistic path costs you:

  1. numeric  : pytorch backend, autograd on the exact statevector (fast, exact)
  2. sampled  : counts from the statevector simulator with finite shots

    python vqe_template.py            # 4 qubits, both paths
    python vqe_template.py 6 3 200    # n_qubits layers steps

Adapt by replacing `hamiltonian_terms()` and `ansatz()`. Keep the exact-
diagonalization check: a VQE number without a reference is an upper bound of
unknown quality.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
import torch

import tyxonq as tq

I2 = sp.identity(2, format="csr")
SX = sp.csr_matrix(np.array([[0, 1], [1, 0]], dtype=float))
SZ = sp.csr_matrix(np.array([[1, 0], [0, -1]], dtype=float))


# --------------------------------------------------------------------------
# Hamiltonian:  H = sum_i Z_i Z_{i+1}  -  h * sum_i X_i     (open chain)
# --------------------------------------------------------------------------
def exact_ground_energy(n: int, h: float = 1.0) -> float:
    def kron(ops):
        m = sp.identity(1, format="csr")
        for o in ops:
            m = sp.kron(m, o, format="csr")
        return m

    H = sum(kron([SZ if j in (i, i + 1) else I2 for j in range(n)]) for i in range(n - 1))
    H = H - h * sum(kron([SX if j == i else I2 for j in range(n)]) for i in range(n))
    return float(sla.eigsh(H.astype(float), k=1, which="SA")[0][0])


def ansatz(params, n: int, layers: int, measure: bool = False):
    """Hardware-efficient ansatz: H layer, then (CNOT chain + RZ + RX) per layer."""
    c = tq.Circuit(n)
    for i in range(n):
        c.h(i)
    k = 0
    for _ in range(layers):
        for i in range(n - 1):
            c.cx(i, i + 1)
        for i in range(n):
            c.rz(i, theta=params[k]); k += 1
        for i in range(n):
            c.rx(i, theta=params[k]); k += 1
    if measure:
        for i in range(n):
            c.measure_z(i)
    return c


# --------------------------------------------------------------------------
# Path 1 — exact statevector + torch autograd
# --------------------------------------------------------------------------
def energy_numeric(params, n: int, layers: int, h: float = 1.0):
    psi = ansatz(params, n, layers).state().reshape([2] * n)

    def apply1(op, q, t):
        t = torch.movedim(t, q, 0).reshape(2, -1)
        return torch.movedim((op @ t).reshape([2] * n), 0, q)

    X = torch.tensor([[0, 1], [1, 0]], dtype=psi.dtype)
    Z = torch.tensor([[1, 0], [0, -1]], dtype=psi.dtype)
    flat = psi.reshape(-1)
    e = 0.0
    for i in range(n - 1):
        e = e + torch.vdot(flat, apply1(Z, i, apply1(Z, i + 1, psi)).reshape(-1)).real
    for i in range(n):
        e = e - h * torch.vdot(flat, apply1(X, i, psi).reshape(-1)).real
    return e


def run_numeric(n: int, layers: int, steps: int, h: float = 1.0, seed: int = 0):
    tq.set_backend("pytorch")
    torch.manual_seed(seed)
    p = torch.tensor(np.random.default_rng(seed).uniform(0, 0.3, 2 * n * layers),
                     requires_grad=True)
    opt = torch.optim.Adam([p], lr=0.08)
    for _ in range(steps):
        opt.zero_grad()
        e = energy_numeric(p, n, layers, h)
        e.backward()
        opt.step()
    return float(e.detach()), p.detach().numpy()


# --------------------------------------------------------------------------
# Path 2 — finite-shot counts (what a device would return)
# --------------------------------------------------------------------------
def energy_sampled(params, n: int, layers: int, shots: int = 8192, h: float = 1.0):
    """ZZ terms from Z-basis counts; X terms need a basis rotation (H before measure)."""
    from tyxonq.postprocessing import metrics

    tq.set_backend("numpy")
    cz = ansatz(params, n, layers, measure=True)
    counts_z = cz.compile().device(provider="simulator",
                                   device="statevector").run(shots=shots)[0]["result"]
    e = sum(metrics.expectation(counts_z, z=[i, i + 1]) for i in range(n - 1))

    cx = ansatz(params, n, layers)          # rotate X -> Z basis
    for i in range(n):
        cx.h(i)
    for i in range(n):
        cx.measure_z(i)
    counts_x = cx.compile().device(provider="simulator",
                                   device="statevector").run(shots=shots)[0]["result"]
    e -= h * sum(metrics.expectation(counts_x, z=[i]) for i in range(n))
    return float(e)


def main() -> None:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    layers = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    steps = int(sys.argv[3]) if len(sys.argv) > 3 else 80
    shots = 8192

    t = time.time()
    e_vqe, params = run_numeric(n, layers, steps)
    t_num = time.time() - t
    e_exact = exact_ground_energy(n)
    e_shot = energy_sampled(params, n, layers, shots=shots)

    print(f"TFIM  n={n}  layers={layers}  steps={steps}")
    print(f"  exact ground state        : {e_exact:+.6f}")
    print(f"  VQE, numeric autograd     : {e_vqe:+.6f}   "
          f"(err {abs(e_vqe - e_exact):.2e}, {t_num:.1f}s)   [exact path]")
    print(f"  same params, {shots} shots : {e_shot:+.6f}   "
          f"(shot noise {abs(e_shot - e_vqe):.2e})   [sampled path]")
    print("\nThe numeric-vs-exact gap is ansatz expressivity (add layers);")
    print("the sampled-vs-numeric gap is shot noise (~1/sqrt(shots)).")
    print("A sampled estimate can land BELOW the exact ground energy — it is an")
    print("unbiased noisy estimator, not a variational bound. Quote an error bar.")


if __name__ == "__main__":
    main()
