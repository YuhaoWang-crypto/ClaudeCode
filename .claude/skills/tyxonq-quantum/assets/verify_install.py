"""TyxonQ install smoke test.

Runs every execution path the skill relies on and checks each against a known
reference value. All checks are CPU-only and finish in well under a minute.

    python verify_install.py

Every line should end in OK. A FAIL means the environment (or a TyxonQ
release) changed — fix that before interpreting any physics.
"""

from __future__ import annotations

import sys
import time

import numpy as np

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'OK  ' if ok else 'FAIL'}  {name}{'  — ' + detail if detail else ''}")
    if not ok:
        FAILURES.append(name)


def main() -> int:
    t0 = time.time()

    import tyxonq as tq

    check("import tyxonq", True, f"v{tq.__version__}")
    tq.set_backend("numpy")

    # --- 1. chain API + the shots contract ---------------------------------
    c = tq.Circuit(2).h(0).cx(0, 1).measure_z(0).measure_z(1)
    res = (
        c.compile()
        .device(provider="simulator", device="statevector")
        .postprocessing(method=None)
        .run(shots=4096)
    )
    counts = res[0]["result"]
    total = sum(counts.values())
    check("bell state counts", total == 4096 and set(counts) <= {"00", "11"},
          f"{counts}")

    # Known quirk: shots given only to .device() are ignored (falls back to 1024).
    quirk = c.compile().device(provider="simulator", device="statevector",
                               shots=4096).run()
    check("gotcha: device(shots=) is ignored", sum(quirk[0]["result"].values()) == 1024,
          "pass shots to .run() instead")

    # --- 2. counts-based expectation ---------------------------------------
    from tyxonq.postprocessing import metrics

    zz = metrics.expectation(counts, z=[0, 1])
    check("<Z0Z1> on Bell state == 1", abs(zz - 1.0) < 1e-9, f"{zz:.6f}")

    # --- 3. exact statevector + operator expectation ------------------------
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    psi = tq.Circuit(2).h(0).cx(0, 1).state()
    amp_ok = abs(abs(psi[0]) - 2 ** -0.5) < 1e-9 and abs(psi[1]) < 1e-9
    check("statevector amplitudes", amp_ok, f"{np.round(psi, 4)}")
    ev = tq.Circuit(2).h(0).cx(0, 1).expectation((Z, [0]), (Z, [1]))
    check("expectation((Z,[0]),(Z,[1]))", abs(float(np.real(ev)) - 1.0) < 1e-9)

    # --- 4. noise (density-matrix path) ------------------------------------
    noisy = tq.Circuit(2).h(0).cx(0, 1).with_noise("depolarizing", p=0.05).run(shots=4000)
    nc = noisy[0]["result"]
    leak = (nc.get("01", 0) + nc.get("10", 0)) / sum(nc.values())
    check("depolarizing noise leaks off-diagonal weight", 0.01 < leak < 0.30,
          f"P(01)+P(10) = {leak:.3f}")

    # --- 5. matrix product state simulator ----------------------------------
    m = tq.Circuit(4)
    for i in range(4):
        m.h(i)
    for i in range(3):
        m.cx(i, i + 1)
    for i in range(4):
        m.measure_z(i)
    mres = m.device(provider="simulator", device="matrix_product_state").run(shots=512)
    check("MPS simulator", sum(mres[0]["result"].values()) == 512)

    # --- 6. QASM export -----------------------------------------------------
    qasm = tq.Circuit(2).h(0).cx(0, 1).to_openqasm()
    check("to_openqasm()", isinstance(qasm, str) and "OPENQASM" in qasm)

    # --- 7. pytorch backend / autograd --------------------------------------
    tq.set_backend("pytorch")
    import torch

    p = torch.tensor([0.3], requires_grad=True)
    s = tq.Circuit(1).rx(0, theta=p[0]).state()
    prob0 = (s.conj() * s).real[0]
    prob0.backward()
    # d/dθ cos²(θ/2) = -sin(θ)/2
    check("torch autograd through circuit",
          abs(float(p.grad[0]) + np.sin(0.3) / 2) < 1e-6,
          f"grad={float(p.grad[0]):.6f}")
    tq.set_backend("numpy")

    # --- 8. quantum chemistry ------------------------------------------------
    from tyxonq.applications.chem import molecule
    from tyxonq.applications.chem.algorithms.uccsd import UCCSD

    e_h2 = UCCSD(molecule.h2).kernel()
    check("UCCSD H2/STO-3G == FCI", abs(e_h2 + 1.137274) < 2e-5, f"{e_h2:.6f} Ha")

    u = UCCSD(molecule.water(), active_space=(4, 4), run_fci=True)
    e_h2o = u.kernel()
    check("UCCSD H2O (4e,4o) vs CASCI", abs(e_h2o - u.e_fci) < 1e-4,
          f"{e_h2o:.6f} vs {u.e_fci:.6f} Ha, {u.n_qubits} qubits")

    from tyxonq.applications.chem.algorithms.hea import HEA

    e_hea = HEA(n_qubits=2, layers=2,
                hamiltonian=[(0.5, [("Z", 0)]), (0.3, [("X", 0), ("X", 1)])]).kernel()
    check("HEA vs analytic ground state", abs(e_hea + np.hypot(0.5, 0.3)) < 1e-4,
          f"{e_hea:.6f}")

    print(f"\n{'ALL CHECKS PASSED' if not FAILURES else 'FAILURES: ' + ', '.join(FAILURES)}"
          f"   ({time.time() - t0:.1f}s)")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
