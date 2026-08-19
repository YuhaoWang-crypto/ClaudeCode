"""Demo 5 — measure the wall: how big an active space can this actually do?

Every claim about "quantum computing for drug design" reduces to one number:
how many spatial orbitals can you correlate. This demo measures it on your
machine instead of asserting it — UCCSD wall time and qubit count vs active
space, on a fixed molecule — and extrapolates to the active spaces real
targets need.

    python demo5_scaling_limits.py           # up to (10,10) = 20 qubits, ~1 min
    python demo5_scaling_limits.py 12        # up to (12,12) = 24 qubits, several min

Statevector simulation cost grows as 2^(2*n_orbitals) in memory alone, so the
table stops being runnable long before it stops being chemically interesting.
That gap is the entire story.

Reference measurement on one idle CPU core set: 16 qubits 1.5 s, 20 qubits
111 s, and a (12,12)/24-qubit run killed unfinished after 40 minutes.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import tyxonq as tq
from pyscf import gto

from tyxonq.applications.chem.algorithms.uccsd import UCCSD

# a fixed, chemically ordinary probe molecule with enough orbitals to slice
GEOM = "O 0.0 0.0 0.117300; H 0.0 0.757200 -0.469200; H 0.0 -0.757200 -0.469200"
BASIS = "cc-pvdz"

# active spaces a few real targets are usually quoted as needing
TARGETS = [
    ("drug-like fragment, valence only", 20),
    ("P450 Compound I (Fe-oxo core)", 30),
    ("[4Fe-4S] cluster", 40),
    ("FeMoco (nitrogenase)", 54),
    ("full enzyme active site (~100 atoms)", 200),
]


def main() -> None:
    tq.set_backend("numpy")
    max_orb = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    print(f"UCCSD scaling on H2O / {BASIS}  (statevector simulator, CPU)\n")
    print(f"{'active space':>14} {'qubits':>7} {'params':>7} {'dim(FCI vec)':>13} {'time':>9}")

    measured = []
    for n in range(2, max_orb + 1, 2):
        m = gto.M(atom=GEOM, basis=BASIS, verbose=0)
        t = time.time()
        u = UCCSD(m, active_space=(n, n))
        u.kernel()
        dt = time.time() - t
        measured.append((u.n_qubits, dt))
        print(f"{f'({n},{n})':>14} {u.n_qubits:>7} {u.n_params:>7} {2 ** u.n_qubits:>13,} {dt:>8.1f}s")

    # empirical growth per added orbital pair, from the last two points
    (q1, t1), (q2, t2) = measured[-2], measured[-1]
    factor = (t2 / t1) if t1 > 0 else float("nan")
    print(f"\nlast step: {q1} -> {q2} qubits cost {factor:.1f}x more wall time")
    print("(wall times are contention-sensitive — run this on an otherwise idle"
          " machine before quoting them)")

    print(f"\n{'target':>38} {'orbitals':>9} {'qubits':>7} {'statevector memory':>20}")
    for name, norb in TARGETS:
        nq = 2 * norb
        gib = 2 ** nq * 16 / 2 ** 30
        mem = f"{gib:.3g} GiB" if gib < 1e6 else f"10^{np.log10(gib):.0f} GiB"
        print(f"{name:>38} {norb:>9} {nq:>7} {mem:>20}")

    print("""
How to read this

  Simulation is capped by memory at roughly 20-24 orbitals in a statevector,
  and by wall time well before that. Real hardware has a different cap — the
  TyxonQ Homebrew_S2 and comparable NISQ devices offer tens of physical
  qubits with no error correction, so circuit depth, not qubit count, is what
  actually fails first for a UCCSD ansatz.

  Neither cap reaches the targets in the second table. Treat "quantum
  chemistry for this enzyme" as a roadmap statement, not something you can run
  today; what you CAN run today is a 10-20 orbital model of the reactive core,
  which is a legitimate scientific object as long as you say that is what it
  is.""")


if __name__ == "__main__":
    main()
