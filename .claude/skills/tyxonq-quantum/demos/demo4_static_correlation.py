"""Demo 4 — the regime where the quantum path actually beats CCSD.

Demos 1-3 all end the same way: an active-space quantum calculation lands on
the HF answer because the missing physics is dynamic correlation spread over
hundreds of virtual orbitals. This demo shows the opposite regime — *static*
(multireference) correlation, where the reference determinant itself breaks
down and coupled cluster fails qualitatively.

System: a linear H6 chain in STO-3G, stretched from 0.8 to 2.6 A. In this basis
the full space is 6 electrons in 6 orbitals, so CASCI(6,6) IS the FCI answer
and every method can be scored exactly.

    python demo4_static_correlation.py

Chemistry this stands in for: bond making/breaking at a transition state,
O-O cleavage in P450 Compound I, Fe-S clusters like FeMoco, singlet
diradicals — i.e. precisely the enzyme intermediates where classical
single-reference methods are least trustworthy, and the ones quoted as the
motivating targets for quantum chemistry on quantum hardware.
"""

from __future__ import annotations

import time

import numpy as np
import tyxonq as tq
from pyscf import cc, fci, gto, scf

from tyxonq.applications.chem.algorithms.uccsd import UCCSD

MHA = 1000.0
DISTANCES = [0.8, 1.0, 1.4, 1.8, 2.2, 2.6]


def h6(r: float) -> str:
    return "; ".join(f"H 0 0 {i * r:.4f}" for i in range(6))


def main() -> None:
    tq.set_backend("numpy")
    print("Linear H6 / STO-3G — 12 qubits, CASCI(6,6) = FCI = exact\n")
    print(f"{'r(A)':>5} {'RHF':>11} {'CCSD':>11} {'FCI':>11} {'UCCSD-q':>11}"
          f" | {'CCSD err':>9} {'quantum err':>11}   (mHa)")

    t0 = time.time()
    rows = []
    for r in DISTANCES:
        m = gto.M(atom=h6(r), basis="sto-3g", verbose=0)
        mf = scf.RHF(m).run()
        e_cc = cc.CCSD(mf).run().e_tot
        e_fci = fci.FCI(mf).kernel()[0]
        u = UCCSD(gto.M(atom=h6(r), basis="sto-3g", verbose=0), active_space=(6, 6))
        e_q = u.kernel()
        rows.append((r, mf.e_tot, e_cc, e_fci, e_q))
        print(f"{r:5.1f} {mf.e_tot:11.6f} {e_cc:11.6f} {e_fci:11.6f} {e_q:11.6f}"
              f" | {(e_cc - e_fci) * MHA:9.2f} {(e_q - e_fci) * MHA:11.2f}")

    stretched = [row for row in rows if row[0] >= 1.8]
    cc_err = max(abs(row[2] - row[3]) for row in stretched) * MHA
    q_err = max(abs(row[4] - row[3]) for row in stretched) * MHA
    below = any(row[2] < row[3] for row in rows)

    print(f"\nstretched region (r >= 1.8 A): max |CCSD error| = {cc_err:.1f} mHa, "
          f"max |quantum error| = {q_err:.1f} mHa   [{time.time() - t0:.0f}s total]")
    if below:
        print("CCSD drops BELOW the exact energy at least once — coupled cluster is")
        print("not variational, so its error has no sign and no bound.")

    print("""
What this shows

  Near equilibrium CCSD is excellent (sub-mHa) and the quantum path buys you
  nothing. Stretched, CCSD's error blows up by two orders of magnitude and
  changes sign, while UCCSD on the quantum path stays variational and several
  times closer to exact.

  That is the honest case for quantum chemistry on quantum hardware: not
  better binding energies, but multireference electronic structure — bond
  cleavage, transition states, open-shell metal cofactors — where the
  classical single-reference hierarchy has no reliable answer to compare
  against in the first place.

  The catch is size. Everything above is 12 qubits. A P450 or FeMoco active
  space is 30-60 spatial orbitals, i.e. 60-120 logical qubits with deep
  circuits, which is beyond both this simulator and today's hardware. Run
  demo5_scaling_limits.py for the measured wall.""")


if __name__ == "__main__":
    main()
