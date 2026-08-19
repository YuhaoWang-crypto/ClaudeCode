"""Demo 2 — intermolecular binding energy, done properly (and its verdict for docking).

Supermolecular interaction energy of the water dimer (S22 geometry, R(O-O) =
2.910 A) with counterpoise correction:

    E_int = E(AB) - E(A in dimer basis) - E(B in dimer basis)

computed at HF, MP2 and CCSD (classical, all orbitals) and on the TyxonQ
quantum path (active-space UCCSD), against the CCSD(T)/CBS benchmark of
-5.02 kcal/mol.

This is the *right* way to ask a quantum computer for a binding energy — and
the result shows exactly why it does not yet help for drug-like complexes.

    python demo2_binding_energy.py
"""

from __future__ import annotations

import time

import numpy as np
import tyxonq as tq
from pyscf import cc, gto, mp, scf

from tyxonq.applications.chem.algorithms.uccsd import UCCSD

HA2KCAL = 627.5094740631
BASIS = "aug-cc-pvdz"
REF_CCSDT_CBS = -5.02          # S22 benchmark, water dimer

# S22 water-dimer geometry (Angstrom)
A = [("O", -1.551007, -0.114520, 0.0),
     ("H", -1.934259, 0.762503, 0.0),
     ("H", -0.599677, 0.040712, 0.0)]
B = [("O", 1.350625, 0.111469, 0.0),
     ("H", 1.680398, -0.373741, -0.758561),
     ("H", 1.680398, -0.373741, 0.758561)]


def spec(atoms, ghost: bool = False) -> str:
    tag = "ghost-" if ghost else ""
    return "\n".join(f"{tag}{s} {x} {y} {z}" for s, x, y, z in atoms)


def classical(atomspec: str):
    m = gto.M(atom=atomspec, basis=BASIS, verbose=0)
    mf = scf.RHF(m).run()
    return mf.e_tot, mp.MP2(mf).run().e_tot, cc.CCSD(mf).run().e_tot


def quantum(atomspec: str, active_space):
    m = gto.M(atom=atomspec, basis=BASIS, verbose=0)
    u = UCCSD(m, active_space=active_space, run_fci=True)
    e = u.kernel()
    assert abs(e - u.e_fci) < 1e-4, "UCCSD did not reach the CASCI limit"
    return e, u.n_qubits


def main() -> None:
    tq.set_backend("numpy")

    dimer = spec(A) + "\n" + spec(B)
    a_cp = spec(A) + "\n" + spec(B, ghost=True)      # A with B's basis functions
    b_cp = spec(A, ghost=True) + "\n" + spec(B)

    m = gto.M(atom=dimer, basis=BASIS, verbose=0)
    roo = np.linalg.norm(m.atom_coord(0) - m.atom_coord(3)) * 0.529177
    print(f"water dimer, R(O-O) = {roo:.3f} A, basis {BASIS} ({m.nao} orbitals)\n")

    print("classical references (all electrons, all orbitals):")
    d = classical(dimer)
    a = classical(a_cp)
    b = classical(b_cp)
    a0 = classical(spec(A))
    b0 = classical(spec(B))
    for i, name in enumerate(("HF", "MP2", "CCSD")):
        cp = (d[i] - a[i] - b[i]) * HA2KCAL
        raw = (d[i] - a0[i] - b0[i]) * HA2KCAL
        print(f"  {name:5s} E_int = {cp:+.3f} kcal/mol (counterpoise)   "
              f"{raw:+.3f} without CP  [BSSE {raw - cp:+.3f}]")

    print("\nquantum path (active-space UCCSD on the same geometries):")
    for as_d, as_m in (((4, 4), (2, 2)), ((8, 8), (4, 4))):
        t = time.time()
        e_d, nq = quantum(dimer, as_d)
        e_a, _ = quantum(a_cp, as_m)
        e_b, _ = quantum(b_cp, as_m)
        e_int = (e_d - e_a - e_b) * HA2KCAL
        print(f"  dimer{as_d} / monomers{as_m}  ({nq} qubits, {time.time() - t:.1f}s)"
              f"   E_int = {e_int:+.3f} kcal/mol")

    print(f"\n  benchmark CCSD(T)/CBS                                E_int = {REF_CCSDT_CBS:+.3f} kcal/mol")

    print("""
Verdict — the number to take away:

  The quantum active-space answer reproduces HF, not CCSD. A (8,8) active
  space is 8 orbitals out of 82; hydrogen-bond and dispersion binding lives in
  the dynamic correlation of ALL the remaining virtual orbitals, which the
  active space discards. Enlarging the active space to cover it means hundreds
  of qubits.

  Two more things dominate the error budget before the quantum method is even
  reached: basis-set incompleteness (aug-cc-pVDZ is ~0.9 kcal/mol from CBS)
  and BSSE. Note the trap in the table above — uncorrected CCSD lands almost
  exactly on the CCSD(T)/CBS benchmark, purely because BSSE over-binding
  cancels basis-set under-binding. Right answer, both errors intact. Any
  binding energy quoted without counterpoise and a basis-set discussion is
  luck, not accuracy.

  So: for protein-ligand affinity, use docking + MM-GBSA/FEP (classical), not
  this. See reference/scope-and-docking.md for the division of labour.""")


if __name__ == "__main__":
    main()
