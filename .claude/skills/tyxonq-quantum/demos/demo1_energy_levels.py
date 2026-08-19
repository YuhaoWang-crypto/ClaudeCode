"""Demo 1 — energy levels: orbital diagram, HOMO-LUMO gap, and a correlated ionization energy.

Three increasingly meaningful notions of "energy level", on H2O/aug-cc-pVDZ:

  1. Koopmans      — HF orbital energies. Free, but a mean-field artifact.
  2. HOMO-LUMO gap — the difference of two orbital energies. Not an observable.
  3. Delta-SCF IP  — E(cation) - E(neutral), an actual measurable ionization
                     energy, computed here on the quantum (UCCSD) path and
                     compared with Delta-HF, Delta-CCSD and experiment.

The point of the demo is the gap between (2) and (3): a HOMO-LUMO number is
cheap and often quoted, but only the total-energy difference is comparable to
experiment — and it exposes how much correlation your active space is missing.

    python demo1_energy_levels.py
"""

from __future__ import annotations

import time

import tyxonq as tq
from pyscf import cc, gto, scf

from tyxonq.applications.chem import ROUCCSD, UCCSD

EV = 27.211386245988
BASIS = "aug-cc-pvdz"
# experimental equilibrium geometry of water
GEOM = "O 0.0 0.0 0.117300; H 0.0 0.757200 -0.469200; H 0.0 -0.757200 -0.469200"
EXP_IP_EV = 12.62          # NIST vertical ionization energy of H2O


def mol(charge: int = 0, spin: int = 0) -> gto.Mole:
    return gto.M(atom=GEOM, basis=BASIS, charge=charge, spin=spin, verbose=0)


def main() -> None:
    tq.set_backend("numpy")
    print(f"H2O / {BASIS}\n")

    # --- 1+2. orbital levels and the HOMO-LUMO gap -------------------------
    u0 = UCCSD(mol(), active_space=(6, 6), init_method="zeros")
    gap = u0.get_homo_lumo_gap()
    # NOTE: v1.2.0 returns keys homo_energy/lumo_energy/gap/homo_idx/lumo_idx/
    # system_type — there is no 'gap_ev' (the upstream example script assumes
    # one and crashes). Convert yourself.
    print("Mean-field levels (Koopmans):")
    print(f"  HOMO  #{gap['homo_idx']}  {gap['homo_energy']:+.4f} Ha  "
          f"({gap['homo_energy'] * EV:+.2f} eV)")
    print(f"  LUMO  #{gap['lumo_idx']}  {gap['lumo_energy']:+.4f} Ha  "
          f"({gap['lumo_energy'] * EV:+.2f} eV)")
    print(f"  HOMO-LUMO gap        {gap['gap'] * EV:.2f} eV   "
          f"[{gap['system_type']}]")
    print("  ^ orbital-energy difference, not an observable; Koopmans IP = "
          f"{-gap['homo_energy'] * EV:.2f} eV\n")

    # --- 3. ionization energy as a total-energy difference ----------------
    print("Ionization energy from total energies (Delta-SCF style):")

    t = time.time()
    e_neutral = UCCSD(mol(), active_space=(6, 6), run_fci=True).kernel()
    e_cation = ROUCCSD(mol(1, 1), active_space=(5, 6), run_fci=True).kernel()
    ip_q = (e_cation - e_neutral) * EV
    t_q = time.time() - t

    mf0 = scf.RHF(mol()).run()
    mfp = scf.ROHF(mol(1, 1)).run()
    ip_hf = (mfp.e_tot - mf0.e_tot) * EV
    ip_cc = (cc.CCSD(mfp).run().e_tot - cc.CCSD(mf0).run().e_tot) * EV

    print(f"  quantum UCCSD/ROUCCSD (12 qubits, {t_q:.1f}s) : {ip_q:6.3f} eV")
    print(f"  Delta-HF   (classical, all orbitals)         : {ip_hf:6.3f} eV")
    print(f"  Delta-CCSD (classical, all orbitals)         : {ip_cc:6.3f} eV")
    print(f"  experiment                                   : {EXP_IP_EV:6.2f} eV")

    print(f"\nRead this honestly: the quantum active-space IP ({ip_q:.2f} eV) lands on top")
    print(f"of Delta-HF ({ip_hf:.2f} eV), not on Delta-CCSD ({ip_cc:.2f} eV). A (6,6)/(5,6)")
    print("active space out of 41 orbitals recovers almost none of the dynamic")
    print("correlation that carries the last 1.3 eV. Active-space size, not the")
    print("quantum algorithm, is the accuracy bottleneck here.")


if __name__ == "__main__":
    main()
