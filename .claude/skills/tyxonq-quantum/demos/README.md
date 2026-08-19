# Demos — docking, binding energies, energy levels, enzyme intermediates

Five self-contained scripts that answer "can TyxonQ do X?" by computing X and
scoring it against a classical reference or experiment. All run on CPU in
minutes. Every number below was produced by these scripts, not quoted from a
paper.

```bash
python demo1_energy_levels.py        #  ~10 s
python demo2_binding_energy.py       #  ~40 s
python demo3_proton_transfer.py      #  ~3 min
python demo4_static_correlation.py   #  ~1 min
python demo5_scaling_limits.py       #  ~1 min (add an argument to push further)
```

## What each one establishes

| Demo | System | Result | Verdict |
|---|---|---|---|
| 1 · energy levels | H₂O, aug-cc-pVDZ | Koopmans IP 13.86 eV; ΔHF 11.15; **quantum 11.13**; ΔCCSD 12.48; exp. 12.62 eV | levels and IPs run fine; a 12-qubit active space lands on HF |
| 2 · binding energy | water dimer (S22), counterpoise | HF −3.57, MP2 −4.37, CCSD −4.16, **quantum −3.57**, CCSD(T)/CBS −5.02 kcal/mol | the quantum answer *is* the HF answer; classical CCSD(T) wins |
| 3 · intermediate / barrier | Zundel H₅O₂⁺ proton transfer | barrier 0.00 kcal/mol at R(O–O)=2.4 Å → 3.44 at 2.6 → ~12 at 2.8 | mechanism-scale energetics work; CCSD still sets the accurate number (1.80 at 2.6 Å) |
| 4 · static correlation | stretched H₆, STO-3G (FCI-exact) | CCSD error up to **−50.5 mHa, below exact**; quantum error ≤ 14.6 mHa, variational | this is the regime that justifies quantum chemistry on quantum hardware |
| 5 · scaling | H₂O active-space sweep | 16 qubits 1.5 s → 20 qubits 111 s (74×); 24 qubits unfinished after 40 min; FeMoco would need 108 qubits | the wall, measured |

## Docking

There is no docking demo because there is no docking. TyxonQ has no pose search
and no scoring function, and a protein-ligand complex is 3-5 orders of
magnitude beyond a tractable active space. Use classical docking for poses and
MM-GBSA/FEP for affinity, then bring a truncated active-site cluster here if —
and only if — the chemistry is multireference. `reference/scope-and-docking.md`
has the routing table and the exact phrasing to use in a report.

## The one-paragraph summary

TyxonQ computes real quantum chemistry, correctly, up to about 20 orbitals.
Everything in demos 1-3 is dominated by dynamic correlation living outside any
such active space, so the quantum path reproduces Hartree-Fock while classical
CCSD/CCSD(T) does better for less. Demo 4 is the exception that defines the
useful scope: when the reference determinant breaks down — bond cleavage,
transition states, open-shell metal cofactors — coupled cluster fails
qualitatively and the variational quantum ansatz does not. That is the case
worth pursuing, and it is also the case that needs 60-120 logical qubits for a
real enzyme cofactor, which no device has today.
