# Scope: what TyxonQ can and cannot do for binding, energy levels and mechanism

Short version, with the measured evidence in `demos/`:

| Question | TyxonQ? | Why |
|---|---|---|
| Dock a small molecule into a protein pocket | **No** | Docking is a search over rigid-body poses + torsions scored by an empirical/force-field function. There is no quantum-chemistry step anywhere in it, and no docking module in the framework. |
| Protein-ligand binding free energy (ΔG, Kd) | **No** | Needs solvation, entropy and conformational sampling over 10⁴-10⁵ atoms. Quantum chemistry gives a gas-phase electronic energy for tens of atoms. |
| Interaction energy of a small complex (H-bond, dimer) | **Yes, but** | Runs today — see demo 2 — but a tractable active space reproduces HF, not CCSD, so classical CCSD(T) is strictly better and cheaper. |
| Orbital levels, HOMO-LUMO gap | **Yes** | `ucc.get_homo_lumo_gap()`; it is HF/Koopmans, i.e. mean-field. |
| Ionization energy / electron affinity | **Yes** | Total-energy differences on the quantum path — demo 1. Accuracy limited by active-space size. |
| Reaction / intermediate energy profile, barriers | **Yes, model-scale** | Demo 3 reproduces the low-barrier-hydrogen-bond mechanism. Model geometries, tens of atoms. |
| Multireference intermediates: bond cleavage, Fe-S clusters, P450 Compound I | **Yes in principle — the actual motivation** | Demo 4 shows the quantum path beating CCSD where CCSD breaks down. Real cofactors need 60-120 logical qubits; today's simulators stop near 20-24 and NISQ hardware fails on depth first. |

## The measured evidence (all reproduced in this repo, CPU, minutes)

- **Binding (demo 2)** — water dimer, counterpoise-corrected, aug-cc-pVDZ:
  HF −3.57, MP2 −4.37, CCSD −4.16, **quantum UCCSD(8,8) −3.57**,
  benchmark CCSD(T)/CBS −5.02 kcal/mol. The quantum number *is* the HF number.
- **Levels (demo 1)** — H₂O ionization energy: Koopmans 13.86, ΔHF 11.15,
  **quantum 11.13**, ΔCCSD 12.48, experiment 12.62 eV.
- **Barriers (demo 3)** — Zundel H₅O₂⁺: barrier 0.0 kcal/mol at R(O–O)=2.4 Å
  (single well) rising to ~8-12 kcal/mol at 2.8 Å. Quantum tracks HF (3.44 vs
  3.48 kcal/mol at 2.6 Å); CCSD halves it to 1.80.
- **Where quantum wins (demo 4)** — stretched H₆: max |CCSD error| 50.5 mHa
  and *below* the exact energy (non-variational failure); max |quantum error|
  14.6 mHa, variational throughout.
- **The wall (demo 5)** — on one idle CPU: 16 qubits 1.5 s, 20 qubits 111 s
  (74× for four more qubits), 24 qubits killed unfinished after 40 min.
  FeMoco's ~54-orbital active space would need 108 qubits.

One consistent conclusion: **for closed-shell, dynamic-correlation-dominated
quantities (binding energies, IPs, ordinary barriers) the quantum path buys
nothing today.** It is worth reaching for only when the reference determinant
itself fails.

## The correct division of labour for a docking + binding project

1. **Pose generation** — docking (AutoDock Vina / smina), or a co-folding model.
   → skills: `sbdd-repro-pipeline` (smina/Vina + QSAR + MD), `boltz-denovo-design`
   (Boltz-2 affinity/co-fold), `rfantibody-epitope-campaign` for biologics.
2. **Binding free energy** — MM-GBSA, alchemical FEP/TI, or APR.
   → skills: `protein-ligand-md`, `cd-pfas-md` (APR/FEP recipes),
   `mlip-surface-binding` (ML-potential binding energies, Kd/Ka).
3. **Electronic-structure refinement of a truncated site** — only where it
   changes a decision: a metal center's spin state, a covalent-inhibitor bond
   forming, a proton-transfer step, a redox potential. Cut a 20-60 atom
   cluster from the pose, cap the bonds, and run correlated quantum chemistry.
   → classical DFT/CCSD(T) first (`materials-compute`, `qe-modal-bader-density`);
   **TyxonQ here**, for the multireference cases where those are unreliable.

TyxonQ enters at step 3 and only at step 3. Sending it a protein is a category
error; sending it the reactive core of a mechanism question is exactly right.

## Honest phrasing for a report

- Do say: "the electronic energy profile of a 9-atom model of the proton-transfer
  step, computed with a 12-qubit active space, reproduces the expected
  donor-acceptor distance dependence."
- Do not say: "we computed the binding energy of the ligand with a quantum
  computer." Neither the system size nor the active space supports it, and the
  simulation ran on a CPU unless you actually submitted to hardware.
- Always state: model size, basis set, active space, whether counterpoise was
  applied, whether the run was simulator or hardware, and what the classical
  reference gives on the same geometry.
