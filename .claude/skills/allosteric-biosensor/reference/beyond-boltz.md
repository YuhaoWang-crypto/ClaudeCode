# Beyond Boltz — computing the dynamic range & the readout

The honest answer to "can an open-source platform compute the ON/OFF dynamic
range (DR) or the colorimetric/electrochemical signal?"

## No single tool computes DR — it is an emergent kinetic/entropic quantity

`DR = kobs(+ligand) / kobs(−ligand)`. It is set by how ligand binding, via a
change in **conformational entropy**, modulates the **reporter enzyme's catalytic
rate**. Structure prediction (Boltz/AlphaFold) gives a *static* model + a binding
confidence — it is structurally necessary-condition evidence, not DR. The readout
is just a linear reporter of turnover, so "predict the signal" ≡ "predict the
kcat modulation":

- **Colorimetric** (nitrocefin, 486 nm): absorbance = ε·[product]; the product
  accumulation rate is kcat·[E]. Beer–Lambert is trivial — the physics is kcat.
- **Electrochemical** (PQQ-GDH bioelectrode): current = n·F·(turnover) × electron-
  transfer efficiency; the ET step is Marcus (reorganization energy λ, coupling).
  Again the switch acts through kcat modulation.

So every route below targets a *piece* of the chain; none returns DR end-to-end.

## Open-source platforms and what each actually gives

| Quantity | Open-source tool | What it yields | Feasibility here |
|---|---|---|---|
| Conformational **entropy / flexibility** ΔS(apo→holo) | **OpenMM**, GROMACS (MD) + quasi-harmonic/Schlitter | the entropic driver the paper links to the switch (RMSF, ΔS, active-site pre-organization) | tractable; needs GPU + µs + ligand FF for convergence |
| Reporter **catalytic barrier** ΔG‡ (kcat of the chemical step) | **QM/MM**: CP2K, Psi4/PySCF+OpenMM, NWChem, xtb; or **EVB** (Q6) | ΔΔG‡ between ON and OFF conformers → the kcat ratio behind DR | rigorous but hard: needs *reliable* ON & OFF ensembles (this switch has no global change; the paper's own MD was inconclusive) |
| Ligand **binding free energy** ΔG_bind (Kd) | **FEP / TI / MM-GBSA** (OpenMM, Amber, Yank) | the receptor affinity the paper measures by titration | tractable; validates "receptor still binds" more rigorously than Boltz iptm |
| Electron-transfer rate (electrochemical signal) | Marcus theory via QM (λ, ΔG, H_ab) | the amperometric current per turnover | orthogonal to the switch; only relevant once kcat is known |
| Colorimetric signal | Beer–Lambert (analytic) | absorbance from [product](t) | trivial; not the bottleneck |

## Why QM/MM is the "right physics" but not the practical answer *here*

QM/MM is the correct tool for an enzyme's chemical-step barrier (β-lactamase
acylation/deacylation of nitrocefin). If you had trustworthy **ON** and **OFF**
active-site structures, ΔΔG‡(ON−OFF) would be a first-principles handle on the
kcat ratio. But this switch works by a **subtle, entropic** mechanism with **no
global conformational change** — so obtaining distinct, converged ON/OFF states
is exactly the missing input, and the paper reports its own MD was inconclusive.
QM/MM would amplify whatever conformational bias you feed it. It is a
downstream step, gated on solving the entropy/ensemble problem first.

## The pragmatic open-source validation ladder (what this repo does / enables)

1. **Boltz apo-vs-holo active-site ordering** (`coupling.py`, free, runs now):
   per-residue pLDDT at the reporter catalytic residues, apo vs holo. Tests the
   *allosteric coupling* (ligand at the receptor ⇒ more-ordered reporter active
   site) that DR depends on. ⚠️ pLDDT is confidence, not physical order.
2. **MD flexibility / entropy** (`md_entropy.py`, OpenMM): the physical route to
   ΔS(apo→holo) and active-site RMSF. The rigorous open-source test of the
   entropic mechanism. Needs GPU + µs + a parameterized ligand to converge; the
   module ships a short CPU smoke run that proves the harness works.
3. **MM-GBSA / FEP** (OpenMM/Amber): ΔG_bind vs the paper's measured Kd.
4. **QM/MM / EVB** (CP2K, Q6): ΔΔG‡ for kcat — only once (2) yields credible
   ON/OFF ensembles.
5. **Bench** — a `kobs(+L)/kobs(−L)` titration is still the sole ground-truth DR.

Rigor labels carry through: MD RMSF / QM-MM barriers are ✅ *measurements on a
model*; equating them to a *measured* DR is ⚠️ until the titration is done.
