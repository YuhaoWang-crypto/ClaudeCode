# From binding energy to an electrochemical sensor: a computational descriptor pipeline

**How to reuse the UMA + PySCF/QE toolchain from this project to screen materials
and molecules for electrochemical / drug-detection sensors.**

For a sensor, the useful quantity is **not Kₐ** — binding only tells you the analyte
sticks. The measurable signal is the *electronic response to adsorption / redox*:
charge transfer, work-function change, DOS change, and the redox potential. This
document maps each of those to a concrete open-source method and demonstrates the
molecular-descriptor stage on the classic **dopamine / ascorbic acid / uric acid
(DA/AA/UA)** sensing benchmark.

---

## 1. Sensor signal ↔ descriptor ↔ method

| Sensor modality (what it reads) | Key descriptor | Engine | In this repo |
|---|---|---|---|
| **Voltammetric (CV/DPV)** peak potential | redox potential **E°(vs SHE)**; HOMO/IP (oxidation), LUMO/EA (reduction) | PySCF + implicit solvent (ddCOSMO/PCM) + thermodynamic cycle | ✅ `06_sensor_descriptors.py` |
| **Chemiresistive / FET** conductance | adsorption **charge transfer ΔQ** (Bader), doping sign (n/p) | UMA screen + QE/GPAW + Henkelman `bader` | ✅ `04_*`, `05_*` |
| **Work-function** (Kelvin probe / FET V_th) | **Δφ** = φ(slab+analyte) − φ(slab) | QE/GPAW slab (vacuum electrostatic potential − E_F) | driver in `05_*` |
| **Electronic-state** response | **DOS/PDOS** shift, Fermi-level shift, in-gap states | QE/GPAW slab PDOS | driver in `05_*` |
| **Selectivity / affinity** | adsorption energy E_ads (target vs interferents) | **UMA** (fast, large panels) | ✅ `01b_*`, `02_*` |
| Intrinsic reactivity | χ, η, ω, dipole | PySCF | ✅ `03_*`, `06_*` |

The three "electronic" descriptors an MLIP can't give (redox E, ΔQ, Δφ/DOS) all come
from open-source DFT — the same engines this project already uses.

## 2. Screening workflow

```
① UMA fast screen  (cheap, scales to hundreds of molecules)
   for each (analyte × electrode/recognition-material): relax + E_ads
   → rank by selectivity (target binds, interferents don't)
        │ top hits
        ▼
② DFT electronic fingerprint  (accurate, on the few hits)
   molecular (PySCF):  HOMO/LUMO, IP/EA, redox E (+solvent), dipole
   surface   (QE/GPAW): ΔQ (Bader), Δφ, DOS/PDOS shift, E_F shift
        ▼
③ Fingerprint table: each analyte → descriptor vector
   compare orthogonality → which material separates target from interferents
        ▼
④ (optional) use descriptors as a virtual cross-reactive sensor array,
   or ML-regress descriptors → measured response (electronic-nose pattern ID)
```

## 3. Redox potential (the electrochemistry core)

For a one-electron oxidation A → A⁺ + e⁻:

```
ΔG_ox = G(A⁺, solvated) − G(A, solvated)          (ddCOSMO/PCM water)
E_ox(vs SHE) = ΔG_ox − 4.44 V                      (4.44 V ≈ absolute SHE)
```

Needs the neutral and the ion each geometry-relaxed + a solvated single point.
Typical accuracy ~0.1–0.3 V — enough to **rank/separate** analytes, which is what
sensor selectivity requires.

## 4. Demonstration — DA / AA / UA fingerprint

`06_sensor_descriptors.py`: UMA(omol) geometries + PySCF B3LYP/6-31G*, ddCOSMO water.

| Analyte | HOMO (eV) | gap (eV) | η (eV) | ω (eV) | μ (D) | IP_gas (eV) | E_ox vs SHE (V) |
|---|---|---|---|---|---|---|---|
| **Dopamine**      | −5.45 | 5.74 | 2.87 | 1.16 | 2.05 | 7.21 | **+0.94** |
| **Uric acid**     | −5.90 | 5.10 | 2.55 | 2.20 | 3.02 | 7.68 | **+1.05** |
| **Ascorbic acid** | −6.04 | 5.69 | 2.84 | 1.79 | 7.19 | 7.84 | **+1.51** |

![DA/AA/UA sensor fingerprint](figures/sensor_fingerprint.png)

Each analyte occupies a **distinct position in descriptor space** (right panel) —
the basis on which a sensor material can discriminate them. The predicted 1e⁻
oxidation peaks (left panel) put dopamine and uric acid close (+0.94/+1.05 V) and
ascorbic acid well separated (+1.51 V).

### Honest caveat — and why it is itself useful

The **experimental** DA/AA/UA peak ordering is not reproduced exactly by this pure
1-electron oxidation, because their real oxidation is **proton-coupled (PCET,
2e⁻/2H⁺ for catechol→quinone and the enediol of ascorbate)** and therefore
**pH-dependent**. This is exactly the kind of physics a screening pipeline should
flag: to get quantitative peak potentials you must compute the **PCET** oxidation
(add proton loss via a thermodynamic cycle with G(H⁺,aq) and a −0.059·pH term),
not the bare radical cation. The bare descriptors (HOMO, η, ω, ΔQ, Δφ) remain valid
*relative* reactivity/selectivity indicators; the redox stage just needs the PCET
correction for absolute volts. Adding it is a small extension of the same cycle.

## 5. Reproduce

```bash
python 06_sensor_descriptors.py   # DA/AA/UA descriptor fingerprint
python plot_sensor.py             # figures/sensor_fingerprint.png
```
Swap the `ANALYTES` dict for your own SMILES to screen any drug/interferent panel.
Add the surface descriptors (ΔQ, Δφ, DOS) with `05_full_systems_inputs.py` on HPC.
