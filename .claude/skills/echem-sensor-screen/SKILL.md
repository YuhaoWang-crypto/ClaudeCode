---
name: echem-sensor-screen
description: >
  Screen small molecules for electrochemical / drug-detection sensors by
  predicting the descriptors that drive a sensor's signal — adsorption energy,
  charge transfer (Bader), work-function change, DOS shift, and especially the
  redox / oxidation potential (including proton-coupled PCET). Use when the user
  wants to predict oxidation/reduction potentials, build a molecular-descriptor
  "fingerprint" for a panel of analytes, screen recognition/electrode materials,
  design a voltammetric (CV/DPV) or chemiresistive/FET sensor, or do
  drug/metabolite detection material selection with fairchem UMA + PySCF.
---

# Electrochemical sensor descriptor screening

Predict, from open-source models, the quantities an electrochemical sensor
actually reads — not just binding affinity. Pipeline: **UMA** (geometries, fast
adsorption-energy screening) + **PySCF** (electronic descriptors, redox
potentials with implicit solvent) + **QE/GPAW** (surface charge transfer, work
function, DOS for the full transducer).

## When to use
- "predict the oxidation/reduction potential of these molecules"
- "screen materials/molecules for a (dopamine / glucose / drug) sensor"
- "which analytes can this electrode distinguish?" / peak-separation questions
- "compute a descriptor fingerprint for this panel of molecules"

## Setup (once)
```bash
bash scripts/setup.sh     # torch(CPU) + fairchem-core + ase + pyscf + rdkit + matplotlib
export HF_TOKEN=...        # access to facebook/UMA (gated)
```

## Core tool — molecular descriptor fingerprint + redox potential
`scripts/sensor_screen.py` takes a JSON panel `{name: SMILES}` and computes, per
analyte: HOMO/LUMO/gap, χ/η/ω, dipole, adiabatic IP, the 1e⁻ oxidation potential,
and the **PCET 1H⁺/1e⁻ oxidation potential vs SHE at a given pH** (enumerates every
O–H/N–H, takes the most stable neutral radical). Outputs a JSON table + a figure
(predicted voltammogram + fingerprint heatmap).

```bash
cd scripts
python sensor_screen.py --analytes ../examples/analytes.json --ph 7 --outdir out
```
Swap `analytes.json` for any target drugs + interferents (SMILES) to screen them.

### Method notes (important for correctness)
- **Redox descriptor must use the right mechanism.** For catechols/phenols/enediols
  /N–H heterocycles (dopamine, ascorbic acid, uric acid, most drugs) oxidation is
  **proton-coupled (PCET)**; the bare 1e⁻ radical-cation route gives the wrong
  ordering. `sensor_screen.py` reports both — trust `E_ox_PCET_*` for peak
  potentials, use the 1e⁻ value only as a rough reactivity proxy.
- Constants: `G(H+,aq) = -11.72 eV`, `E_abs(SHE) = 4.44 V` (using 4.28 V shifts all
  values by +0.16 V uniformly — order/separation unchanged). Nernst: −0.05916·pH
  per 1H⁺/1e⁻.
- Accuracy ~0.1–0.3 V absolute; the **relative ordering / peak separation** (what
  sensor selectivity needs) is the reliable output. Validated on DA/AA/UA:
  predicted order AA < DA < UA matches experiment, values within ~0.1 V.

## Sensor signal → descriptor → engine
| Sensor modality | Descriptor | Engine |
|---|---|---|
| Voltammetric peak (CV/DPV) | redox E° (PCET, pH) | PySCF + ddCOSMO |
| Chemiresistive / FET | charge transfer ΔQ (Bader) | UMA + QE/GPAW + `bader` |
| Work-function (Kelvin/FET) | Δφ | QE/GPAW slab |
| Electronic-state | DOS/PDOS, E_F shift | QE/GPAW slab |
| Selectivity | adsorption energy | UMA (fast screen) |

For the **surface** descriptors (ΔQ, Δφ, DOS) on a real electrode, use the
companion skill `fairchem-dft-repro` (its `04_*`/`05_*` scripts do Δρ + Bader and
generate Quantum ESPRESSO inputs).

## Workflow to recommend to the user
1. UMA fast-screen adsorption energy over the (analyte × material) grid → rank by
   selectivity.
2. `sensor_screen.py` on the top hits → molecular fingerprint + PCET redox.
3. Add surface ΔQ/Δφ/DOS (fairchem-dft-repro skill, needs GPU/HPC for full size).
4. Compare fingerprints → pick the material that best separates target from
   interferents; optionally treat descriptors as a virtual cross-reactive array.

## Outputs
`out/sensor_screen.json` (descriptor table + oxidation ordering) and
`out/sensor_screen.png` (voltammogram + fingerprint heatmap).
