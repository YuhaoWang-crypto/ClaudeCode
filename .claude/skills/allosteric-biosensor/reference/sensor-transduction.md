# Coupling a binder to a signal — transducer engineering

Learned from **Zhang, Ke, Zhi, Jin … Cao, "De novo design of ligand binding and sensing
with a physics based generative approach"** (bioRxiv 2026.07.13.738243, Westlake /
ProBuilder). Implemented in `biosensor_pipeline/transducer.py`.

## The lesson this paper adds
Designing a *binder* is only half a sensor. The hard, usually-skipped half is the
**transducer coupling** — turning a binding event into signal. Most de-novo binders
"emphasize static, high-affinity binding and seldom undergo ligand-induced conformational
change", so they don't work as sensors. This paper designs binders that are *sensor-ready*
and gives concrete coupling mechanisms + engineering knobs.

## Why our earlier de-novo / redesign runs failed (and the fix)
The paper contrasts two paradigms:
- **Dock a pre-made scaffold** into the ligand (the old way) — "scaffold libraries rarely
  contain structures perfectly suited for every ligand." Our Boltz de-novo/redesign runs
  against the flexible vitamin-D secosteroid are this paradigm, and they failed exactly as
  predicted (100% sample failure / non-convergence).
- **Grow the backbone around the ligand** (their ProBuilder way): enumerate favorable
  side-chain contacts with a **rotamer interaction field (RIF)** → use RIF residue backbone
  positions as **anchors** → fold outward from a "folding root" by Monte-Carlo fragment
  assembly, scoring **folding (RPX) and binding (RIF) simultaneously**, on **helical-bundle**
  topology (4–6 helices, ~120 aa), upweighting privileged H-bonds. This is the tractable
  route for arbitrary small molecules, and the topology **splits cleanly for sensors**.

Takeaway for our generative route: restrict to **helical bundles**, anchor design on
**RIF-derived key contacts** (privileged H-bond as a design constraint from the start, not a
post-hoc score), and score fold+bind jointly. Tools: RFdiffusionAA/LigandMPNN or ProBuilder;
our Boltz co-fold is the validator, not the generator.

## Three coupling mechanisms (`transducer.py`)
| mechanism | how | function |
|---|---|---|
| **A. Split-bundle reassembly** | split a helical bundle into 2 fragments, each fused to **cpGFP**; ligand binding drives reassociation → fluorescence | `build_split_bundle_sensor`, `verify_split_bundle` |
| **B. Ligand-induced folding** | a marginally-stable binder folds only when ligand-bound; folding read by **enzyme** or **FRET** pair | `build_induced_folding_sensor` |
| **C. Metal coordination** | metal (Zn tetrahedral H/H/D/E; Ni octahedral) completes a coordination sphere → induced folding | `metal_coordination_plan` |

## The knobs that actually make a split sensor work (the valuable part)
Split sensors fail from **ligand-independent self-association → high background**. The paper's
fixes, encoded as deterministic ops:
- **Fragment duplication** (`fragment_duplication`) — duplicate a segment to raise the energy
  barrier for self-association → low background, ligand-gated activation. *The key move.*
- **Interface-weakening mutations** (`interface_weakening`) — buried-hydrophobic → polar at the
  fragment interface, reducing spontaneous reassembly.
- **Terminal truncation** (`truncate_terminus`) — trim a sub-element to tune responsiveness.

`recommend_transducer(analyte_class, binder_topology, readout)` picks the mechanism; the knobs
are the dynamic-range dial (complementing our site×linker and Kd-matching knobs).

## Honesty
- ✅ every construction/knob here is deterministic and unit-tested.
- ⚠️ whether a fusion switches, its background, and its dynamic range are **hypotheses** until
  a wet-lab fluorescence/enzyme titration — the same ground truth as the rest of the skill.
