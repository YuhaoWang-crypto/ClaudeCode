# <Study/Project name> — open-source reproduction

## 1. Overview
What the original study computed, the systems, and the deliverables to reproduce.

## 2. Methods
- fairchem UMA checkpoint + task(s) used; ASE optimizer + F_max.
- PySCF level (functional/basis) for orbitals/density; solvent model if any.
- Bader tool version.

## 3. Feasibility map (MLIP vs electronic-structure)
| Deliverable | Original method | MLIP? | Reproduced with |
|---|---|---|---|
| adsorption/binding energy | … | ✅ | UMA |
| geometry optimization | … | ✅ | UMA |
| FMO / charge density / Bader / DOS | … | ❌ | PySCF / QE |

## 4. Results (compare to original, per deliverable)
State the UMA/PySCF number next to the original, note sign/magnitude/trend
agreement, and explain expected offsets (functional differences).

## 5. Caveats
Reference-state discipline, functional offsets, cluster-model choices,
approximate ΔG, grid discretisation for Bader, size limits on CPU.

## 6. Scorecard
One row per deliverable → ✅ reproduced / ✅ pipeline-proven (needs HPC) / ❌.

## 7. Reproduce
The exact commands (setup.sh, then each script) and where outputs land.
