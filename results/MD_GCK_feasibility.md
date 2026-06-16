# GCK Zn-switch MD — in-environment feasibility check

What this is: the cheap, CPU-only sanity checks we *can* run in this sandbox
(OpenMM 8.5.2 + pdbfixer + Amber14, **no GPU, no AmberTools**) before handing the
real work to Colab. It is **not** the substrate-binding ΔΔG — see
`docs/MD_zn_switch_protocol.md` for the full MM-GBSA workflow.

## Models (`data/models/`)
GCK closed/active conformer `1V4S:A` + glucose (substrate, in the active site),
distal His-pair `G193H / F195H` (top hit from `results/GCK.md`, region 186–199,
~24 Å from the active centre).

| model | atoms (protein) | role |
|-------|-----------------|------|
| `GCK_WT_glc.pdb`        | 6952 | control |
| `GCK_mutHis_glc.pdb`    | 6959 | apo / no metal (MM-GBSA state 1) |
| `GCK_mutHis_Zn_glc.pdb` | 6959 + Zn | metal-bound (MM-GBSA state 2) |

## 1. Amber14 + GBn2 single-point energy (protein only)
| model | single-point E (kJ/mol) |
|-------|-------------------------|
| WT       | ~1.0 × 10⁴ |
| mutHis   | ~0.9 × 10⁴ |

Both are ~10⁴ kJ/mol = the mild strain expected of an unminimised, freshly
protonated crystal model — **orders of magnitude below the 10⁵–10⁶ of a model
with atomic clashes**. ΔE(mutant − WT) ≈ −1 × 10³ kJ/mol. Conclusion: grafting
the His-pair introduces **no catastrophic clash**; the model is geometrically
feasible. (Exact value drifts run-to-run because pdbfixer places added H/atoms
non-deterministically — these are pre-minimisation single points, not minima.)

## 2. Closest contacts at the engineered site
| model | His-sidechain ↔ rest of protein | Zn ↔ protein |
|-------|-------------------------------|--------------|
| WT       | 3.24 Å | — |
| mutHis   | 2.52 Å | — |
| mutHis+Zn| 2.52 Å | 1.86 Å |

His-sidechain↔rest 2.52 Å is a normal van der Waals contact (no hard clash <2.0 Å).
Zn↔protein 1.86 Å is near a Zn–N coordination distance — a reasonable *starting*
point; the proper coordination geometry is set in the restrained minimisation.

## What could NOT be done here (needs GPU/Colab)
- **Full energy minimisation / MD**: a 450-residue protein with GBn2 implicit
  solvent on CPU did not converge within 15 min (killed). MD must run on GPU.
- **Substrate-binding ΔΔG (the actual question)**: glucose needs GAFF2/AM1-BCC
  (antechamber) and Zn needs ZAFF/12-6 params — neither toolchain is installable
  here. Run on the making-it-rain `Protein_ligand.ipynb` per
  `docs/MD_zn_switch_protocol.md`.

## Reproduce
```bash
python src/build_md_models.py    # writes the 3 PDBs into data/models/
python src/md_check.py           # single-point energies + contact report
```
