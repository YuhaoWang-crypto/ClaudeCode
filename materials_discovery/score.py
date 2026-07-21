"""
score.py -- the pluggable energy scorer (GNoME's GNN + DFT stages).

Three back-ends, same signature `formation_energy_per_atom(composition) -> eV/atom`:

  1. surrogate_energy   -- ⚠️ transparent offline model (Pauling-ionicity). NOT DFT.
                           Lets the whole loop run/test here without weights or cloud.
  2. uma_energy         -- ✅ Meta FAIR UMA (fairchem-core) relax+energy. The modern
                           analogue of GNoME's trained GNN. Needs weights/GPU;
                           routed through the sibling `uma-crystal-mof` skill.
  3. qe_modal_energy    -- ✅ Quantum ESPRESSO SCF on Modal GPU -> DFT formation
                           energy: the active-learning verification label. Routed
                           through the sibling `qe-modal-bader-density` / `mlp-modal`
                           skills.

The discovery loop (discover.py) uses (1) or (2) as the cheap prescreen and (3) to
verify the top candidates -- exactly GNoME's cheap-predict / expensive-label split.
"""
from __future__ import annotations
from .composition import Composition, ELECTRONEG as PAULING

_ANIONS = {"O", "S", "Se", "F", "Cl", "Br", "I", "N"}


def surrogate_energy(comp: Composition) -> float:
    """⚠️ Offline surrogate formation energy (eV/atom). Ionicity heuristic:
    more electronegativity contrast between cations and anions -> more stable.
    Monotonic and deterministic; good enough to exercise the hull logic, NOT DFT."""
    frac = comp.fractional()
    cations = {el: f for el, f in frac.items() if el not in _ANIONS}
    anions = {el: f for el, f in frac.items() if el in _ANIONS}
    if not cations or not anions:
        return 0.0  # elemental / no ionic bonding in this crude model
    # mean electronegativity difference weighted by fractions
    chi_c = sum(PAULING[el] * f for el, f in cations.items()) / sum(cations.values())
    chi_a = sum(PAULING[el] * f for el, f in anions.items()) / sum(anions.values())
    dchi = abs(chi_a - chi_c)
    anion_frac = sum(anions.values())
    # Pauling-like ionic stabilization; scaled to a realistic eV/atom range
    return -0.55 * (dchi ** 2) * (2.0 * anion_frac)


def uma_energy(comp: Composition, structure_cif: str | None = None) -> float:
    """✅ UMA (fairchem-core) relaxed formation energy per atom.

    Requires fairchem-core + UMA weights (gated) and ideally a GPU -- use the
    `uma-crystal-mof` skill's environment. Implementation sketch:

        from fairchem.core import pretrained_mlip, FAIRChemCalculator
        from ase.io import read
        atoms = read(structure_cif)
        atoms.calc = FAIRChemCalculator(pretrained_mlip.get_predict_unit("uma-s-1"),
                                        task_name="omat")
        # relax with ASE, then subtract elemental references to get formation energy
    """
    raise NotImplementedError(
        "uma_energy requires fairchem-core + UMA weights; run inside the "
        "uma-crystal-mof skill environment. Use surrogate_energy for the offline demo.")


def qe_modal_energy(comp: Composition, structure_cif: str | None = None) -> float:
    """✅ DFT (Quantum ESPRESSO pw.x SCF) formation energy per atom, run on Modal GPU.

    This is GNoME's verification label. Route through the `qe-modal-bader-density`
    or `mlp-modal` skill (builds the QE image, stages SSSP pseudopotentials, runs
    pw.x). Compute E(compound)/atom minus the elemental reference energies to get a
    formation energy comparable to the hull references.
    """
    raise NotImplementedError(
        "qe_modal_energy runs pw.x on Modal; use the qe-modal-bader-density skill. "
        "Offline, the demo uses surrogate_energy.")


SCORERS = {"surrogate": surrogate_energy, "uma": uma_energy, "qe_modal": qe_modal_energy}
