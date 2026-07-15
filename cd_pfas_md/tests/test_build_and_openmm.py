"""Tests that need optional heavy deps (OpenMM, RDKit); skipped if absent.

These validate the two pieces of hand-written physics/chemistry that would
otherwise only surface on Modal:
  * the APR restraint force + unit conversion (against the real OpenMM API)
  * the modified-host graft (stoichiometry + formal charge)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT.parent))

openmm = pytest.importorskip("openmm", reason="OpenMM not installed")


def test_apr_restraint_matches_analytic_energy():
    """0.5*k*(d-r0)^2 in real OpenMM == analytic, proving groups + unit conversion."""
    import numpy as np
    from openmm import unit as u
    from cd_pfas_md.src import run_apr

    system = openmm.System()
    for _ in range(5):
        system.addParticle(12.0 * u.amu)
    L = 5.0
    system.setDefaultPeriodicBoxVectors(
        openmm.Vec3(L, 0, 0) * u.nanometer,
        openmm.Vec3(0, L, 0) * u.nanometer,
        openmm.Vec3(0, 0, L) * u.nanometer,
    )
    pos = np.array([[0, 0, 0], [0.2, 0, 0], [0.1, 0.2, 0],
                    [9.9, 9.9, 9.9], [1.0, 0, 0]]) * u.nanometer
    manifest = {
        "restraint_k_dist": 5.0,
        "attach_percent": [0.0, 100.0],
        "pull_distances_ang": [4.0, 6.0, 8.0, 10.0],
        "host_axis_index": [0, 1, 2],
        "guest_anchor_index": 4,
    }
    k, r0 = run_apr._window_k_r0(manifest, "p003")      # k=5, r0=10 A
    run_apr._apply_window_restraints(system, manifest, "p003")

    ctx = openmm.Context(system, openmm.VerletIntegrator(1.0 * u.femtosecond),
                         openmm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(pos)
    E = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)

    host = np.array([[0, 0, 0], [2, 0, 0], [1, 2, 0]])  # in A
    d = float(np.linalg.norm(np.array([10.0, 0, 0]) - host.mean(axis=0)))
    E_analytic = 0.5 * k * (d - r0) ** 2 * 4.184        # kcal -> kJ
    assert abs(E - E_analytic) / abs(E_analytic) < 1e-4


def test_restraint_requires_resolved_anchors():
    from cd_pfas_md.src import run_apr
    system = openmm.System()
    system.addParticle(12.0)
    with pytest.raises(KeyError):
        run_apr._apply_window_restraints(system, {"restraint_k_dist": 5.0,
                                                  "attach_percent": [0.0],
                                                  "pull_distances_ang": [4.0]}, "a000")


# -------------------------- RDKit graft builder ---------------------------
rdkit = pytest.importorskip("rdkit", reason="RDKit not installed")


@pytest.mark.parametrize("mod_id,exp_charge,exp_formula", [
    ("mono_6_trimethylammonium", +1, "C45H78NO34+"),
    ("fluorous_tagged", 0, "C54H73F21O35"),
    ("fluorous_cationic", +1, "C57H81F21NO34+"),
])
def test_modified_host_stoichiometry(tmp_path, mod_id, exp_charge, exp_formula):
    from cd_pfas_md.src import build_modified_host as bmh, utils
    bcd = utils.resolve("data/hosts/bcd.sdf")
    info = bmh.build_modification(mod_id, bcd, tmp_path)
    assert info["formal_charge"] == exp_charge
    assert info["formula"] == exp_formula
    # structure must round-trip
    from rdkit import Chem
    assert Chem.MolFromMolFile(info["sdf"], removeHs=False) is not None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
