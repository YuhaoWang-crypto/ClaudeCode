"""Build the APR (attach-pull-release) window structures + restraints with pAPRika.

APR computes an ABSOLUTE binding free energy by three legs:
  attach  : turn on host/guest restraints in the bound state
  pull    : translate the guest along an axis out of the cavity
  release : analytic correction returning the restrained guest to standard state

Reference: Henriksen, Fenley & Gilson, JCTC 2015; the pAPRika β-CD tutorials
(this is exactly the class of system pAPRika was built and benchmarked for).

This module only PREPARES the windows (topology + coordinates + restraint json).
The actual MD is run by run_apr.py and analyzed by analyze_apr.py.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.build_apr")


def _tleap_solvate(host_mol2, host_frcmod, guest_mol2, guest_frcmod,
                   cfg, work: Path) -> tuple[Path, Path]:
    """Combine host+guest, add dummy anchor atoms, solvate + neutralize with tleap.

    Returns (prmtop, inpcrd) for the bound complex in a box.
    """
    utils.require_tool("tleap")
    ff = cfg["forcefield"]
    water = ff["water"]
    water_leaprc = {"tip3p": "leaprc.water.tip3p", "opc": "leaprc.water.opc"}[water]
    salt = ff["ions"]

    prmtop = work / "complex.prmtop"
    inpcrd = work / "complex.inpcrd"
    # Dummy atoms (DM1/DM2/DM3) are the fixed APR frame; pAPRika adds these, but
    # we scaffold the tleap input so the combined system is box-ready.
    leap = f"""
source leaprc.gaff2
source {water_leaprc}
loadamberparams {Path(host_frcmod).name}
loadamberparams {Path(guest_frcmod).name}
HST = loadmol2 {Path(host_mol2).name}
GST = loadmol2 {Path(guest_mol2).name}
model = combine {{ HST GST }}
solvatebox model TIP3PBOX 14.0 iso
addionsrand model {salt['negative']} 0
addionsrand model {salt['positive']} 0
{'addionsrand model '+salt['positive']+' 0 '+salt['negative']+' 0' if salt.get('neutralize') else ''}
saveamberparm model {prmtop.name} {inpcrd.name}
quit
"""
    leap_in = work / "tleap.in"
    leap_in.write_text(leap)
    utils.run(["tleap", "-f", leap_in.name], cwd=work)
    if not prmtop.exists():
        raise RuntimeError(f"tleap did not produce {prmtop}; check {work/'leap.log'}")
    return prmtop, inpcrd


def build(config_path: str, params_dir: str, work_dir: str, guest_key: str) -> Path:
    """Set up all APR windows for one host·guest pair. Returns the windows dir."""
    cfg = utils.load_config(config_path)
    params = utils.resolve(params_dir)
    work = utils.resolve(work_dir) / f"apr_{guest_key}"
    work.mkdir(parents=True, exist_ok=True)

    host = cfg["host"]
    hp = params / "host"
    gp = params / guest_key
    host_mol2 = hp / f"{host['name']}.gaff2.mol2"
    host_frcmod = hp / f"{host['name']}.frcmod"
    guest_mol2 = gp / f"{cfg['guests'][guest_key]['name']}.gaff2.mol2"
    guest_frcmod = gp / f"{cfg['guests'][guest_key]['name']}.frcmod"
    for f in (host_mol2, host_frcmod, guest_mol2, guest_frcmod):
        if not f.exists():
            raise FileNotFoundError(f"missing parameter file {f}; run parameterize.py first")

    prmtop, inpcrd = _tleap_solvate(
        host_mol2, host_frcmod, guest_mol2, guest_frcmod, cfg, work
    )

    # --- pAPRika restraint definition -------------------------------------
    try:
        from paprika.restraints import DAT_restraint, static_DAT_restraint
        from paprika.restraints.restraints import create_window_list
    except ImportError as exc:  # keep the module importable without paprika
        raise ImportError(
            "paprika is required for build_apr; `conda env create -f environment.yml`"
        ) from exc

    apr = cfg["apr"]
    attach = [f * 100.0 for f in apr["attach_fractions"]]   # pAPRika uses % fractions
    pull = apr["pull_distances_ang"]

    # Anchor selections: 3 host atoms (define axis) + guest anchor atom(s).
    # These MUST be edited to your actual atom names after parameterization.
    host_axis = host["anchors"]["dummy_axis_atoms"]
    guest_anchor = f":{cfg['guests'][guest_key]['name'][:3].upper()}@1"  # placeholder

    windows = create_window_list_safe(attach, pull)
    manifest = {
        "guest": guest_key,
        "prmtop": str(prmtop),
        "inpcrd": str(inpcrd),
        "attach_percent": attach,
        "pull_distances_ang": pull,
        "host_axis_atoms": host_axis,
        "guest_anchor": guest_anchor,
        "restraint_k_dist": apr["restraint_k_dist"],
        "restraint_k_angle": apr["restraint_k_angle"],
        "windows": windows,
        "temperature_K": cfg["thermo"]["temperature_K"],
    }
    (work / "apr_manifest.json").write_text(json.dumps(manifest, indent=2))
    for w in windows:
        (work / "windows" / w).mkdir(parents=True, exist_ok=True)

    log.info(
        "APR build for %s: %d attach + %d pull windows -> %s",
        guest_key, len(attach), len(pull), work,
    )
    log.warning(
        "EDIT the anchor selections in apr_manifest.json (host_axis_atoms / "
        "guest_anchor) to your real Amber atom names before running MD."
    )
    return work


def create_window_list_safe(attach_pct, pull_dist) -> list[str]:
    """Deterministic window directory names: a000..aNNN then p000..pMMM."""
    names = [f"a{ i:03d}" for i in range(len(attach_pct))]
    names += [f"p{i:03d}" for i in range(len(pull_dist))]
    return names


def main() -> None:
    ap = argparse.ArgumentParser(description="Build APR windows for a host-guest pair.")
    ap.add_argument("--config", default="config/system.yaml")
    ap.add_argument("--params", default="work/params")
    ap.add_argument("--work", default="work")
    ap.add_argument("--guest", required=True, help="guest key from system.yaml (e.g. dye, pfoa)")
    args = ap.parse_args()
    build(args.config, args.params, args.work, args.guest)


if __name__ == "__main__":
    main()
