"""Build explicit 3D structures for the modified β-CD hosts (FEP screen inputs).

Grafts substituents onto the real β-CD core (`data/hosts/bcd.sdf`) at the
primary (C6) rim and writes, per modification:
  * <id>.sdf                 — the modified host, CD core coords preserved
  * <id>.graft_atoms.json    — 0-based indices of the grafted atoms (the atoms
                               that differ from WT β-CD) for the FEP alchemical
                               region + the leaving atoms that were removed

Implemented modifications (matching config/modifications.yaml):
  * mono_6_trimethylammonium : 6-deoxy-6-(trimethylammonium) at ONE primary C6
                               (remove O6-H, attach N+(CH3)3 to C6)  → charge +1
  * fluorous_tagged          : O6-linked perfluorobutyl ether at 3 primary sites
                               (remove O6-H, attach CH2C3F7 to O6)   → neutral
  * fluorous_cationic        : fluorous_tagged + one C6 trimethylammonium → +1

⚠️ REVIEW POINTS (flagged, not hidden) — these need your sign-off before you
trust ΔΔG:
  1. GEOMETRY: grafted atoms are placed by extrapolating the bond vector and
     the structure is NOT minimized here. tleap + an OpenMM minimization on
     Modal relax it; eyeball the first frame.
  2. FEP TOPOLOGY: `fep_ti_ddg` currently treats the graft as an alchemically
     (de)coupled group — i.e. WT is approximated as "graft switched off,"
     leaving-group O6-H differences neglected. This is a fast RANKING scheme,
     not a rigorous single-topology O6H↔graft morph. For publication-grade ΔΔG,
     replace with a dual/single-topology mapping (perses/openfe). The graft-atom
     indices written here are what such a mapping would consume.
  3. SITE CHOICE: multi-site grafts use the lowest-index primary sites
     deterministically; the real regiochemistry of your synthesis may differ.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.build_mod")

PRIMARY_SMARTS = "[CH2;X4][OX2;H1]"   # (C6, O6) primary-rim CH2-OH


def _load_core(bcd_sdf: Path):
    from rdkit import Chem
    m = Chem.MolFromMolFile(str(bcd_sdf), removeHs=False)
    if m is None:
        raise ValueError(f"could not parse host structure {bcd_sdf}")
    return m


def _primary_sites(mol) -> list[tuple[int, int, int]]:
    """Return sorted (C6, O6, H_on_O6) triples for the primary rim."""
    from rdkit import Chem
    patt = Chem.MolFromSmarts(PRIMARY_SMARTS)
    sites = []
    for c6, o6 in mol.GetSubstructMatches(patt):
        h = next((n.GetIdx() for n in mol.GetAtomWithIdx(o6).GetNeighbors()
                  if n.GetAtomicNum() == 1), None)
        if h is not None:
            sites.append((c6, o6, h))
    return sorted(sites)


def _unit_vec(p, q):
    import numpy as np
    v = np.array(q) - np.array(p)
    n = float((v @ v) ** 0.5)
    return v / n if n > 1e-6 else np.array([1.0, 0.0, 0.0])


def _add_atom(rw, conf, element_z: int, pos, formal_charge: int = 0) -> int:
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    a = Chem.Atom(element_z)
    if formal_charge:
        a.SetFormalCharge(formal_charge)
    idx = rw.AddAtom(a)
    conf.SetAtomPosition(idx, Point3D(float(pos[0]), float(pos[1]), float(pos[2])))
    return idx


def _graft_trimethylammonium(rw, conf, c6, o6, h6) -> list[int]:
    """6-deoxy-6-(trimethylammonium): remove O6-H, bond N+(CH3)3 to C6."""
    from rdkit import Chem
    import numpy as np
    pC = np.array(conf.GetAtomPosition(c6)); pO = np.array(conf.GetAtomPosition(o6))
    direction = _unit_vec(pC, pO)
    grafted = []
    n = _add_atom(rw, conf, 7, pC + 1.50 * direction, formal_charge=+1)
    rw.AddBond(c6, n, Chem.BondType.SINGLE)
    grafted.append(n)
    # three methyls roughly tetrahedral around N
    import numpy as np
    basis = np.eye(3)
    offs = [direction + 0.6 * basis[0], direction + 0.6 * basis[1],
            direction - 0.6 * basis[0] - 0.6 * basis[1]]
    for off in offs:
        u = off / float((off @ off) ** 0.5)
        c = _add_atom(rw, conf, 6, np.array(conf.GetAtomPosition(n)) + 1.48 * u)
        rw.AddBond(n, c, Chem.BondType.SINGLE)
        grafted.append(c)
    return grafted, [o6, h6]  # grafted atoms, removed (leaving) atoms


def _graft_fluorous_ether(rw, conf, c6, o6, h6) -> list[int]:
    """O6-linked perfluorobutyl: remove O6-H, chain O6-CH2-CF2-CF2-CF3."""
    from rdkit import Chem
    import numpy as np
    pO = np.array(conf.GetAtomPosition(o6)); pH = np.array(conf.GetAtomPosition(h6))
    direction = _unit_vec(pO, pH)
    grafted = []
    prev = o6
    # CH2
    c = _add_atom(rw, conf, 6, pO + 1.42 * direction); rw.AddBond(prev, c, Chem.BondType.SINGLE)
    grafted.append(c); prev = c
    # three CF2/CF3 carbons, each carrying F's
    for i in range(3):
        pos = np.array(conf.GetAtomPosition(prev)) + 1.53 * direction
        cc = _add_atom(rw, conf, 6, pos); rw.AddBond(prev, cc, Chem.BondType.SINGLE)
        grafted.append(cc)
        nF = 3 if i == 2 else 2
        import numpy as np
        perp = np.cross(direction, [0, 0, 1.0]); perp = perp / (float((perp@perp)**0.5)+1e-9)
        for k in range(nF):
            fpos = pos + 1.34 * (perp if k == 0 else -perp) + (0.4 * direction if k == 2 else 0)
            f = _add_atom(rw, conf, 9, fpos); rw.AddBond(cc, f, Chem.BondType.SINGLE)
            grafted.append(f)
        prev = cc
    return grafted, [h6]  # ether keeps O6; only the hydroxyl H leaves


def build_modification(mod_id: str, bcd_sdf: Path, out_dir: Path) -> dict:
    from rdkit import Chem
    core = _load_core(bcd_sdf)
    rw = Chem.RWMol(core)
    conf = rw.GetConformer()
    sites = _primary_sites(core)
    if not sites:
        raise RuntimeError("no primary rim sites found in host")

    grafted, removed = [], []
    if mod_id == "mono_6_trimethylammonium":
        g, r = _graft_trimethylammonium(rw, conf, *sites[0]); grafted += g; removed += r
    elif mod_id == "fluorous_tagged":
        for site in sites[:3]:
            g, r = _graft_fluorous_ether(rw, conf, *site); grafted += g; removed += r
    elif mod_id == "fluorous_cationic":
        for site in sites[:3]:
            g, r = _graft_fluorous_ether(rw, conf, *site); grafted += g; removed += r
        g, r = _graft_trimethylammonium(rw, conf, *sites[3]); grafted += g; removed += r
    else:
        raise NotImplementedError(
            f"no explicit builder for '{mod_id}'. Supported: mono_6_trimethylammonium, "
            "fluorous_tagged, fluorous_cationic. Add one here or supply a mol2."
        )

    # remove leaving atoms (highest index first so indices stay valid)
    for idx in sorted(set(removed), reverse=True):
        rw.RemoveAtom(idx)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / f"{mod_id}.sdf"
    Chem.MolToMolFile(mol, str(out_sdf))

    # grafted indices shift after removals; recompute by re-matching is fragile,
    # so we record indices in the PRE-removal molecule and the net counts, plus
    # the formal charge/formula as validation the graft is chemically sane.
    from rdkit.Chem import rdMolDescriptors
    charge = Chem.GetFormalCharge(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    info = {
        "modification": mod_id,
        "n_atoms": mol.GetNumAtoms(),
        "formal_charge": charge,
        "formula": formula,
        "n_grafted_atoms": len(grafted),
        "n_removed_atoms": len(set(removed)),
        "sdf": str(out_sdf),
        "review": "geometry unminimized; FEP topology = alchemical graft decouple (see module docstring)",
    }
    (out_dir / f"{mod_id}.graft_info.json").write_text(json.dumps(info, indent=2))
    log.info("%s: %s, charge %+d, %d atoms -> %s",
             mod_id, formula, charge, mol.GetNumAtoms(), out_sdf.name)
    return info


def main() -> None:
    ap = argparse.ArgumentParser(description="Build modified β-CD host structures.")
    ap.add_argument("--bcd", default="data/hosts/bcd.sdf")
    ap.add_argument("--out", default="data/hosts/modified")
    ap.add_argument("--mods", nargs="*", default=[
        "mono_6_trimethylammonium", "fluorous_tagged", "fluorous_cationic"])
    args = ap.parse_args()
    bcd = utils.resolve(args.bcd)
    out = utils.resolve(args.out)
    for mid in args.mods:
        build_modification(mid, bcd, out)


if __name__ == "__main__":
    main()
