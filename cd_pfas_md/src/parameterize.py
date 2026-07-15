"""Force-field parameterization for the host and guests.

Guests (dye, PFAS) -> GAFF2 + AM1-BCC (or RESP) charges via antechamber/parmchk2.
Host (β-CD)         -> GAFF2, or an override parameter set if you adopt q4md-CD.

PFAS-SPECIFIC HANDLING
----------------------
Perfluoro chains are the single most error-prone part of this system:
  * Net charge must be the deprotonated anion (-1). We assert on it.
  * AM1-BCC underpolarizes C-F; for publication-grade numbers prefer RESP at
    HF/6-31G* (set charge_method: resp). The BCC route is fine for ranking.
  * Standard GAFF2 torsions around -CF2-CF2- can be too floppy. With
    fluorine_policy: strict we run a second parmchk2 pass at -s 2 (penalty
    reporting) and emit any parameter with a nonzero penalty to a review file
    so you can substitute values from a fluorocarbon-specific set if needed.

Nothing here needs a GPU; it is quick CPU work via AmberTools.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from . import utils

log = utils.get_logger("cd_pfas.parameterize")


def smiles_to_mol2(smiles: str, name: str, net_charge: int, out_dir: Path) -> Path:
    """Embed a 3D conformer from SMILES and write a Tripos mol2 for antechamber."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES for {name!r}: {smiles}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE) != 0:
        raise RuntimeError(f"conformer embedding failed for {name}")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)

    formal = Chem.GetFormalCharge(mol)
    if formal != net_charge:
        raise ValueError(
            f"{name}: SMILES formal charge {formal} != declared net_charge {net_charge}. "
            "Fix the protonation state in the SMILES (PFAS/dye head groups must be "
            "the deprotonated anion) before parameterizing."
        )

    out = out_dir / f"{name}.mol2"
    Chem.MolToMolFile(mol, str(out.with_suffix(".mol")))
    # Convert .mol -> .mol2 through antechamber (keeps atom typing consistent).
    utils.require_tool("antechamber")
    utils.run([
        "antechamber", "-i", out.with_suffix(".mol"), "-fi", "mdl",
        "-o", out, "-fo", "mol2", "-nc", str(net_charge), "-dr", "no",
    ], cwd=out_dir)
    log.info("%s: wrote %s (formal charge %d OK)", name, out.name, net_charge)
    return out


def structure_to_mol2(structure: str, name: str, net_charge: int, out_dir: Path) -> Path:
    """Convert an existing 3D structure (sdf/pdb/mol2) to a mol2 for antechamber.

    Used for the β-CD host, whose real 3D geometry comes from the PDB Chemical
    Component Dictionary (ligand BCD, C42H70O35) rather than RDKit embedding —
    macrocycle conformer generation is unreliable, a crystallographic-ideal
    structure is not.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    src = utils.resolve(structure)
    if not src.exists():
        raise FileNotFoundError(f"structure file not found: {src}")
    fmt = {".sdf": "sdf", ".mdl": "mdl", ".pdb": "pdb", ".mol2": "mol2",
           ".mol": "mdl"}.get(src.suffix.lower())
    if fmt is None:
        raise ValueError(f"unsupported structure format {src.suffix} for {name}")
    out = out_dir / f"{name}.mol2"
    utils.require_tool("antechamber")
    utils.run([
        "antechamber", "-i", src, "-fi", fmt,
        "-o", out, "-fo", "mol2", "-nc", str(net_charge), "-dr", "no",
    ], cwd=out_dir)
    log.info("%s: converted %s (%s) -> %s", name, src.name, fmt, out.name)
    return out


def parameterize_guest(guest: utils.Guest, out_dir: Path) -> dict[str, Path]:
    """Run antechamber (charges + GAFF2 atom types) and parmchk2 (missing params).

    Returns paths to the generated {mol2, frcmod}.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    utils.require_tool("antechamber")
    utils.require_tool("parmchk2")

    if guest.structure:
        src_mol2 = structure_to_mol2(guest.structure, guest.name, guest.net_charge, out_dir)
    elif guest.mol2:
        src_mol2 = utils.resolve(guest.mol2)
    elif guest.smiles:
        src_mol2 = smiles_to_mol2(guest.smiles, guest.name, guest.net_charge, out_dir)
    else:
        raise ValueError(f"guest {guest.key}: provide structure, mol2, or smiles")

    charge_flag = {"bcc": "bcc", "resp": "resp"}[guest.charge_method]
    if charge_flag == "resp":
        log.warning(
            "%s: RESP requested — antechamber will need a QM engine (Gaussian or, "
            "via pyRED/psi4, an external RESP fit). Ensure it is configured.",
            guest.name,
        )

    charged_mol2 = out_dir / f"{guest.name}.gaff2.mol2"
    utils.run([
        "antechamber",
        "-i", src_mol2, "-fi", "mol2",
        "-o", charged_mol2, "-fo", "mol2",
        "-c", charge_flag, "-at", "gaff2",
        "-nc", str(guest.net_charge), "-dr", "no",
    ], cwd=out_dir)

    frcmod = out_dir / f"{guest.name}.frcmod"
    args = ["parmchk2", "-i", charged_mol2, "-f", "mol2", "-o", frcmod, "-s", "gaff2"]
    utils.run(args, cwd=out_dir)

    if guest.fluorine_policy == "strict" and _has_fluorine(guest):
        _fluorine_review(charged_mol2, out_dir, guest.name)

    log.info("%s: parameterized -> %s , %s", guest.name, charged_mol2.name, frcmod.name)
    return {"mol2": charged_mol2, "frcmod": frcmod}


def _has_fluorine(guest: utils.Guest) -> bool:
    return bool(guest.smiles and "F" in guest.smiles) or guest.name.lower().startswith("pf")


def _fluorine_review(mol2: Path, out_dir: Path, name: str) -> None:
    """Second parmchk2 pass that reports parameter-assignment penalties.

    Any bonded parameter parmchk2 had to GUESS (nonzero penalty) is flagged. For
    -CF2- torsions and C-F stretches these guesses are the usual source of error;
    review the file and, if needed, splice in values from a fluorocarbon set.
    """
    review = out_dir / f"{name}.fluorine_review.frcmod"
    utils.run([
        "parmchk2", "-i", mol2, "-f", "mol2", "-o", review,
        "-s", "gaff2", "-a", "Y",  # -a Y: print ALL params with penalty scores
    ], cwd=out_dir)
    penalties = [
        ln for ln in review.read_text().splitlines()
        if "penalty score" in ln.lower() and " 0.0" not in ln
    ]
    if penalties:
        log.warning(
            "%s: %d bonded parameters were GUESSED by parmchk2 (see %s). "
            "Review CF2/CF3 torsions before trusting absolute ΔG.",
            name, len(penalties), review.name,
        )
    else:
        log.info("%s: no guessed fluorine parameters — GAFF2 covered the atom types.", name)


def parameterize_host(host_cfg: dict, temperature_k: float, out_dir: Path) -> dict[str, Path]:
    """Parameterize β-CD (or a modified host). Neutral host by default."""
    guest_like = utils.Guest.from_config("host", host_cfg, temperature_k)
    guest_like.charge_method = host_cfg.get("charge_method", "bcc")
    guest_like.name = host_cfg["name"]
    return parameterize_guest(guest_like, out_dir)


def main() -> None:
    ap = argparse.ArgumentParser(description="Parameterize host + guests from config.")
    ap.add_argument("--config", default="config/system.yaml")
    ap.add_argument("--out", default="work/params")
    args = ap.parse_args()

    cfg = utils.load_config(args.config)
    out_dir = utils.resolve(args.out)
    T = cfg["thermo"]["temperature_K"]

    parameterize_host(cfg["host"], T, out_dir / "host")
    for key, gcfg in cfg["guests"].items():
        guest = utils.Guest.from_config(key, gcfg, T)
        parameterize_guest(guest, out_dir / key)
    log.info("all parameterization complete -> %s", out_dir)


if __name__ == "__main__":
    main()
