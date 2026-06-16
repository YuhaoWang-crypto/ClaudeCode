"""
build_md_models.py
------------------
Build the input PDB models for a distal His-pair / Zn-switch MD experiment,
ready for the "making-it-rain" Protein_ligand notebook (Amber/OpenMM) or for an
in-environment OpenMM minimisation.

Demo target: GCK (human glucokinase), closed/active conformer 1V4S:A, which has
glucose (GLC, the substrate) already bound in the active site. We graft a distal
His-pair (top hit from results/GCK.md) and place a Zn2+ between them.

Outputs (into data/models/):
  GCK_WT_glc.pdb          WT protein + glucose                 (control)
  GCK_mutHis_glc.pdb      His-pair mutant + glucose            (apo / no metal)
  GCK_mutHis_Zn_glc.pdb   His-pair mutant + glucose + Zn2+     (metal-bound)

The two latter models are the apo vs Zn-bound pair: MM-GBSA on glucose in each,
then ddG = dG(+Zn) - dG(no Zn) estimates whether the distal metal switch
perturbs substrate binding (= a proxy for activity perturbation).
"""
from __future__ import annotations
import os
import numpy as np
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Modeller, element
from openmm import Vec3

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(ROOT, "data", "models")

SRC_PDB = os.path.join(RAW, "1V4S.pdb")
CHAIN = "A"
MUTATIONS = ["GLY-193-HIS", "PHE-195-HIS"]   # top His-pair from results/GCK.md
LIG_RESNAME = "GLC"                            # glucose = substrate in the active site


def _protein_fixer(apply_mut: bool) -> PDBFixer:
    """Protein-only PDBFixer: strip hetero, (optionally) mutate, rebuild atoms+H."""
    fixer = PDBFixer(filename=SRC_PDB)
    # keep only the protein chain we care about
    chains = [c.id for c in fixer.topology.chains()]
    fixer.removeChains([i for i, cid in enumerate(chains) if cid != CHAIN])
    if apply_mut:
        fixer.applyMutations(MUTATIONS, CHAIN)
    fixer.findMissingResidues()
    # don't try to build long missing terminal loops -> keep only internal gaps
    fixer.missingResidues = {}
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer


def _glucose_atoms():
    """Pull the glucose (GLC) atoms + coords straight out of the crystal file."""
    import prody
    prody.confProDy(verbosity="none")
    st = prody.parsePDB(SRC_PDB)
    glc = st.select(f"chain {CHAIN} and resname {LIG_RESNAME}")
    if glc is None:
        raise SystemExit("glucose not found in source structure")
    return glc


def _hetatm(serial, name, resname, chain, resnum, xyz, element):
    """Column-correct PDB HETATM line (resName 18-20, chainID 22, element 77-78)."""
    # atom name: 4-char field cols 13-16; <4-char names conventionally start col 14
    nm = name if len(name) >= 4 else " " + name
    return (f"HETATM{serial:5d} {nm:<4s} {resname:>3s} {chain:1s}{resnum:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          "
            f"{element:>2s}\n")


def _write_combined(fixer, glc, zn_xyz, path):
    """Write protein (fixer topology/positions) + glucose [+ Zn] to one PDB."""
    with open(path, "w") as fh:
        # protein first, via OpenMM writer into a buffer we post-process
        import io
        buf = io.StringIO()
        PDBFile.writeFile(fixer.topology, fixer.positions, buf, keepIds=True)
        protein_lines = [ln for ln in buf.getvalue().splitlines()
                         if ln.startswith(("ATOM", "TER"))]
        fh.write("\n".join(protein_lines) + "\n")
        # glucose as HETATM, preserving crystal coords
        ser = 90000
        coords = glc.getCoords()
        names = glc.getNames()
        for nm, xyz in zip(names, coords):
            ser += 1
            el = nm[0]
            fh.write(_hetatm(ser, nm, LIG_RESNAME, CHAIN, 9999, xyz, el))
        if zn_xyz is not None:
            ser += 1
            fh.write(_hetatm(ser, "ZN", "ZN", CHAIN, 9998, zn_xyz, "ZN"))
        fh.write("END\n")


def _zn_position(fixer):
    """Place Zn2+ ~equidistant from the two engineered His ring nitrogens."""
    mut_resnums = {int(m.split("-")[1]) for m in MUTATIONS}
    ns = {}  # resnum -> {atomname: xyz}
    pos = fixer.positions
    for atom in fixer.topology.atoms():
        r = atom.residue
        if r.chain.id == CHAIN and r.name == "HIS" and int(r.id) in mut_resnums:
            if atom.name in ("ND1", "NE2"):
                p = pos[atom.index].value_in_unit(pos.unit)
                ns.setdefault(int(r.id), {})[atom.name] = np.array([p.x, p.y, p.z]) * 10.0  # nm->A
    rids = sorted(ns)
    if len(rids) < 2:
        raise SystemExit(f"could not locate both His ring nitrogens: {ns.keys()}")
    a, b = ns[rids[0]], ns[rids[1]]
    # choose the closest N-N pair between the two His
    best = None
    for na, va in a.items():
        for nb, vb in b.items():
            d = np.linalg.norm(va - vb)
            if best is None or d < best[0]:
                best = (d, va, vb, na, nb)
    d, va, vb, na, nb = best
    zn = (va + vb) / 2.0
    return zn, d, (rids[0], na), (rids[1], nb)


def main():
    os.makedirs(OUT, exist_ok=True)
    glc = _glucose_atoms()

    # WT control
    wt = _protein_fixer(apply_mut=False)
    _write_combined(wt, glc, None, os.path.join(OUT, "GCK_WT_glc.pdb"))

    # mutant, no metal (apo)
    mut = _protein_fixer(apply_mut=True)
    _write_combined(mut, glc, None, os.path.join(OUT, "GCK_mutHis_glc.pdb"))

    # mutant + Zn
    zn, dnn, sa, sb = _zn_position(mut)
    _write_combined(mut, glc, zn, os.path.join(OUT, "GCK_mutHis_Zn_glc.pdb"))

    print("built models in", OUT)
    print(f"  His pair: {MUTATIONS}")
    print(f"  Zn placed between {sa} and {sb}; that N-N distance = {dnn:.2f} A")
    print(f"  Zn xyz (A): {zn.round(2).tolist()}")
    for f in ("GCK_WT_glc.pdb", "GCK_mutHis_glc.pdb", "GCK_mutHis_Zn_glc.pdb"):
        p = os.path.join(OUT, f)
        n = sum(1 for ln in open(p) if ln.startswith(("ATOM", "HETATM")))
        print(f"  {f}: {n} atoms")


if __name__ == "__main__":
    main()
