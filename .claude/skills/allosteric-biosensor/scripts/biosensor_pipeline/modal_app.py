"""
modal_app.py -- run the apo/holo MD entropy protocol on a cloud GPU via Modal.

The production MD (md_entropy.run_leg / compare_apo_holo) needs a GPU and the
conda-forge small-molecule-FF stack (openff-toolkit + openmmforcefields). Modal
(https://modal.com) provides both: this script defines a GPU container image
with those deps and runs the two legs (apo, holo) in parallel, then reports
ΔS(holo−apo) and the active-site RMSF change.

Setup (once):
    pip install modal && modal token new

Run (from the repo root; needs a clean protein PDB per leg + the ligand SDF):
    modal run biosensor_pipeline/modal_app.py \
        --apo-pdb  biosensor_out/vitd_apo.pdb \
        --holo-pdb biosensor_out/vitd_holo_protein.pdb \
        --ligand-sdf biosensor_out/vitd_ligand.sdf \
        --smiles "CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C" \
        --prod-ns 50

Prepare inputs from a Boltz model with `cif_to_pdb()` below (strips the ligand
for the protein PDB; the holo ligand comes from Boltz's sample_0_ligands.sdf).

HONESTY: ΔS and RMSF are ✅ measurements on the trajectories. Mapping ΔS to a
dynamic-range number is ⚠️ until calibrated against a wet-lab kobs assay. ΔS<0
(holo more ordered) supports the paper's entropic-switch mechanism.
"""

from __future__ import annotations

try:
    import modal
except ImportError:  # let the file import for docs/tests without modal installed
    modal = None


IMAGE = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "openmm=8.1", "openmmforcefields", "openff-toolkit", "mdtraj",
        "numpy", "biotite", channels=["conda-forge"],
    )
    .add_local_python_source("biosensor_pipeline")
    if modal is not None else None
)

app = modal.App("biosensor-md", image=IMAGE) if modal is not None else None


# --- local helper: Boltz CIF -> clean protein-only PDB (run locally) ----------
def cif_to_pdb(cif_path: str, out_pdb: str, drop_ligand: bool = True) -> str:
    import warnings; warnings.filterwarnings("ignore")
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure.io.pdb as biopdb
    import biotite.structure as struc
    from .md_entropy import _add_cterm_oxt
    cif = (pdbx.CIFFile.read(cif_path) if cif_path.lower().endswith((".cif", ".mmcif"))
           else pdbx.BinaryCIFFile.read(cif_path))
    arr = pdbx.get_structure(cif, model=1)
    if drop_ligand:
        arr = arr[struc.filter_amino_acids(arr)]
    chain = max(set(arr.chain_id), key=lambda c: (arr.chain_id == c).sum())
    arr = arr[(arr.chain_id == chain) & (arr.element != "H")]
    arr = _add_cterm_oxt(arr, struc)
    f = biopdb.PDBFile(); f.set_structure(arr); f.write(out_pdb)
    return out_pdb


if modal is not None:

    @app.function(gpu="A10G", timeout=3600)
    def leg_remote(protein_text: str, ligand_sdf_text: str = None, ligand_smiles: str = None,
                   solvent: str = "explicit", equil_ps: float = 200.0, prod_ns: float = 20.0,
                   catalytic_resids: list = None) -> dict:
        """Run one MD leg on the GPU and return entropy + RMSF."""
        import os, tempfile
        from biosensor_pipeline.md_entropy import run_leg
        d = tempfile.mkdtemp()
        pp = os.path.join(d, "protein.pdb")
        open(pp, "w").write(protein_text)
        lig = None
        if ligand_sdf_text:
            lig = os.path.join(d, "lig.sdf")
            open(lig, "w").write(ligand_sdf_text)
        return run_leg(pp, ligand_sdf=lig, ligand_smiles=ligand_smiles, solvent=solvent,
                       equil_ps=equil_ps, prod_ns=prod_ns, catalytic_resids=catalytic_resids)

    @app.local_entrypoint()
    def main(apo_pdb: str, holo_pdb: str, ligand_sdf: str = None, smiles: str = None,
             solvent: str = "explicit", equil_ps: float = 200.0, prod_ns: float = 20.0):
        # TEM-1 catalytic residues in the chimera are before the insertion, so their
        # chimera res_id == reporter 0-idx + 1 (see scoring.map_reporter_index).
        from biosensor_pipeline.systems import TEM1
        cat = [v + 1 for v in TEM1.catalytic.values()]
        apo_txt = open(apo_pdb).read()
        holo_txt = open(holo_pdb).read()
        lig_txt = open(ligand_sdf).read() if ligand_sdf else None
        # two GPU legs in parallel
        a = leg_remote.spawn(apo_txt, None, None, solvent, equil_ps, prod_ns, cat)
        h = leg_remote.spawn(holo_txt, lig_txt, smiles, solvent, equil_ps, prod_ns, cat)
        apo, holo = a.get(), h.get()
        dS = round(holo["S_JmolK"] - apo["S_JmolK"], 1)
        print("=== apo/holo MD entropy (Modal GPU) ===")
        print(f"  apo : S={apo['S_JmolK']} J/mol/K  active-site RMSF={apo.get('active_site_rmsf_A')} A")
        print(f"  holo: S={holo['S_JmolK']} J/mol/K  active-site RMSF={holo.get('active_site_rmsf_A')} A")
        print(f"  ΔS(holo-apo) = {dS} J/mol/K  -> entropy reduced on binding: {dS < 0}")
        print("  ✅ measurements; ⚠️ ΔS->DR needs wet-lab calibration. ΔS<0 supports the entropic switch.")
