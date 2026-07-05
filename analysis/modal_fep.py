import modal, json

app = modal.App("is621-fep")
vol = modal.Volume.from_name("is621-fep-vol", create_if_missing=True)
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("perses", "openmm", "openmmtools", "openff-toolkit",
                           "openmmforcefields", "pymbar", "ambertools", "pdbfixer",
                           "mdtraj", "cudatoolkit=11.8", channels=["conda-forge"]))

FF = ["amber14/protein.ff14SB.xml", "amber14/DNA.OL15.xml",
      "amber14/RNA.OL3.xml", "amber14/tip3p.xml"]


@app.function(image=img, gpu="A10G", timeout=14400, volumes={"/vol": vol})
def run_leg(pdb_text: str, resid: str, proposed: str, label: str,
            n_states: int = 12, n_iter: int = 500, n_steps: int = 250):
    import os, traceback, numpy as np
    res = {"label": label, "resid": resid, "proposed": proposed,
           "n_states": n_states, "n_iter": n_iter}
    try:
        os.makedirs("/work", exist_ok=True); os.chdir("/work")
        # add hydrogens (standard resnames) via pdbfixer
        from pdbfixer import PDBFixer
        import openmm.app as app_mm, openmm.unit as unit, openmm as mm
        open("in.pdb", "w").write(pdb_text)
        fx = PDBFixer(filename="in.pdb")
        fx.findMissingResidues(); fx.missingResidues = {}
        fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.0)
        app_mm.PDBFile.writeFile(fx.topology, fx.positions, open("h.pdb", "w"))

        from perses.app.relative_point_mutation_setup import PointMutationExecutor
        pme = PointMutationExecutor(
            protein_filename="h.pdb", mutation_chain_id="A", mutation_residue_id=resid,
            proposed_residue=proposed, forcefield_files=FF, is_solvated=False,
            conduct_endstate_validation=False,
            transform_waters_into_ions_for_charge_changes=True)
        htf = pme.get_apo_htf()

        from perses.samplers.multistate import HybridRepexSampler
        from openmmtools.multistate import MultiStateReporter, MultiStateSamplerAnalyzer
        from openmmtools import mcmc
        from perses.annihilation.lambda_protocol import LambdaProtocol
        move = mcmc.LangevinDynamicsMove(timestep=4.0*unit.femtoseconds,
               collision_rate=1.0/unit.picoseconds, n_steps=n_steps, reassign_velocities=False)
        hss = HybridRepexSampler(mcmc_moves=move, hybrid_factory=htf, online_analysis_interval=None)
        reporter = MultiStateReporter(f"/work/{label}.nc", checkpoint_interval=50)
        hss.setup(n_states=n_states, temperature=300*unit.kelvin, storage_file=reporter,
                  lambda_protocol=LambdaProtocol(functions="default"))
        hss.extend(n_iter)
        ana = MultiStateSamplerAnalyzer(reporter)
        f_ij, df_ij = ana.get_free_energy()
        kT_kcal = (ana.kT).value_in_unit(unit.kilocalorie_per_mole)
        res["dG_kcal"] = float(f_ij[0, -1] * kT_kcal)
        res["dG_err_kcal"] = float(df_ij[0, -1] * kT_kcal)
        res["ok"] = True
    except Exception:
        res["ok"] = False; res["traceback"] = traceback.format_exc()[-3000:]
    with open(f"/vol/{label}.json", "w") as f:
        json.dump(res, f)
    vol.commit()
    return {k: v for k, v in res.items() if k != "traceback"} | (
        {"tb_tail": res.get("traceback", "")[-600:]} if not res.get("ok") else {})


@app.local_entrypoint()
def main():
    import os
    apo = open("mdqm/fep_apo.pdb").read()
    if os.environ.get("FEP_VALIDATE") == "1":
        r = run_leg.remote(apo, "231", "TYR", "F231Y_validate", 12, 5)
        print("VALIDATION RESULT:", r); return
    rna = open("mdqm/fep_complex_rna.pdb").read()
    h1 = run_leg.spawn(apo, "231", "TYR", "F231Y_apo", 12, 800)
    h2 = run_leg.spawn(rna, "231", "TYR", "F231Y_complexRNA", 12, 800)
    print("spawned F231Y_apo:", h1.object_id, " F231Y_complexRNA:", h2.object_id)
