import modal

app = modal.App("is621-fep-probe")
# openmmtools + pymbar alchemical stack (lighter than full perses); test it imports & builds
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("openmm", "openmmtools", "pymbar", "ambertools", "parmed",
                           "mdtraj", "numpy", "cudatoolkit=11.8", channels=["conda-forge"]))


@app.function(image=img, gpu="A10G", timeout=900)
def probe():
    out = {}
    import openmm, openmmtools, pymbar, parmed, mdtraj, numpy as np
    out["openmm"] = openmm.__version__
    out["openmmtools"] = openmmtools.__version__
    out["pymbar"] = pymbar.__version__
    # GPU platform?
    from openmm import Platform
    plats = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    out["platforms"] = plats
    # can we build an alchemical factory + state (core FEP machinery)?
    from openmmtools import alchemy, testsystems
    ts = testsystems.AlanineDipeptideVacuum()
    factory = alchemy.AbsoluteAlchemicalFactory()
    region = alchemy.AlchemicalRegion(alchemical_atoms=list(range(6)))
    alch_sys = factory.create_alchemical_system(ts.system, region)
    state = alchemy.AlchemicalState.from_system(alch_sys)
    state.lambda_electrostatics = 0.5
    out["alchemy_ok"] = True
    # pymbar MBAR smoke
    from pymbar import MBAR
    out["mbar_import"] = True
    return out


@app.local_entrypoint()
def main():
    import json
    print(json.dumps(probe.remote(), indent=2, default=str))
