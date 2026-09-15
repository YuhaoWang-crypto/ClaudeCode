"""
Modal GPU runner for E1 (build + equilibrate + produce) and the E4 NEMD
viscosity stage.  The local CPU does ~17 ns/day for these boxes; an A10G does
several hundred, so production goes here.  Results land in the Modal volume
`elyte-work` and are pulled back with:

    modal volume get elyte-work / elyte_work/

Usage:
    modal run electrolyte_pipeline/modal_run.py                       # all systems
    modal run electrolyte_pipeline/modal_run.py --labels LiTFSI_r10   # one
    modal run electrolyte_pipeline/modal_run.py --stage nemd          # viscosity NEMD
"""
from __future__ import annotations

import os

import modal

image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "openmm=8.*", "cuda-version=12.4", "openff-toolkit", "openff-interchange",
        "openff-nagl", "openff-nagl-models", "packmol", "numpy", "pandas", "scipy",
        channels=["conda-forge"],
    )
    .add_local_python_source("electrolyte_pipeline")
)
vol = modal.Volume.from_name("elyte-work", create_if_missing=True)
app = modal.App("elyte-md", image=image)
VOL = "/vol"


@app.function(gpu="A10G", timeout=6 * 3600, volumes={VOL: vol})
def run_system(label: str, equil_ns: float, prod_ns: float, scale: float | None, suffix: str) -> dict:
    from electrolyte_pipeline import e0_systems, e1_build, e1_run
    ff = e0_systems.ForceFieldSpec()
    if scale is not None:
        ff.ion_charge_scale = scale
    name = label + suffix
    outdir = os.path.join(VOL, name)
    e1_build.build(label, ff=ff, outdir=outdir)
    rec = e1_run.run(name, workdir=VOL, equil_ns=equil_ns, prod_ns=prod_ns)
    vol.commit()
    return rec


@app.function(gpu="A10G", timeout=6 * 3600, volumes={VOL: vol})
def run_nemd(label: str, amplitudes: list[float], ns: float) -> dict:
    from electrolyte_pipeline import e4_nemd_viscosity
    out = e4_nemd_viscosity.run(label, workdir=VOL, amplitudes=amplitudes, ns=ns)
    vol.commit()
    return out


@app.function(gpu="A10G", timeout=6 * 3600, volumes={VOL: vol})
def run_interface(equil_ns: float, prod_ns: float) -> dict:
    from electrolyte_pipeline import e7_interface
    out = e7_interface.run(workdir=VOL, equil_ns=equil_ns, prod_ns=prod_ns)
    vol.commit()
    return out


@app.local_entrypoint()
def main(labels: str = "", stage: str = "md", equil_ns: float = 2.0, prod_ns: float = 10.0,
         scale: float = -1.0, suffix: str = "", nemd_ns: float = 2.0,
         amplitudes: str = "0.01,0.02"):
    from electrolyte_pipeline.e0_systems import SERIES
    labs = [s for s in labels.split(",") if s] or [c.label for c in SERIES]
    sc = None if scale < 0 else scale
    if stage == "md":
        args = [(l, equil_ns, prod_ns, sc, suffix) for l in labs]
        for rec in run_system.starmap(args):
            print(rec["label"], rec["label_note"], "rho=%.4f" % rec["prod_stationarity"]["density"]["mean"],
                  "wall %.0fs" % rec["wall_s"])
    elif stage == "interface":
        print(run_interface.remote(equil_ns, prod_ns))
    elif stage == "nemd":
        amps = [float(a) for a in amplitudes.split(",")]
        for out in run_nemd.starmap([(l, amps, nemd_ns) for l in labs]):
            print(out)
