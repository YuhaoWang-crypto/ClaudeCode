#!/usr/bin/env python3
"""
Run Quantum ESPRESSO DFT on Modal (https://modal.com) -- serverless cloud CPU.

WHY MODAL FOR DFT
-----------------
Plane-wave DFT (QE, VASP) is dominated by dense linear algebra + MPI and is
overwhelmingly *CPU*-bound, not GPU-bound. So for the STAGE-1 DFT / AIMD
labelling you want many CPU cores, NOT a GPU. Modal is a good fit because:
  * you pay per-second for a fat multi-core box only while the SCF runs;
  * the whole environment (QE + pseudopotentials) is defined in code, so runs
    are reproducible and you never touch a cluster module system.
Reserve GPUs (see 03_mlp/modal_train_deepmd.py) for MLP *training* and for
LAMMPS+DeePMD *inference*, which are the GPU-bound parts of this project.

PREREQUISITES
-------------
    pip install modal
    modal token new          # one-time browser auth
Then:
    modal run modal_run_dft.py --infile qe/scf_graphene.in

This script is a runnable TEMPLATE. It will execute on Modal once you have an
account; it does not run inside a restricted sandbox with no Modal token.
"""
import pathlib
import sys

try:
    import modal
except ImportError:  # let the file be importable/inspectable without modal
    modal = None
    print("[note] `pip install modal` to actually launch cloud runs.",
          file=sys.stderr)

# ---- Container image: Quantum ESPRESSO from conda-forge + MPI ----------------
if modal is not None:
    qe_image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("wget", "libgomp1")
        .run_commands(
            # micromamba install of a prebuilt, MPI-enabled QE
            "wget -qO /tmp/mm.sh https://micro.mamba.pm/install.sh || true",
            "conda config --add channels conda-forge || true",
        )
        .micromamba_install(
            "qe=7.2", "openmpi", "numpy", "ase",
            channels=["conda-forge"],
        )
    )

    app = modal.App("desal-dft")

    # A Volume persists pseudopotentials + outputs across runs.
    vol = modal.Volume.from_name("desal-dft-data", create_if_missing=True)

    @app.function(image=qe_image, cpu=16.0, memory=32768,
                  timeout=60 * 60 * 4, volumes={"/data": vol})
    def run_scf(infile_text: str, prefix: str = "job"):
        """Run pw.x on `infile_text` using all allocated cores. Returns stdout."""
        import subprocess
        import os

        workdir = f"/data/{prefix}"
        os.makedirs(workdir, exist_ok=True)
        in_path = os.path.join(workdir, "scf.in")
        with open(in_path, "w") as fh:
            fh.write(infile_text)

        ncore = int(os.environ.get("MODAL_CPU", "16"))
        cmd = ["mpirun", "-np", str(ncore), "pw.x", "-in", in_path]
        print("Running:", " ".join(cmd))
        proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
        out_path = os.path.join(workdir, "scf.out")
        with open(out_path, "w") as fh:
            fh.write(proc.stdout)
        vol.commit()
        # pull the final energy line for a quick sanity check
        energy = [l for l in proc.stdout.splitlines() if "!" in l and "total energy" in l]
        return {
            "returncode": proc.returncode,
            "energy_line": energy[-1] if energy else "(SCF did not converge?)",
            "stderr_tail": proc.stderr[-2000:],
        }

    @app.local_entrypoint()
    def main(infile: str = "qe/scf_graphene.in", prefix: str = "graphene"):
        text = pathlib.Path(infile).read_text()
        result = run_scf.remote(text, prefix=prefix)
        print("=" * 60)
        print("returncode:", result["returncode"])
        print("energy    :", result["energy_line"])
        if result["returncode"] != 0:
            print("stderr:", result["stderr_tail"])
