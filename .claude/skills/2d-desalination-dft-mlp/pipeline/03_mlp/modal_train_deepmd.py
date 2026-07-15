#!/usr/bin/env python3
"""
Train a DeePMD-kit machine-learning potential on Modal GPU.

WHY GPU HERE (and not for the DFT step)
---------------------------------------
Neural-network potential *training* is GPU-bound: it is exactly the kind of
dense tensor workload a single A10G/A100 accelerates 10-50x over CPU. The same
goes for LAMMPS+DeePMD *inference* during the ns-scale NEMD runs. This is the
complement to 02_dft/modal_run_dft.py (which uses CPU for the plane-wave DFT).

So the natural Modal split for this project is:
    DFT / AIMD labelling      -> CPU boxes   (modal_run_dft.py)
    MLP training              -> 1 GPU       (this file)
    NEMD production (LAMMPS)  -> 1 GPU       (analogous, swap the command)

PREREQUISITES
-------------
    pip install modal
    modal token new
    # upload your DFT-labelled training data to the Modal Volume first, e.g.
    modal volume put desal-mlp-data ./data /data
Then:
    modal run modal_train_deepmd.py --input deepmd_input.json

This is a runnable template; it launches once you have a Modal account + data.
"""
import pathlib
import sys

try:
    import modal
except ImportError:
    modal = None
    print("[note] `pip install modal` to actually launch cloud runs.",
          file=sys.stderr)

if modal is not None:
    # DeePMD-kit ships an official CUDA image on PyPI (deepmd-kit[gpu]).
    dp_image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("git")
        .pip_install(
            "deepmd-kit[gpu]",   # pulls the TensorFlow/PyTorch CUDA backend
            "dpdata",            # DFT -> DeePMD data conversion
            "tensorflow",
        )
    )

    app = modal.App("desal-mlp")
    vol = modal.Volume.from_name("desal-mlp-data", create_if_missing=True)

    @app.function(image=dp_image, gpu="A10G",
                  timeout=60 * 60 * 12, volumes={"/data": vol})
    def train(input_json_text: str, run_name: str = "mlp_run"):
        """Run `dp train` then freeze + compress the model. Returns metrics."""
        import subprocess
        import os
        import json

        workdir = f"/data/{run_name}"
        os.makedirs(workdir, exist_ok=True)
        # rewrite data paths to the mounted volume
        cfg = json.loads(input_json_text)
        with open(os.path.join(workdir, "input.json"), "w") as fh:
            json.dump(cfg, fh, indent=2)

        def _run(cmd):
            print("Running:", " ".join(cmd))
            return subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)

        train_proc = _run(["dp", "train", "input.json"])
        _run(["dp", "freeze", "-o", "graph.pb"])
        _run(["dp", "compress", "-i", "graph.pb", "-o", "graph-compress.pb"])
        vol.commit()

        # learning-curve tail = final train/valid RMSE on energy & force
        lcurve = os.path.join(workdir, "lcurve.out")
        tail = ""
        if os.path.exists(lcurve):
            with open(lcurve) as fh:
                tail = "".join(fh.readlines()[-5:])
        return {
            "returncode": train_proc.returncode,
            "lcurve_tail": tail,
            "stderr_tail": train_proc.stderr[-2000:],
        }

    @app.local_entrypoint()
    def main(input: str = "deepmd_input.json", run_name: str = "mlp_run"):
        text = pathlib.Path(input).read_text()
        result = train.remote(text, run_name=run_name)
        print("=" * 60)
        print("returncode:", result["returncode"])
        print("final learning curve (step  rmse_e  rmse_f ...):")
        print(result["lcurve_tail"] or "(no lcurve produced)")
        if result["returncode"] != 0:
            print("stderr:", result["stderr_tail"])
