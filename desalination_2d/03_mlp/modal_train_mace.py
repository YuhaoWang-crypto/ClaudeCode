#!/usr/bin/env python3
"""
Train a MACE potential on Modal GPU -- the equivariant-MLP counterpart to
modal_train_deepmd.py, for a DeePMD-vs-MACE(-vs-NequIP) bake-off on identical
data.

Like DeePMD training, MACE training is GPU-bound. The only real differences
from the DeePMD launcher are the pip image (mace-torch) and the train command
(mace_run_train reading an extended-XYZ dataset).

PREREQUISITES
-------------
    pip install modal
    modal token new
    # convert + upload data (extended-XYZ, not npy):
    python prepare_deepmd_data.py --to-extxyz ../data/train/pore_crossing \
           --out ../data/train/pore_crossing.xyz
    modal volume put desal-mlp-data ../data /data
Then:
    modal run modal_train_mace.py --config mace_config.yaml
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
    mace_image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("git")
        .pip_install("mace-torch", "ase", "torch")
    )
    app = modal.App("desal-mace")
    vol = modal.Volume.from_name("desal-mlp-data", create_if_missing=True)

    @app.function(image=mace_image, gpu="A10G",
                  timeout=60 * 60 * 12, volumes={"/data": vol})
    def train(config_text: str, run_name: str = "mace_run"):
        """Run mace_run_train with the given config; return final error table."""
        import subprocess
        import os
        import re

        workdir = f"/data/{run_name}"
        os.makedirs(workdir, exist_ok=True)
        cfg_path = os.path.join(workdir, "mace_config.yaml")
        # point the config at the mounted data volume
        cfg = config_text.replace("../data", "/data")
        with open(cfg_path, "w") as fh:
            fh.write(cfg)

        cmd = ["mace_run_train", "--config", cfg_path]
        print("Running:", " ".join(cmd))
        proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
        vol.commit()

        # scrape the last reported validation force/energy RMSE
        tail = proc.stdout[-4000:]
        rmse_lines = [l for l in tail.splitlines()
                      if re.search(r"RMSE|rmse|Epoch", l)][-8:]
        return {
            "returncode": proc.returncode,
            "metrics_tail": "\n".join(rmse_lines),
            "stderr_tail": proc.stderr[-2000:],
        }

    @app.local_entrypoint()
    def main(config: str = "mace_config.yaml", run_name: str = "mace_run"):
        text = pathlib.Path(config).read_text()
        result = train.remote(text, run_name=run_name)
        print("=" * 60)
        print("returncode:", result["returncode"])
        print("metrics tail:")
        print(result["metrics_tail"] or "(no metrics parsed)")
        if result["returncode"] != 0:
            print("stderr:", result["stderr_tail"])
