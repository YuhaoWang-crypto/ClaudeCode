"""Minimal end-to-end proof: build an image + run on a real GPU on Modal."""
import modal

app = modal.App("claude-gpu-probe")

# "install" step: a container image with a scientific package baked in
image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("numpy==2.4.6")
)

@app.function(gpu="T4", image=image, timeout=120)
def probe():
    import subprocess, numpy as np
    smi = subprocess.run(["nvidia-smi",
                          "--query-gpu=name,memory.total,driver_version",
                          "--format=csv,noheader"],
                         capture_output=True, text=True).stdout.strip()
    # a trivial numeric op to confirm the userland stack runs
    x = float(np.linalg.eigvalsh(np.array([[2.0, 1.0], [1.0, 2.0]]))[0])
    return {"gpu": smi, "numpy": np.__version__, "eig_min": x}

@app.local_entrypoint()
def main():
    print("RESULT:", probe.remote())
