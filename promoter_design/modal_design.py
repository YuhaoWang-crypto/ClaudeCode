"""
Deployable Modal app: run the stimulus-responsive promoter design pipeline on
Modal using the proto-tools image (Evo2 generation + AlphaGenome scoring).

    modal deploy modal_design.py                 # register app (no compute cost)
    modal run modal_design.py::assemble          # CPU: build FASTA cassettes
    modal run modal_design.py::design \           # GPU: real Evo2+AlphaGenome run
        --stimulus interferon_typeII --target THP1
    modal run modal_design.py::design \           # stimulus-AND-cell-type gate
        --stimulus hypoxia --target HepG2 --lineage hepatocyte

Notes
- The image is the SAME one built/verified by modal_setup.py (proto-language +
  git proto-tools). Building it here reuses Modal's cached layers.
- `design` needs a GPU + the 'huggingface' secret (gated Evo2 / AlphaGenome).
  AlphaGenome also requires you to have accepted its license on your HF account.
- `assemble` is CPU-only and free-ish; use it to sanity-check the mount.
"""
import modal

app = modal.App("proto-promoter-design")

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential")
    .pip_install("proto-language")
    .run_commands("pip install git+https://github.com/evo-design/proto-tools.git")
    # Mount our design code so it is importable inside the container.
    .add_local_dir(".", remote_path="/root/pd", ignore=["designs/*", "__pycache__/*"])
)

GPU = "A100"          # Evo2 needs a large-VRAM GPU; adjust to your model size.
hf = modal.Secret.from_name("huggingface")


@app.function(image=image, timeout=1200)
def assemble(copies: int = 3, min_promoter: str = "minCMV"):
    """CPU: build all cassettes (single-stimulus + composite/AND) and return FASTA text."""
    import sys, subprocess, os
    sys.path.insert(0, "/root/pd")
    os.makedirs("/root/pd/designs", exist_ok=True)
    subprocess.run([sys.executable, "build_constructs.py", "--copies", str(copies),
                    "--min_promoter", min_promoter], cwd="/root/pd", check=True)
    subprocess.run([sys.executable, "dual_and_designs.py"], cwd="/root/pd", check=True)
    import os, glob
    out = {}
    for f in glob.glob("/root/pd/designs/*.fasta"):
        with open(f) as fh:
            out[os.path.basename(f)] = fh.read()
    return out


@app.function(image=image, gpu=GPU, secrets=[hf], timeout=3600)
def design(stimulus: str = "interferon_typeII", target: str = "THP1",
           lineage: str | None = None, num_steps: int = 12, num_results: int = 4):
    """GPU: run the proto_language optimiser (Evo2 generator + AlphaGenome
    contrastive cell-type scoring) for one (stimulus, target cell [, lineage])."""
    import sys
    sys.path.insert(0, "/root/pd")
    import design_pipeline as dp
    if not dp._READY:
        return {"error": f"proto stack not importable on GPU box: {dp._IMPORT_ERR}"}
    opt = dp.design(stimulus, target, lineage=lineage,
                    num_steps=num_steps, num_results=num_results, use_evo2=True)
    seqs = [getattr(s, "sequence", "") for s in getattr(opt, "segments", [])]
    return {"stimulus": stimulus, "target": target, "lineage": lineage,
            "designed_spacers_or_segments": seqs}


@app.function(image=image, gpu=GPU, secrets=[hf], timeout=1800)
def smoke():
    """GPU: minimal proof that Evo2 + AlphaGenome are callable in this image."""
    import subprocess
    print(subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
                          "--format=csv,noheader"], capture_output=True, text=True).stdout)
    import proto_tools
    from proto_tools.tools.causal_models import evo2 as _evo2          # noqa
    from proto_tools.tools.sequence_scoring import alphagenome as _ag  # noqa
    return {"proto_tools": proto_tools.__version__,
            "evo2_module": _evo2.__name__, "alphagenome_module": _ag.__name__}


@app.local_entrypoint()
def main(stimulus: str = "interferon_typeII", target: str = "THP1",
         lineage: str = ""):
    r = design.remote(stimulus, target, lineage or None)
    print(r)
