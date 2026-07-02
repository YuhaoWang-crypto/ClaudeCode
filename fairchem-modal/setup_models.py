"""
Warm the FAIR chemistry model cache on Modal
============================================

These models are NOT pip packages -- they are gated checkpoints pulled from
HuggingFace via fairchem. "Installing" them = confirming your HF token can
load each, and pre-downloading the weights into the persistent Modal volume
`hf-cache` so every later run (from any machine with your Modal creds) is an
instant cache hit.

  modal run setup_models.py::list_models     # what's available in this fairchem
  modal run setup_models.py::warm            # download all requested checkpoints
"""

import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch", "fairchem-core>=2.0.0", "huggingface_hub")
)
app = modal.App("fairchem-setup", image=image)
hf_secret = modal.Secret.from_name("huggingface")
model_cache = modal.Volume.from_name("hf-cache", create_if_missing=True)
CACHE_DIR = "/root/.cache/huggingface"


@app.function(secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=600)
def _list():
    from fairchem.core.calculate import pretrained_mlip as pm
    names = None
    for attr in ("available_models", "MODEL_REGISTRY", "_MODEL_REGISTRY", "models"):
        if hasattr(pm, attr):
            obj = getattr(pm, attr)
            try:
                names = sorted(obj.keys()) if hasattr(obj, "keys") else sorted(obj)
            except Exception:
                names = list(obj)
            break
    if names is None:  # fallback: provoke the KeyError that lists them
        try:
            pm.get_predict_unit("___does_not_exist___")
        except Exception as e:
            names = str(e)
    return names


@app.local_entrypoint()
def list_models():
    print("Available fairchem model tags:")
    print(_list.remote())


@app.function(secrets=[hf_secret], volumes={CACHE_DIR: model_cache}, timeout=3600)
def _warm(tags: list):
    from fairchem.core import pretrained_mlip
    report = []
    for tag in tags:
        try:
            pretrained_mlip.get_predict_unit(tag, device="cpu")  # forces download
            report.append({"model": tag, "status": "cached OK"})
        except Exception as e:
            report.append({"model": tag, "status": f"FAILED: {type(e).__name__}: {str(e)[:160]}"})
    model_cache.commit()
    return report


@app.local_entrypoint()
def warm(tags: str = ""):
    import json
    # default: one representative checkpoint per model family the user has access to
    default = ["uma-s-1p2", "uma-m-1p1", "esen-sm-conserving-all-omol",
               "esen-sm-full-odac25", "esen-sm-conserving-all-oc25"]
    tag_list = [t.strip() for t in tags.split(",") if t.strip()] or default
    print(f"Warming {len(tag_list)} checkpoints into hf-cache ...")
    print(json.dumps(_warm.remote(tag_list), indent=2))
