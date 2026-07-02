"""
Install evo-design/proto-tools (+ proto-language) on Modal and verify it.

The image build runs on Modal's cloud (which can reach GitHub), so the git
install succeeds there even when the local sandbox proxy blocks it.

    modal run modal_setup.py                 # CPU: build image + import check
    modal run modal_setup.py::gpu_check      # GPU: import check on a GPU box

proto-tools pulls gated models (AlphaGenome, ESM3, ...) that need an accepted
license + HF token at RUN time; import/listing does not. Provide HF_TOKEN via a
Modal secret named 'huggingface' when you start running real inference.
"""
import modal

app = modal.App("proto-tools-setup")

# Build proto-tools from source. git needed for the VCS install.
proto_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential")
    .pip_install("proto-language")  # language layer (PyPI)
    .run_commands(
        "pip install --no-build-isolation "
        "git+https://github.com/evo-design/proto-tools.git || "
        "pip install git+https://github.com/evo-design/proto-tools.git"
    )
)


def _report():
    import importlib, pkgutil, json
    out = {}
    try:
        import proto_tools
        out["proto_tools_version"] = getattr(proto_tools, "__version__", "?")
        out["proto_tools_path"] = list(getattr(proto_tools, "__path__", []))
    except Exception as e:
        out["import_error"] = f"{type(e).__name__}: {e}"
        return out
    # Enumerate the scoring/generation tools we care about.
    wanted = ["alphagenome", "evo2", "enformer", "borzoi", "esm3", "esm2"]
    found = {}
    try:
        import proto_tools.tools as T
        for _, name, _ in pkgutil.walk_packages(T.__path__, "proto_tools.tools."):
            for w in wanted:
                if name.endswith(w) or name.endswith(w + "." + w):
                    found.setdefault(w, []).append(name)
    except Exception as e:
        out["walk_error"] = f"{type(e).__name__}: {e}"
    out["found_tools"] = found
    # Does proto_language import cleanly now that proto_tools is present?
    try:
        import proto_language
        out["proto_language_ok"] = True
        out["proto_language_version"] = getattr(proto_language, "__version__", "?")
    except Exception as e:
        out["proto_language_ok"] = False
        out["proto_language_error"] = f"{type(e).__name__}: {e}"
    return out


@app.function(image=proto_image, timeout=2400)
def check():
    import json
    r = _report()
    print("PROTO_TOOLS_REPORT " + json.dumps(r, indent=2))
    return r


@app.function(image=proto_image, gpu="T4", timeout=2400)
def gpu_check():
    import json, subprocess
    print(subprocess.run(["nvidia-smi"], capture_output=True, text=True).stdout[:600])
    r = _report()
    print("PROTO_TOOLS_REPORT " + json.dumps(r, indent=2))
    return r


@app.local_entrypoint()
def main():
    r = check.remote()
    print("\n=== summary ===")
    print("import_error:", r.get("import_error", "none"))
    print("proto_language_ok:", r.get("proto_language_ok"))
    print("found_tools keys:", list((r.get("found_tools") or {}).keys()))
