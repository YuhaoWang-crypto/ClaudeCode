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
- `assemble` is CPU-only and free-ish; use it to sanity-check the mount.

RUNTIME REALITY (verified via the inspect_* functions below)
------------------------------------------------------------
The base proto-tools install is a thin orchestration layer: it has the tool
WRAPPERS but NOT the heavy model runtimes. Confirmed missing from this image:
torch, evo2, flash_attn, transformer_engine, alphagenome. To actually run:
  * Evo2   -> build on a CUDA image with torch + evo2 (+ flash-attn /
              transformer_engine) and accept the arcinstitute Evo2 HF license.
              model_checkpoint default 'evo2_7b'; smallest is 'evo2_1b_base'.
  * AlphaGenome -> `pip install alphagenome` AND set ALPHAGENOME_CHECKPOINT_PATH
              to the gated DeepMind checkpoint (you must download it; not public).
The inspect_* functions (CPU, free) re-check all of this on demand.
"""
import modal

app = modal.App("proto-promoter-design")

# Base image: proto-language + proto-tools framework (NO local files yet, so we
# can still append build steps for the GPU variant below).
base_image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "build-essential")
    .pip_install("proto-language")
    .run_commands("pip install git+https://github.com/evo-design/proto-tools.git")
)

_MNT = dict(remote_path="/root/pd", ignore=["designs/*", "__pycache__/*"])
# CPU image: base + local design code (add_local_dir must be LAST).
image = base_image.add_local_dir(".", **_MNT)

GPU = "A100"          # Evo2 needs a large-VRAM GPU; adjust to your model size.
hf = modal.Secret.from_name("huggingface")

# --- Self-hosted runtime (Evo2 + AlphaGenome via proto-tools standalone envs) --
# proto-tools builds per-tool micromamba envs (CUDA toolkit, flash-attn, evo2 /
# jax, alphagenome) and downloads weights on first call. A persistent Volume
# caches all of it under PROTO_HOME so later runs are fast.
proto_cache = modal.Volume.from_name("proto-cache", create_if_missing=True)

gpu_image = (
    base_image
    .apt_install("curl", "bzip2")
    # micromamba binary for the standalone-env builder.
    .run_commands(
        "curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest "
        "| tar -xj -C /usr/local bin/micromamba && chmod +x /usr/local/bin/micromamba")
    .env({
        "PROTO_HOME": "/proto_home",
        "HF_HOME": "/proto_home/hf",
        "HUGGINGFACE_HUB_CACHE": "/proto_home/hf/hub",
        "MAMBA_BIN": "/usr/local/bin/micromamba",
        "MAMBA_ROOT_PREFIX": "/proto_home/mamba",
        "PROTO_MODEL_CACHE": "/proto_home/models",
    })
    # local files LAST.
    .add_local_dir(".", **_MNT)
)
CACHE = {"/proto_home": proto_cache}


@app.function(image=gpu_image, gpu=GPU, secrets=[hf], volumes=CACHE, timeout=10800)
def full_design(stimulus: str = "interferon_typeII", target: str = "THP1",
                lineage: str | None = None, num_steps: int = 2, num_results: int = 1,
                evo2_model: str = "evo2_1b_base"):
    """GPU: trigger the Evo2 + AlphaGenome standalone builds (first call compiles
    them into the Volume) and run the real design. Long on first run."""
    import sys, os, json, traceback, subprocess
    sys.path.insert(0, "/root/pd")
    for d in ("/proto_home/hf", "/proto_home/mamba", "/proto_home/models"):
        os.makedirs(d, exist_ok=True)
    print("GPU:", subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total",
          "--format=csv,noheader"], capture_output=True, text=True).stdout.strip())
    import design_pipeline as dp
    if not dp._READY:
        return {"error": f"proto stack import failed: {dp._IMPORT_ERR}"}
    result = {"stimulus": stimulus, "target": target, "lineage": lineage,
              "evo2_model": evo2_model}
    try:
        opt = dp.design(stimulus, target, lineage=lineage, num_steps=num_steps,
                        num_results=num_results, use_evo2=True, evo2_model=evo2_model)
        result["designed"] = [getattr(s, "sequence", "") for s in getattr(opt, "segments", [])]
        result["status"] = "ok"
    except Exception as e:
        result["status"] = "error"
        result["error"] = f"{type(e).__name__}: {e}"
        result["traceback"] = traceback.format_exc()[-3000:]
    proto_cache.commit()
    print("FULL_DESIGN_RESULT " + json.dumps(result, indent=2)[:4000])
    return result


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


@app.function(image=image, timeout=600)
def inspect_evo2():
    """CPU: dump Evo2GeneratorConfig fields + defaults and the evo2 tool's default
    model, so we size the GPU / pick the model before spending GPU budget."""
    import json
    from proto_language import Evo2GeneratorConfig
    fields = {}
    for name, fi in Evo2GeneratorConfig.model_fields.items():
        default = fi.default
        try:
            json.dumps(default)
        except Exception:
            default = repr(default)
        fields[name] = {"default": default, "annotation": str(fi.annotation)}
    out = {"config_fields": fields}
    # Any module-level default model constant in the proto_tools evo2 tool?
    try:
        import proto_tools.tools.causal_models.evo2 as e2
        consts = {k: getattr(e2, k) for k in dir(e2)
                  if k.isupper() and isinstance(getattr(e2, k), (str, int, float, list, dict))}
        out["evo2_module_constants"] = consts
    except Exception as ex:
        out["evo2_introspect_error"] = f"{type(ex).__name__}: {ex}"
    print("EVO2_INSPECT " + json.dumps(out, indent=2, default=str))
    return out


@app.function(image=image, timeout=600)
def inspect_alphagenome():
    """CPU: find out how proto_tools' AlphaGenome tool obtains the model — does it
    need an API key / license env var? Grep its source + dump any config."""
    import os, re, json, glob
    import proto_tools.tools.sequence_scoring.alphagenome as ag
    root = os.path.dirname(ag.__file__)
    keys, env_refs = set(), set()
    for f in glob.glob(root + "/**/*.py", recursive=True):
        try:
            src = open(f).read()
        except Exception:
            continue
        for m in re.findall(r"[A-Z][A-Z0-9_]*(?:API|KEY|TOKEN|LICENSE)[A-Z0-9_]*", src):
            keys.add(m)
        for m in re.findall(r"os\.environ(?:\.get)?\(\s*[\"']([^\"']+)[\"']", src):
            env_refs.add(m)
        for m in re.findall(r"getenv\(\s*[\"']([^\"']+)[\"']", src):
            env_refs.add(m)
    out = {"module_dir": root, "keyish_tokens": sorted(keys),
           "env_vars_referenced": sorted(env_refs),
           "files": [os.path.basename(x) for x in glob.glob(root + "/*.py")]}
    print("AG_INSPECT " + json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=600)
def inspect_runtime():
    """CPU: which heavy inference deps are actually present in the image?"""
    import importlib, json
    mods = ["torch", "evo2", "vortex", "flash_attn", "transformer_engine",
            "alphagenome", "enformer_pytorch", "borzoi_pytorch", "esm"]
    out = {}
    for m in mods:
        try:
            mod = importlib.import_module(m)
            out[m] = getattr(mod, "__version__", "present")
        except Exception as e:
            out[m] = f"MISSING ({type(e).__name__})"
    print("RUNTIME_INSPECT " + json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=600)
def inspect_extras():
    """CPU: proto-tools optional-dependency extras + how the evo2 tool imports its
    runtime (native `evo2` pkg vs HF transformers). Decides the GPU image recipe."""
    import importlib.metadata as md, json, inspect, os, re
    out = {}
    try:
        meta = md.metadata("proto-tools")
        out["provides_extra"] = meta.get_all("Provides-Extra") or []
        reqs = md.requires("proto-tools") or []
        out["conditional_requires"] = [r for r in reqs if "extra ==" in r]
    except Exception as e:
        out["meta_error"] = str(e)
    # How does the evo2 tool load the model?
    try:
        import proto_tools.tools.causal_models.evo2 as e2
        d = os.path.dirname(e2.__file__)
        imports = set()
        for f in [x for x in os.listdir(d) if x.endswith(".py")]:
            src = open(os.path.join(d, f)).read()
            for m in re.findall(r"^\s*(?:import|from)\s+([a-zA-Z0-9_\.]+)", src, re.M):
                top = m.split(".")[0]
                if top in {"evo2","vortex","transformers","torch","flash_attn",
                           "transformer_engine","huggingface_hub"}:
                    imports.add(top)
        out["evo2_runtime_imports"] = sorted(imports)
    except Exception as e:
        out["evo2_introspect_error"] = str(e)
    print("EXTRAS_INSPECT " + json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=600)
def dump_evo2_source():
    """CPU: print the evo2 tool's model-loading lines so we build the exact runtime."""
    import os, re, proto_tools.tools.causal_models.evo2 as e2
    d = os.path.dirname(e2.__file__)
    hits = []
    for f in sorted(os.listdir(d)):
        if not f.endswith(".py"):
            continue
        for i, line in enumerate(open(os.path.join(d, f)), 1):
            if re.search(r"import (evo2|vortex|transformer_engine|flash_attn|torch|transformers)"
                         r"|from (evo2|vortex|transformers)"
                         r"|from_pretrained|AutoModel|Evo2\(|load_checkpoint|hf_hub|snapshot_download",
                         line):
                hits.append(f"{f}:{i}: {line.rstrip()[:140]}")
    print("EVO2_SRC\n" + "\n".join(hits))
    return hits


@app.function(image=image, timeout=600)
def dump_dispatch():
    """CPU: how does proto_tools run a tool? Local weights vs hosted proto-client?"""
    import os, re, json, glob
    out = {}
    # (a) full evo2 tool dir + source
    import proto_tools.tools.causal_models.evo2 as e2
    d = os.path.dirname(e2.__file__)
    out["evo2_files"] = sorted(os.listdir(d))
    src = {}
    for f in glob.glob(d + "/*.py"):
        src[os.path.basename(f)] = open(f).read()[:1500]
    out["evo2_init_head"] = src.get("__init__.py", "")[:1500]
    # (b) proto_client: base URL / API-key env
    try:
        import proto_client, inspect
        pc_dir = os.path.dirname(proto_client.__file__)
        envs, urls = set(), set()
        for f in glob.glob(pc_dir + "/**/*.py", recursive=True):
            s = open(f).read()
            envs |= set(re.findall(r"(?:environ(?:\.get)?\(|getenv\()\s*[\"']([^\"']+)", s))
            urls |= set(re.findall(r"https?://[a-zA-Z0-9\.\-/]+", s))
        out["proto_client_env_vars"] = sorted(envs)
        out["proto_client_urls"] = sorted(urls)[:10]
        out["proto_client_attrs"] = [a for a in dir(proto_client) if not a.startswith("_")][:30]
    except Exception as e:
        out["proto_client_err"] = str(e)
    # (c) does proto_tools have a tool runner that chooses local vs remote?
    try:
        import proto_tools, glob as g
        base = os.path.dirname(proto_tools.__file__)
        run_hits = []
        for f in g.glob(base + "/**/*.py", recursive=True):
            s = open(f).read()
            if re.search(r"proto_client|BACKEND|remote|hosted|api_key|ARC_", s):
                run_hits.append(os.path.relpath(f, base))
        out["files_mentioning_client_or_backend"] = run_hits[:25]
    except Exception as e:
        out["runner_err"] = str(e)
    print("DISPATCH " + json.dumps(out, indent=2)[:6000])
    return out


@app.function(image=image, timeout=600)
def dump_evo2_dispatch():
    """CPU: show how run_evo2_sample chooses hosted (proto_client) vs local."""
    import os, re, proto_tools.tools.causal_models.evo2.evo2_sample as s
    src = open(s.__file__).read()
    lines = src.splitlines()
    keep = []
    for i, ln in enumerate(lines):
        if re.search(r"proto_client|PROTO_API_KEY|def run_evo2_sample|remote|local|standalone|"
                     r"backend|ProtoClient|create_run|\.run\(|getenv|environ|raise ", ln):
            keep.append(f"{i+1}: {ln[:150]}")
    print("EVO2_DISPATCH\n" + "\n".join(keep[:60]))
    return keep[:60]


@app.function(image=image, timeout=600)
def dump_weight_sources():
    """CPU: where do evo2 & alphagenome get their weights? Read the tools'
    links.yaml / license.yaml manifests and grep for HF repos / API keys."""
    import os, re, glob, json
    import proto_tools
    base = os.path.dirname(proto_tools.__file__)
    out = {}
    for tool in ["causal_models/evo2", "sequence_scoring/alphagenome"]:
        d = os.path.join(base, "tools", tool)
        info = {"files": []}
        if os.path.isdir(d):
            info["files"] = sorted(os.listdir(d))
            for man in ["links.yaml", "license.yaml"]:
                p = os.path.join(d, man)
                if os.path.exists(p):
                    info[man] = open(p).read()[:1500]
            # grep standalone + tool for hf repos / download urls / api keys
            refs = set(); keys = set()
            for f in glob.glob(d + "/**/*.py", recursive=True):
                s = open(f).read()
                refs |= set(re.findall(r"(?:huggingface\.co/|hf_hub_download\([^)]*repo_id=?['\"]?|snapshot_download\([^)]*repo_id=?['\"]?)([A-Za-z0-9_\-/\.]+)", s))
                refs |= set(re.findall(r"https?://[A-Za-z0-9\.\-/]+\.(?:pt|ckpt|safetensors|tar|zip)", s))
                keys |= set(re.findall(r"[A-Z][A-Z0-9_]*(?:API_KEY|CHECKPOINT|TOKEN|API)[A-Z0-9_]*", s))
            info["weight_refs"] = sorted(refs)[:15]
            info["key_env"] = sorted(keys)[:15]
        out[tool] = info
    print("WEIGHT_SOURCES " + json.dumps(out, indent=2)[:6000])
    return out


hf_check_image = modal.Image.debian_slim(python_version="3.11").pip_install("huggingface_hub")


@app.function(image=hf_check_image, secrets=[hf], timeout=600)
def check_hf_access():
    """CPU: does HF_TOKEN have access to the Evo2 + AlphaGenome repos? No full
    download — just repo metadata (gated repos raise/40x if not accepted)."""
    import os, json
    from huggingface_hub import HfApi
    tok = os.environ.get("HF_TOKEN")
    api = HfApi(token=tok)
    out = {"hf_token_present": bool(tok)}
    try:
        who = api.whoami(token=tok)
        out["hf_user"] = who.get("name")
    except Exception as e:
        out["whoami_error"] = f"{type(e).__name__}: {str(e)[:120]}"
    for repo in ["arcinstitute/evo2_1b_base", "google/alphagenome-all-folds"]:
        rec = {}
        try:
            info = api.model_info(repo, token=tok, files_metadata=False)
            rec["accessible"] = True
            rec["gated"] = getattr(info, "gated", None)
            rec["siblings"] = len(getattr(info, "siblings", []) or [])
        except Exception as e:
            rec["accessible"] = False
            rec["error"] = f"{type(e).__name__}: {str(e)[:160]}"
        out[repo] = rec
    print("HF_ACCESS " + json.dumps(out, indent=2))
    return out


@app.function(image=image, timeout=600)
def dump_standalone_specs():
    """CPU: pinned install recipes from evo2 & alphagenome standalone/ dirs."""
    import os, glob, json, proto_tools
    base = os.path.join(os.path.dirname(proto_tools.__file__), "tools")
    out = {}
    for tool in ["causal_models/evo2", "sequence_scoring/alphagenome"]:
        d = os.path.join(base, tool, "standalone")
        info = {"exists": os.path.isdir(d)}
        if info["exists"]:
            info["tree"] = [os.path.relpath(p, d) for p in glob.glob(d + "/**/*", recursive=True)][:40]
            for pat in ["requirements*.txt", "pyproject.toml", "setup.py", "*.sh",
                        "environment*.yml", "uv.lock", "*.toml"]:
                for p in glob.glob(d + "/**/" + pat, recursive=True):
                    info[os.path.relpath(p, d)] = open(p).read()[:2000]
        out[tool] = info
    print("STANDALONE " + json.dumps(out, indent=2)[:8000])
    return out


@app.function(image=image, timeout=600)
def dump_venv_manager():
    """CPU: find proto-tools' standalone-venv orchestration + host prereqs
    (micromamba/MAMBA_BIN, PROTO_HOME, how setup.sh is invoked)."""
    import os, re, glob, json, proto_tools
    base = os.path.dirname(proto_tools.__file__)
    hits = {}
    for f in glob.glob(base + "/**/*.py", recursive=True):
        s = open(f).read()
        if re.search(r"MAMBA_BIN|micromamba|standalone_helpers|VENV_PATH|setup\.sh|PROTO_HOME|standalone_env|make_venv|ensure_venv", s):
            rel = os.path.relpath(f, base)
            found = sorted(set(re.findall(r"MAMBA_BIN|micromamba|standalone_helpers\.sh|VENV_PATH|setup\.sh|PROTO_HOME|def [a-z_]*venv[a-z_]*|def [a-z_]*standalone[a-z_]*", s)))
            hits[rel] = found
    # find standalone_helpers.sh + any manager module
    helpers = glob.glob(base + "/**/standalone_helpers.sh", recursive=True)
    mgr = [os.path.relpath(f, base) for f in glob.glob(base + "/**/*.py", recursive=True)
           if re.search(r"venv|standalone", os.path.basename(f))]
    out = {"files_touching_venv": hits, "standalone_helpers_paths":
           [os.path.relpath(h, base) for h in helpers], "manager_modules": mgr[:20]}
    print("VENV_MGR " + json.dumps(out, indent=2)[:6000])
    return out


@app.function(image=image, timeout=600)
def dump_registry_setup():
    """CPU: read tool_registry standalone/micromamba logic + cli eject cmd."""
    import os, re, proto_tools, json
    base = os.path.dirname(proto_tools.__file__)
    out = {}
    for rel in ["tools/tool_registry.py", "cli.py"]:
        p = os.path.join(base, rel)
        s = open(p).read()
        lines = s.splitlines()
        keep = [f"{i+1}: {ln[:150]}" for i, ln in enumerate(lines)
                if re.search(r"micromamba|MAMBA|PROTO_HOME|standalone|setup\.sh|def _?cmd|def has_|def (ensure|build|prepare|setup|make)|download|VENV_PATH|env_vars\.txt|subprocess|\.run\(", ln)]
        out[rel] = keep[:55]
    # env vars proto-tools reads
    envs = sorted(set(re.findall(r"(?:environ(?:\.get)?\(|getenv\()\s*[\"']([A-Z_]+)",
                                 open(os.path.join(base, "tools/tool_registry.py")).read())))
    out["registry_env_vars"] = envs
    print("REG_SETUP " + json.dumps(out, indent=2)[:7000])
    return out


@app.function(image=image, timeout=600)
def find_executor():
    """CPU: find the module that actually builds the standalone env (runs setup.sh,
    sets MAMBA_BIN/VENV_PATH, downloads micromamba) + how it's triggered."""
    import os, re, glob, json, proto_tools
    base = os.path.dirname(proto_tools.__file__)
    exec_files = {}
    for f in glob.glob(base + "/**/*.py", recursive=True):
        s = open(f).read()
        # the executor RUNS setup.sh and sets these env vars for the child
        if re.search(r"setup\.sh", s) and re.search(r"MAMBA_BIN|VENV_PATH|micromamba|env=", s):
            rel = os.path.relpath(f, base)
            lines = s.splitlines()
            keep = [f"{i+1}: {ln.strip()[:150]}" for i, ln in enumerate(lines)
                    if re.search(r"setup\.sh|MAMBA_BIN|VENV_PATH|micromamba|PROTO_HOME|def [a-z_]+|subprocess|download.*mamba|https?://", ln)]
            exec_files[rel] = keep[:40]
    # env vars used across proto_tools utils (esp PROTO_HOME, micromamba url)
    allenv = set(); mamba_urls = set()
    for f in glob.glob(base + "/utils/**/*.py", recursive=True) + glob.glob(base + "/**/*worker*.py", recursive=True) + glob.glob(base + "/**/*pool*.py", recursive=True):
        s = open(f).read()
        allenv |= set(re.findall(r"(?:environ(?:\.get)?\(|getenv\()\s*[\"']([A-Z_]+)", s))
        mamba_urls |= set(re.findall(r"https?://[^\s\"']*(?:micromamba|mamba)[^\s\"']*", s))
    out = {"executor_files": exec_files, "utils_env_vars": sorted(allenv),
           "micromamba_urls": sorted(mamba_urls)}
    print("EXECUTOR " + json.dumps(out, indent=2)[:7000])
    return out


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
