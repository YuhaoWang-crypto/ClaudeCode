import modal

app = modal.App("is621-perses-api")
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("perses", "openmm", "openmmtools", "openff-toolkit",
                           "openmmforcefields", "pymbar", "ambertools",
                           "cudatoolkit=11.8", channels=["conda-forge"]))


@app.function(image=img, timeout=600)
def introspect():
    import inspect
    from perses.app.relative_point_mutation_setup import PointMutationExecutor
    out = {}
    out["init_sig"] = str(inspect.signature(PointMutationExecutor.__init__))
    out["doc"] = (PointMutationExecutor.__init__.__doc__ or "")[:2500]
    out["methods"] = [m for m in dir(PointMutationExecutor) if not m.startswith("_")]
    # also look for the sampler/analysis helpers
    try:
        import perses.app.setup_relative_calculation as src
        out["setup_module_fns"] = [f for f in dir(src) if not f.startswith("_")][:40]
    except Exception as e:
        out["setup_err"] = str(e)[:200]
    return out


@app.local_entrypoint()
def main():
    import json
    r = introspect.remote()
    print("INIT SIGNATURE:\n", r["init_sig"])
    print("\nATTRS/METHODS:\n", r["methods"])
    print("\nDOC:\n", r["doc"])
    print("\nSETUP FNS:\n", r.get("setup_module_fns"))
