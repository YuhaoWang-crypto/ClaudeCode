import modal

app = modal.App("is621-perses-probe")
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("perses", "openmm", "openmmtools", "openff-toolkit",
                           "openmmforcefields", "pymbar", "ambertools",
                           "cudatoolkit=11.8", channels=["conda-forge"]))


@app.function(image=img, gpu="A10G", timeout=1200)
def probe():
    out = {}
    try:
        import perses, openmm
        out["perses"] = getattr(perses, "__version__", "?")
        out["openmm"] = openmm.__version__
        from perses.app.relative_point_mutation_setup import PointMutationExecutor
        out["point_mutation_executor_import"] = True
    except Exception as e:
        out["import_error"] = repr(e)[:400]
    return out


@app.local_entrypoint()
def main():
    import json
    print(json.dumps(probe.remote(), indent=2, default=str))
