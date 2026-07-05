import modal

app = modal.App("is621-perses-api2")
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("perses", "openmm", "openmmtools", "openff-toolkit",
                           "openmmforcefields", "pymbar", "ambertools",
                           "cudatoolkit=11.8", channels=["conda-forge"]))


@app.function(image=img, timeout=600)
def introspect():
    import inspect
    from perses.app.relative_point_mutation_setup import PointMutationExecutor as PME
    out = {}
    out["full_sig"] = str(inspect.signature(PME.__init__))
    src = inspect.getsource(PME.__init__)
    # lines that assign self.<attr> = ... (to learn what objects are produced)
    out["self_assigns"] = [ln.strip() for ln in src.splitlines() if "self." in ln and "=" in ln][:60]
    return out


@app.local_entrypoint()
def main():
    r = introspect.remote()
    print("FULL SIGNATURE:\n", r["full_sig"])
    print("\nSELF ASSIGNMENTS (attributes produced):")
    for a in r["self_assigns"]:
        print("  ", a)
