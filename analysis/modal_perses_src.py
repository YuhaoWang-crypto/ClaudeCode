import modal
app = modal.App("is621-perses-src")
img = (modal.Image.micromamba(python_version="3.10")
       .micromamba_install("perses", "openmm", "openmmtools", "openff-toolkit",
                           "openmmforcefields", "pymbar", "ambertools",
                           "cudatoolkit=11.8", channels=["conda-forge"]))
@app.function(image=img, timeout=300)
def src():
    p = "/opt/conda/lib/python3.10/site-packages/perses/app/relative_point_mutation_setup.py"
    lines = open(p).read().splitlines()
    return "\n".join(f"{i+1}: {lines[i]}" for i in range(220, 275))
@app.local_entrypoint()
def main():
    print(src.remote())
