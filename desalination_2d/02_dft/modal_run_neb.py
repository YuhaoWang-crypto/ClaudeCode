#!/usr/bin/env python3
"""
Run a Quantum ESPRESSO CI-NEB (Milestone 2) on Modal cloud CPU, and return the
energy barrier E_a directly.

Like plane-wave SCF, neb.x is CPU/MPI-bound -- it drives many pw.x evaluations
along the path -- so this uses a fat CPU box, not a GPU. It shares the exact
same QE image as modal_run_dft.py.

PREREQUISITES
-------------
    pip install modal
    modal token new
    # put your pseudopotentials on the Volume once:
    modal volume put desal-dft-data ./pseudo /pseudo
Then generate the input and launch:
    python build_neb_endpoints.py --membrane ../figures/graphene_pore.xyz --ion Na
    modal run modal_run_neb.py --infile qe/neb_na_ready.in --ion Na

The barrier comes back parsed: forward E_a, reverse barrier, reaction energy.
This is your first physical selectivity number -- compare across Na+, Cl-, H2O.
"""
import pathlib
import sys

try:
    import modal
except ImportError:
    modal = None
    print("[note] `pip install modal` to actually launch cloud runs.",
          file=sys.stderr)


def _barrier_from_dat(text):
    """Parse a neb .dat (reaction_coord energy error) string -> barrier dict.
    Kept inline so the Modal container needs no local imports."""
    coords, energies = [], []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        try:
            coords.append(float(parts[0]))
            energies.append(float(parts[1]))
        except (ValueError, IndexError):
            continue
    if len(energies) < 2:
        return {"error": "could not parse >=2 images"}
    imax = max(range(len(energies)), key=lambda i: energies[i])
    return {
        "forward_barrier_eV": round(energies[imax] - energies[0], 4),
        "reverse_barrier_eV": round(energies[imax] - energies[-1], 4),
        "reaction_energy_eV": round(energies[-1] - energies[0], 4),
        "saddle_image": imax + 1,
        "n_images": len(energies),
    }


if modal is not None:
    qe_image = (
        modal.Image.debian_slim(python_version="3.11")
        .apt_install("wget", "libgomp1")
        .micromamba_install("qe=7.2", "openmpi", "numpy",
                            channels=["conda-forge"])
    )
    app = modal.App("desal-neb")
    vol = modal.Volume.from_name("desal-dft-data", create_if_missing=True)

    @app.function(image=qe_image, cpu=16.0, memory=32768,
                  timeout=60 * 60 * 12, volumes={"/data": vol})
    def run_neb(infile_text: str, prefix: str = "neb"):
        """Run neb.x and return stdout tail + parsed barrier from <prefix>.dat."""
        import subprocess
        import os
        import glob

        workdir = f"/data/{prefix}"
        os.makedirs(workdir, exist_ok=True)
        # symlink shared pseudos into the run dir
        if os.path.isdir("/data/pseudo") and not os.path.exists(
                os.path.join(workdir, "pseudo")):
            os.symlink("/data/pseudo", os.path.join(workdir, "pseudo"))

        in_path = os.path.join(workdir, "neb.in")
        with open(in_path, "w") as fh:
            fh.write(infile_text)

        ncore = int(os.environ.get("MODAL_CPU", "16"))
        # neb.x parallelises over images with -ni; keep it simple with plain MPI
        cmd = ["mpirun", "-np", str(ncore), "neb.x", "-inp", in_path]
        print("Running:", " ".join(cmd))
        proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
        with open(os.path.join(workdir, "neb.out"), "w") as fh:
            fh.write(proc.stdout)
        vol.commit()

        dat_files = glob.glob(os.path.join(workdir, "*.dat"))
        barrier = {"error": "no .dat produced"}
        if dat_files:
            with open(dat_files[0]) as fh:
                barrier = _barrier_from_dat(fh.read())
        return {
            "returncode": proc.returncode,
            "barrier": barrier,
            "stderr_tail": proc.stderr[-2000:],
        }

    @app.local_entrypoint()
    def main(infile: str = "qe/neb_na_ready.in", ion: str = "Na"):
        text = pathlib.Path(infile).read_text()
        result = run_neb.remote(text, prefix=f"neb_{ion.lower()}")
        print("=" * 60)
        print("returncode:", result["returncode"])
        b = result["barrier"]
        if "forward_barrier_eV" in b:
            ea = b["forward_barrier_eV"]
            print(f"{ion}: forward barrier E_a = {ea} eV "
                  f"({ea * 96.485:.1f} kJ/mol)")
            print(f"     reverse = {b['reverse_barrier_eV']} eV, "
                  f"reaction dE = {b['reaction_energy_eV']} eV, "
                  f"saddle image {b['saddle_image']}/{b['n_images']}")
        else:
            print("barrier parse:", b)
            print("stderr:", result["stderr_tail"])
