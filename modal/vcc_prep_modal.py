"""Run the official validator at full scale, on a machine that can hold it.

``vcc prep`` keeps the whole submission as an in-memory CSR.  Measured on the
15 GB container this project runs in:

    40,000 cells   238 M stored values    6.9 GB peak   passes
   120,000 cells   715 M stored values   14.0 GB peak   SIGKILL
   360,000 cells  2,078 M stored values   est. 40-60 GB  not attemptable

The file itself is well inside the scorer's own 4.75e9 cap, so this is a local
validator limit rather than a spec problem -- but it is the last unticked check
on the submission, and it needs about 64 GB.

Nothing large is uploaded.  The submission is *rebuilt* on Modal from public
data plus the official controls bundle, which is the same code path that
produced it locally, so a pass here is a pass for the real artifact.

    modal run modal/vcc_prep_modal.py --stage build     # ~40 min, writes to a Volume
    modal run modal/vcc_prep_modal.py --stage prep      # ~15 min at 64 GB
    modal run modal/vcc_prep_modal.py --stage fetch     # pull the .vcc back

No credential is sent to Modal.  The controls bundle and the source screen are
uploaded once with ``modal volume put``; the token stays where it already is.
"""

from __future__ import annotations

import modal

REPO = "https://github.com/yuhaowang-crypto/claudecode.git"
BRANCH = "claude/virtual-cell-model-prediction-qu31ry"
GWPS = "https://ndownloader.figshare.com/files/35774443"

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("git", "curl", "unzip")
    .pip_install("numpy", "scipy", "pandas", "anndata", "h5py", "scikit-learn",
                 "vcc-cli", "psutil")
    .run_commands(f"git clone --depth 1 -b {BRANCH} {REPO} /repo")
)

app = modal.App("vcc-prep")
vol = modal.Volume.from_name("vcc-submission", create_if_missing=True)
DATA = "/vol"


def _redact(s: str) -> str:
    import re
    return re.sub(r"vcc_pat_[A-Za-z0-9]+", "***", s)


@app.function(image=image, volumes={DATA: vol}, timeout=60 * 90, cpu=8,
              memory=32768)
def build() -> str:
    """Rebuild the submission on the Volume, from public data + the controls."""
    import os
    import subprocess
    from pathlib import Path

    for p in ("vcc_data/replogle", "vcc_official", "vcc_submission"):
        Path(f"{DATA}/{p}").mkdir(parents=True, exist_ok=True)
    for src, dst in (("vcc_data", "/home/user/vcc_data"),
                     ("vcc_official", "/home/user/vcc_official"),
                     ("vcc_submission", "/home/user/vcc_submission")):
        Path("/home/user").mkdir(exist_ok=True)
        if not Path(dst).exists():
            os.symlink(f"{DATA}/{src}", dst)

    gwps = Path(f"{DATA}/vcc_data/replogle/K562_gwps_raw_bulk.h5ad")
    if not gwps.exists():
        subprocess.check_call(["curl", "-sSL", "--retry", "4", "-o", str(gwps), GWPS])
    for need in ("vcc_official/gene_names.csv", "vcc_official/context_A.h5ad"):
        if not Path(f"{DATA}/{need}").exists():
            raise SystemExit(
                f"missing {need} on the volume -- upload the controls bundle "
                f"first:\n  modal volume put vcc-submission <local> /{need}")

    out = subprocess.run(
        "cd /repo && PYTHONPATH=/repo python -m virtualcell.vcc2026 --build "
        "--out /vol/vcc_submission/prediction.h5ad",
        shell=True, capture_output=True, text=True)
    vol.commit()
    return _redact(out.stdout[-3000:] + out.stderr[-2000:])


@app.function(image=image, volumes={DATA: vol}, timeout=60 * 60, cpu=8,
              memory=65536)
def prep() -> str:
    """The check the 15 GB container could not run."""
    import resource
    import subprocess

    out = subprocess.run(
        "cd /vol/vcc_submission && vcc prep prediction.h5ad "
        "-g /vol/vcc_official/gene_names.csv "
        "--perts /vol/vcc_official/pert_counts.csv -o prediction.vcc -f",
        shell=True, capture_output=True, text=True)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6
    vol.commit()

    from pathlib import Path
    vcc = Path("/vol/vcc_submission/prediction.vcc")
    size = f"{vcc.stat().st_size / 1e6:.0f} MB" if vcc.exists() else "NOT WRITTEN"
    return (_redact(out.stdout[-4000:] + out.stderr[-2000:])
            + f"\n\nexit={out.returncode}  peak child RSS={peak:.1f} GB"
            + f"\nprediction.vcc: {size}")


@app.function(image=image, volumes={DATA: vol}, timeout=60 * 30, memory=8192)
def manifest() -> dict:
    """What is on the Volume, so nothing has to be guessed."""
    from pathlib import Path
    return {str(p.relative_to(DATA)): p.stat().st_size
            for p in Path(DATA).rglob("*") if p.is_file()}


@app.local_entrypoint()
def main(stage: str = "prep"):
    if stage == "build":
        print(build.remote())
    elif stage == "prep":
        print(prep.remote())
    elif stage == "manifest":
        for k, v in sorted(manifest.remote().items()):
            print(f"{v / 1e6:>10.1f} MB  {k}")
    else:
        raise SystemExit("stage must be build, prep or manifest")
