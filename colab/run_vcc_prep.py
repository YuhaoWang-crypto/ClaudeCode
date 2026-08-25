"""Rebuild the 2026 submission and run the official validator, on a big-RAM VM.

Runtime: **High-RAM** (Runtime -> Change runtime type -> High-RAM). No GPU
needed; this is a memory job, not a compute one. ~50 GB is enough; the standard
12.7 GB runtime will be killed exactly as the 15 GB container was.

Everything is rebuilt from public data plus your own controls bundle, so nothing
large has to be moved between machines. Roughly:

    downloads   ~1.1 GB   (genome-wide K562 screen + your controls bundle)
    build       ~40 min   (360,000 cells x 18,533 genes, 2.08e9 stored values)
    prep        ~15 min   (peak RSS 40-60 GB)

Put your VCC key in Colab's Secrets panel as VCC_TOKEN. It is read from there,
never printed, and never written to disk in plaintext by this script.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = Path("/content/vcc_data")
OFFICIAL = Path("/content/vcc_official")
OUT = Path("/content/vcc_submission")
GWPS = "https://ndownloader.figshare.com/files/35774443"   # K562 gwps, 375 MB
ATLAS = {  # the four-line atlas, needed for the per-gene prior
    "K562_essential_raw_bulk.h5ad": "35773070",
    "rpe1_raw_bulk.h5ad": "35775581",
}


def sh(cmd: str, **kw) -> int:
    """Run a shell command, redacting anything that looks like a token."""
    print(f"$ {re.sub(r'vcc_pat_[A-Za-z0-9]+', '***', cmd)}", flush=True)
    return subprocess.call(cmd, shell=True, **kw)


def token() -> str:
    try:
        from google.colab import userdata            # type: ignore
        t = userdata.get("VCC_TOKEN")
        if t:
            return t.strip()
    except Exception:
        pass
    t = os.environ.get("VCC_TOKEN", "").strip()
    if not t:
        sys.exit("Set VCC_TOKEN in Colab's Secrets panel (key icon, left bar). "
                 "Do not paste it into a cell.")
    return t


def main() -> None:
    import psutil

    gb = psutil.virtual_memory().total / 1e9
    print(f"RAM available: {gb:.0f} GB")
    if gb < 30:
        print("\n!! This runtime is too small. Runtime -> Change runtime type "
              "-> High-RAM.\n   One context alone was SIGKILLed at 14 GB.")
        sys.exit(1)

    sh("pip -q install anndata scanpy h5py scipy scikit-learn pandas vcc-cli")
    DATA.mkdir(parents=True, exist_ok=True)
    OFFICIAL.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. the source screen: Replogle genome-wide K562, 375 MB pseudobulk
    gwps = DATA / "replogle" / "K562_gwps_raw_bulk.h5ad"
    gwps.parent.mkdir(parents=True, exist_ok=True)
    if not gwps.exists():
        sh(f"curl -sSL --retry 4 -o {gwps} {GWPS}")
    for name, fid in ATLAS.items():
        f = gwps.parent / name
        if not f.exists():
            sh(f"curl -sSL --retry 4 -o {f} "
               f"https://ndownloader.figshare.com/files/{fid}")

    # 2. the official controls bundle, via your own key
    if not (OFFICIAL / "gene_names.csv").exists():
        p = subprocess.run("vcc login --token-stdin", shell=True,
                           input=token(), text=True, capture_output=True)
        if p.returncode != 0:
            sys.exit(re.sub(r"vcc_pat_[A-Za-z0-9]+", "***", p.stderr)[-800:])
        sh(f"cd {OFFICIAL} && vcc datasets download controls -o controls.zip "
           f"&& unzip -oq controls.zip")

    # 3. build, with the paths this repo's modules expect
    # The modules resolve data by absolute path, so link rather than edit them.
    env = dict(os.environ, PYTHONPATH=str(REPO))
    sh("mkdir -p /home/user")
    sh(f"ln -sfn {DATA} /home/user/vcc_data 2>/dev/null || true")
    sh(f"ln -sfn {OFFICIAL} /home/user/vcc_official 2>/dev/null || true")
    sh(f"ln -sfn {OUT} /home/user/vcc_submission 2>/dev/null || true")

    rc = sh(f"cd {REPO} && PYTHONPATH={REPO} python -m virtualcell.vcc2026 "
            f"--build --out {OUT}/prediction.h5ad", env=env)
    if rc != 0:
        sys.exit("build failed")

    # 4. the check this project could not run locally
    rc = sh(f"cd {OUT} && vcc prep prediction.h5ad "
            f"-g {OFFICIAL}/gene_names.csv --perts {OFFICIAL}/pert_counts.csv "
            f"-o prediction.vcc -f 2>&1 | tail -40")

    peak = 0.0
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6
    except Exception:
        pass

    vcc = OUT / "prediction.vcc"
    print("\n=== PASTE BACK ===")
    print(f"runtime RAM      : {gb:.0f} GB")
    print(f"peak child RSS   : {peak:.1f} GB")
    print(f"vcc prep exit    : {rc}")
    print(f"prediction.vcc   : "
          f"{vcc.stat().st_size / 1e6:.0f} MB" if vcc.exists() else
          "prediction.vcc   : NOT WRITTEN")
    print("=== END ===")
    if vcc.exists():
        print("\nTo submit, from your own terminal (not from here):")
        print("  vcc submit prediction.vcc -m 'ContextTransfer, gwps K562'")


if __name__ == "__main__":
    main()
