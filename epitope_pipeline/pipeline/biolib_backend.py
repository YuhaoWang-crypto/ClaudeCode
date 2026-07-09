"""BioLib cloud backend (pybiolib) for DTU tools published on biolib.com.

Some DTU tools are published on BioLib and run in the BioLib cloud with no
local install and no license download — you only need a free BioLib account
(``biolib login`` or the ``BIOLIB_TOKEN`` env var). Confirmed apps:

    DTU/BepiPred-3     linear/discontinuous B-cell epitopes   (FASTA in)
    DTU/DiscoTope-3    conformational B-cell epitopes         (PDB in)
    DTU/NetSurfP-3     surface accessibility / secondary str  (FASTA in)
    DTU/DeepTMHMM      transmembrane topology                 (FASTA in)

The licensed MHC binders (NetMHCpan / NetMHCIIpan / NetCTLpan) are NOT on
BioLib and must use the local-CLI backend.

This module is imported lazily so the CLI-only path never needs pybiolib.
"""

from __future__ import annotations

import os
from typing import List, Optional


class BiolibNotAvailable(RuntimeError):
    pass


def _import_biolib():
    try:
        import biolib  # noqa
        return biolib
    except ImportError as exc:  # pragma: no cover
        raise BiolibNotAvailable(
            "pybiolib not installed. `pip install pybiolib` and set "
            "BIOLIB_TOKEN or run `biolib login`.") from exc


def is_authenticated() -> bool:
    """True if a BIOLIB_TOKEN / signed-in session is present.

    Loading an app works anonymously, but running a job needs identity for
    most apps; we treat presence of a token as 'can run'. Apps flagged
    requires_user_identity=False may also run anonymously.
    """
    return bool(os.environ.get("BIOLIB_TOKEN") or os.environ.get("BIOLIB_API_TOKEN"))


def run_app(app_uri: str, args: str, outdir: str,
            timeout: int = 1800) -> str:
    """Run a BioLib app and save its output files to `outdir`.

    Returns the absolute output directory. Raises on failure so the runner can
    downgrade the model to an 'error' result without crashing the batch.
    """
    biolib = _import_biolib()
    os.makedirs(outdir, exist_ok=True)
    app = biolib.load(app_uri)
    job = app.cli(args=args, blocking=True)
    code = job.get_exit_code()
    if code != 0:
        stderr = b""
        try:
            stderr = job.get_stderr()
        except Exception:
            pass
        raise RuntimeError(
            f"{app_uri} exited {code}: {stderr.decode(errors='ignore')[:300]}")
    job.save_files(output_dir=outdir)
    return os.path.abspath(outdir)
