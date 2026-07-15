"""Run the β-CD / PFAS binding pipeline on Modal GPUs.

    pip install modal && modal setup           # one-time
    modal run cd_pfas_md/modal_app.py                       # smoke test (dye APR, tiny sims)
    modal run cd_pfas_md/modal_app.py --mode prod --guest dye
    modal run cd_pfas_md/modal_app.py::screen --mode prod --guest pfoa
    modal volume get cd-pfas-work /apr_dye ./pulled          # fetch results

Architecture
------------
* One Modal Volume ("cd-pfas-work") holds everything under /work, so parameter
  files, per-window trajectories, and result JSON persist across containers and
  can be pulled to your laptop.
* Parameterization + build + analysis run on CPU containers (AmberTools/pAPRika).
* Each APR window and each FEP λ-leg is its OWN GPU container, fanned out with
  Modal `.starmap`, so N windows cost wall-clock ≈ 1 window, not N.
* `--mode smoke` shrinks every simulation to seconds so you can validate the
  whole graph on a T4 before spending real GPU-hours in `--mode prod`.

The MD engine is OpenMM with the CUDA platform (conda-forge openmm ships CUDA
builds; Modal supplies the driver). See MODAL.md for GPU sizing + cost.
"""
from __future__ import annotations

import modal

# --------------------------------------------------------------------------
# Container image: micromamba conda env identical to environment.yml, plus the
# package source mounted at /root/cd_pfas_md so `import cd_pfas_md` works.
# --------------------------------------------------------------------------
image = (
    modal.Image.micromamba(python_version="3.10")
    .micromamba_install(
        "openmm=8.*",
        "openmmtools=0.23.*",
        "pymbar=4.*",
        "paprika=1.2.*",
        "ambertools=23.*",
        "openff-toolkit=0.14.*",
        "rdkit=2023.*",
        "mdanalysis",
        "numpy",
        "pandas",
        "pyyaml",
        "scipy",
        channels=["conda-forge"],
    )
    .env({"PYTHONPATH": "/root"})
    # ship the package (configs + src) into the image
    .add_local_dir(".", remote_path="/root/cd_pfas_md")
)

app = modal.App("cd-pfas-md", image=image)
vol = modal.Volume.from_name("cd-pfas-work", create_if_missing=True)

WORK = "/work"                 # Volume mount point == pipeline work dir
CONFIG = "/root/cd_pfas_md/config/system.yaml"
MODS = "/root/cd_pfas_md/config/modifications.yaml"

# smoke-test overrides: make every leg finish in seconds to validate plumbing.
SMOKE = {
    "minimize_steps": 200, "thermalize_ps": 2, "equilibrate_ps": 4,
    "production_ps": 10, "report_interval_ps": 2,
}
SMOKE_FEP = {"per_window_ps": 10, "lambda_windows": 6, "iterations": 5}
GPU_TYPE = "A10G"              # T4 (cheap) | A10G (default) | A100 (large systems)


def _load_cfg(mode: str):
    from cd_pfas_md.src import utils
    cfg = utils.load_config(CONFIG)
    if mode == "smoke":
        cfg["simulation"].update(SMOKE)
    return cfg


# ==========================================================================
# APR legs
# ==========================================================================
@app.function(cpu=4.0, volumes={WORK: vol}, timeout=3600)
def parameterize_all(mode: str = "smoke") -> str:
    from pathlib import Path
    from cd_pfas_md.src import parameterize
    cfg = _load_cfg(mode)
    T = cfg["thermo"]["temperature_K"]
    out = Path(WORK) / "params"
    parameterize.parameterize_host(cfg["host"], T, out / "host")
    from cd_pfas_md.src import utils
    for key, gcfg in cfg["guests"].items():
        parameterize.parameterize_guest(utils.Guest.from_config(key, gcfg, T), out / key)
    vol.commit()
    return str(out)


@app.function(cpu=4.0, volumes={WORK: vol}, timeout=3600)
def build_apr_stage(guest: str, mode: str = "smoke") -> list[str]:
    from cd_pfas_md.src import build_apr
    vol.reload()
    work = build_apr.build(CONFIG, f"{WORK}/params", WORK, guest)
    # resolve heuristic anchor atom indices into the manifest so windows can run
    from cd_pfas_md.src import anchors
    anchors.resolve_into_manifest(work)
    vol.commit()
    import json
    manifest = json.loads((work / "apr_manifest.json").read_text())
    return manifest["windows"]


@app.function(gpu=GPU_TYPE, volumes={WORK: vol}, timeout=7200)
def run_apr_window(guest: str, window: str, mode: str = "smoke") -> str:
    from pathlib import Path
    from cd_pfas_md.src import run_apr
    vol.reload()
    cfg = _load_cfg(mode)
    wdir = run_apr.run_window(Path(WORK) / f"apr_{guest}", window, cfg)
    vol.commit()
    return f"{window}:done"


@app.function(cpu=4.0, volumes={WORK: vol}, timeout=3600)
def analyze_apr_stage(guest: str, mode: str = "smoke") -> dict:
    import json
    from pathlib import Path
    from cd_pfas_md.src import analyze_apr, utils
    vol.reload()
    cfg = _load_cfg(mode)
    work = Path(WORK) / f"apr_{guest}"
    result = analyze_apr.compute_binding_free_energy(work, cfg)
    g = utils.Guest.from_config(guest, cfg["guests"][guest], cfg["thermo"]["temperature_K"])
    report = analyze_apr.calibrate(result, g)
    (work / "binding_result.json").write_text(json.dumps(report, indent=2))
    vol.commit()
    return report


# ==========================================================================
# FEP/TI ΔΔG screen (host WT->MOD, complex + apo legs)
# ==========================================================================
@app.function(gpu=GPU_TYPE, volumes={WORK: vol}, timeout=10800)
def fep_modification(mod_id: str, guest: str, mode: str = "smoke") -> dict:
    from dataclasses import asdict
    from pathlib import Path
    from cd_pfas_md.src import fep_ti_ddg, utils
    vol.reload()
    cfg_sys = _load_cfg(mode)
    cfg_screen = utils.load_config(MODS)
    if mode == "smoke":
        cfg_screen.update(SMOKE_FEP)
    mod = next(m for m in cfg_screen["modifications"] if m["id"] == mod_id)
    res = fep_ti_ddg.screen_modification(
        mod, guest, cfg_sys, cfg_screen, Path(f"{WORK}/params"), Path(f"{WORK}/ddg")
    )
    vol.commit()
    return asdict(res)


# ==========================================================================
# Prefilter (CPU, seconds) — cheap triage, no structures needed
# ==========================================================================
@app.function(cpu=2.0)
def prefilter() -> dict:
    from cd_pfas_md.src import prefilter
    return prefilter.run(CONFIG, MODS)


# ==========================================================================
# Orchestrators
# ==========================================================================
@app.local_entrypoint()
def main(guest: str = "dye", mode: str = "smoke"):
    """Full APR calibration for one guest: parameterize -> build -> windows -> analyze."""
    from cd_pfas_md.src import build_apr, utils
    cfg = utils.load_config("cd_pfas_md/config/system.yaml")
    windows = build_apr.create_window_list_safe(
        cfg["apr"]["attach_fractions"], cfg["apr"]["pull_distances_ang"]
    )
    print(f"[modal] parameterizing (mode={mode}) ...")
    parameterize_all.remote(mode)
    print(f"[modal] building {len(windows)} APR windows for guest={guest} ...")
    build_apr_stage.remote(guest, mode)
    print(f"[modal] running {len(windows)} windows on {GPU_TYPE} in parallel ...")
    done = list(run_apr_window.starmap([(guest, w, mode) for w in windows]))
    print(f"[modal] {len(done)} windows finished; analyzing ...")
    report = analyze_apr_stage.remote(guest, mode)
    print("\n=== APR binding result ===")
    for k, v in report.items():
        print(f"  {k}: {v}")
    print("\nPull full outputs with:  modal volume get cd-pfas-work /apr_%s ./pulled" % guest)


@app.local_entrypoint()
def screen(guest: str = "pfoa", mode: str = "smoke"):
    """FEP/TI ΔΔG screen: fan every modification across GPU containers, rank."""
    from cd_pfas_md.src import utils
    mods = utils.load_config("cd_pfas_md/config/modifications.yaml")
    ids = [m["id"] for m in mods["modifications"]]
    print(f"[modal] ensuring parameters ...")
    parameterize_all.remote(mode)
    print(f"[modal] screening {len(ids)} modifications vs {guest} on {GPU_TYPE} ...")
    results = list(fep_modification.starmap([(mid, guest, mode) for mid in ids]))
    results.sort(key=lambda r: (r["ddG_bind_kcal"] != r["ddG_bind_kcal"], r["ddG_bind_kcal"]))
    print(f"\n=== ΔΔG screen vs {guest} (negative = tighter than β-CD) ===")
    for i, r in enumerate(results, 1):
        v = r["ddG_bind_kcal"]
        s = f"{v:+.2f} ± {r['ddG_sem']:.2f}" if v == v else "  n/a (skipped/needs structure)"
        print(f"  {i:>2}. {r['modification_id']:<26} {s:>22}  {r['note'][:44]}")


@app.local_entrypoint()
def triage():
    """Just the CPU prefilter — instant, no GPU, no structures required."""
    from cd_pfas_md.src import prefilter as pf
    report = prefilter.remote()          # runs in a Modal CPU container
    pf.print_report(report)
