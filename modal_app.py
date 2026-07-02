"""Modal app: run the mdscreen MD/analysis pipeline on serverless GPUs.

Usage (after `pip install modal` and `modal token new`):

    # single protein-ligand run
    modal run modal_app.py::complex --pdb 3poz.pdb --ligand "CCO" --ns 2

    # protein-only MD
    modal run modal_app.py::protein --pdb 3poz.pdb --ns 5

    # small-molecule-only MD
    modal run modal_app.py::ligand --ligand ligand.sdf --ns 5

    # full activity/allostery screen from a YAML config
    modal run modal_app.py::screen --config examples/screen_config.yaml

Outputs (trajectories, plots, CSVs, JSON) are written to a persistent Modal
Volume named ``mdscreen-outputs`` and can be pulled locally with:

    modal volume get mdscreen-outputs <run_id> ./results
"""
from __future__ import annotations

import os
from typing import Optional

import modal

# ---------------------------------------------------------------------------
# Image: OpenMM + OpenFF + AmberTools + analysis stack, built on conda-forge.
# ---------------------------------------------------------------------------
image = (
    modal.Image.micromamba(python_version="3.11")
    .micromamba_install(
        "openmm=8.1",
        "openmmforcefields>=0.12",
        "openff-toolkit>=0.15",
        "openff-interchange",
        "ambertools>=23",
        "pdbfixer",
        "rdkit",
        "mdanalysis>=2.7",
        "numpy",
        "scipy",
        "matplotlib",
        "networkx",
        "pyyaml",
        channels=["conda-forge"],
    )
    # CUDA runtime so OpenMM's CUDA platform is available on the GPU worker.
    .env({"OPENMM_CPU_THREADS": "4"})
    .add_local_python_source("mdscreen")
)

app = modal.App("mdscreen", image=image)
volume = modal.Volume.from_name("mdscreen-outputs", create_if_missing=True)

GPU = os.environ.get("MDSCREEN_GPU", "A10G")          # override: L4, A100, H100
TIMEOUT = 60 * 60 * 6                                   # 6h per job

OUT_ROOT = "/outputs"


# ---------------------------------------------------------------------------
# helpers shared by every remote function
# ---------------------------------------------------------------------------
def _stage_ref(ref, suffix: str, workdir: str) -> str:
    """Materialise an input inside the container.

    ``ref`` is either a 4-letter PDB id, a SMILES string, or (bytes, filename).
    Returns a value suitable to pass to the pipeline (path or id/SMILES).
    """
    if isinstance(ref, (tuple, list)) and len(ref) == 2:
        data, fname = ref
        path = os.path.join(workdir, fname)
        mode = "wb" if isinstance(data, (bytes, bytearray)) else "w"
        with open(path, mode) as fh:
            fh.write(data)
        return path
    return ref  # PDB id or SMILES: passed through untouched


def _sim_cfg(cfg_dict: dict):
    from mdscreen.config import SimConfig
    return SimConfig.from_dict(cfg_dict or {})


# ---------------------------------------------------------------------------
# remote GPU functions
# ---------------------------------------------------------------------------
@app.function(gpu=GPU, timeout=TIMEOUT, volumes={OUT_ROOT: volume})
def complex_md(run_id: str, name: str, pdb_ref, ligand_ref, cfg_dict: dict,
               compute_binding: bool = True) -> dict:
    import tempfile
    from mdscreen import pipeline
    work = tempfile.mkdtemp()
    pdb = _stage_ref(pdb_ref, ".pdb", work)
    ligand = _stage_ref(ligand_ref, ".sdf", work)
    outdir = os.path.join(OUT_ROOT, run_id, name)
    summary = pipeline.run_complex_md(pdb, ligand, _sim_cfg(cfg_dict), outdir,
                                      compute_binding=compute_binding)
    volume.commit()
    return summary


@app.function(gpu=GPU, timeout=TIMEOUT, volumes={OUT_ROOT: volume})
def protein_md(run_id: str, pdb_ref, cfg_dict: dict, name: str = "apo") -> dict:
    import tempfile
    from mdscreen import pipeline
    work = tempfile.mkdtemp()
    pdb = _stage_ref(pdb_ref, ".pdb", work)
    outdir = os.path.join(OUT_ROOT, run_id, name)
    summary = pipeline.run_protein_md(pdb, _sim_cfg(cfg_dict), outdir)
    volume.commit()
    return summary


@app.function(gpu=GPU, timeout=TIMEOUT, volumes={OUT_ROOT: volume})
def ligand_md(run_id: str, ligand_ref, cfg_dict: dict, name: str = "lig",
              solvate: bool = True) -> dict:
    import tempfile
    from mdscreen import pipeline
    work = tempfile.mkdtemp()
    ligand = _stage_ref(ligand_ref, ".sdf", work)
    outdir = os.path.join(OUT_ROOT, run_id, name)
    summary = pipeline.run_ligand_md(ligand, _sim_cfg(cfg_dict), outdir,
                                     solvate=solvate)
    volume.commit()
    return summary


@app.function(gpu=GPU, timeout=TIMEOUT, volumes={OUT_ROOT: volume})
def allostery_compare(run_id: str, apo_name: str, holo_name: str,
                      site_resids: list) -> dict:
    from mdscreen import allostery
    base = os.path.join(OUT_ROOT, run_id)
    apo = os.path.join(base, apo_name)
    holo = os.path.join(base, holo_name)
    out = os.path.join(holo, "allostery")
    result = allostery.compare_apo_holo(
        os.path.join(apo, "apo_topology.pdb"),
        os.path.join(apo, "apo_production.dcd"),
        os.path.join(holo, "complex_topology.pdb"),
        os.path.join(holo, "complex_production.dcd"),
        out, binding_site_resids=site_resids)
    volume.commit()
    return result


# ---------------------------------------------------------------------------
# local orchestration (runs on your machine, dispatches to GPUs)
# ---------------------------------------------------------------------------
def _ref_from_arg(arg: str):
    """Convert a CLI argument into a remote-safe reference (upload local files)."""
    from mdscreen.utils import is_smiles
    if arg and os.path.exists(arg):
        with open(arg, "rb") as fh:
            return (fh.read(), os.path.basename(arg))
    return arg  # PDB id or SMILES


def _run_id() -> str:
    # deterministic-ish id without Date.now (Modal provides an app run context)
    import uuid
    return "run_" + uuid.uuid4().hex[:10]


@app.local_entrypoint()
def complex(pdb: str, ligand: str, ns: float = 2.0, name: str = "cpx",
            gpu: Optional[str] = None):
    from mdscreen.config import SimConfig
    cfg = SimConfig(production_ns=ns).to_dict()
    rid = _run_id()
    print(f"[mdscreen] run_id={rid}  gpu={gpu or GPU}")
    summary = complex_md.remote(rid, name, _ref_from_arg(pdb),
                                _ref_from_arg(ligand), cfg)
    _print_summary(summary)
    print(f"\nPull results:  modal volume get mdscreen-outputs {rid} ./results")


@app.local_entrypoint()
def protein(pdb: str, ns: float = 5.0):
    from mdscreen.config import SimConfig
    cfg = SimConfig(production_ns=ns).to_dict()
    rid = _run_id()
    print(f"[mdscreen] run_id={rid}")
    summary = protein_md.remote(rid, _ref_from_arg(pdb), cfg)
    _print_summary(summary)
    print(f"\nPull results:  modal volume get mdscreen-outputs {rid} ./results")


@app.local_entrypoint()
def ligand(ligand: str, ns: float = 5.0, solvate: bool = True):
    from mdscreen.config import SimConfig
    cfg = SimConfig(production_ns=ns).to_dict()
    rid = _run_id()
    print(f"[mdscreen] run_id={rid}")
    summary = ligand_md.remote(rid, _ref_from_arg(ligand), cfg, solvate=solvate)
    _print_summary(summary)
    print(f"\nPull results:  modal volume get mdscreen-outputs {rid} ./results")


@app.local_entrypoint()
def screen(config: str):
    """Parallel activity + allostery screen. Fans ligands out across GPUs."""
    import yaml
    from mdscreen.config import ScreenConfig
    from mdscreen.pipeline import classify_ligand
    from mdscreen.utils import write_json

    with open(config) as fh:
        sc = ScreenConfig.from_dict(yaml.safe_load(fh))
    rid = _run_id()
    cfg = sc.sim.to_dict()
    names = sc.ligand_names or [f"L{i+1}" for i in range(len(sc.ligands))]
    pdb_ref = _ref_from_arg(sc.receptor)
    print(f"[mdscreen] screen run_id={rid}  ligands={len(sc.ligands)}  gpu={GPU}")

    # apo reference (one job) for allostery
    apo_summary = None
    if sc.run_apo_reference and sc.compute_allostery:
        apo_summary = protein_md.remote(rid, pdb_ref, cfg, name="apo")

    # complex jobs in parallel across GPUs
    args = [(rid, name, pdb_ref, _ref_from_arg(lig), cfg, sc.compute_binding)
            for lig, name in zip(sc.ligands, names)]
    summaries = list(complex_md.starmap(args))

    results = []
    for name, lig, summary in zip(names, sc.ligands, summaries):
        verdict = classify_ligand(summary, sc)
        if sc.compute_allostery and apo_summary is not None:
            site = [c["resid"] for c in summary.get("contact_residues", [])]
            try:
                verdict["allostery"] = allostery_compare.remote(rid, "apo", name, site)
            except Exception as exc:
                print(f"  allostery failed for {name}: {exc}")
        results.append({"name": name, "ligand": lig, "verdict": verdict})

    results.sort(key=lambda r: (r["verdict"].get("dg_bind_kcal")
                                if r["verdict"].get("dg_bind_kcal") is not None else 1e9))
    for i, r in enumerate(results, 1):
        r["rank"] = i

    print("\n==================== SCREEN RANKING ====================")
    print(f"{'rank':<5}{'name':<10}{'ΔG(kcal/mol)':<14}{'pose':<8}{'verdict'}")
    for r in results:
        v = r["verdict"]
        dg = v.get("dg_bind_kcal")
        print(f"{r['rank']:<5}{r['name']:<10}"
              f"{(f'{dg:.2f}' if dg is not None else 'n/a'):<14}"
              f"{('stable' if v.get('stable_pose') else 'drift'):<8}"
              f"{v.get('verdict')}")
    print(f"\nPull results:  modal volume get mdscreen-outputs {rid} ./results")
    # also stash the report in the volume
    _persist_report.remote(rid, {"receptor": sc.receptor, "ranked_results": results})


@app.function(volumes={OUT_ROOT: volume})
def _persist_report(run_id: str, report: dict):
    from mdscreen.utils import write_json
    write_json(os.path.join(OUT_ROOT, run_id, "screen_report.json"), report)
    volume.commit()


def _print_summary(summary: dict):
    print("\n---------------- summary ----------------")
    for k, v in summary.items():
        if isinstance(v, (int, float, str)):
            print(f"  {k:<28} {v}")
    mm = summary.get("mmgbsa")
    if isinstance(mm, dict) and "dg_bind_kcal" in mm:
        print(f"  {'MM/GBSA ΔG (kcal/mol)':<28} {mm['dg_bind_kcal']:.2f} "
              f"± {mm.get('dg_std_kcal', 0):.2f}")
