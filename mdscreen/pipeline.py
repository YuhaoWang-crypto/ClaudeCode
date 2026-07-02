"""High-level workflows that chain prepare -> simulate -> analyse -> score.

These functions are engine-agnostic: they run identically on a local machine
(if the MD stack is installed) or inside a Modal GPU container. ``modal_app.py``
simply calls them remotely.
"""
from __future__ import annotations

import os
from typing import List, Optional

from .config import SimConfig, ScreenConfig
from .utils import log, ensure_dir, write_json


# ---------------------------------------------------------------------------
# Single workflows
# ---------------------------------------------------------------------------
def run_protein_md(pdb: str, cfg: SimConfig, outdir: str) -> dict:
    from . import prepare, simulate, analyze
    ensure_dir(outdir)
    prepared = prepare.build_protein(pdb, cfg)
    sim = simulate.run(prepared, cfg, outdir, tag="apo")
    summary = analyze.analyze(sim.topology_pdb, sim.trajectory_dcd, outdir,
                              ligand_resname=None, tag="apo")
    summary["artifacts"] = sim.__dict__
    return summary


def run_ligand_md(ligand: str, cfg: SimConfig, outdir: str,
                  name: str = "LIG", solvate: bool = True) -> dict:
    from . import prepare, simulate, analyze
    ensure_dir(outdir)
    prepared = prepare.build_ligand_only(ligand, cfg, name, solvate=solvate)
    sim = simulate.run(prepared, cfg, outdir, tag="lig")
    summary = analyze.analyze(sim.topology_pdb, sim.trajectory_dcd, outdir,
                              ligand_resname=sim.ligand_resname, tag="lig")
    summary["artifacts"] = sim.__dict__
    return summary


def run_complex_md(pdb: str, ligand: str, cfg: SimConfig, outdir: str,
                   name: str = "LIG", compute_binding: bool = True) -> dict:
    from . import prepare, simulate, analyze, binding
    ensure_dir(outdir)
    prepared = prepare.build_complex(pdb, ligand, cfg, name)
    sim = simulate.run(prepared, cfg, outdir, tag="complex")
    summary = analyze.analyze(sim.topology_pdb, sim.trajectory_dcd, outdir,
                              ligand_resname=sim.ligand_resname, tag="complex")
    if compute_binding:
        try:
            b = binding.compute_mmgbsa(
                sim.topology_pdb, sim.trajectory_dcd, pdb, ligand, cfg, outdir,
                ligand_name=name)
            summary["mmgbsa"] = b.__dict__
        except Exception as exc:  # pragma: no cover
            log.warning("MM/GBSA failed: %s", exc)
            summary["mmgbsa"] = {"error": str(exc)}
    summary["artifacts"] = sim.__dict__
    return summary


# ---------------------------------------------------------------------------
# Screening campaign
# ---------------------------------------------------------------------------
def classify_ligand(summary: dict, sc: ScreenConfig) -> dict:
    """Turn raw metrics into an activity verdict."""
    dg = (summary.get("mmgbsa") or {}).get("dg_bind_kcal")
    pose_rmsd = summary.get("ligand_rmsd_mean_nm")
    contacts = summary.get("ligand_contacts_mean", 0.0)

    stable_pose = pose_rmsd is not None and pose_rmsd <= sc.pose_rmsd_stable_nm
    verdict = "inconclusive"
    confidence = "low"
    if dg is not None:
        if dg <= sc.strong_dg_threshold and stable_pose:
            verdict, confidence = "likely active (strong binder, stable pose)", "high"
        elif dg <= sc.active_dg_threshold and stable_pose:
            verdict, confidence = "likely active", "medium"
        elif dg <= sc.active_dg_threshold and not stable_pose:
            verdict, confidence = "borderline (favourable ΔG but drifting pose)", "low"
        else:
            verdict, confidence = "likely inactive (weak binding)", "medium"
    return {
        "dg_bind_kcal": dg,
        "pose_rmsd_mean_nm": pose_rmsd,
        "stable_pose": stable_pose,
        "contacts_mean": contacts,
        "verdict": verdict,
        "confidence": confidence,
    }


def run_screen(sc: ScreenConfig, outdir: str,
               runner=None) -> dict:
    """Run a full activity/allostery screen over a ligand series.

    ``runner`` is an optional callable ``(fn_name, kwargs) -> dict`` used to
    dispatch each per-ligand job (e.g. to Modal). If None, jobs run in-process.
    """
    from . import allostery
    ensure_dir(outdir)
    results = []
    names = sc.ligand_names or [f"L{i+1}" for i in range(len(sc.ligands))]

    # Optional apo reference for allostery.
    apo_artifacts = None
    if sc.run_apo_reference and (sc.compute_allostery):
        apo_dir = os.path.join(outdir, "apo")
        apo_summary = _dispatch(runner, "run_protein_md",
                                dict(pdb=sc.receptor, cfg=sc.sim.to_dict(),
                                     outdir=apo_dir))
        apo_artifacts = apo_summary.get("artifacts")

    for lig, name in zip(sc.ligands, names):
        lig_dir = os.path.join(outdir, name)
        log.info("=== Screening ligand %s ===", name)
        summary = _dispatch(runner, "run_complex_md",
                            dict(pdb=sc.receptor, ligand=lig, cfg=sc.sim.to_dict(),
                                 outdir=lig_dir, name="LIG",
                                 compute_binding=sc.compute_binding))
        verdict = classify_ligand(summary, sc)

        # Allostery: compare holo vs apo.
        if sc.compute_allostery and apo_artifacts:
            holo = summary.get("artifacts", {})
            site = [c["resid"] for c in summary.get("contact_residues", [])]
            try:
                allo = allostery.compare_apo_holo(
                    apo_artifacts["topology_pdb"], apo_artifacts["trajectory_dcd"],
                    holo["topology_pdb"], holo["trajectory_dcd"],
                    os.path.join(lig_dir, "allostery"), binding_site_resids=site)
                verdict["allostery"] = allo
            except Exception as exc:  # pragma: no cover
                log.warning("Allostery comparison failed for %s: %s", name, exc)

        results.append({"name": name, "ligand": lig, "verdict": verdict,
                        "metrics": {k: v for k, v in summary.items()
                                    if isinstance(v, (int, float, str))}})

    # Rank by ΔG (most negative first).
    def _key(r):
        dg = r["verdict"].get("dg_bind_kcal")
        return dg if dg is not None else 1e9
    results.sort(key=_key)
    for rank, r in enumerate(results, 1):
        r["rank"] = rank

    report = {"receptor": sc.receptor, "n_ligands": len(sc.ligands),
              "ranked_results": results}
    write_json(os.path.join(outdir, "screen_report.json"), report)
    _write_screen_table(outdir, results)
    log.info("Screen complete: %d ligands ranked", len(results))
    return report


def _dispatch(runner, fn_name: str, kwargs: dict) -> dict:
    if runner is not None:
        return runner(fn_name, kwargs)
    # in-process: reconstruct SimConfig and call directly
    from . import pipeline as _self
    kwargs = dict(kwargs)
    if "cfg" in kwargs and isinstance(kwargs["cfg"], dict):
        kwargs["cfg"] = SimConfig.from_dict(kwargs["cfg"])
    return getattr(_self, fn_name)(**kwargs)


def _write_screen_table(outdir: str, results: List[dict]) -> None:
    import csv
    path = os.path.join(outdir, "screen_ranking.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "name", "dg_bind_kcal", "pose_rmsd_nm",
                    "stable_pose", "verdict", "confidence", "allosteric_verdict"])
        for r in results:
            v = r["verdict"]
            allo = (v.get("allostery") or {}).get("allosteric_verdict", "")
            w.writerow([r.get("rank"), r["name"], v.get("dg_bind_kcal"),
                        v.get("pose_rmsd_mean_nm"), v.get("stable_pose"),
                        v.get("verdict"), v.get("confidence"), allo])
