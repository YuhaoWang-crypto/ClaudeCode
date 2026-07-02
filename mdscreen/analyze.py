"""Trajectory analysis with MDAnalysis.

Computes the standard stability/interaction observables that let you judge
whether a pose is stable and a ligand is engaging the protein:

  * backbone RMSD (protein stability / convergence)
  * per-residue RMSF (flexibility, hinge/allosteric hotspots)
  * radius of gyration (compactness / (un)folding)
  * ligand RMSD relative to frame 0 (pose stability)
  * protein-ligand contacts + H-bond count (engagement)

Outputs CSVs, PNG plots and a JSON summary of scalar metrics.
"""
from __future__ import annotations

import os
from typing import Optional

import numpy as np

from .utils import log, ensure_dir, write_json


def _universe(topology_pdb: str, trajectory_dcd: str):
    import MDAnalysis as mda
    return mda.Universe(topology_pdb, trajectory_dcd)


def analyze(topology_pdb: str, trajectory_dcd: str, outdir: str,
            ligand_resname: Optional[str] = None, tag: str = "sim") -> dict:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align

    ensure_dir(outdir)
    u = _universe(topology_pdb, trajectory_dcd)
    n_frames = len(u.trajectory)
    log.info("Analysing %d frames", n_frames)

    protein = u.select_atoms("protein")
    summary: dict = {"tag": tag, "n_frames": n_frames}

    # --- backbone RMSD ----------------------------------------------------
    if len(protein) > 0:
        rmsd_sel = "backbone"
        R = rms.RMSD(u, u, select=rmsd_sel, ref_frame=0)
        R.run()
        rmsd_ang = R.results.rmsd[:, 2]           # column 2 = RMSD in Å
        _save_series(outdir, f"{tag}_rmsd_backbone", rmsd_ang, "RMSD (Å)")
        summary["rmsd_backbone_mean_nm"] = float(np.mean(rmsd_ang) / 10.0)
        summary["rmsd_backbone_final_nm"] = float(rmsd_ang[-1] / 10.0)

        # --- per-residue RMSF (align first) -------------------------------
        rmsf = _compute_rmsf(u)
        if rmsf is not None:
            resids, values = rmsf
            _save_xy(outdir, f"{tag}_rmsf", resids, values, "Residue", "RMSF (Å)")
            summary["rmsf_mean_ang"] = float(np.mean(values))
            summary["rmsf_max_ang"] = float(np.max(values))
            # flexible residues (>1.5x mean) are candidate hinge/allosteric sites
            thr = np.mean(values) + np.std(values)
            summary["flexible_residues"] = [int(r) for r, v in zip(resids, values) if v > thr][:30]

    # --- radius of gyration ----------------------------------------------
    rg = []
    for _ in u.trajectory:
        rg.append(protein.radius_of_gyration() if len(protein) else np.nan)
    rg = np.array(rg)
    if np.isfinite(rg).any():
        _save_series(outdir, f"{tag}_rgyr", rg, "Rg (Å)")
        summary["rgyr_mean_ang"] = float(np.nanmean(rg))

    # --- ligand pose stability + engagement ------------------------------
    if ligand_resname:
        lig = u.select_atoms(f"resname {ligand_resname}")
        if len(lig) > 0:
            lig_rmsd = _ligand_rmsd(u, ligand_resname)
            _save_series(outdir, f"{tag}_ligand_rmsd", lig_rmsd, "Ligand RMSD (Å)")
            summary["ligand_rmsd_mean_nm"] = float(np.mean(lig_rmsd) / 10.0)
            summary["ligand_rmsd_final_nm"] = float(lig_rmsd[-1] / 10.0)

            contacts, hbonds = _contacts_and_hbonds(u, ligand_resname)
            summary["ligand_contacts_mean"] = float(np.mean(contacts)) if contacts else 0.0
            summary["ligand_hbonds_mean"] = float(hbonds)
            summary["contact_residues"] = _contact_residues(u, ligand_resname)

    write_json(os.path.join(outdir, f"{tag}_analysis.json"), summary)
    log.info("Analysis summary: %s", {k: v for k, v in summary.items()
                                      if isinstance(v, (int, float))})
    return summary


# ---------------------------------------------------------------------------
def _compute_rmsf(u):
    from MDAnalysis.analysis import align, rms
    import MDAnalysis as mda

    ca = u.select_atoms("protein and name CA")
    if len(ca) == 0:
        return None
    # align trajectory to average structure in memory
    average = align.AverageStructure(u, u, select="protein and name CA",
                                     ref_frame=0).run()
    ref = average.results.universe
    align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()
    rmsf = rms.RMSF(ca).run()
    return ca.resids, rmsf.results.rmsf


def _ligand_rmsd(u, ligand_resname):
    from MDAnalysis.analysis import rms
    R = rms.RMSD(u, u, select="protein and backbone",
                 groupselections=[f"resname {ligand_resname}"], ref_frame=0)
    R.run()
    # column 3 is the RMSD of the first groupselection (ligand)
    return R.results.rmsd[:, 3]


def _contacts_and_hbonds(u, ligand_resname):
    """Mean number of protein heavy atoms within 4 Å of the ligand, plus a
    crude H-bond count (donor/acceptor pairs within 3.5 Å)."""
    contacts = []
    lig = u.select_atoms(f"resname {ligand_resname} and not name H*")
    for _ in u.trajectory:
        near = u.select_atoms(
            f"protein and not name H* and around 4.0 (resname {ligand_resname})")
        contacts.append(len(near))
    # H-bonds via MDAnalysis HBA if available
    hbonds = 0.0
    try:
        from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis as HBA
        hba = HBA(universe=u,
                  between=[f"resname {ligand_resname}", "protein"])
        hba.guess_hydrogens = True
        hba.run()
        counts = hba.count_by_time()
        hbonds = float(np.mean(counts)) if len(counts) else 0.0
    except Exception as exc:  # pragma: no cover
        log.warning("H-bond analysis skipped: %s", exc)
    return contacts, hbonds


def _contact_residues(u, ligand_resname, cutoff=4.0, freq=0.5):
    """Residues in contact with the ligand in >= `freq` of frames = binding site."""
    from collections import Counter
    counter = Counter()
    n = 0
    for _ in u.trajectory:
        n += 1
        sel = u.select_atoms(
            f"protein and around {cutoff} (resname {ligand_resname})")
        for res in {(a.resname, int(a.resid)) for a in sel}:
            counter[res] += 1
    return [{"resname": r[0], "resid": r[1], "frequency": round(c / n, 3)}
            for r, c in counter.most_common(25) if c / n >= freq]


# ---------------------------------------------------------------------------
def _save_series(outdir, name, y, ylabel):
    x = np.arange(len(y))
    _save_xy(outdir, name, x, y, "Frame", ylabel)


def _save_xy(outdir, name, x, y, xlabel, ylabel):
    import csv
    csv_path = os.path.join(outdir, f"{name}.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([xlabel, ylabel])
        for xi, yi in zip(x, y):
            w.writerow([xi, yi])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(x, y, lw=1.2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(name)
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=120)
        plt.close(fig)
    except Exception as exc:  # pragma: no cover
        log.warning("Plot for %s skipped: %s", name, exc)
