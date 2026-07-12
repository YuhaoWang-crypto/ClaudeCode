"""
Demo #4 — Molten-salt property workflow + DeePMD-kit training-set packaging.

Motivation (from the chat context)
-----------------------------------
Client (郭硕): molten-salt system, 50-100 atoms, AIMD sampling at 6 temperature
points, DFT functional, wants a DeePMD training set and ultimately macroscopic
properties (density, viscosity, diffusion). The service provider (九章算) needs
to show the end-to-end pipeline: sample structures -> label energy/forces ->
package for DeePMD -> train ML potential -> run long MD -> extract properties.

WHAT THIS IS
------------
A runnable, end-to-end demonstration of that WORKFLOW using a cheap, transparent
Lennard-Jones liquid as a stand-in for the real molten salt:
    1. run NVT molecular dynamics at several "temperature points" (ASE),
    2. extract macroscopic/structural observables per temperature
       (radial distribution function g(r) and self-diffusion coefficient D),
    3. package the labelled frames (energy + forces) into the exact
       DeePMD-kit `deepmd/npy` directory format that the client asked for.

WHAT THIS IS *NOT*
------------------
The forces/energies here come from a Lennard-Jones potential, NOT DFT. The real
job needs CP2K/VASP AIMD (PBE or hybrid), 50-100-atom molten-salt cells, proper
Born-Mayer-Huggins/Coulomb interactions, and DeePMD/NEP training on HPC. This
script proves the plumbing and the deliverable format; swap the ASE LJ
calculator for parsed CP2K/VASP outputs and the same packaging code applies.
"""

import os
import numpy as np
from ase import units
from ase.build import bulk
from ase.calculators.lj import LennardJones
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

HERE = os.path.dirname(os.path.abspath(__file__))

# LJ parameters loosely scaled to an "Ar-like" liquid; only relative trends matter
SIGMA = 3.4       # Angstrom
EPS = 0.0103      # eV  (~120 K)
RC = 3.0 * SIGMA

# six "temperature points" (like the client's 6-temperature AIMD plan)
TEMPERATURES = [300, 500, 700, 900, 1100, 1300]  # Kelvin


def build_liquid(n_cell=3):
    """~4*n_cell^3 atoms of an fcc solid we will melt into a liquid."""
    atoms = bulk("Ar", "fcc", a=5.26, cubic=True).repeat((n_cell, n_cell, n_cell))
    atoms.calc = LennardJones(sigma=SIGMA, epsilon=EPS, rc=RC)
    return atoms


def rdf(atoms_positions, box, rmax, nbins=120):
    """Simple periodic radial distribution function g(r) from a set of frames."""
    n = atoms_positions[0].shape[0]
    vol = np.prod(box)
    rho = n / vol
    hist = np.zeros(nbins)
    edges = np.linspace(0, rmax, nbins + 1)
    for pos in atoms_positions:
        d = pos[:, None, :] - pos[None, :, :]
        d -= box * np.round(d / box)        # minimum image
        r = np.sqrt((d ** 2).sum(-1))
        r = r[np.triu_indices(n, k=1)]
        h, _ = np.histogram(r, bins=edges)
        hist += h
    hist /= len(atoms_positions)
    shell = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (0.5 * n * rho * shell)
    centers = 0.5 * (edges[1:] + edges[:-1])
    return centers, g


def run_temperature(T, n_equil=1500, n_prod=2500, dt_fs=5.0):
    atoms = build_liquid()
    MaxwellBoltzmannDistribution(atoms, temperature_K=T)
    dyn = Langevin(atoms, dt_fs * units.fs, temperature_K=T, friction=0.02)

    # equilibrate
    dyn.run(n_equil)

    # production: collect frames (positions, energies, forces) for RDF + DeePMD
    frames_pos, energies, forces = [], [], []
    msd_ref = atoms.get_positions().copy()
    box = atoms.get_cell().lengths()
    msd_t, msd_v = [], []

    def sample():
        frames_pos.append(atoms.get_positions().copy())
        energies.append(atoms.get_potential_energy())
        forces.append(atoms.get_forces().copy())
        disp = atoms.get_positions() - msd_ref
        msd_v.append((disp ** 2).sum(1).mean())
        msd_t.append(dyn.get_time() / units.fs)

    dyn.attach(sample, interval=25)
    dyn.run(n_prod)

    # self-diffusion from MSD slope: MSD = 6 D t
    msd_t = np.array(msd_t) - msd_t[0]
    msd_v = np.array(msd_v)
    half = len(msd_t) // 2
    slope = np.polyfit(msd_t[half:], msd_v[half:], 1)[0]  # A^2 / fs
    D = slope / 6.0 * 1e-1  # A^2/fs -> 1e-5 cm^2/s  (1 A^2/fs = 1e-1 cm^2/s ... see note)
    r, g = rdf(frames_pos, box, rmax=RC)

    return dict(T=T, r=r, g=g, D=slope / 6.0, frames=frames_pos,
                energies=np.array(energies), forces=np.array(forces),
                box=box, symbols=atoms.get_chemical_symbols())


def pack_deepmd(results, outdir):
    """Write frames into DeePMD-kit deepmd/npy format using dpdata."""
    import dpdata
    os.makedirs(outdir, exist_ok=True)
    written = []
    for res in results:
        nframes = len(res["frames"])
        natoms = res["frames"][0].shape[0]
        cell = np.zeros((nframes, 3, 3))
        for i in range(3):
            cell[:, i, i] = res["box"][i]
        data = {
            "atom_names": ["Ar"],
            "atom_numbs": [natoms],
            "atom_types": np.zeros(natoms, dtype=int),
            "cells": cell,
            "coords": np.array(res["frames"]),
            "energies": res["energies"],
            "forces": res["forces"],
            "orig": np.zeros(3),
        }
        ls = dpdata.LabeledSystem(data=data, fmt="dict")
        d = os.path.join(outdir, f"T{res['T']}K")
        ls.to("deepmd/npy", d)
        written.append((d, nframes, natoms))
    return written


def main():
    print("Running LJ-liquid MD at 6 temperature points (proxy for molten-salt AIMD)...")
    results = []
    for T in TEMPERATURES:
        res = run_temperature(T)
        results.append(res)
        print(f"  T={T:>5} K   D={res['D']:.4e} A^2/fs   "
              f"first g(r) peak at r={res['r'][np.argmax(res['g'])]:.2f} A")

    # ---- package DeePMD-kit training set (the client's requested deliverable) ----
    outdir = os.path.join(HERE, "deepmd_dataset")
    written = pack_deepmd(results, outdir)
    print(f"\nDeePMD-kit dataset written under: {outdir}")
    for d, nf, na in written:
        rel = os.path.relpath(d, HERE)
        print(f"  {rel}/  ({nf} frames x {na} atoms)  "
              f"-> set.000/{{coord,box,energy,force}}.npy + type.raw")

    # ---- figure: g(r) and Arrhenius-like D(T) ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(results)))
    for res, c in zip(results, cmap):
        ax[0].plot(res["r"], res["g"], color=c, label=f"{res['T']} K")
    ax[0].set_xlabel("r (Å)"); ax[0].set_ylabel("g(r)")
    ax[0].set_title("Radial distribution function vs temperature")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

    Ts = np.array([r["T"] for r in results])
    Ds = np.array([r["D"] for r in results])
    ax[1].semilogy(1000.0 / Ts, Ds, "o-", color="#d1495b")
    ax[1].set_xlabel("1000 / T  (1/K)")
    ax[1].set_ylabel("self-diffusion D  (Å²/fs)")
    ax[1].set_title("Arrhenius plot of diffusion (transport property)")
    ax[1].grid(alpha=0.3, which="both")
    fig.suptitle("Demo #4 · Molten-salt property + DeePMD-packaging workflow\n"
                 "LJ-liquid proxy · ASE NVT MD  (swap LJ→CP2K/VASP AIMD for production)",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    out = os.path.join(HERE, "moltensalt_demo.png")
    fig.savefig(out, dpi=150)
    print(f"\nsaved figure -> {out}")


if __name__ == "__main__":
    main()
