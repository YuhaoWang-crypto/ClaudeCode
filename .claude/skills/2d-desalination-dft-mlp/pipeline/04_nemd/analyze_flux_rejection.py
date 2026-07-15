#!/usr/bin/env python3
"""
Post-processing for pressure-driven (NEMD) desalination trajectories.

Two observables define a reverse-osmosis membrane:

  * Water flux      -- net number of water molecules that cross the membrane
                       plane per unit time per unit area, usually reported in
                       L / (cm^2 . day . MPa) or in molecules/ns.
  * Salt rejection  -- R = 1 - c_permeate / c_feed, the fraction of ions the
                       membrane keeps out of the permeate side.

This module is engine-agnostic: it consumes plain numpy arrays of z-positions
over time (which you can extract from a LAMMPS dump, an ASE trajectory, or a
GROMACS .xtc). Keeping the physics in importable functions lets the toy
validation in 05_toy_validation exercise exactly this code path against a
known ground truth.

Cross-counting convention
-------------------------
A molecule is counted as "permeated" the first time its z goes from the feed
side (z < z_membrane) to the permeate side (z > z_membrane). Back-crossings
are subtracted so `net_crossings` is a true net flux.
"""
from __future__ import annotations
import argparse
import numpy as np


def count_crossings(z_series, z_membrane):
    """Count net feed->permeate crossings for a set of particles over time.

    Parameters
    ----------
    z_series : ndarray, shape (n_frames, n_particles)
        z-coordinate of each particle at each frame.
    z_membrane : float
        z of the membrane plane.

    Returns
    -------
    net : int
        (feed->permeate) minus (permeate->feed) crossing events, summed
        over all particles and frames.
    cumulative : ndarray, shape (n_frames,)
        Running net number of particles on the permeate side over time.
    """
    z_series = np.asarray(z_series, dtype=float)
    # two-valued side: 0 = feed (z <= membrane), 1 = permeate (z > membrane).
    # A boolean (not np.sign) avoids double-counting a particle that lands
    # exactly on the membrane plane, where sign() would give a spurious 0.
    side = (z_series > z_membrane).astype(int)
    dsig = np.diff(side, axis=0)       # +1 => feed->permeate, -1 => back
    forward = (dsig > 0).sum(axis=1)   # per-frame forward crossings
    backward = (dsig < 0).sum(axis=1)  # per-frame backward crossings
    net_per_frame = forward - backward
    cumulative = np.concatenate([[0], np.cumsum(net_per_frame)])
    return int(cumulative[-1]), cumulative


def water_flux(z_water, z_membrane, times_ns, area_A2):
    """Return (net_permeated, flux_molecules_per_ns, cumulative_curve).

    flux is the slope of net permeated water vs time, i.e. steady-state
    throughput. `area_A2` (Angstrom^2) lets callers convert to per-area units.
    """
    net, cumulative = count_crossings(z_water, z_membrane)
    times_ns = np.asarray(times_ns, dtype=float)
    if times_ns[-1] > times_ns[0]:
        # robust slope via least squares over the (assumed) steady portion
        slope = np.polyfit(times_ns, cumulative, 1)[0]
    else:
        slope = 0.0
    flux_per_ns = float(slope)
    flux_per_ns_per_nm2 = flux_per_ns / (area_A2 / 100.0)  # A^2 -> nm^2
    return {
        "net_permeated": net,
        "flux_molecules_per_ns": flux_per_ns,
        "flux_molecules_per_ns_per_nm2": flux_per_ns_per_nm2,
        "cumulative": cumulative,
    }


def salt_rejection(z_ions_final, z_membrane, n_ions_feed_initial):
    """Rejection from final ion positions.

    R = 1 - (ions that reached permeate) / (ions initially in feed).
    Returns R in [0, 1] (clipped).
    """
    z_ions_final = np.asarray(z_ions_final, dtype=float)
    n_permeated = int((z_ions_final > z_membrane).sum())
    if n_ions_feed_initial <= 0:
        return 1.0, 0
    R = 1.0 - n_permeated / n_ions_feed_initial
    return float(np.clip(R, 0.0, 1.0)), n_permeated


def flux_to_permeance(flux_molecules_per_ns, area_A2, pressure_MPa):
    """Convert molecular flux to an experimental permeance
    in L/(cm^2 . day . MPa) so results are comparable to the RO literature.

    1 water molecule = 2.99e-23 L (18 g/mol / 1 g/cm^3 -> 18 cm^3/mol).
    """
    if pressure_MPa <= 0:
        return float("nan")
    vol_per_molecule_L = 18.015 / 6.022e23 / 1000.0  # cm^3 -> L
    molecules_per_day = flux_molecules_per_ns * 1e9 * 86400.0  # /ns -> /day
    vol_L_per_day = molecules_per_day * vol_per_molecule_L
    area_cm2 = area_A2 * 1e-16
    return vol_L_per_day / area_cm2 / pressure_MPa


def summarize(z_water, z_ions, z_membrane, times_ns, area_A2,
              n_ions_feed_initial, pressure_MPa):
    """One-call summary returning a dict of all headline numbers."""
    fx = water_flux(z_water, z_membrane, times_ns, area_A2)
    R, n_ion_perm = salt_rejection(z_ions[-1], z_membrane, n_ions_feed_initial)
    permeance = flux_to_permeance(fx["flux_molecules_per_ns"], area_A2,
                                  pressure_MPa)
    return {
        "net_water_permeated": fx["net_permeated"],
        "water_flux_per_ns": round(fx["flux_molecules_per_ns"], 4),
        "water_flux_per_ns_per_nm2": round(
            fx["flux_molecules_per_ns_per_nm2"], 4),
        "water_permeance_L_cm2_day_MPa": round(permeance, 3),
        "salt_rejection": round(R, 4),
        "ions_permeated": n_ion_perm,
        "_cumulative_water": fx["cumulative"],
    }


def load_lammps_dump(path, water_types, ion_types):
    """Parse a LAMMPS 'custom id type x y z' dump into per-species z-series.

    Returns (z_water, z_ions, times) where z_water has shape (n_frames,
    n_water_atoms). Frame index is used as the time axis (scale by your
    dump interval * timestep to get real ns). Only the z column is retained.

    This is intentionally dependency-free (no MDAnalysis) so the pipeline
    stays lightweight; for large trajectories prefer MDAnalysis/ovito.
    """
    water_types = set(water_types)
    ion_types = set(ion_types)
    frames_w, frames_i, times = [], [], []
    with open(path) as fh:
        lines = fh.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            times.append(int(lines[i + 1]))
            i += 2
            continue
        if lines[i].startswith("ITEM: NUMBER OF ATOMS"):
            natoms = int(lines[i + 1])
            i += 2
            continue
        if lines[i].startswith("ITEM: ATOMS"):
            cols = lines[i].split()[2:]  # e.g. id type x y z
            ci_type, ci_z = cols.index("type"), cols.index("z")
            zw, zi = [], []
            for j in range(i + 1, i + 1 + natoms):
                parts = lines[j].split()
                t = int(parts[ci_type])
                z = float(parts[ci_z])
                if t in water_types:
                    zw.append(z)
                elif t in ion_types:
                    zi.append(z)
            frames_w.append(zw)
            frames_i.append(zi)
            i += 1 + natoms
            continue
        i += 1

    return (np.array(frames_w), np.array(frames_i),
            np.arange(len(frames_w), dtype=float))


def _demo():
    """Tiny self-check on synthetic data (run: python analyze_flux_rejection.py)."""
    rng = np.random.default_rng(0)
    n_frames, n_water = 200, 50
    # 20 waters drift steadily across the membrane, others stay in feed
    z = np.full((n_frames, n_water), -5.0)
    for i in range(20):
        start = rng.integers(0, n_frames - 20)
        z[start:, i] = np.linspace(-5.0, 5.0, n_frames - start)
    z_ions = np.full((n_frames, 10), -5.0)  # all ions rejected
    times = np.linspace(0, 2.0, n_frames)
    out = summarize(z, z_ions, 0.0, times, area_A2=150.0,
                    n_ions_feed_initial=10, pressure_MPa=100.0)
    out.pop("_cumulative_water")
    print("Self-check summary:", out)
    assert out["net_water_permeated"] == 20, out
    assert out["salt_rejection"] == 1.0, out
    print("OK: analyzer recovers 20 permeated waters and 100% rejection.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--demo", action="store_true",
                    help="Run a built-in synthetic self-check")
    args = ap.parse_args()
    if args.demo:
        _demo()
    else:
        print("Import this module, or run with --demo. "
              "See 05_toy_validation for a full pipeline test.")
