"""
2D FDTD (TM: E_z, H_x, H_y) on a Yee grid — MEEP-equivalent, pure NumPy.

Used to demonstrate backscattering-immune transport along a topological domain
wall of the real silicon-rod Wu-Hu crystal: a wave launched on the wall follows
a sharp bend and passes a defect with little reflection, unlike a trivial guide.

Units: c = 1, a = 1, mu = 1; frequencies in c/a (= omega a / 2 pi c). A graded
conductivity "sponge" absorbs at the four borders (adequate for a demo; a true
UPML would be the production upgrade).
"""
from __future__ import annotations
import numpy as np

A = 1.0
SQ3 = np.sqrt(3.0)
EPS_ROD, EPS_BG, R_ROD = 11.7, 1.0, 0.12 * A
R_TOPO, R_TRIV = 0.36 * A, 0.30 * A


def _cluster(cx, cy, R):
    return [(cx + R * np.cos(np.pi / 3 * k), cy + R * np.sin(np.pi / 3 * k))
            for k in range(6)]


def build_domain(Ncx, Ncy, interface, Nx=16, drop_rod=None):
    """
    Rectangular domain of Ncx x Ncy orthorhombic cells (each a x sqrt3 a, two
    clusters). `interface(cx, cy) -> True` marks the TOPOLOGICAL side.
    `drop_rod=(x,y)` optionally removes the nearest rod (a point defect).
    Returns eps (Nxg, Nyg), h.
    """
    h = A / Nx
    cell_h = SQ3 * A
    Nyc = int(round(cell_h / h))
    Nxg, Nyg = Ncx * Nx, Ncy * Nyc
    Lx, Ly = Nxg * h, Nyg * h
    xs = (np.arange(Nxg) + 0.5) * h
    ys = (np.arange(Nyg) + 0.5) * h
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    eps = np.full((Nxg, Nyg), EPS_BG)

    rods = []
    for ci in range(Ncx):
        for cj in range(Ncy):
            for (ox, oy) in ((0.0, 0.0), (A / 2, cell_h / 2)):
                cx, cy = ci * A + ox, cj * cell_h + oy
                R = R_TOPO if interface(cx, cy) else R_TRIV
                rods += _cluster(cx, cy, R)
    if drop_rod is not None:
        dx0, dy0 = drop_rod
        rods.sort(key=lambda p: (p[0] - dx0) ** 2 + (p[1] - dy0) ** 2)
        rods = rods[1:]  # remove nearest rod
    for (rx, ry) in rods:
        eps[(X - rx) ** 2 + (Y - ry) ** 2 < R_ROD ** 2] = EPS_ROD
    return eps, h, Lx, Ly


def _sponge(Nxg, Nyg, width, smax=2.0):
    """Graded absorbing conductivity at the four borders."""
    sx = np.zeros(Nxg)
    sy = np.zeros(Nyg)
    r = np.arange(width)
    ramp = smax * ((width - r) / width) ** 3
    sx[:width] += ramp[::-1]
    sx[-width:] += ramp
    sy[:width] += ramp[::-1]
    sy[-width:] += ramp
    return sx[:, None] + sy[None, :]


def run_fdtd(eps, h, freq, nsteps, src_xy, mon_xy=None, mon_lines=None,
             snapshots=None, pulse=False, sponge_w=20):
    """
    Leapfrog TM FDTD. Returns dict with optional monitor time series and snapshots.

    src_xy, mon_xy : (x, y) in units of a. freq in c/a.
    pulse=False -> CW sine (for movies); True -> Gaussian pulse (for spectra).
    """
    Nxg, Nyg = eps.shape
    dt = 0.5 * h  # CFL for 2D
    Ez = np.zeros((Nxg, Nyg))
    Hx = np.zeros((Nxg, Nyg))
    Hy = np.zeros((Nxg, Nyg))
    sig = _sponge(Nxg, Nyg, sponge_w)
    damp = np.exp(-sig * dt)

    si = int(src_xy[0] / h); sj = int(src_xy[1] / h)
    mons = {}
    if mon_xy:
        for name, (mx, my) in mon_xy.items():
            mons[name] = (int(mx / h), int(my / h))
    rec = {k: np.zeros(nsteps) for k in mons}
    lines = {}
    if mon_lines:
        for name, mx in mon_lines.items():
            lines[name] = int(mx / h)
    reclines = {k: np.zeros((nsteps, Nyg)) for k in lines}
    frames = []
    t0, tw = 30 * dt / dt * dt, 12.0  # pulse centre/width (in time units)
    w = 2 * np.pi * freq

    for n in range(nsteps):
        t = n * dt
        # H update
        Hx[:, :-1] -= (dt / h) * (Ez[:, 1:] - Ez[:, :-1])
        Hy[:-1, :] += (dt / h) * (Ez[1:, :] - Ez[:-1, :])
        Hx *= damp; Hy *= damp
        # E update
        curl = np.zeros_like(Ez)
        curl[1:, 1:] = (Hy[1:, 1:] - Hy[:-1, 1:]) / h - (Hx[1:, 1:] - Hx[1:, :-1]) / h
        Ez += dt / eps * curl
        Ez *= damp
        # soft source
        if pulse:
            amp = np.exp(-((t - 3 * tw) ** 2) / (2 * (tw / 3) ** 2)) * np.sin(w * t)
        else:
            ramp = min(1.0, t / (20.0))
            amp = ramp * np.sin(w * t)
        Ez[si, sj] += amp
        for k, (mx, my) in mons.items():
            rec[k][n] = Ez[mx, my]
        for k, xi in lines.items():
            reclines[k][n] = Ez[xi, :]
        if snapshots and n in snapshots:
            frames.append((n, Ez.copy()))
    return {"monitors": rec, "lines": reclines, "frames": frames,
            "dt": dt, "eps": eps, "h": h}
