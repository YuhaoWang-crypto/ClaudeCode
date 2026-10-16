"""Coarse-grained Brownian dynamics of a surface-tethered DNA20-MB probe (OpenMM, CPU).

Model (all values are simplifying assumptions, not sequence-specific measurements):
  anchor (fixed, z = 0 = SAM top) -0.3- C6 bead -0.3- DNA1 -0.6- ... -0.6- DNA20 -0.3- L1 -0.3- L2 -0.3- MB
  harmonic bonds (100 kT/nm^2), cosine bending U = kappa (1 + cos theta) with kappa = l_p / b,
  excluded volume for pairs >= 3 bonds apart and a repulsive plane at z = 0, both as truncated
  harmonic repulsion U = (k_rep/2)(sigma - r)^2 for r < sigma, k_rep = 1000 kT/nm^2, sigma = sum of
  radii (pairs) or bead radius (wall).  (Deviation from the protocol's WCA: WCA on 0.15 nm linker
  beads makes explicit Brownian steps unstable for every allowed dt >= 0.1 ps.)
  Radii: DNA 0.25 nm, MB 0.35 nm, linker/C6 beads 0.15 nm (assumed).
  No surface attraction, no base pairing/stacking, no electrostatics, no hydrodynamic interaction
  (free-draining Stokes friction 6 pi eta a, water at 298.15 K, no near-wall correction).
Time is coarse-grained Brownian time with the stated friction; it is NOT experimentally validated.

Output per run (npz), arrays shaped (frames, replicas): t_ps, z_mb, rg_dna, r_end (anchor-MB), n_contact (MB to DNA beads 1..18
within r_MB + r_DNA + 0.3 nm), and every 100th frame of coordinates for audit.
"""
import argparse
import json
import os
import time

import numpy as np
import openmm as mm
import openmm.unit as u

KB = 0.008314462618  # kJ/mol/K
AMU_KG = 1.66053906660e-27


def build_topology(n_dna=20, n_link=3, b_dna=0.6, b_link=0.3, n_c6=2,
                   r_dna=0.25, r_mb=0.35, r_link=0.15, lp_dna=1.0, lp_link=1.0):
    kinds, radii, bonds, blen = [], [], [], []
    kinds.append("anchor"); radii.append(r_link)
    for _ in range(n_c6 - 1):
        kinds.append("c6"); radii.append(r_link)
    for _ in range(n_dna):
        kinds.append("dna"); radii.append(r_dna)
    for _ in range(n_link - 1):
        kinds.append("link"); radii.append(r_link)
    kinds.append("mb"); radii.append(r_mb)
    n = len(kinds)
    for i in range(n - 1):
        a, b = kinds[i], kinds[i + 1]
        bonds.append((i, i + 1))
        blen.append(b_dna if (a == "dna" and b == "dna") else b_link)
    angles = []
    for i in range(1, n - 1):
        centre_dna = kinds[i] == "dna" and kinds[i - 1] == "dna" and kinds[i + 1] == "dna"
        kappa = lp_dna / b_dna if centre_dna else lp_link / b_link
        angles.append((i - 1, i, i + 1, kappa))
    return dict(kinds=kinds, radii=np.array(radii), bonds=bonds, blen=np.array(blen), angles=angles,
                i_mb=n - 1, dna_idx=[i for i, k in enumerate(kinds) if k == "dna"])


def make_system(top, n_rep=1, T=298.15, eta_pa_s=0.890e-3, k_bond_kT=100.0, gamma_ps=1.0, k_rep_kT=1000.0):
    """n_rep independent, non-interacting copies of the probe in one System (no cross-replica forces);
    each copy is an independent Brownian trajectory with its own noise."""
    kT = KB * T
    n = len(top["kinds"])
    sysm = mm.System()
    for r in range(n_rep):
        for i, (k, a) in enumerate(zip(top["kinds"], top["radii"])):
            if k == "anchor":
                sysm.addParticle(0.0)            # massless = fixed in OpenMM
            else:
                zeta = 6 * np.pi * eta_pa_s * a * 1e-9           # kg/s
                m_amu = zeta / (gamma_ps * 1e12 * AMU_KG)        # m*gamma = zeta
                sysm.addParticle(m_amu)
    hb = mm.HarmonicBondForce()
    ang = mm.CustomAngleForce("kap*(1+cos(theta))")
    ang.addPerAngleParameter("kap")
    for r in range(n_rep):
        o = r * n
        for (i, j), l0 in zip(top["bonds"], top["blen"]):
            hb.addBond(o + i, o + j, float(l0), float(k_bond_kT * kT))
        for i, j, k, kap in top["angles"]:
            ang.addAngle(o + i, o + j, o + k, [float(kap * kT)])
    sysm.addForce(hb)
    sysm.addForce(ang)
    wca = mm.CustomBondForce("step(sig-r)*0.5*krep*(sig-r)^2")
    wca.addGlobalParameter("krep", k_rep_kT * kT)
    wca.addPerBondParameter("sig")
    bonded = {}
    for i, j in top["bonds"]:
        bonded.setdefault(i, set()).add(j)
        bonded.setdefault(j, set()).add(i)
    for i in range(n):
        for j in range(i + 1, n):
            if j in bonded.get(i, ()) or any(j in bonded.get(m, ()) for m in bonded.get(i, ())):
                continue
            if top["kinds"][i] == "anchor" or top["kinds"][j] == "anchor":
                continue
            for r in range(n_rep):
                wca.addBond(r * n + i, r * n + j, [float(top["radii"][i] + top["radii"][j])])
    sysm.addForce(wca)
    wall = mm.CustomExternalForce("step(s-z)*0.5*krep*(s-z)^2")
    wall.addGlobalParameter("krep", k_rep_kT * kT)
    wall.addPerParticleParameter("s")
    for r in range(n_rep):
        for i, (k, a) in enumerate(zip(top["kinds"], top["radii"])):
            if k != "anchor":
                wall.addParticle(r * n + i, [float(a)])
    sysm.addForce(wall)
    return sysm


def initial_positions(top, mode, rng):
    n = len(top["kinds"])
    pos = np.zeros((n, 3))
    if mode == "extended":
        for i in range(1, n):
            pos[i] = pos[i - 1] + np.array([0, 0, top["blen"][i - 1]])
        pos[:, :2] += rng.normal(0, 0.01, (n, 2))
        pos[0] = 0
    else:   # collapsed: self-avoiding-ish random walk kept close to the surface
        for i in range(1, n):
            for _ in range(1000):
                d = rng.normal(size=3)
                d /= np.linalg.norm(d)
                p = pos[i - 1] + top["blen"][i - 1] * d
                ok = p[2] > top["radii"][i] * 1.2 and p[2] < 2.0
                ok = ok and all(np.linalg.norm(p - pos[j]) > 0.9 * (top["radii"][i] + top["radii"][j])
                                for j in range(i - 1))
                if ok:
                    pos[i] = p
                    break
            else:
                pos[i] = pos[i - 1] + np.array([0, 0, top["blen"][i - 1]])
    return pos


def descriptors(xyz, top):
    imb = top["i_mb"]
    dna = top["dna_idx"]
    z_mb = xyz[:, imb, 2]
    d = xyz[:, dna, :]
    rg = np.sqrt(((d - d.mean(1, keepdims=True)) ** 2).sum(2).mean(1))
    r_end = np.linalg.norm(xyz[:, imb, :] - xyz[:, 0, :], axis=1)
    cut = top["radii"][imb] + top["radii"][dna[0]] + 0.3
    dd = np.linalg.norm(xyz[:, dna[:-2], :] - xyz[:, imb:imb + 1, :], axis=2)
    ncon = (dd < cut).sum(1)
    return z_mb, rg, r_end, ncon


def run(out, seed, mode, dt_ps, t_ns, t_eq_ns, save_ps, T=298.15, threads=1, n_rep=1, **topkw):
    """mode: 'extended', 'collapsed' or 'mixed' (alternating start per replica)."""
    rng = np.random.default_rng(seed)
    top = build_topology(**topkw)
    n = len(top["kinds"])
    sysm = make_system(top, n_rep=n_rep, T=T)
    integ = mm.BrownianIntegrator(T * u.kelvin, 1.0 / u.picosecond, dt_ps * u.picosecond)
    integ.setRandomNumberSeed(int(seed) * 7919 + 13)
    plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(sysm, integ, plat, {"Threads": str(threads)})
    modes = [mode if mode != "mixed" else ("extended" if r % 2 == 0 else "collapsed") for r in range(n_rep)]
    ctx.setPositions(np.vstack([initial_positions(top, m, rng) for m in modes]) * u.nanometer)
    mm.LocalEnergyMinimizer.minimize(ctx, 1.0, 200)
    nsave = max(1, int(round(save_ps / dt_ps)))
    t0 = time.time()
    integ.step(int(round(t_eq_ns * 1000 / dt_ps)))
    nfr = int(round(t_ns * 1000 / save_ps))
    desc = {k: np.empty((nfr, n_rep), np.float32) for k in ("z", "rg", "r_end", "nc")}
    xyz_sub = []
    chunk = 500
    buf = np.empty((chunk, n_rep * n, 3), np.float32)
    k0 = 0
    for k in range(nfr):
        integ.step(nsave)
        buf[k - k0] = ctx.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)
        if k - k0 == chunk - 1 or k == nfr - 1:
            m = k - k0 + 1
            x = buf[:m].reshape(m, n_rep, n, 3)
            for r in range(n_rep):
                z_mb, rg, r_end, ncon = descriptors(x[:, r], top)
                desc["z"][k0:k + 1, r] = z_mb
                desc["rg"][k0:k + 1, r] = rg
                desc["r_end"][k0:k + 1, r] = r_end
                desc["nc"][k0:k + 1, r] = ncon
            xyz_sub.append(x[::100, : min(n_rep, 4)].copy())
            k0 = k + 1
    wall = time.time() - t0
    meta = dict(seed=seed, mode=mode, n_rep=n_rep, dt_ps=dt_ps, t_ns=t_ns, t_eq_ns=t_eq_ns, save_ps=save_ps,
                T=T, wall_s=wall, topology={k: v for k, v in topkw.items()},
                n_particles=n, radii=top["radii"].tolist(), kinds=top["kinds"], replica_modes=modes)
    np.savez_compressed(out, t_ps=np.arange(1, nfr + 1) * save_ps, z_mb=desc["z"], rg_dna=desc["rg"],
                        r_end=desc["r_end"], n_contact=desc["nc"].astype(np.int16),
                        xyz_sub=np.concatenate(xyz_sub), meta=json.dumps(meta))
    return meta


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mode", default="mixed", choices=["extended", "collapsed", "mixed"])
    ap.add_argument("--n_rep", type=int, default=1)
    ap.add_argument("--dt", type=float, default=0.5)
    ap.add_argument("--t_ns", type=float, default=100.0)
    ap.add_argument("--t_eq_ns", type=float, default=20.0)
    ap.add_argument("--save_ps", type=float, default=10.0)
    ap.add_argument("--n_link", type=int, default=3)
    ap.add_argument("--lp_link", type=float, default=1.0)
    ap.add_argument("--n_dna", type=int, default=20, help="number of nucleotides (1 bead each)")
    ap.add_argument("--r_reporter", type=float, default=0.35, help="reporter radius (nm): MB ~0.35, ferrocene ~0.3")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    m = run(a.out, a.seed, a.mode, a.dt, a.t_ns, a.t_eq_ns, a.save_ps, n_rep=a.n_rep,
            n_link=a.n_link, lp_link=a.lp_link, n_dna=a.n_dna, r_mb=a.r_reporter)
    print(json.dumps({k: m[k] for k in ("seed", "mode", "dt_ps", "t_ns", "wall_s")}))
