"""
E7 — interface structure: beyond the bulk RDF, look along the surface normal   [digest p9]

Model (stated explicitly, digest "把界面模型写清楚"):
  electrode   3-layer AB graphite basal plane, orthorhombic cell (a=2.46 Å), atoms FIXED
              (mass 0), uncharged, LJ only (OPLS aromatic C: σ=3.55 Å, ε=0.07 kcal/mol),
              i.e. an UNCHARGED, non-polarisable surface (cf. Fig.14a which is also uncharged);
  electrolyte LiTFSI/DME at solvent:Li = 10 (same FF/charges as E1), filled at the E1 bulk
              density in a slab; periodic in x,y,z so the liquid sees the graphite top face
              on one side and the bottom face on the other (two equivalent interfaces);
  ensemble    NVT (a barostat would move the fixed slab); the liquid density in the
              centre of the slab is checked against the E1 bulk value and reported;
  coordinate  z measured from the plane of the top carbon layer (surface-aligned, digest A).

Outputs: number-density profiles ρ_s(z) for Li, O(DME), N(TFSI) and DME centre of mass,
DME orientation P(cos θ) near vs far from the surface (θ between the O–O vector and ẑ),
and W_Li(z) = −kT ln(ρ_Li(z)/ρ_Li,bulk) as the free-energy view of the same profile.
"""
from __future__ import annotations

import json
import os
import time

import numpy as np

from electrolyte_pipeline.e0_systems import WORK, FIG, SPECIES, ForceFieldSpec, DynamicsSpec

LABEL = "IF_LiTFSI_r10"
A_CC = 1.42                       # Å
NX, NY, NLAYERS = 14, 8, 3        # 14 * 2.46 = 34.44 Å ; 8 * 4.26 = 34.09 Å
D_LAYER = 3.35
GAP = 3.0                         # Å between top carbon plane and packmol region
LIQ_Z = 40.0                      # Å of liquid
C_SIGMA_NM, C_EPS_KJ = 0.355, 0.07 * 4.184
COUNTS = {"DME": 251, "Li": 25, "TFSI": 25}   # r10 composition. Fill calibrated against the E1 bulk density:
                                              # 236/24/24 -> centre 0.921 g/cm3 (-9 %), 258/26/26 -> 1.062 (+4.5 %),
                                              # interpolated to 251/25/25 for the production run.


def graphite_xyz():
    a = np.sqrt(3) * A_CC                          # 2.46
    b = 3 * A_CC                                   # 4.26
    basis = np.array([[0, 0], [a / 2, b / 6], [a / 2, b / 2], [0, 2 * b / 3]])   # rectangular 4-atom cell
    shiftB = np.array([a / 2, b / 6])              # AB stacking shift
    pts = []
    for L in range(NLAYERS):
        z = L * D_LAYER
        sh = shiftB * (L % 2)
        for i in range(NX):
            for j in range(NY):
                for bx, by in basis:
                    pts.append([(bx + i * a + sh[0]) % (NX * a), (by + j * b + sh[1]) % (NY * b), z])
    return np.array(pts), NX * a, NY * b


def build(outdir: str):
    from openff.toolkit import ForceField
    from openff.interchange import Interchange
    from openff.interchange.components._packmol import pack_box
    from openff.units import unit
    import openmm
    from openmm import app
    from electrolyte_pipeline.e1_build import _species_molecule
    os.makedirs(outdir, exist_ok=True)
    ff, dyn = ForceFieldSpec(), DynamicsSpec()
    carb, Lx, Ly = graphite_xyz()
    z_top = (NLAYERS - 1) * D_LAYER
    names = list(COUNTS)
    mols = [_species_molecule(n, ff) for n in names]
    box = np.diag([Lx, Ly, LIQ_Z]) * unit.angstrom
    top = pack_box(molecules=mols, number_of_copies=[COUNTS[n] for n in names], box_vectors=box,
                   tolerance=2.0 * unit.angstrom)
    Lz = z_top + GAP + LIQ_Z + GAP + D_LAYER * 0 + 0.0     # slab (0..z_top) + gap + liquid + gap, then periodic image of slab bottom
    Lz = z_top + 2 * GAP + LIQ_Z
    top.box_vectors = np.diag([Lx, Ly, Lz]) * unit.angstrom
    sage = ForceField(ff.name)
    ic = Interchange.from_smirnoff(sage, top, charge_from_molecules=mols, allow_nonintegral_charges=True)
    system = ic.to_openmm(combine_nonbonded_forces=True)
    omm_top = ic.to_openmm_topology()
    pos = ic.positions.to_openmm().value_in_unit(openmm.unit.angstrom)
    pos = np.array(pos) + np.array([0, 0, z_top + GAP])    # lift the liquid above the slab
    nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
    nb.setNonbondedMethod(openmm.NonbondedForce.PME); nb.setCutoffDistance(ff.cutoff_nm * openmm.unit.nanometer)
    nb.setUseSwitchingFunction(True); nb.setSwitchingDistance(ff.switch_nm * openmm.unit.nanometer)
    nb.setUseDispersionCorrection(False)                     # inhomogeneous system
    # add fixed carbon atoms
    chain = omm_top.addChain("G")
    n_liq = system.getNumParticles()
    per_layer = len(carb) // NLAYERS
    for k in range(len(carb)):
        if k % per_layer == 0:                               # one residue per layer: atom names stay unique (<=4 chars)
            res = omm_top.addResidue("GRA", chain)
        system.addParticle(0.0)                              # mass 0 -> fixed
        nb.addParticle(0.0, C_SIGMA_NM, C_EPS_KJ)
        omm_top.addAtom(f"C{k % per_layer}", app.element.carbon, res)
    allpos = np.vstack([pos, carb]) * openmm.unit.angstrom
    omm_top.setPeriodicBoxVectors(np.diag([Lx, Ly, Lz]) * openmm.unit.angstrom)
    system.setDefaultPeriodicBoxVectors(*(np.diag([Lx, Ly, Lz]) * openmm.unit.angstrom))
    with open(os.path.join(outdir, "system.xml"), "w") as fh:
        fh.write(openmm.XmlSerializer.serialize(system))
    with open(os.path.join(outdir, "initial.pdb"), "w") as fh:
        app.PDBFile.writeFile(omm_top, allpos, fh)
    atoms = []
    res_species = []
    for n in names:
        res_species += [n] * COUNTS[n]
    for r in omm_top.residues():
        for a in r.atoms():
            atoms.append({"index": a.index, "element": a.element.symbol, "molecule": r.index,
                          "species": res_species[r.index] if r.index < len(res_species) else "GRA"})
    json.dump(atoms, open(os.path.join(outdir, "atoms.json"), "w"))
    rec = {"composition": {"label": LABEL, "counts": {**COUNTS, "C_graphite": int(len(carb))}},
           "electrode": {"type": "graphite basal plane, AB, 3 layers, fixed, uncharged",
                         "cell_A": [Lx, Ly, Lz], "z_top_A": z_top, "gap_A": GAP,
                         "C_LJ": {"sigma_nm": C_SIGMA_NM, "eps_kJ_mol": C_EPS_KJ}},
           "interactions": {**ff.__dict__, "dispersion_correction": False},
           "dynamics": {**dyn.__dict__, "ensemble": "NVT"},
           "species": {s: SPECIES[s] for s in COUNTS}}
    json.dump(rec, open(os.path.join(outdir, "system_record.json"), "w"), indent=2)
    print(f"[E7 build] {system.getNumParticles()} particles ({len(carb)} fixed C), cell {Lx:.2f} x {Ly:.2f} x {Lz:.2f} Å")
    return outdir


def run(workdir: str = WORK, equil_ns: float = 1.0, prod_ns: float = 4.0) -> dict:
    import openmm
    from openmm import app, unit
    from electrolyte_pipeline.e1_run import _platform, _reporter, convergence_check
    dyn = DynamicsSpec()
    d = os.path.join(workdir, LABEL)
    build(d)
    system = openmm.XmlSerializer.deserialize(open(os.path.join(d, "system.xml")).read())
    pdb = app.PDBFile(os.path.join(d, "initial.pdb"))
    integ = openmm.LangevinMiddleIntegrator(dyn.temperature_K * unit.kelvin, dyn.friction_per_ps / unit.picosecond,
                                            dyn.timestep_fs * unit.femtosecond)
    plat, props = _platform()
    sim = app.Simulation(pdb.topology, system, integ, plat, props)
    sim.context.setPositions(pdb.positions)
    sim.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
    t0 = time.time()
    sim.minimizeEnergy(tolerance=10 * unit.kilojoule_per_mole / unit.nanometer)
    sim.context.setVelocitiesToTemperature(dyn.temperature_K * unit.kelvin, 3)
    spp = int(1000 / dyn.timestep_fs)
    sim.reporters = [_reporter(os.path.join(d, "equil.csv"), int(dyn.thermo_ps * spp))]
    sim.step(int(equil_ns * 1000 * spp))
    sim.reporters = [_reporter(os.path.join(d, "prod.csv"), int(dyn.thermo_ps * spp)),
                     app.DCDReporter(os.path.join(d, "prod.dcd"), int(dyn.frame_ps * spp), enforcePeriodicBox=True)]
    sim.step(int(prod_ns * 1000 * spp))
    st = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    with open(os.path.join(d, "final.pdb"), "w") as fh:
        app.PDBFile.writeFile(sim.topology, st.getPositions(), fh)
    conv = convergence_check(os.path.join(d, "equil.csv"))
    rec = {"label": LABEL, "platform": plat.getName(), "equil_ns": equil_ns, "prod_ns": prod_ns,
           "frame_ps": dyn.frame_ps, "equil_convergence": conv, "wall_s": round(time.time() - t0, 1)}
    json.dump(rec, open(os.path.join(d, "run_record.json"), "w"), indent=2)
    print(f"[E7 run] {plat.getName()}: {equil_ns}+{prod_ns} ns, wall {rec['wall_s']} s, "
          f"Epot drift {conv['epot']['drift']:.0f} kJ/mol")
    return rec


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------
def analyse(workdir: str = WORK, dz: float = 0.2, bulk_ref: str = "LiTFSI_r10") -> dict:
    from electrolyte_pipeline.traj import Traj
    tr = Traj(LABEL, workdir)
    Lz = tr.record["electrode"]["cell_A"][2]
    z_top = tr.record["electrode"]["z_top_A"]
    area = tr.record["electrode"]["cell_A"][0] * tr.record["electrode"]["cell_A"][1]
    groups = {"Li": tr.sel(species="Li"), "O(DME)": tr.sel(species="DME", element="O"),
              "N(TFSI)": tr.sel(species="TFSI", element="N"), "F(TFSI)": tr.sel(species="TFSI", element="F")}
    dme_idx = tr.sel(species="DME"); dme_mol = tr.molid[dme_idx]
    dme_O = tr.sel(species="DME", element="O")
    edges = np.arange(0, Lz - z_top + dz, dz)                     # z' = z - z_top  (0 at the top carbon plane)
    zc = 0.5 * (edges[1:] + edges[:-1])
    hist = {k: np.zeros(len(zc)) for k in groups}
    hist["DME(COM)"] = np.zeros(len(zc))
    cos_near, cos_far = [], []
    nfr = 0
    for ts in tr.frames():
        pos = tr.u.atoms.positions.astype(float)
        zz = (pos[:, 2] - z_top) % Lz
        for k, idx in groups.items():
            hist[k] += np.histogram(zz[idx], bins=edges)[0]
        # DME COM and orientation (O–O vector vs z)
        m = tr.mass[dme_idx]
        um, inv = np.unique(dme_mol, return_inverse=True)
        p = pos[dme_idx].copy()
        # make molecules whole in z (they are whole already in x,y with enforcePeriodicBox=True? not guaranteed) -> unwrap within molecule
        first = np.zeros((len(um), 3))
        for i, mol in enumerate(um):
            sel = np.where(inv == i)[0]
            ref = p[sel[0]]
            dd = p[sel] - ref
            dd -= np.array(tr.box(ts)) * np.round(dd / np.array(tr.box(ts)))
            p[sel] = ref + dd
        com = np.zeros((len(um), 3)); np.add.at(com, inv, p * m[:, None]); com /= np.bincount(inv, weights=m)[:, None]
        zcom = (com[:, 2] - z_top) % Lz
        hist["DME(COM)"] += np.histogram(zcom, bins=edges)[0]
        oo = p[np.isin(dme_idx, dme_O)].reshape(-1, 2, 3)
        v = oo[:, 1] - oo[:, 0]; v /= np.linalg.norm(v, axis=1)[:, None]
        c = np.abs(v[:, 2])
        near = (zcom < 6.0) | (zcom > Lz - z_top - 6.0)
        cos_near += c[near].tolist(); cos_far += c[(zcom > 15) & (zcom < Lz - z_top - 15)].tolist()
        nfr += 1
    prof = {k: h / (nfr * area * dz) for k, h in hist.items()}      # number / Å³
    # bulk reference = centre of the liquid slab
    centre = (zc > 15) & (zc < Lz - z_top - 15)
    bulk = {k: float(p[centre].mean()) for k, p in prof.items()}
    e1 = os.path.join(workdir, bulk_ref, "e4_properties.json")
    rho_bulk_e1 = json.load(open(e1))["density"]["mean"] if os.path.exists(e1) else None
    # mass density in the centre
    mass_hist = np.zeros(len(zc))
    for ts in tr.frames(stride=5):
        zz = (tr.u.atoms.positions[:, 2].astype(float) - z_top) % Lz
        liquid = tr.species != "GRA"
        mass_hist += np.histogram(zz[liquid], bins=edges, weights=tr.mass[liquid])[0]
    rho_mass = mass_hist / (len(list(tr.frames(stride=5))) * area * dz) * 1.66054   # amu/Å³ -> g/cm³
    kT = 0.0019872041 * 298.15
    with np.errstate(divide="ignore"):
        W_li = -kT * np.log(prof["Li"] / bulk["Li"])
    first_peaks = {k: float(zc[np.argmax(p * (zc < 8))]) for k, p in prof.items()}
    res = {"label": LABEL, "z_reference": "top carbon plane (z'=0); z' increases into the liquid",
           "n_frames": nfr, "dz_A": dz, "z": zc.tolist(), "profiles_per_A3": {k: p.tolist() for k, p in prof.items()},
           "bulk_centre_density_per_A3": bulk, "first_peak_z_A": first_peaks,
           "centre_mass_density_g_cm3": float(rho_mass[centre].mean()), "E1_bulk_density_g_cm3": rho_bulk_e1,
           "W_Li_kcal": np.where(np.isfinite(W_li), W_li, None).tolist(),
           "DME_orientation": {"mean_|cos|_near_surface(<6A)": float(np.mean(cos_near)),
                               "mean_|cos|_bulk": float(np.mean(cos_far)),
                               "note": "|cos θ| of the O–O vector vs surface normal; 0.5 = isotropic"},
           "note": "⚠️ uncharged, non-polarisable, rigid graphite; no potential control; one composition"}
    json.dump(res, open(os.path.join(tr.dir, "e7_interface.json"), "w"))
    _plot(res)
    return res


def _plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    z = np.array(res["z"])
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for k, p in res["profiles_per_A3"].items():
        p = np.array(p) * 1000                     # nm^-3
        ax[0].plot(z, p, label=f"{k} (1st peak {res['first_peak_z_A'][k]:.1f} Å)")
    ax[0].set_xlim(0, 22); ax[0].set_xlabel("z' from top C plane (Å)"); ax[0].set_ylabel("number density (nm⁻³)")
    ax[0].legend(fontsize=7); ax[0].set_title("A  layering at an uncharged graphite face", fontsize=9)
    for k in ("Li", "N(TFSI)"):
        p = np.array(res["profiles_per_A3"][k]) / res["bulk_centre_density_per_A3"][k]
        ax[1].plot(z, p, label=k)
    ax[1].axhline(1, color="gray", lw=0.5); ax[1].set_xlim(0, 22); ax[1].set_xlabel("z' (Å)"); ax[1].set_ylabel("ρ(z)/ρ_bulk")
    ax[1].legend(fontsize=8); ax[1].set_title("B  ions relative to the slab centre", fontsize=9)
    W = np.array([np.nan if w is None else w for w in res["W_Li_kcal"]])
    ax[2].plot(z, W, "k-"); ax[2].set_xlim(0, 22); ax[2].set_ylim(-2, 3)
    ax[2].set_xlabel("z' (Å)"); ax[2].set_ylabel("W_Li(z) = −kT ln ρ/ρ_bulk (kcal/mol)")
    o = res["DME_orientation"]
    ax[2].set_title(f"D  Li PMF along z; DME |cosθ| near {o['mean_|cos|_near_surface(<6A)']:.2f} vs bulk {o['mean_|cos|_bulk']:.2f}", fontsize=9)
    fig.suptitle(f"E7  graphite/LiTFSI-DME interface (NVT, rigid uncharged slab; centre ρ={res['centre_mass_density_g_cm3']:.3f} g/cm³ vs bulk E1 {res['E1_bulk_density_g_cm3']})", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "e7_interface.png"), dpi=140); plt.close(fig)


def report(workdir: str = WORK):
    if not os.path.exists(os.path.join(workdir, LABEL, "prod.dcd")):
        print("E7  interface: no trajectory (run `modal run electrolyte_pipeline/modal_run.py --stage interface`)")
        return None
    r = analyse(workdir)
    print("E7  graphite(0001)/LiTFSI–DME interface, z' from the top carbon plane")
    for k, v in r["first_peak_z_A"].items():
        print(f"    first peak {k:9s} at z' = {v:.1f} Å")
    print(f"    centre mass density {r['centre_mass_density_g_cm3']:.3f} g/cm³ (E1 bulk {r['E1_bulk_density_g_cm3']})")
    o = r["DME_orientation"]
    print(f"    DME O–O |cosθ| vs normal: near surface {o['mean_|cos|_near_surface(<6A)']:.2f}, bulk {o['mean_|cos|_bulk']:.2f} (0.5 isotropic)")
    print(f"    {r['note']}")
    return r


if __name__ == "__main__":
    import sys
    if "--run" in sys.argv:
        run(equil_ns=0.01, prod_ns=0.01)
    report()
