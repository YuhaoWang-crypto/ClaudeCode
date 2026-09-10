"""How fast does the report's MD actually run on THIS box? Measure, don't guess.

Chain breaks in the 7CBC crystal are split into separate chains so each gets
its own charged termini; that changes the electrostatics slightly and nothing
about the timing, which is all this script is for.
"""
import time, string
import numpy as np
import openmm as mm, openmm.app as app
from openmm import unit

PDB = "/home/user/ClaudeCode/sparkcage/structure/pdb/7cbc.pdb"
AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())

# heavy atoms each residue must have; anything short of this is a disordered
# side chain in the crystal and gets truncated to alanine, which is what
# pdbfixer would do and which does not change a timing measurement
HEAVY = {
    "ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2",
    "ASP": "CB CG OD1 OD2", "CYS": "CB SG", "GLN": "CB CG CD OE1 NE2",
    "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ",
    "MET": "CB CG SD CE", "PHE": "CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "CB CG CD",
    "SER": "CB OG", "THR": "CB OG1 CG2", "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2",
}


def segmented_pdb(lo, hi, out):
    res, order = {}, []
    for ln in open(PDB):
        if not ln.startswith("ATOM") or ln[21] != "A" or ln[16] not in " A":
            continue
        if ln[17:20].strip() not in AA:
            continue
        n = int(ln[22:26])
        if not (lo <= n <= hi):
            continue
        if n not in res:
            res[n] = []; order.append(n)
        res[n].append(ln)
    for n, lines in res.items():
        name = lines[0][17:20].strip()
        have = {l[12:16].strip() for l in lines}
        if not set(HEAVY[name].split()) <= have:
            keep = {"N", "CA", "C", "O", "CB", "OXT"}
            res[n] = [l[:17] + "ALA" + l[20:] for l in lines
                      if l[12:16].strip() in keep]
    # split where the numbering jumps or the peptide bond is long
    segs, cur = [], [order[0]]
    for a, b in zip(order, order[1:]):
        c = [l for l in res[a] if l[12:16].strip() == "C"]
        n_ = [l for l in res[b] if l[12:16].strip() == "N"]
        gap = b - a != 1
        if not gap and c and n_:
            p = np.array([float(c[0][30+8*i:38+8*i]) for i in range(3)])
            q = np.array([float(n_[0][30+8*i:38+8*i]) for i in range(3)])
            gap = np.linalg.norm(p - q) > 2.0
        if gap:
            segs.append(cur); cur = []
        cur.append(b)
    segs.append(cur)
    def xyz(ln):
        return np.array([float(ln[30+8*i:38+8*i]) for i in range(3)])

    with open(out, "w") as fh:
        i = 0
        for ch, seg in zip(string.ascii_uppercase, segs):
            for n in seg:
                lines = list(res[n])
                if n == seg[-1] and not any(l[12:16].strip() == "OXT" for l in lines):
                    # crystal C-termini and chain-break ends carry no OXT; place one
                    # by reflecting the carbonyl O through the CA-C axis
                    d = {l[12:16].strip(): l for l in lines}
                    if {"CA", "C", "O"} <= set(d):
                        ca, c, o = xyz(d["CA"]), xyz(d["C"]), xyz(d["O"])
                        u = (c - ca) / np.linalg.norm(c - ca)
                        v = o - c
                        oxt = c + (2*np.dot(v, u)*u - v)
                        t = d["O"]
                        lines.append(t[:12] + " OXT" + t[16:30]
                                     + "".join(f"{q:8.3f}" for q in oxt) + t[54:])
                for ln in lines:
                    i += 1
                    fh.write(ln[:6] + f"{i:5d}" + ln[11:21] + ch + ln[22:])
            fh.write("TER\n")
        fh.write("END\n")
    return len(segs), sum(len(s) for s in segs)


def build_and_time(lo, hi, label, steps=2000):
    nseg, nres = segmented_pdb(lo, hi, "/tmp/seg.pdb")
    pdb = app.PDBFile("/tmp/seg.pdb")
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod = app.Modeller(pdb.topology, pdb.positions)
    mod.addHydrogens(ff)
    system = ff.createSystem(mod.topology, nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.2*unit.nanometer,
                             constraints=app.HBonds, hydrogenMass=3*unit.amu)
    integ = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 4*unit.femtosecond)
    sim = app.Simulation(mod.topology, system, integ, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy(maxIterations=200)
    sim.context.setVelocitiesToTemperature(310*unit.kelvin)
    sim.step(200)
    t0 = time.time(); sim.step(steps); dt = time.time() - t0
    nspd = (steps*4e-6) / (dt/86400.0)
    print(f"  {label:26s} {nres:3d} res / {nseg} seg / {mod.topology.getNumAtoms():5d} atoms"
          f"   {nspd:7.2f} ns/day")
    return nspd


if __name__ == "__main__":
    print("OpenMM 8.6, CPU platform, 4 cores, amber14 + GBn2 implicit, 4 fs + HMR, 1.2 nm cutoff\n")
    lat = build_and_time(245, 319, "latch only (245-319)")
    full = build_and_time(1, 400, "full 7CBC chain A")
    print(f"""
  What the report's MD would cost on THIS machine
    M2  MB-latch, 4 ns                 {4.0/lat*24:7.1f} h
    M2  latch-graphene, 1 ns           {1.0/lat*24:7.1f} h   (graphene sheet not counted)
    M3  umbrella, 16 windows x 0.6 ns  {9.6/full*24:7.1f} h
    M3  steered MD, 1.5 ns             {1.5/full*24:7.1f} h
    M5  Ala scan, 15 x 2 minimisations       minutes, not hours
""")
