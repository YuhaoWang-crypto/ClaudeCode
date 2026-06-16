"""
make_colab_notebook.py
----------------------
Emit a self-contained Colab notebook that runs the GCK distal-Zn-switch
substrate-binding MM-GBSA for the two states (no-Zn vs +Zn) built by
src/build_md_models.py. The user opens it in Colab, picks a GPU runtime, and
Runs All.

Design choices to keep it robust (it is generated here but can only be *run* on
Colab with a GPU + AmberTools, so we minimise fragile steps):
  * ligand (glucose) -> GAFF2/AM1-BCC via antechamber
  * Zn -> 12-6 divalent ion params + harmonic His(N)-Zn distance restraints,
    with the restrained nitrogens auto-picked as the 2 protein N closest to Zn
    (no reliance on tleap residue renumbering)
  * MM-GBSA computed directly in OpenMM with GBn2 by re-evaluating
    complex / receptor / ligand energies over the trajectory (no MMPBSA.py setup)
"""
import json, os

REPO = "https://github.com/YuhaoWang-crypto/ClaudeCode.git"
BRANCH = "claude/determined-cray-fv87ut"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "colab", "GCK_Zn_switch_MMGBSA.ipynb")


def md(*lines):
    return {"cell_type": "markdown", "metadata": {}, "source": list(lines)}


def code(src):
    return {"cell_type": "code", "metadata": {}, "execution_count": None,
            "outputs": [], "source": src.strip("\n").splitlines(keepends=True)}


cells = [
md("# GCK distal Zn-switch → glucose-binding ΔΔG (MM-GBSA)\n",
   "\n",
   "Tests whether a distal engineered His-pair Zn site (G193H/F195H, ~24 Å from the\n",
   "active site) perturbs **glucose binding** in human glucokinase — a computational\n",
   "proxy for activity perturbation.\n",
   "\n",
   "**ΔΔG = ΔG_bind(+Zn) − ΔG_bind(no Zn)**.  |ΔΔG| large ⇒ the distal metal is\n",
   "allosterically coupled to the active site.\n",
   "\n",
   "### Before you run\n",
   "1. **Runtime ▸ Change runtime type ▸ GPU**.\n",
   "2. **Runtime ▸ Run all**. Cell 1 installs conda and auto-restarts the kernel —\n",
   "   that is expected; just let it continue.\n",
   "\n",
   "Honest caveats: short MD (set for speed) → treat ΔΔG as a *first estimate*;\n",
   "use ≥3 replicas and longer production for a real number. Zn is a 12-6 ion with\n",
   "restraints (approximate); a ZAFF bonded model is more physical. Direction\n",
   "(activate/inhibit) must ultimately be confirmed by a kcat/Km ± Zn assay."),

code("""
# 1. conda on Colab (auto-restarts the kernel; this is normal)
!pip -q install condacolab
import condacolab
condacolab.install()
"""),

code("""
# 2. MD / parametrisation toolchain
!mamba install -q -y -c conda-forge ambertools openmm openmmforcefields openff-toolkit parmed mdanalysis
import openmm, parmed
print("openmm", openmm.version.version, "| parmed", parmed.__version__)
"""),

code(f"""
# 3. fetch the pre-built models from the branch
import os
if not os.path.isdir("ClaudeCode"):
    !git clone -q --branch {BRANCH} {REPO}
MODELS = "ClaudeCode/data/models"
# If the repo is private and the clone failed, uncomment to upload the 3 PDBs:
# from google.colab import files; up = files.upload()
#   then set MODELS = "." after uploading GCK_*_glc.pdb
print(sorted(os.listdir(MODELS)))
"""),

code("""
# 4. split each model into protein / glucose / (Zn)
import os
WORK = "work"; os.makedirs(WORK, exist_ok=True)

def split_pdb(src, tag):
    prot, glc, zn = [], [], []
    for ln in open(src):
        if ln.startswith(("ATOM", "HETATM", "TER")):
            res = ln[17:20].strip()
            if ln.startswith("ATOM") or res not in ("GLC", "ZN"):
                prot.append(ln)
            elif res == "GLC":
                glc.append(ln)
            elif res == "ZN":
                zn.append(ln)
    pp = f"{WORK}/{tag}_protein.pdb"; open(pp, "w").write("".join(prot) + "END\\n")
    gp = f"{WORK}/{tag}_glc.pdb";     open(gp, "w").write("".join(glc) + "END\\n")
    zp = None
    if zn:
        zp = f"{WORK}/{tag}_zn.pdb";   open(zp, "w").write("".join(zn) + "END\\n")
    return pp, gp, zp

states = {}
states["noZn"] = split_pdb(f"{MODELS}/GCK_mutHis_glc.pdb", "noZn")
states["Zn"]   = split_pdb(f"{MODELS}/GCK_mutHis_Zn_glc.pdb", "Zn")
print(states)
"""),

code("""
# 5. parametrise glucose once (GAFF2 + AM1-BCC); clean it with pdb4amber first
!pdb4amber -i work/Zn_glc.pdb -o work/glc_clean.pdb >/dev/null 2>&1
!cd work && antechamber -i glc_clean.pdb -fi pdb -o glc.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 -s 0 >/dev/null
!cd work && parmchk2 -i glc.mol2 -f mol2 -o glc.frcmod -s gaff2
print("glucose params:", os.path.exists("work/glc.mol2"), os.path.exists("work/glc.frcmod"))
"""),

code("""
# 6. tleap: build solvated Amber topology for a state (with/without Zn)
def build(tag, with_zn):
    prot, glc, zn = states[tag]
    !pdb4amber -i {prot} -o work/{tag}_prot_clean.pdb -y >/dev/null 2>&1
    zn_lines = ""
    if with_zn:
        zn_lines = (
            "loadamberparams frcmod.ions234lm_126_tip3p\\n"
            f"ZN = loadpdb {zn}\\n")
        combine = "complex = combine { prot GLC ZN }"
    else:
        combine = "complex = combine { prot GLC }"
    leap = f'''source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams work/glc.frcmod
GLC = loadmol2 work/glc.mol2
prot = loadpdb work/{tag}_prot_clean.pdb
{zn_lines}{combine}
solvateBox complex TIP3PBOX 12.0
addIonsRand complex Na+ 0 Cl- 0
saveamberparm complex work/{tag}.prmtop work/{tag}.inpcrd
savepdb complex work/{tag}_solv.pdb
quit
'''
    open(f"work/{tag}.leap", "w").write(leap)
    !tleap -f work/{tag}.leap > work/{tag}.leap.log 2>&1
    ok = os.path.exists(f"work/{tag}.prmtop")
    print(tag, "prmtop:", ok)
    if not ok:
        !tail -25 work/{tag}.leap.log

build("noZn", with_zn=False)
build("Zn",   with_zn=True)
"""),

code('''
# 7. minimise + (restrained, for +Zn) equilibrate + short production MD in OpenMM
import numpy as np, parmed
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloBarostat, CustomBondForce, Platform
from openmm.app import PME, HBonds, DCDReporter, StateDataReporter

NS_PROD = 1.0           # production length (ns) -- raise for a real number
PLATFORM = Platform.getPlatformByName("CUDA")

def zn_restraints(system, parm):
    """Harmonic Zn-N(His) restraints on the 2 protein N closest to Zn."""
    xyz = parm.coordinates
    zn_idx = [a.idx for a in parm.atoms if a.residue.name == "ZN"]
    if not zn_idx:
        return None
    zi = zn_idx[0]; zpos = xyz[zi]
    cand = [(np.linalg.norm(xyz[a.idx] - zpos), a.idx) for a in parm.atoms
            if a.element_name == "N" and a.residue.name in ("HIS", "HIE", "HID", "HIP")]
    cand.sort()
    pair = [i for _, i in cand[:2]]
    f = CustomBondForce("0.5*k*(r-r0)^2")
    f.addPerBondParameter("k"); f.addPerBondParameter("r0")
    for ni in pair:
        f.addBond(zi, ni, [200000.0, 0.21])   # 0.21 nm, strong
    system.addForce(f)
    print("  Zn restrained to N atoms", pair)
    return f

def run_md(tag):
    parm = parmed.load_file(f"work/{tag}.prmtop", f"work/{tag}.inpcrd")
    system = parm.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*unit.nanometer,
                               constraints=HBonds)
    zn_restraints(system, parm)
    system.addForce(MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
    integ = LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    sim = app.Simulation(parm.topology, system, integ, PLATFORM)
    sim.context.setPositions(parm.positions)
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(300*unit.kelvin)
    sim.step(50000)                                   # 100 ps equilibration
    sim.reporters.append(DCDReporter(f"work/{tag}.dcd", 5000))
    sim.reporters.append(StateDataReporter(f"work/{tag}.csv", 5000, step=True,
                         potentialEnergy=True, temperature=True))
    sim.step(int(NS_PROD*500000))                     # production
    print(tag, "MD done ->", f"work/{tag}.dcd")

run_md("noZn")
run_md("Zn")
'''),

code('''
# 8. MM-GBSA (single-trajectory) directly in OpenMM with GBn2:
#    dG_bind = <E_complex - E_receptor - E_ligand> over the trajectory (no entropy)
import numpy as np, parmed, mdtraj as mdt
from openmm import unit, VerletIntegrator, Platform
from openmm.app import Simulation, NoCutoff, GBn2, HBonds
CPU = Platform.getPlatformByName("CPU")

def gb_energy(struct):
    sysm = struct.createSystem(implicitSolvent=GBn2, nonbondedMethod=NoCutoff,
                               constraints=HBonds)
    sim = Simulation(struct.topology, sysm, VerletIntegrator(1*unit.femtosecond), CPU)
    def E(positions):
        sim.context.setPositions(positions)
        return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilocalorie_per_mole)
    return E

def mmgbsa(tag, stride=5):
    parm = parmed.load_file(f"work/{tag}.prmtop")
    # strip solvent/ions -> dry complex; receptor = complex w/o GLC; ligand = GLC
    dry = parm["!:WAT,Na+,Cl-,K+"]
    rec = dry["!:GLC"]
    lig = dry[":GLC"]
    Ecx, Erc, Elg = gb_energy(dry), gb_energy(rec), gb_energy(lig)
    traj = mdt.load(f"work/{tag}.dcd", top=f"work/{tag}.prmtop")
    # atom index maps from full topology to each subsystem
    keep = {a.idx for a in dry.atoms}
    rec_keep = {a.idx for a in rec.atoms}
    lig_keep = {a.idx for a in lig.atoms}
    # build boolean masks over full atom order
    full = list(range(parm.ptr("NATOM")))
    dvals = []
    for fr in range(0, traj.n_frames, stride):
        xyz = traj.xyz[fr]  # nm, full system order
        cx = xyz[sorted(keep)] * unit.nanometer
        rc = xyz[sorted(rec_keep)] * unit.nanometer
        lg = xyz[sorted(lig_keep)] * unit.nanometer
        dvals.append(Ecx(cx) - Erc(rc) - Elg(lg))
    dvals = np.array(dvals)
    print(f"  {tag}: dG_bind = {dvals.mean():.2f} +/- {dvals.std()/np.sqrt(len(dvals)):.2f} "
          f"kcal/mol (n={len(dvals)})")
    return dvals.mean(), dvals.std()/np.sqrt(len(dvals))

g_no, e_no = mmgbsa("noZn")
g_zn, e_zn = mmgbsa("Zn")
'''),

code('''
# 9. result
ddg = g_zn - g_no
err = (e_no**2 + e_zn**2) ** 0.5
print("="*56)
print(f"dG_bind(no Zn) = {g_no:7.2f} kcal/mol")
print(f"dG_bind(+Zn)   = {g_zn:7.2f} kcal/mol")
print(f"ddG (+Zn - noZn) = {ddg:6.2f} +/- {err:.2f} kcal/mol")
print("="*56)
print("|ddG| within ~1 kcal/mol of zero (and of the error) => the distal Zn")
print("switch does NOT measurably perturb glucose binding in this short run.")
print("A clear non-zero ddG => coupling; sign = tighter(-)/weaker(+) binding.")
print("Confirm with longer MD + replicas, and ultimately a kcat/Km +/- Zn assay.")
'''),
]

nb = {
    "cells": cells,
    "metadata": {
        "accelerator": "GPU",
        "colab": {"provenance": [], "toc_visible": True},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {"name": "python"},
    },
    "nbformat": 4,
    "nbformat_minor": 0,
}

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as fh:
    json.dump(nb, fh, indent=1)
print("wrote", OUT, "with", len(cells), "cells")
