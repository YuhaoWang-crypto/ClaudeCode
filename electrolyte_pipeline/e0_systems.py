"""
E0 — model definition  [digest p4: "建模信息至少写到这个程度"]

Everything that later figures depend on is defined here once:
composition (counts, not just "1 M"), units and their conversion basis,
boundary conditions, interaction model (force field, charges, scaling),
and dynamics settings.  `system_record()` dumps all of it as JSON next to
every trajectory so that each plot is traceable (digest p4, "建议保留四类文件").
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK = os.path.join(ROOT, "elyte_work")          # trajectories, logs (git-ignored)
FIG = os.path.join(ROOT, "figures", "electrolyte")
os.makedirs(WORK, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# --------------------------------------------------------------------------
# Species.  SMILES are the topology source; NAGL (GNN AM1-BCC) supplies charges.
# --------------------------------------------------------------------------
SPECIES = {
    "DME":  {"smiles": "COCCOC", "charge": 0, "mw": 90.121,
             "note": "1,2-dimethoxyethane; two ether O can chelate Li+ (bidentate)"},
    "EC":   {"smiles": "O=C1OCCO1", "charge": 0, "mw": 88.062,
             "note": "ethylene carbonate; carbonyl O is the usual Li+ site"},
    "TFSI": {"smiles": "[N-](S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F", "charge": -1,
             "mw": 280.146, "note": "bis(trifluoromethanesulfonyl)imide anion"},
    "Li":   {"smiles": "[Li+]", "charge": +1, "mw": 6.941, "note": "lithium cation"},
    "Na":   {"smiles": "[Na+]", "charge": +1, "mw": 22.990, "note": "sodium cation"},
}

# --------------------------------------------------------------------------
# Interaction model — one place, so that a change is a recorded change.
# --------------------------------------------------------------------------
@dataclass
class ForceFieldSpec:
    name: str = "openff-2.2.1.offxml"          # OpenFF Sage 2.2.1 (SMIRNOFF)
    charge_model: str = "openff-gnn-am1bcc-0.1.0-rc.3.pt"   # NAGL, no AmberTools needed
    ion_charge_scale: float = 0.8              # scaled charges for non-polarisable FF
    cutoff_nm: float = 1.0
    switch_nm: float = 0.9
    pme: bool = True
    constraints: str = "h-bonds"
    note: str = ("Generic SMIRNOFF force field; NOT tuned for LiTFSI/DME. "
                 "Ion charges scaled by 0.8 (common ECC-style practice for "
                 "transport in non-polarisable models). Sensitivity to this "
                 "choice is checked explicitly in E4.")


@dataclass
class DynamicsSpec:
    temperature_K: float = 298.15
    pressure_bar: float = 1.0
    timestep_fs: float = 2.0
    friction_per_ps: float = 1.0
    barostat_interval: int = 25
    minimise_tol_kj_mol_nm: float = 10.0
    equil_npt_ns: float = 2.0
    prod_ns: float = 10.0
    frame_ps: float = 2.0          # trajectory frame spacing
    thermo_ps: float = 0.2         # energy/volume log spacing


@dataclass
class Composition:
    label: str
    counts: dict                   # species -> number of molecules / ions
    box_nm: float                  # initial cubic packmol box edge (NPT relaxes it)
    target_note: str = ""          # e.g. "≈1 M, solvent:Li = 10"

    @property
    def n_atoms_est(self) -> int:
        na = {"DME": 16, "EC": 10, "TFSI": 15, "Li": 1, "Na": 1}
        return sum(na[s] * n for s, n in self.counts.items())

    @property
    def total_charge(self) -> int:
        return sum(SPECIES[s]["charge"] * n for s, n in self.counts.items())


def molarity_from_box(counts: dict, box_volume_nm3: float, cation="Li") -> float:
    """mol/L from a box volume — the digest asks for the conversion basis: here
    it is N_cation / (N_A * V_box).  1 nm^3 = 1e-24 L."""
    n = counts.get(cation, 0)
    return n / (6.02214076e23 * box_volume_nm3 * 1e-24)


def molality(counts: dict, solvent="DME", cation="Li") -> float:
    """mol salt per kg solvent (independent of density)."""
    m_solv_kg = counts[solvent] * SPECIES[solvent]["mw"] / 1000 / 6.02214076e23 * 6.02214076e23 / 1000
    # -> counts*mw in g/mol; per molecule: mw/N_A g; total kg = counts*mw/N_A/1000
    kg = counts[solvent] * SPECIES[solvent]["mw"] / 6.02214076e23 / 1000
    return counts[cation] / 6.02214076e23 / kg


# --------------------------------------------------------------------------
# The concentration series.  Salt:solvent RATIO is the primary variable
# (digest p5/p6 use "solvent:Li" as well); molarity is *derived* from the
# equilibrated NPT density and reported alongside.
# --------------------------------------------------------------------------
N_DME = 200
SERIES = [
    Composition("DME_pure",   {"DME": N_DME},                              3.4, "0 M, solvent only"),
    Composition("LiTFSI_r20", {"DME": N_DME, "Li": 10,  "TFSI": 10},       3.5, "solvent:Li = 20  (≈0.45 M)"),
    Composition("LiTFSI_r10", {"DME": N_DME, "Li": 20,  "TFSI": 20},       3.6, "solvent:Li = 10  (≈0.85 M)"),
    Composition("LiTFSI_r5",  {"DME": N_DME, "Li": 40,  "TFSI": 40},       3.8, "solvent:Li = 5   (≈1.5 M)"),
    Composition("LiTFSI_r3",  {"DME": N_DME, "Li": 67,  "TFSI": 67},       4.0, "solvent:Li = 3   (≈2.2 M)"),
    # Li/Na comparison at one ratio (digest p8 Fig.25)
    Composition("NaTFSI_r10", {"DME": N_DME, "Na": 20,  "TFSI": 20},       3.6, "solvent:Na = 10  (≈0.85 M)"),
]
SERIES_BY_LABEL = {c.label: c for c in SERIES}

# --------------------------------------------------------------------------
# Experimental reference values for the comparison in E4.
# Values are *approximate literature numbers* entered by hand; the digest's
# own rule (p8 table, "原始实验来源与定义是否一致") applies: check the source
# before quoting them in a paper.  Left empty where I am not confident.
# --------------------------------------------------------------------------
EXPERIMENT = {
    "DME_pure": {
        "density_g_cm3": 0.8637, "viscosity_mPa_s": 0.42, "dielectric": 7.2,
        "T_K": 298.15,
        "source": "CRC Handbook / Riddick–Bunger 'Organic Solvents' (pure DME, 25 °C)",
    },
    # LiTFSI/DME concentration series — read from Fig.25 (open circles) of the
    # review, which in turn cites refs 458–460.  Reading a log axis by eye is
    # ±20 %; these are guide values only.
    "LiTFSI_r20": {"viscosity_mPa_s": 0.6, "source": "Chem.Rev. Fig.25 open circles, ~0.45 M (eye-read)"},
    "LiTFSI_r10": {"viscosity_mPa_s": 0.85, "source": "Chem.Rev. Fig.25 open circles, ~0.9 M (eye-read)"},
    "LiTFSI_r5":  {"viscosity_mPa_s": 2.0, "source": "Chem.Rev. Fig.25 open circles, ~1.5 M (eye-read)"},
    "LiTFSI_r3":  {"viscosity_mPa_s": 6.0, "source": "Chem.Rev. Fig.25 open circles, ~2.2 M (eye-read)"},
}


def system_record(comp: Composition, ff: ForceFieldSpec, dyn: DynamicsSpec,
                  extra: dict | None = None) -> dict:
    rec = {
        "composition": {"label": comp.label, "counts": comp.counts,
                        "target_note": comp.target_note,
                        "total_charge": comp.total_charge,
                        "n_atoms_est": comp.n_atoms_est,
                        "molality_mol_kg": (molality(comp.counts,
                                                     cation="Na" if "Na" in comp.counts else "Li")
                                            if ("Li" in comp.counts or "Na" in comp.counts) else 0.0)},
        "boundary": {"periodic": True, "initial_box_nm": comp.box_nm,
                     "ensemble": "NPT (equil + production)"},
        "interactions": asdict(ff),
        "dynamics": asdict(dyn),
        "species": {s: SPECIES[s] for s in comp.counts},
    }
    if extra:
        rec.update(extra)
    return rec


def write_record(path: str, rec: dict):
    with open(path, "w") as f:
        json.dump(rec, f, indent=2)


if __name__ == "__main__":
    ff, dyn = ForceFieldSpec(), DynamicsSpec()
    print(f"{'label':12s} {'DME':>4s} {'Li/Na':>5s} {'TFSI':>4s} {'atoms':>6s} {'m (mol/kg)':>10s}  note")
    for c in SERIES:
        rec = system_record(c, ff, dyn)
        print(f"{c.label:12s} {c.counts.get('DME',0):4d} "
              f"{c.counts.get('Li', c.counts.get('Na',0)):5d} {c.counts.get('TFSI',0):4d} "
              f"{c.n_atoms_est:6d} {rec['composition']['molality_mol_kg']:10.3f}  {c.target_note}")
