"""
Demo #2 — First-principles "initial screening" of metal centres for a
biomimetic (peroxidase-like) single-atom catalyst.

Motivation (from the chat context)
-----------------------------------
The reference work builds an Mn-organic-complex (MnOC) scaffold, swaps in
different metal centres to make an M-MnOC series, and screens ~ten metals
computationally before settling on a doped design. The user's own system is a
Cr-doped single-atom material targeting peroxidase-like (H2O2-activating)
activity. Step 1 of such a study is an ELECTRONIC-STRUCTURE / CHARGE-
DISTRIBUTION pre-screen of candidate metals.

WHAT THIS IS
------------
A genuine DFT calculation (PySCF, PBE/def2-SVP, unrestricted Kohn-Sham) run
end-to-end in this environment on a deliberately MINIMAL active-site model:
a single metal centre bound to one hydroxo group, [M-OH].

The metal-oxygen unit is the mechanistically relevant fragment: peroxidase
chemistry runs through a high-valent metal-oxo/hydroxo species (the "Compound I"
analogue), so a metal's ability to host and polarise an M-O bond is a
reasonable first descriptor.

For each metal we scan spin states and keep the lowest-energy one, then read off
screening descriptors:
    * ground-state spin multiplicity
    * HOMO / LUMO / HOMO-LUMO gap  (electronic stability proxy)
    * Mulliken charge on the metal  (Lewis acidity / charge distribution)
    * Mulliken spin population on the metal  (radical character)

WHAT THIS IS *NOT*
------------------
This is NOT the real periodic Cr-Co3O4 / Ir-MnOC material. A production study
needs: the full periodic slab in VASP/QE, geometry optimisation, the real
coordination environment, solvent, and reaction-barrier (NEB) calculations.
This molecular fragment is a fast, transparent proxy for the *workflow* only.
Treat the numbers as relative screening signals, not publication values.
"""

import numpy as np
from pyscf import gto, dft
from pyscf.scf import addons

HARTREE2EV = 27.211386245988

# candidate metal centres (relevant to Cr-doping and the M-MnOC screen)
METALS = ["Cr", "Mn", "Fe", "Co", "Ni", "Cu"]

# fixed minimal-model geometry: M-O ~1.80 A, O-H 0.97 A, M-O-H angle 110 deg
# (same geometry for every metal => fair relative comparison)
def geom(metal):
    import math
    a = math.radians(110.0)
    return [
        [metal, (0.0, 0.0, 0.0)],
        ["O",   (1.80, 0.0, 0.0)],
        ["H",   (1.80 + 0.97 * math.cos(a), 0.97 * math.sin(a), 0.0)],
    ]

Z = {"Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29}


def candidate_spins(metal):
    """Return spin (2S = n_alpha - n_beta) values consistent with electron parity.

    Neutral [M-OH]: total electrons = Z(M) + 8(O) + 1(H) = Z + 9.
    """
    n_elec = Z[metal] + 9
    if n_elec % 2 == 0:          # even electrons -> spin 0,2,4 (mult 1,3,5)
        return [0, 2, 4]
    else:                        # odd electrons  -> spin 1,3,5 (mult 2,4,6)
        return [1, 3, 5]


def run_one(metal, spin):
    mol = gto.M(atom=geom(metal), basis="def2-svp", spin=spin, charge=0,
                verbose=0)
    mf = dft.UKS(mol)
    mf.xc = "pbe"
    mf = addons.remove_linear_dep_(mf)
    mf.level_shift = 0.3          # help open-shell TM SCF converge
    mf.max_cycle = 200
    e = mf.kernel()
    if not mf.converged:
        mf = mf.newton()         # second-order fallback
        e = mf.kernel()
    if not mf.converged:
        return None
    return mf, e


def descriptors(mf):
    mol = mf.mol
    # frontier orbitals across both spin channels
    homo = -1e9
    lumo = 1e9
    for spin_ch in (0, 1):
        mo_e = mf.mo_energy[spin_ch]
        occ = mf.mo_occ[spin_ch]
        occ_e = mo_e[occ > 0]
        vir_e = mo_e[occ == 0]
        if len(occ_e):
            homo = max(homo, occ_e.max())
        if len(vir_e):
            lumo = min(lumo, vir_e.min())
    gap = (lumo - homo) * HARTREE2EV

    # Mulliken charges and spins (metal is atom 0)
    pop, chg = mf.mulliken_pop(verbose=0)
    q_metal = chg[0]
    dm = mf.make_rdm1()
    spin_dm = dm[0] - dm[1]
    s = mf.mol.aoslice_by_atom()
    ao0, ao1 = s[0][2], s[0][3]
    ovlp = mf.get_ovlp()
    spin_pop_metal = np.einsum("ij,ji->", spin_dm[ao0:ao1, :], ovlp[:, ao0:ao1])

    return dict(homo_eV=homo * HARTREE2EV, lumo_eV=lumo * HARTREE2EV,
                gap_eV=gap, q_metal=q_metal, spin_metal=spin_pop_metal)


def main():
    rows = []
    for m in METALS:
        best = None
        for sp in candidate_spins(m):
            try:
                res = run_one(m, sp)
            except Exception as exc:
                print(f"  [{m}] spin={sp} failed: {type(exc).__name__}")
                res = None
            if res is None:
                continue
            mf, e = res
            if best is None or e < best[0]:
                best = (e, sp, mf)
        if best is None:
            print(f"[{m}] no converged spin state — skipping")
            continue
        e, sp, mf = best
        d = descriptors(mf)
        d.update(metal=m, mult=sp + 1, energy=e)
        rows.append(d)
        print(f"[{m}] ground mult={sp+1}  E={e:.4f} Ha  "
              f"gap={d['gap_eV']:.2f} eV  q(M)={d['q_metal']:+.2f}  "
              f"spin(M)={d['spin_metal']:+.2f}")

    # simple screening summary table
    print("\n" + "=" * 74)
    print(f"{'metal':<6}{'mult':>5}{'HOMO/eV':>10}{'LUMO/eV':>10}"
          f"{'gap/eV':>9}{'q(M)':>8}{'spin(M)':>9}")
    print("-" * 74)
    for r in sorted(rows, key=lambda x: -x["gap_eV"]):
        print(f"{r['metal']:<6}{r['mult']:>5}{r['homo_eV']:>10.2f}"
              f"{r['lumo_eV']:>10.2f}{r['gap_eV']:>9.2f}"
              f"{r['q_metal']:>8.2f}{r['spin_metal']:>9.2f}")
    print("=" * 74)
    print("Ranked by HOMO-LUMO gap (electronic-stability proxy). Larger gap ~ more")
    print("robust closed-shell-like centre; larger q(M) ~ stronger Lewis acidity;")
    print("larger spin(M) ~ more radical / redox-active metal-oxo character.")
    return rows


if __name__ == "__main__":
    rows = main()

    # bar-chart summary of the screening descriptors
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt, os
        rows_s = sorted(rows, key=lambda x: -x["gap_eV"])
        labels = [r["metal"] for r in rows_s]
        fig, ax = plt.subplots(1, 3, figsize=(12, 4))
        ax[0].bar(labels, [r["gap_eV"] for r in rows_s], color="#2c6fbb")
        ax[0].set_title("HOMO–LUMO gap (stability proxy)")
        ax[0].set_ylabel("eV")
        ax[1].bar(labels, [r["q_metal"] for r in rows_s], color="#e08e0b")
        ax[1].set_title("Mulliken charge on metal (Lewis acidity)")
        ax[1].set_ylabel("|e|")
        ax[2].bar(labels, [abs(r["spin_metal"]) for r in rows_s], color="#d1495b")
        ax[2].set_title("Spin population on metal (radical character)")
        ax[2].set_ylabel(r"$\mu_B$")
        for a in ax:
            a.grid(alpha=0.3, axis="y")
        fig.suptitle("Demo #2 · DFT electronic-structure pre-screen of metal centres\n"
                     "minimal [M–OH] active-site model · PBE/def2-SVP · UKS  "
                     "(workflow proxy, not the real periodic material)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.9))
        out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "metal_center_screen.png")
        fig.savefig(out, dpi=150)
        print(f"\nsaved figure -> {out}")
    except Exception as exc:
        print(f"(plot skipped: {exc})")
