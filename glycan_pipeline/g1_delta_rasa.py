"""
Module G1 - DeltaRASA glycan shielding, per glycosylation site.

This is the computational core of Lunn-Halbert et al.: the scalar that decides
whether a designed N-glycan will stay oligomannose. Everything else in the
design pipeline (RFdiffusionAA -> LigandMPNN -> RF3/AF3) produces candidates;
DeltaRASA is what ranks them.

What this module does
---------------------
1. Proves its own SASA engine against ``Bio.PDB.SASA.ShrakeRupley``.
2. Shows the quadrature convergence of DeltaRASA in ``n_points``, because
   DeltaRASA is a small difference of two similar RASA values and so is more
   sensitive to sphere resolution than either term.
3. Computes per-site DeltaRASA for every Asn-linked glycan in a structure,
   under two definitions of "context":
      - ``protein``               protein atoms only  (the paper's definition)
      - ``protein+other_glycans`` protein plus every OTHER glycan tree
   The difference between the two is the subject of module G2.

Rigor
-----
(rigorous) The SASA, RASA and DeltaRASA numbers: deterministic geometry,
verified against Biopython to 1e-14 A^2, converged in n_points.

(caution) Calling a *crystallographic* tree "Man6" reports what was MODELLED,
which is a lower bound on the real glycan - disordered arms are simply absent
from the coordinates. Per-site composition from a structure is never a
substitute for LC-MS/MS. The paper computes DeltaRASA on a *modelled* Man9 or
Man6 grafted onto every design precisely so the glycan is identical across
designs; crystallographic trees are not standardised that way.
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .glycans import find_glycan_trees, load_model
from .pdbcache import fetch
from .sasa import N_POINTS, delta_rasa, radii_for_elements, validate_against_biopython

FIG_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures"
)

#: Demo structures. 6VXX and 7PFP are two of the structures the paper itself
#: harvested glycan-bearing asparagines from as RFdiffusionAA inputs; 5FYL is
#: a fully glycan-resolved HIV-1 Env trimer, the textbook case of steric
#: restriction producing oligomannose glycans in nature.
DEMO_STRUCTURES = {
    "5FYL": "HIV-1 BG505 SOSIP.664 Env trimer (glycan shield, mannose patch)",
    "7PFP": "RFdiffusionAA glycan-donor input used in the paper",
    "6VXX": "SARS-CoV-2 spike, closed state (paper's glycan-donor input)",
}

PROTEIN_SKIP = {"HOH", "WAT", "DOD"}


def protein_atoms(model, exclude_residues=None):
    """Heavy protein atoms of ``model`` (standard amino-acid residues only)."""
    exclude = exclude_residues or set()
    xyz, els = [], []
    for chain in model:
        for res in chain:
            if res.id[0] != " ":                       # hetero / water
                continue
            if res.get_resname().strip().upper() in PROTEIN_SKIP:
                continue
            if id(res) in exclude:
                continue
            for atom in res:
                if (atom.element or "").strip().upper() == "H":
                    continue
                if atom.get_altloc() not in (" ", "A"):
                    continue
                xyz.append(atom.get_coord())
                els.append(atom.element)
    return np.array(xyz, dtype=np.float64), radii_for_elements(els)


def site_delta_rasa(model, trees=None, context="protein", n_points=N_POINTS):
    """Per-site DeltaRASA for every Asn-linked glycan tree in ``model``.

    ``context``:
      ``protein``                occluders are protein atoms only (paper)
      ``protein+other_glycans``  occluders also include every other tree
    """
    if context not in ("protein", "protein+other_glycans"):
        raise ValueError(f"unknown context: {context}")

    trees = trees if trees is not None else find_glycan_trees(model)
    prot_xyz, prot_r = protein_atoms(model)

    rows = []
    for tree in trees:
        if tree.asn_resseq is None:
            continue
        g_xyz, g_r, _ = tree.coords_radii()
        if len(g_xyz) == 0:
            continue

        occ_xyz, occ_r = prot_xyz, prot_r
        if context == "protein+other_glycans":
            extra_xyz, extra_r = [], []
            for other in trees:
                if other is tree:
                    continue
                o_xyz, o_r, _ = other.coords_radii()
                if len(o_xyz):
                    extra_xyz.append(o_xyz)
                    extra_r.append(o_r)
            if extra_xyz:
                occ_xyz = np.vstack([prot_xyz] + extra_xyz)
                occ_r = np.concatenate([prot_r] + extra_r)

        res = delta_rasa(g_xyz, g_r, occ_xyz, occ_r, n_points=n_points)
        rows.append({
            "site": tree.site,
            "chain": tree.asn_chain,
            "resseq": tree.asn_resseq,
            "n_sugar_residues": tree.n_residues,
            "modelled_species": tree.species(),
            "n_atoms": res["n_atoms"],
            "mean_delta_rasa": res["mean_delta_rasa"],
            "max_delta_rasa": res["max_delta_rasa"],
            "n_contact_atoms": res["n_contact_atoms"],
            "frac_contact_atoms": res["n_contact_atoms"] / max(res["n_atoms"], 1),
            "context": context,
        })
    rows.sort(key=lambda r: -r["mean_delta_rasa"])
    return rows


def convergence_check(model, trees, n_points_list=(100, 240, 480, 960, 1920)):
    """DeltaRASA of one representative tree vs sphere resolution.

    Justifies the default ``N_POINTS``: Biopython's own default of 100 points
    leaves a visible quadrature bias in a quantity built from a difference.
    """
    tree = max((t for t in trees if t.asn_resseq is not None),
               key=lambda t: t.n_residues)
    prot_xyz, prot_r = protein_atoms(model)
    g_xyz, g_r, _ = tree.coords_radii()

    out = []
    for n in n_points_list:
        res = delta_rasa(g_xyz, g_r, prot_xyz, prot_r, n_points=n)
        out.append({
            "n_points": n,
            "mean_delta_rasa": res["mean_delta_rasa"],
            "n_contact_atoms": res["n_contact_atoms"],
        })
    ref = out[-1]["mean_delta_rasa"]
    for row in out:
        row["abs_dev_from_finest"] = abs(row["mean_delta_rasa"] - ref)
    return {"site": tree.site, "n_atoms": len(g_xyz), "rows": out}


def analyse(pdb_id: str = "5FYL", n_points: int = N_POINTS):
    path = fetch(pdb_id)
    model = load_model(path)
    trees = find_glycan_trees(model)
    linked = [t for t in trees if t.asn_resseq is not None]

    return {
        "pdb_id": pdb_id,
        "description": DEMO_STRUCTURES.get(pdb_id, ""),
        "n_trees": len(trees),
        "n_linked": len(linked),
        "validation": validate_against_biopython(n_points=100),
        "convergence": convergence_check(model, linked),
        "protein_context": site_delta_rasa(model, linked, "protein", n_points),
        "both_context": site_delta_rasa(
            model, linked, "protein+other_glycans", n_points),
    }


def make_figure(r, path=None):
    path = path or os.path.join(FIG_DIR, "g1_delta_rasa.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)

    prot = r["protein_context"]
    both = {row["site"]: row for row in r["both_context"]}
    sites = [row["site"] for row in prot]
    x = np.arange(len(sites))
    v_prot = np.array([row["mean_delta_rasa"] for row in prot])
    v_both = np.array([both[s]["mean_delta_rasa"] for s in sites])

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2),
                             gridspec_kw={"width_ratios": [3, 1]})

    ax = axes[0]
    ax.bar(x - 0.21, v_prot, width=0.42, color="#4C72B0",
           label="context = protein only (paper definition)")
    ax.bar(x + 0.21, v_both, width=0.42, color="#55A868",
           label="context = protein + neighbouring glycans")
    ax.set_xticks(x)
    ax.set_xticklabels(sites, rotation=90, fontsize=7)
    ax.set_ylabel(r"mean per-atom $\Delta$RASA (%)")
    ax.set_title(f"{r['pdb_id']}: glycan shielding per site\n{r['description']}",
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis="y", alpha=0.3)

    ax = axes[1]
    rows = r["convergence"]["rows"]
    n = [row["n_points"] for row in rows]
    dev = [row["abs_dev_from_finest"] for row in rows]
    ax.plot(n, dev, "o-", color="#C44E52")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.axvline(N_POINTS, ls="--", c="gray", lw=1)
    ax.text(N_POINTS, max(dev) if max(dev) > 0 else 1e-2, " default",
            fontsize=7, color="gray", va="top")
    ax.set_xlabel("sphere points per atom")
    ax.set_ylabel(r"|$\Delta$RASA $-$ finest| (%)")
    ax.set_title(f"quadrature convergence\n({r['convergence']['site']})",
                 fontsize=9)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def report(pdb_id: str = "5FYL"):
    r = analyse(pdb_id)
    print("=" * 68)
    print(f"MODULE G1 - DeltaRASA glycan shielding ({r['pdb_id']})")
    print("=" * 68)

    v = r["validation"]
    print("SASA engine self-check vs Bio.PDB.SASA.ShrakeRupley:")
    print(f"    max per-atom deviation : {v['max_abs_atom_diff_A2']:.2e} A^2"
          f"   -> {'AGREES' if v['agrees'] else 'MISMATCH'}")

    print(f"\nquadrature convergence on {r['convergence']['site']} "
          f"({r['convergence']['n_atoms']} atoms):")
    for row in r["convergence"]["rows"]:
        print(f"    n_points={row['n_points']:>5}  "
              f"mean dRASA={row['mean_delta_rasa']:6.3f}%  "
              f"dev={row['abs_dev_from_finest']:.4f}%")

    print(f"\n{r['n_linked']} Asn-linked glycan trees "
          f"(of {r['n_trees']} total sugar trees)")
    print("\nper-site mean DeltaRASA, context = protein only (paper definition)")
    print(f"    {'site':>10} {'modelled':>10} {'atoms':>6} "
          f"{'meandRASA':>10} {'maxdRASA':>9} {'contact':>8}")
    for row in r["protein_context"][:15]:
        print(f"    {row['site']:>10} {row['modelled_species']:>10} "
              f"{row['n_atoms']:>6} {row['mean_delta_rasa']:>9.2f}% "
              f"{row['max_delta_rasa']:>8.2f}% "
              f"{row['n_contact_atoms']:>3}/{row['n_atoms']:<4}")

    p = make_figure(r)
    print(f"\nfigure written to: {p}")
    return r


if __name__ == "__main__":
    report()
