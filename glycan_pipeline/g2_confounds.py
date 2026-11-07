"""
Module G2 - Two confounds that will wreck a naive DeltaRASA screen.

DeltaRASA is simple enough that it is tempting to point it at whatever
coordinates you have and rank the results. Doing that on experimental
structures produces rankings driven by artefacts rather than by shielding.
This module measures both artefacts on a real glycoprotein so the size of the
problem is a number rather than a worry.

Confound 1 - glycan size
    Mean per-atom DeltaRASA averages over the whole tree. A one-residue GlcNAc
    stub sits flat against the protein, so nearly every atom of it is occluded
    and the mean is high. A full Man8 projects most of its atoms into bulk
    solvent, which dilutes the mean no matter how buried its core is. The
    score therefore partly measures how much glycan was modelled.

    This is exactly why the paper does not score experimental glycans: it
    grafts the SAME GlycoShape-modelled Man9 (or Man6, for AF3-server size
    limits) onto every design, so tree size is constant across the comparison.
    Reproduce that discipline or the ranking is not about shielding.

Confound 2 - missing coordinates
    Glycan arms are mobile and frequently unresolved, so a deposited structure
    under-reports glycan-glycan crowding. Since crowding is the mechanism
    behind the HIV-1 "intrinsic mannose patch", a DeltaRASA screen run on
    crystallographic coordinates cannot see the very effect that motivated
    the paper. Quantified here as the gap between modelled tree size and an
    intact Man9GlcNAc2, and as the (small, lower-bound) neighbour-glycan term.

Note also that the paper's own definition deliberately excludes neighbouring
glycans from the context: its designs carry a single glycan, so DeltaRASA
there means protein-pocket burial specifically, not crowding. Both contexts
are computed here so the distinction is explicit.

Rigor
-----
(rigorous) Every correlation and DeltaRASA value below is computed from
deposited coordinates.
(caution) The neighbour-glycan term is a LOWER BOUND, for exactly the reason
confound 2 describes. Do not read a small value as "no crowding here".
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .g1_delta_rasa import FIG_DIR, site_delta_rasa
from .glycans import find_glycan_trees, load_model
from .pdbcache import fetch
from .util import spearman

#: Residue count of an intact Man9GlcNAc2 (2 GlcNAc + 9 Man), the species the
#: paper models onto its designs.
MAN9_RESIDUES = 11

#: Residues kept when standardising size: the shared GlcNAc-GlcNAc-Man core.
CORE_RESIDUES = 3


def analyse(pdb_id: str = "5FYL", assembly: int | None = 1):
    """Measure both confounds on one structure."""
    model = load_model(fetch(pdb_id, assembly=assembly))
    trees = [t for t in find_glycan_trees(model) if t.asn_resseq is not None]

    full = site_delta_rasa(model, trees, "protein")
    crowd = {r["site"]: r for r in
             site_delta_rasa(model, trees, "protein+other_glycans")}

    # Size-standardised control: every tree cut back to its shared core.
    core_trees = [t.truncate(CORE_RESIDUES) for t in trees]
    core = {r["site"]: r for r in site_delta_rasa(model, core_trees, "protein")}

    rows = []
    for r in full:
        site = r["site"]
        c = core.get(site)
        w = crowd.get(site)
        rows.append({
            "site": site,
            "modelled_species": r["modelled_species"],
            "n_sugar_residues": r["n_sugar_residues"],
            "n_atoms": r["n_atoms"],
            "completeness_vs_man9": r["n_sugar_residues"] / MAN9_RESIDUES,
            "drasa_full": r["mean_delta_rasa"],
            "drasa_core3": c["mean_delta_rasa"] if c else float("nan"),
            "core3_n_residues": c["n_sugar_residues"] if c else 0,
            "crowding_term": (w["mean_delta_rasa"] - r["mean_delta_rasa"]
                              if w else float("nan")),
        })

    # Confound 1, measured.
    size = [x["n_atoms"] for x in rows]
    d_full = [x["drasa_full"] for x in rows]
    rho_size_full = spearman(size, d_full)

    matched = [x for x in rows if x["core3_n_residues"] == CORE_RESIDUES]
    rho_size_core = spearman([x["n_atoms"] for x in rows],
                             [x["drasa_core3"] for x in rows])
    rho_full_vs_core = spearman([x["drasa_full"] for x in matched],
                                [x["drasa_core3"] for x in matched])

    # Confound 2, measured.
    completeness = np.array([x["completeness_vs_man9"] for x in rows])
    crowding = np.array([x["crowding_term"] for x in rows])

    # Free precision check: in an assembly the same site recurs once per
    # protomer. Those copies are independent trips through parsing, tree
    # building, neighbour search and quadrature, so their spread is an
    # end-to-end numerical error bar on the whole pipeline.
    by_site: dict[tuple, list[float]] = {}
    for x in rows:
        base = x["site"].split(":")[0].split("-")[0]
        by_site.setdefault((base, x["site"].split(":")[1]), []).append(
            x["drasa_full"])
    replicated = {k: v for k, v in by_site.items() if len(v) > 1}
    spreads = [max(v) - min(v) for v in replicated.values()]

    return {
        "n_replicated_sites": len(replicated),
        "replicate_copies": max((len(v) for v in replicated.values()), default=0),
        "max_replicate_spread": float(max(spreads)) if spreads else float("nan"),
        "median_replicate_spread": float(np.median(spreads)) if spreads
        else float("nan"),
        "pdb_id": pdb_id,
        "assembly": assembly,
        "n_sites": len(rows),
        "rows": rows,
        "rho_size_vs_drasa_full": rho_size_full,
        "rho_size_vs_drasa_core3": rho_size_core,
        "rho_full_vs_core3": rho_full_vs_core,
        "n_size_matched": len(matched),
        "median_completeness_vs_man9": float(np.median(completeness)),
        "frac_sites_below_half_man9": float((completeness < 0.5).mean()),
        "mean_crowding_term": float(np.nanmean(crowding)),
        "max_crowding_term": float(np.nanmax(crowding)),
        "frac_sites_zero_crowding": float((crowding <= 1e-9).mean()),
    }


def make_figure(r, path=None):
    path = path or os.path.join(FIG_DIR, "g2_confounds.png")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    rows = r["rows"]

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    ax = axes[0]
    ax.scatter([x["n_atoms"] for x in rows], [x["drasa_full"] for x in rows],
               s=26, alpha=0.75, color="#4C72B0")
    ax.set_xlabel("glycan atoms modelled")
    ax.set_ylabel(r"mean per-atom $\Delta$RASA (%)")
    ax.set_title("Confound 1: size drives the score\n"
                 rf"Spearman $\rho$ = {r['rho_size_vs_drasa_full']:.2f}",
                 fontsize=10)
    ax.grid(alpha=0.3)

    ax = axes[1]
    matched = [x for x in rows if x["core3_n_residues"] == CORE_RESIDUES]
    ax.scatter([x["drasa_full"] for x in matched],
               [x["drasa_core3"] for x in matched],
               s=26, alpha=0.75, color="#C44E52")
    lim = max([x["drasa_full"] for x in matched]
              + [x["drasa_core3"] for x in matched]) * 1.1
    ax.plot([0, lim], [0, lim], ls="--", c="gray", lw=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(r"$\Delta$RASA, whole modelled tree (%)")
    ax.set_ylabel(r"$\Delta$RASA, size-standardised core (%)")
    ax.set_title("Standardising size changes the ranking\n"
                 rf"Spearman $\rho$ = {r['rho_full_vs_core3']:.2f} "
                 f"(n={len(matched)})", fontsize=10)
    ax.grid(alpha=0.3)

    ax = axes[2]
    comp = sorted(x["completeness_vs_man9"] for x in rows)
    ax.plot(np.arange(1, len(comp) + 1), comp, "o-", ms=4, color="#55A868")
    ax.axhline(1.0, ls="--", c="gray", lw=1)
    ax.text(1, 1.02, "intact Man9GlcNAc2", fontsize=7, color="gray")
    ax.set_xlabel("site (sorted)")
    ax.set_ylabel("modelled residues / 11")
    ax.set_title("Confound 2: most glycans are\nonly partly modelled",
                 fontsize=10)
    ax.grid(alpha=0.3)

    fig.suptitle(f"{r['pdb_id']} assembly {r['assembly']} - "
                 r"artefacts in a naive $\Delta$RASA screen", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def report(pdb_id: str = "5FYL"):
    r = analyse(pdb_id)
    print("=" * 68)
    print(f"MODULE G2 - DeltaRASA confounds ({r['pdb_id']} assembly {r['assembly']})")
    print("=" * 68)
    print(f"glycosylation sites analysed : {r['n_sites']}")

    print("\nPRECISION - symmetry-copy replicates (end-to-end error bar)")
    print(f"    {r['n_replicated_sites']} sites present in "
          f"{r['replicate_copies']} protomer copies")
    print(f"    median / max spread in dRASA across copies : "
          f"{r['median_replicate_spread']:.3f}% / "
          f"{r['max_replicate_spread']:.3f}%")
    print("    -> differences smaller than this are not meaningful.")

    print("\nCONFOUND 1 - glycan size")
    print(f"    Spearman(atoms modelled, dRASA)      : "
          f"{r['rho_size_vs_drasa_full']:+.3f}")
    print(f"    after standardising to a 3-sugar core: "
          f"{r['rho_size_vs_drasa_core3']:+.3f}")
    print(f"    rank agreement full vs standardised  : "
          f"{r['rho_full_vs_core3']:+.3f}  (n={r['n_size_matched']})")
    print("    -> graft ONE identical glycan onto every design before ranking,")
    print("       which is what the paper does with GlycoShape Man9/Man6.")

    print("\nCONFOUND 2 - missing coordinates")
    print(f"    median modelled size vs Man9GlcNAc2  : "
          f"{r['median_completeness_vs_man9']:.2f}")
    print(f"    sites under half an intact Man9      : "
          f"{r['frac_sites_below_half_man9']:.0%}")
    print(f"    neighbour-glycan term, mean / max    : "
          f"{r['mean_crowding_term']:.2f}% / {r['max_crowding_term']:.2f}%")
    print(f"    sites where it is exactly zero       : "
          f"{r['frac_sites_zero_crowding']:.0%}")
    print("    -> that term is a LOWER BOUND; unresolved arms cannot occlude.")

    print("\n    largest neighbour-glycan contributions:")
    for x in sorted(r["rows"], key=lambda z: -z["crowding_term"])[:5]:
        print(f"        {x['site']:>10}  {x['modelled_species']:>9}  "
              f"protein {x['drasa_full']:5.2f}%  +glycans "
              f"{x['crowding_term']:+5.2f}%")

    p = make_figure(r)
    print(f"\nfigure written to: {p}")
    return r


if __name__ == "__main__":
    report()
