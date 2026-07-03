"""Fragment-aware comparison of DrugCLIP GLP1R hits vs measured ChEMBL actives.

Whole-molecule Tanimoto failed because DrugCLIP returns fragments (MW~257) and the
actives are drug-sized (MW~607). A fragment can still be *contained in* a larger
active. So here we ask, per DrugCLIP fragment:

  1. Substructure: is the (neutralised) fragment a substructure of >=1 measured
     GLP1R active?  -> count of actives it appears in.
  2. Scaffold: does the fragment's Bemis-Murcko scaffold match an active's
     scaffold (exact) or sit inside one (substructure)?

Enrichment null = DrugCLIP's fragments for an unrelated target (EGFR). If DrugCLIP
is target-specific, GLP1R fragments should be substructures of GLP1R actives more
often than EGFR fragments are. Fisher exact test.

To avoid trivial matches (a lone benzene matches everything) we require "specific"
fragments: >= MIN_HEAVY heavy atoms and >= MIN_RINGS ring.

Outputs
-------
  data/targets/GLP1R/GLP1R_fragment_consensus.csv   fragments contained in actives
  results/figures/glp1r_fragment_match.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize
from scipy.stats import fisher_exact

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[1]
GLP1R = ROOT / "data" / "targets" / "GLP1R"
FIG = ROOT / "results" / "figures"; FIG.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREY = "#0072B2", "#E69F00", "#8A8A8A"

MIN_HEAVY = 10   # fragment specificity floor
MIN_RINGS = 1
_UNCHARGE = rdMolStandardize.Uncharger()


def neutral_mol(smiles: str):
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    try:
        return _UNCHARGE.uncharge(m)
    except Exception:
        return m


def murcko(m):
    try:
        return MurckoScaffold.GetScaffoldForMol(m)
    except Exception:
        return None


def is_specific(m) -> bool:
    return m is not None and m.GetNumHeavyAtoms() >= MIN_HEAVY and \
        Chem.rdMolDescriptors.CalcNumRings(m) >= MIN_RINGS


def count_substruct_hits(frag, actives) -> int:
    """How many active molecules contain this fragment as a substructure."""
    if frag is None or frag.GetNumHeavyAtoms() == 0:
        return 0
    return sum(a.HasSubstructMatch(frag) for a in actives if a is not None)


def analyse(frag_smiles, active_mols, active_scaffolds):
    frags = [neutral_mol(s) for s in frag_smiles]
    rows = []
    for s, m in zip(frag_smiles, frags):
        spec = is_specific(m)
        n_hits = count_substruct_hits(m, active_mols) if spec else 0
        sc = murcko(m) if m else None
        sc_smiles = Chem.MolToSmiles(sc) if sc and sc.GetNumHeavyAtoms() else ""
        scaffold_exact = sc_smiles in active_scaffolds if sc_smiles else False
        rows.append({"smiles": s, "specific": spec, "n_active_hits": n_hits,
                     "in_active": n_hits > 0, "scaffold": sc_smiles,
                     "scaffold_in_actives": scaffold_exact})
    return pd.DataFrame(rows)


def main() -> None:
    chembl = pd.read_csv(GLP1R / "GLP1R_ligands.csv")
    actives = [neutral_mol(s) for s in chembl[chembl["class"] == "small_molecule"]["smiles"]]
    actives = [a for a in actives if a is not None]
    active_scaffolds = set()
    for a in actives:
        sc = murcko(a)
        if sc and sc.GetNumHeavyAtoms():
            active_scaffolds.add(Chem.MolToSmiles(sc))
    print(f"measured small-molecule actives: {len(actives)} "
          f"({len(active_scaffolds)} unique Murcko scaffolds)")

    drug = pd.read_csv(GLP1R / "GLP1R_drugclip_predicted.csv").drop_duplicates("ligand_id")
    decoy = pd.read_csv(ROOT / "data" / "targets" / "EGFR_decoy" /
                        "EGFR_decoy_drugclip_predicted.csv").drop_duplicates("ligand_id")

    gd = analyse(drug["smiles"], actives, active_scaffolds)
    dd = analyse(decoy["smiles"], actives, active_scaffolds)
    gd = pd.concat([drug.reset_index(drop=True), gd], axis=1)

    # restrict enrichment to specific fragments
    g_spec, d_spec = gd[gd["specific"]], dd[dd["specific"]]
    g_hit = int(g_spec["in_active"].sum()); g_tot = len(g_spec)
    d_hit = int(d_spec["in_active"].sum()); d_tot = len(d_spec)
    odds, p = fisher_exact([[g_hit, g_tot - g_hit], [d_hit, d_tot - d_hit]],
                           alternative="greater")
    g_rate = g_hit / max(g_tot, 1); d_rate = d_hit / max(d_tot, 1)

    print(f"\nSPECIFIC fragments (>= {MIN_HEAVY} heavy atoms, >= {MIN_RINGS} ring):")
    print(f"  GLP1R : {g_hit}/{g_tot} ({g_rate:.1%}) are substructures of >=1 active")
    print(f"  EGFR  : {d_hit}/{d_tot} ({d_rate:.1%}) (decoy)")
    print(f"  Fisher exact (GLP1R>decoy): odds={odds:.2f}, p={p:.3g}")
    print(f"  scaffold-exact matches: GLP1R={int(g_spec['scaffold_in_actives'].sum())}, "
          f"EGFR={int(d_spec['scaffold_in_actives'].sum())}")

    # consensus = specific fragments contained in >=1 active, ranked
    cons = gd[gd["in_active"]].sort_values(
        ["n_active_hits", "drugclip_score"], ascending=False)
    keep = ["ligand_id", "smiles", "source", "drugclip_score", "docking_score",
            "n_active_hits", "scaffold_in_actives", "nearby_residues"]
    cons[keep].to_csv(GLP1R / "GLP1R_fragment_consensus.csv", index=False)
    print(f"\nfragment-corroborated consensus: {len(cons)} fragments "
          f"-> {GLP1R/'GLP1R_fragment_consensus.csv'}")
    if len(cons):
        print(cons[["ligand_id", "source", "drugclip_score", "n_active_hits"]].head(10).to_string(index=False))

    # --- figure -----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 130})
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 5.0))

    # Panel A: substructure containment rate, GLP1R vs decoy
    bars = ax[0].bar(["GLP1R\nDrugCLIP", "EGFR decoy"], [g_rate, d_rate],
                     color=[BLUE, GREY], width=0.6)
    ymax = max(g_rate, d_rate)
    ax[0].set_ylim(0, ymax * 1.25 if ymax > 0 else 1)
    for b, hit, tot in zip(bars, [g_hit, d_hit], [g_tot, d_tot]):
        ax[0].text(b.get_x() + b.get_width() / 2, b.get_height() + ymax * 0.03,
                   f"{hit}/{tot}", ha="center", fontsize=9)
    ax[0].set_ylabel("Fraction contained in a measured active")
    ax[0].set_title(f"A  Fragment-in-active containment\nFisher p={p:.2g}, OR={odds:.2f} (n.s.)")

    # Panel B: for corroborated fragments, how many actives each sits inside
    if len(cons):
        c = cons.head(12).iloc[::-1]
        ypos = np.arange(len(c))
        cols = [BLUE if s == "REAL" else ORANGE for s in c["source"]]
        ax[1].barh(ypos, c["n_active_hits"], color=cols, height=0.68)
        ax[1].set_yticks(ypos); ax[1].set_yticklabels(c["ligand_id"], fontsize=7.5)
        ax[1].set_xlabel("# measured actives containing this fragment")
        ax[1].set_title("B  Fragment-corroborated hits")
    else:
        ax[1].text(0.5, 0.5, "no fragment contained\nin any measured active",
                   ha="center", va="center", transform=ax[1].transAxes)
        ax[1].set_title("B  Fragment-corroborated hits")
    fig.tight_layout()
    out = FIG / "glp1r_fragment_match.png"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote figure -> {out}")


if __name__ == "__main__":
    main()
