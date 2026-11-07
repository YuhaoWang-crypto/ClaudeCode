"""
Glycan bookkeeping: parse structures, find N-glycan trees, map them to their
sequon asparagine, and classify composition along the N-glycan processing
continuum.

Two separate jobs live here:

1. **Structural** (rigorous): pull sugar residues out of a PDB/mmCIF file,
   group them into covalently connected trees, and link each tree to the Asn
   it hangs off. This is pure geometry and is deterministic.

2. **Compositional** (one part rigorous, one part reconstructed): turn a
   composition into a named species (Man9, FHybrid, HexNAc4(F), Core, ...) and
   then onto the paper's ordinal 1-20 "underprocessing" ladder. The naming
   rules are stated in the paper's Methods; the numeric ladder is NOT, and is
   reconstructed here. See ``UNDERPROCESSING_LADDER`` for the full argument.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

from .sasa import radii_for_elements

__all__ = [
    "SUGAR_CLASS",
    "GlycanTree",
    "load_model",
    "find_glycan_trees",
    "classify_composition",
    "underprocessing_score",
    "UNDERPROCESSING_LADDER",
    "MAX_UNDERPROCESSING_SCORE",
]

# ---------------------------------------------------------------------------
# PDB chemical-component codes for the monosaccharides that occur in N-glycans.
# Mapped onto the compositional buckets the paper's LC-MS/MS classification
# uses: HexNAc / Hex / Fuc / NeuAc / Xyl, plus Glc tracked separately because
# a retained terminal glucose is what makes Glc1Man9 the least-processed
# species on the ladder.
# ---------------------------------------------------------------------------
SUGAR_CLASS: dict[str, str] = {
    # N-acetylhexosamines
    "NAG": "HexNAc", "NDG": "HexNAc", "NGA": "HexNAc", "A2G": "HexNAc",
    "GCS": "HexNAc", "BGN": "HexNAc",
    # mannose
    "MAN": "Hex", "BMA": "Hex", "MMA": "Hex",
    # galactose
    "GAL": "Hex", "GLA": "Hex", "GXL": "Hex",
    # glucose (tracked separately as well, see GLUCOSE below)
    "GLC": "Hex", "BGC": "Hex",
    # fucose
    "FUC": "Fuc", "FUL": "Fuc",
    # sialic acids
    "SIA": "NeuAc", "NGC": "NeuAc", "SLB": "NeuAc",
    # xylose (plant/insect glycans)
    "XYS": "Xyl", "XYP": "Xyl",
}

#: Codes that specifically mean glucose, needed to spot Glc1Man9.
GLUCOSE = {"GLC", "BGC"}

#: Codes that specifically mean mannose, needed to call oligomannose species.
MANNOSE = {"MAN", "BMA", "MMA"}

# Bond-length cutoffs (Angstrom). A glycosidic C-O bond is ~1.4 A; the Asn
# ND2-C1 amide bond is ~1.45 A. Cutoffs are deliberately loose to tolerate
# modest coordinate error in low-resolution structures, but well under the
# ~2.4 A at which unrelated atoms start to appear.
GLYCOSIDIC_CUTOFF = 1.9
ND2_C1_CUTOFF = 2.1


@dataclass
class GlycanTree:
    """One covalently connected N-glycan, with its attachment site."""

    residues: list = field(repr=False)
    asn_chain: str | None = None
    asn_resseq: int | None = None
    asn_icode: str = " "

    @property
    def site(self) -> str:
        if self.asn_resseq is None:
            return "unlinked"
        return f"{self.asn_chain}:N{self.asn_resseq}{self.asn_icode.strip()}"

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def atoms(self, heavy_only: bool = True) -> list:
        out = []
        for res in self.residues:
            for atom in res:
                if heavy_only and (atom.element or "").strip().upper() == "H":
                    continue
                if atom.get_altloc() not in (" ", "A"):
                    continue
                out.append(atom)
        return out

    def coords_radii(self, heavy_only: bool = True):
        atoms = self.atoms(heavy_only)
        xyz = np.array([a.get_coord() for a in atoms], dtype=np.float64)
        radii = radii_for_elements([a.element for a in atoms])
        return xyz, radii, atoms

    def adjacency(self) -> dict[int, set]:
        """Glycosidic bond graph over ``self.residues``, keyed by list index."""
        n = len(self.residues)
        adj: dict[int, set] = {i: set() for i in range(n)}
        for i in range(n):
            for j in range(i + 1, n):
                if _linked(self.residues[i], self.residues[j]):
                    adj[i].add(j)
                    adj[j].add(i)
        return adj

    def root_index(self) -> int:
        """Index of the reducing-end residue (the one bonded to the Asn).

        Falls back to the residue whose C1 is nearest the Asn ND2, and finally
        to index 0, so the method is total even for unlinked trees.
        """
        if self.asn_resseq is None:
            return 0
        for i, res in enumerate(self.residues):
            if "C1" not in res:
                continue
            parent = res.get_parent().get_parent()
            for chain in parent:
                if chain.id != self.asn_chain:
                    continue
                for cand in chain:
                    if cand.id[1] != self.asn_resseq or "ND2" not in cand:
                        continue
                    d = np.linalg.norm(
                        res["C1"].get_coord() - cand["ND2"].get_coord())
                    if d < ND2_C1_CUTOFF:
                        return i
        return 0

    def bfs_order(self) -> list[int]:
        """Residue indices ordered by bond distance from the reducing end."""
        adj = self.adjacency()
        root = self.root_index()
        order, seen = [], {root}
        queue = [root]
        while queue:
            i = queue.pop(0)
            order.append(i)
            for j in sorted(adj[i]):
                if j not in seen:
                    seen.add(j)
                    queue.append(j)
        for i in range(len(self.residues)):        # disconnected leftovers
            if i not in seen:
                order.append(i)
        return order

    def truncate(self, max_residues: int) -> "GlycanTree":
        """A copy keeping only the ``max_residues`` residues nearest the Asn.

        Used to standardise glycan size across sites. Mean per-atom DeltaRASA
        is size-confounded (see module G2), so comparing a 1-residue stub with
        an 11-residue Man9 is not a like-for-like comparison; truncating every
        tree to a shared core is the cheap controlled alternative to grafting
        an identical modelled glycan onto every site.
        """
        keep = self.bfs_order()[:max_residues]
        return GlycanTree(
            residues=[self.residues[i] for i in keep],
            asn_chain=self.asn_chain,
            asn_resseq=self.asn_resseq,
            asn_icode=self.asn_icode,
        )

    def composition(self) -> dict:
        """Counts per compositional bucket, plus Hex split into Man/Glc/other."""
        comp = {"HexNAc": 0, "Hex": 0, "Fuc": 0, "NeuAc": 0, "Xyl": 0,
                "Man": 0, "Glc": 0}
        for res in self.residues:
            name = res.get_resname().strip().upper()
            bucket = SUGAR_CLASS.get(name)
            if bucket is None:
                continue
            comp[bucket] += 1
            if name in MANNOSE:
                comp["Man"] += 1
            elif name in GLUCOSE:
                comp["Glc"] += 1
        return comp

    def species(self) -> str:
        return classify_composition(**self.composition())[0]

    def ordinal(self) -> int:
        return classify_composition(**self.composition())[1]


def sequons_in_structure(model) -> dict[tuple, str]:
    """Every N-X-S/T sequon (X != Pro) in ``model``, keyed by (chain, resseq).

    Walks author residue numbering rather than a flattened sequence so that
    chain breaks cannot create a sequon that does not exist, and so the keys
    line up with the Asn numbering used by ``find_glycan_trees``.
    """
    from Bio.PDB.Polypeptide import protein_letters_3to1 as three_to_one

    out: dict[tuple, str] = {}
    for chain in model:
        res = {r.id[1]: r.get_resname().strip().upper()
               for r in chain if r.id[0] == " "}
        for num, name in res.items():
            if name != "ASN":
                continue
            x, z = res.get(num + 1), res.get(num + 2)
            if x is None or z is None:
                continue
            x1 = three_to_one.get(x, "X")
            z1 = three_to_one.get(z, "X")
            if x1 == "P" or z1 not in ("S", "T"):
                continue
            out[(chain.id, num)] = "N" + x1 + z1
    return out


def load_model(path: str, model_index: int = 0):
    """Parse a .cif/.pdb file and return one model. Format from the extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".cif", ".mmcif"):
        from Bio.PDB.MMCIFParser import MMCIFParser
        parser = MMCIFParser(QUIET=True)
    else:
        from Bio.PDB.PDBParser import PDBParser
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", path)
    return structure[model_index]


def _heavy_atoms(residue):
    out = []
    for atom in residue:
        if (atom.element or "").strip().upper() == "H":
            continue
        if atom.get_altloc() not in (" ", "A"):
            continue
        out.append(atom)
    return out


def _linked(res_a, res_b) -> bool:
    """True if two sugar residues share a glycosidic bond."""
    a_atoms = _heavy_atoms(res_a)
    b_atoms = _heavy_atoms(res_b)
    if not a_atoms or not b_atoms:
        return False
    a_xyz = np.array([x.get_coord() for x in a_atoms])
    b_xyz = np.array([x.get_coord() for x in b_atoms])
    a_el = [(x.element or "").strip().upper() for x in a_atoms]
    b_el = [(x.element or "").strip().upper() for x in b_atoms]
    d = np.linalg.norm(a_xyz[:, None, :] - b_xyz[None, :, :], axis=2)
    close = np.argwhere(d < GLYCOSIDIC_CUTOFF)
    for i, j in close:
        pair = {a_el[i], b_el[j]}
        if pair == {"C", "O"}:          # glycosidic bond is always C-O
            return True
    return False


def find_glycan_trees(model, min_residues: int = 1) -> list[GlycanTree]:
    """Group every sugar residue in ``model`` into covalent trees.

    Trees are built by union-find over glycosidic C-O contacts, so the result
    does not depend on how the depositor split glycans across chains (mmCIF
    branched entities get their own chain IDs; legacy PDB files do not).
    Each tree is then linked to an Asn via its reducing-end C1...ND2 bond.
    """
    sugars, asn_nd2 = [], []
    for chain in model:
        for res in chain:
            name = res.get_resname().strip().upper()
            if name in SUGAR_CLASS:
                sugars.append(res)
            elif name == "ASN":
                if "ND2" in res:
                    asn_nd2.append((res, res["ND2"].get_coord()))

    n = len(sugars)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    # Only test residue pairs whose centroids are close enough to bond.
    cents = np.array([
        np.mean([a.get_coord() for a in _heavy_atoms(r)] or [[0, 0, 0]], axis=0)
        for r in sugars
    ])
    for i in range(n):
        if not _heavy_atoms(sugars[i]):
            continue
        for j in range(i + 1, n):
            if np.linalg.norm(cents[i] - cents[j]) > 12.0:
                continue
            if _linked(sugars[i], sugars[j]):
                union(i, j)

    groups: dict[int, list] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(sugars[i])

    nd2_xyz = np.array([c for _, c in asn_nd2]) if asn_nd2 else np.zeros((0, 3))

    trees = []
    for members in groups.values():
        if len(members) < min_residues:
            continue
        tree = GlycanTree(residues=members)
        # Reducing-end link: a sugar C1 bonded to an Asn ND2.
        if len(nd2_xyz):
            best = None
            for res in members:
                if "C1" not in res:
                    continue
                c1 = res["C1"].get_coord()
                d = np.linalg.norm(nd2_xyz - c1, axis=1)
                k = int(np.argmin(d))
                if d[k] < ND2_C1_CUTOFF and (best is None or d[k] < best[0]):
                    best = (float(d[k]), k)
            if best is not None:
                asn_res = asn_nd2[best[1]][0]
                tree.asn_chain = asn_res.get_parent().id
                tree.asn_resseq = asn_res.id[1]
                tree.asn_icode = asn_res.id[2]
        trees.append(tree)

    trees.sort(key=lambda t: (t.asn_chain or "~", t.asn_resseq or 0))
    return trees


# ---------------------------------------------------------------------------
# Composition -> species -> ordinal underprocessing score
# ---------------------------------------------------------------------------
#
# RIGOR LABEL FOR THE LADDER BELOW: (warning) reconstructed, not transcribed.
#
# What the paper *states*:
#   - Methods: "Glycans were classified by composition: HexNAc(2)Hex(9-3) as
#     Man9-Man3, fucosylated forms as FM, complex-type glycans by HexNAc count
#     and fucosylation, and truncated species smaller than Man3 as core
#     glycans."                                              -> naming rules
#   - Fig. 3a legend: glycans get "an ordinal underprocessing score (1-20)
#     based on their position along the N-glycan processing continuum. Lower
#     scores correspond to increasingly processed core/complex glycans, while
#     higher scores correspond to oligomannose glycans. ... Highest score (20)
#     represents the least-processed Glc1Man9."              -> endpoints only
#   - Fig. 3b: x-axis "Total Underprocessing Score" spans 0-2000, selection
#     threshold 1500.
#
# What the paper does NOT give anywhere (main text, figures or Extended Data):
#   the species -> integer lookup. So it is rebuilt here from the categorical
#   axis the paper uses for its own site-specific plots (Extended Data Fig. 5b),
#   read in processing order:
#       M9Glc | M9-M4 | M3/FM3 | FM | Hybrid | FHybrid | HexNAc3(F) |
#       HexNAc4 | HexNAc4(F) | HexNAc5 | HexNAc5(F) | HexNAc6(F) | Core
#   Expanding the collapsed "M9-M4" bin into its six species gives exactly
#   twenty rungs, and placing Glc1Man9 at the top reproduces the stated
#   endpoint. Two independent consistency checks fall out of that and both
#   hold:
#       * max total score = 100% x 20 = 2000, matching the Fig. 3b axis;
#       * the 1500 cutoff becomes "abundance-weighted mean rung >= 15", i.e.
#         Man5-or-less-processed on average, which is what O1-O4 are described
#         as carrying.
#   That is good corroboration, not proof. Swap in the published table if it
#   is ever released; ``underprocessing_score`` takes a ``ladder`` override so
#   nothing downstream has to change.

#: species label -> ordinal rung (20 = least processed).
UNDERPROCESSING_LADDER: dict[str, int] = {
    "Glc1Man9": 20,
    "Man9": 19,
    "Man8": 18,
    "Man7": 17,
    "Man6": 16,
    "Man5": 15,
    "Man4": 14,
    "Man3/FM3": 13,
    "FM": 12,
    "Hybrid": 11,
    "FHybrid": 10,
    "HexNAc3": 9,
    "HexNAc3(F)": 8,
    "HexNAc4": 7,
    "HexNAc4(F)": 6,
    "HexNAc5": 5,
    "HexNAc5(F)": 4,
    "HexNAc6": 3,
    "HexNAc6(F)": 2,
    "Core": 1,
}

#: 100% abundance at the top rung.
MAX_UNDERPROCESSING_SCORE = 100.0 * max(UNDERPROCESSING_LADDER.values())


def classify_composition(
    HexNAc: int = 0,
    Hex: int = 0,
    Fuc: int = 0,
    NeuAc: int = 0,
    Xyl: int = 0,
    Man: int | None = None,
    Glc: int = 0,
    ladder: dict[str, int] | None = None,
) -> tuple[str, int]:
    """Map a composition onto ``(species_label, ordinal_rung)``.

    Follows the paper's stated naming rules; the rung comes from
    ``UNDERPROCESSING_LADDER`` (see the rigor note above).
    """
    ladder = ladder or UNDERPROCESSING_LADDER
    fucosylated = Fuc > 0

    # Least-processed species: glucosylated Man9 straight off the OST product.
    if Glc >= 1 and HexNAc == 2 and Hex - Glc >= 9:
        return "Glc1Man9", ladder["Glc1Man9"]

    if HexNAc == 2 and NeuAc == 0:
        if Hex < 3:
            return "Core", ladder["Core"]
        if Hex == 3:
            return "Man3/FM3", ladder["Man3/FM3"]
        if 4 <= Hex <= 9:
            if fucosylated:
                return "FM", ladder["FM"]
            return f"Man{Hex}", ladder[f"Man{Hex}"]
        # >Man9 without a glucose: treat as the top oligomannose rung.
        return "Man9", ladder["Man9"]

    if HexNAc < 2:
        return "Core", ladder["Core"]

    # Hybrid: one GlcNAc on the 3-arm, mannoses retained on the 6-arm.
    if HexNAc == 3 and Hex >= 5 and NeuAc == 0:
        label = "FHybrid" if fucosylated else "Hybrid"
        return label, ladder[label]

    # Complex: graded by antennarity and core fucosylation.
    n = min(max(HexNAc, 3), 6)
    label = f"HexNAc{n}" + ("(F)" if fucosylated else "")
    if label not in ladder:                     # HexNAc6 non-fucosylated etc.
        label = f"HexNAc{n}"
    return label, ladder[label]


def underprocessing_score(
    abundances: dict[str, float],
    occupancy: float | None = None,
    ladder: dict[str, int] | None = None,
    renormalise: bool = True,
) -> dict:
    """Abundance-weighted underprocessing score for one glycosylation site.

    ``abundances`` maps species label -> relative abundance in percent, over
    the *occupied* fraction of the site (the paper plots occupancy separately
    on the y-axis of Fig. 3b). Unknown labels raise, so typos cannot silently
    score zero.

    Returns the total score (0-2000), the abundance-weighted mean rung, and
    the paper's selection verdict (score > 1500 and occupancy > 50%).
    """
    ladder = ladder or UNDERPROCESSING_LADDER
    unknown = sorted(set(abundances) - set(ladder))
    if unknown:
        raise KeyError(f"species not on the ladder: {unknown}")

    total_abundance = float(sum(abundances.values()))
    if renormalise and total_abundance > 0:
        weights = {k: 100.0 * v / total_abundance for k, v in abundances.items()}
    else:
        weights = dict(abundances)

    score = sum(weights[k] * ladder[k] for k in weights)
    mean_rung = score / 100.0 if weights else 0.0

    selected = None
    if occupancy is not None:
        selected = bool(score > 1500.0 and occupancy > 50.0)

    return {
        "total_underprocessing_score": score,
        "max_possible": MAX_UNDERPROCESSING_SCORE,
        "mean_rung": mean_rung,
        "input_abundance_sum": total_abundance,
        "occupancy": occupancy,
        "passes_paper_filter": selected,
    }
