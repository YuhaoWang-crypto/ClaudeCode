#!/usr/bin/env python3
"""
Stage 0 of a LYTAC / EndoTag-style campaign: build a validated, reusable epitope
specification for each lysosome-targeting receptor (LTR).

Target-independent. The output depends only on the receptor, so one run serves any
protein-of-interest you later decide to degrade.

For each receptor this script:

  1. Reads the authoritative auth <-> UniProt residue mapping out of the mmCIF
     header (_struct_ref / _struct_ref_seq) instead of trusting published residue
     numbers. Published LTR epitopes are quoted in several incompatible frames
     (mature protein, precursor, species ortholog); guessing silently designs
     against the wrong surface.
  2. Writes a cleaned single-chain coordinate file -- the exact file a design run
     should upload -- and indexes everything 0-based against *that* file, which is
     the contract the Boltz design API expects.
  3. Locates every native-ligand footprint geometrically (peptide/protein ligand
     chains plus relevant HET groups) so they can be passed as
     `non_binding_residues`. Orthogonality to the native ligand is the design rule
     EndoTag leans on hardest, so it should be a hard constraint, not a hope.
  4. Nominates candidate epitope patches: solvent-exposed, spatially clustered,
     and a safe distance from every native-ligand footprint.

Nomination is not validation. Patches are hypotheses to feed a design run; only a
design + co-fold pass (and ultimately experiment) says whether an interface is real.

Usage:
    python3 ltr_epitope_spec.py --struct-dir <dir-of-cifs> --out-dir specs
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path

import gemmi
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import cKDTree

# Tien et al. 2013, theoretical maximum per-residue SASA (A^2), used to turn the
# absolute Shrake-Rupley number into a relative exposure that is comparable
# between residue types.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# Residues a de novo binder interface is typically built against. Patches made
# only of charged/polar surface rarely yield a designable hydrophobic core.
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "TYR", "PRO"}

EXPOSURE_CUTOFF = 0.25      # relative SASA above which a residue counts as surface
AVOID_RADIUS = 8.0          # A from a native-ligand heavy atom -> do not design here
PATCH_RADIUS = 11.0         # A; approximate footprint of a minibinder interface
MIN_PATCH_SIZE = 8          # residues; smaller patches do not support an interface

# Construct termini and the edges of unmodelled loops score as highly exposed but
# are artefacts: these receptors are fragments of larger proteins (an LTR
# ectodomain continues into a transmembrane region), so a terminus is where the
# crystallographer cut, not a surface a binder can use. Residues flanking a chain
# break sit at the edge of disorder and are equally unreliable.
TERMINAL_EXCLUDE = 5        # residues trimmed from each end of the construct
GAP_FLANK_EXCLUDE = 2       # residues excluded either side of a chain break


@dataclass
class Receptor:
    """One LTR design target: which structure, chain, and native ligands it has."""

    key: str
    pdb_id: str
    chain: str
    description: str
    # Chains that are native ligands (peptide/protein) whose footprint must be avoided.
    ligand_chains: list[str] = field(default_factory=list)
    # HET residue names that mark a native ligand or cofactor site to avoid.
    ligand_het: list[str] = field(default_factory=list)
    # Sub-regions to report separately, in auth numbering: {"D6": (773, 932)}.
    subdomains: dict[str, tuple[int, int]] = field(default_factory=dict)
    # Published epitope residues, with the frame they were quoted in, for cross-checking.
    published_epitope: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)


RECEPTORS = [
    Receptor(
        key="sortilin",
        pdb_id="3F6K",
        chain="A",
        description="Sortilin / SORT1 Vps10p domain (human). CNS-biased LTR; "
                    "best single-domain EndoTag reported (21 nM).",
        ligand_chains=["N"],
        ligand_het=[],
        published_epitope={
            "source": "EndoTag, Nature 2024 -- 'Site 1'",
            "residues_as_published": [92, 93, 546, 559, 561],
            "expected_names": ["PHE", "VAL", "THR", "THR", "THR"],
            "frame": "mature protein (propeptide removed)",
        },
        notes=[
            "4PO7 resolves a SECOND neurotensin site (chains N and P); both are "
            "merged into the avoid set by --extra-ligand-structs.",
            "6X3L carries compound UMJ at the progranulin-interaction site, a third "
            "native-ligand surface worth avoiding.",
        ],
    ),
    Receptor(
        key="asgpr",
        pdb_id="5JQ1",
        chain="A",
        description="ASGPR / ASGR1 carbohydrate-recognition domain (human). "
                    "Liver-restricted LTR; the classic tri-GalNAc receptor.",
        ligand_chains=[],
        ligand_het=["ZPF", "CA"],
        published_epitope={
            "source": "EndoTag, Nature 2024",
            "residues_as_published": None,
            "frame": "not enumerated in the paper text; ASmb1 reported at 2.7 uM",
        },
        notes=[
            "ZPF is the compact galactosamine mimic -> marks the glycan site that a "
            "protein binder must avoid, or it will compete with GalNAc ligands.",
            "CA is the structural Ca2+ of the C-type lectin fold; its site is part of "
            "the glycan-binding machinery.",
            "5JQ1 contains two CRD copies (chains A and B); ASGPR is a trimer in vivo, "
            "so inter-CRD spacing sets the 2C/3C multivalent linker span.",
        ],
    ),
    Receptor(
        key="igf2r_bovine",
        pdb_id="6UM2",
        chain="A",
        description="IGF2R / CI-M6PR full ectodomain WITH IGF2 bound (Bos taurus). "
                    "Ubiquitous LTR. The only structure holding D6 and D11 in one frame.",
        ligand_chains=["B"],
        ligand_het=[],
        subdomains={"D6": (773, 932), "D11": (1523, 1657)},
        published_epitope={
            "source": "EndoTag, Nature 2024 -- D6 minibinder template",
            "residues_as_published": None,
            "frame": "6UM2 coordinates (bovine)",
        },
        notes=[
            "SPECIES TRAP: this is bovine IGF2R (P08169). The paper designed D6 on "
            "this structure but D11 on human 1GP0 (P11717). A D6 binder optimised on "
            "bovine coordinates is not guaranteed to transfer to human.",
            "IGF2 bridges D11 and D6, so IGF2 itself is the natural bivalent "
            "cross-domain ligand -- which is why single-domain engagement fails to "
            "drive uptake.",
        ],
    ),
    Receptor(
        key="igf2r_d11_human",
        pdb_id="1GP0",
        chain="A",
        description="IGF2R domain 11, human (P11717), 1.4 A. Correct-species "
                    "template for the D11 arm.",
        ligand_chains=[],
        ligand_het=[],
        notes=["No native ligand in this crystal; the IGF2 footprint must be "
               "transferred from 6UM2 by sequence mapping."],
    ),
]


def read_uniprot_mapping(cif_path: Path) -> dict:
    """Pull the auth <-> UniProt alignment straight out of the mmCIF header.

    This is the authoritative frame. Deriving it by matching residue identities
    against a UniProt sequence works only when the structure happens to be the
    same species and construct, which is exactly the assumption that silently
    fails.
    """
    block = gemmi.cif.read(str(cif_path)).sole_block()

    refs = {}
    for row in block.find("_struct_ref.", ["id", "db_name", "pdbx_db_accession"]):
        refs[row.str(0)] = {"db": row.str(1), "accession": row.str(2)}

    per_chain: dict[str, list[dict]] = {}
    cols = ["ref_id", "pdbx_strand_id", "pdbx_auth_seq_align_beg",
            "pdbx_auth_seq_align_end", "db_align_beg", "db_align_end"]
    for row in block.find("_struct_ref_seq.", cols):
        ref = refs.get(row.str(0), {})
        auth_beg, db_beg = int(row.str(2)), int(row.str(4))
        per_chain.setdefault(row.str(1), []).append({
            "accession": ref.get("accession"),
            "auth_range": [auth_beg, int(row.str(3))],
            "uniprot_range": [db_beg, int(row.str(5))],
            "auth_minus_uniprot": auth_beg - db_beg,
        })
    return per_chain


def polymer_residues(structure: gemmi.Structure, chain_name: str) -> list:
    """Ordered polymer residues of one chain, waters and HET groups dropped."""
    chain = structure[0][chain_name]
    return [r for r in chain
            if r.het_flag != "H" and r.name != "HOH" and r.name in MAX_ASA]


def write_clean_chain(structure: gemmi.Structure, chain_name: str, out_path: Path) -> list:
    """Write a single-chain, polymer-only coordinate file and return its residues.

    Everything downstream is indexed 0-based against this file, so the indices a
    design run receives cannot drift from the coordinates it is given.
    """
    clean = gemmi.Structure()
    clean.spacegroup_hm = "P 1"
    clean.cell = gemmi.UnitCell()
    model = gemmi.Model("1")
    new_chain = gemmi.Chain(chain_name)
    kept = polymer_residues(structure, chain_name)
    for residue in kept:
        new_chain.add_residue(residue)
    model.add_chain(new_chain)
    clean.add_model(model)
    clean.setup_entities()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    clean.make_mmcif_document().write_file(str(out_path))
    return kept


def relative_sasa(cif_path: Path, chain_name: str, residues: list) -> dict[int, float]:
    """Relative solvent accessibility per residue, computed on the isolated chain.

    Keyed by auth seq id. The isolated chain is the right context: a binder sees
    this chain's surface, and crystallographic neighbours would bury surface that
    is genuinely available in vivo.
    """
    parser = MMCIFParser(QUIET=True)
    model = parser.get_structure("s", str(cif_path))[0]

    class _Keep:
        def accept_model(self, m): return True
        def accept_chain(self, c): return c.id == chain_name
        def accept_residue(self, r): return r.get_resname() in MAX_ASA
        def accept_atom(self, a): return a.element != "H"

    wanted = {r.seqid.num for r in residues}
    chain = model[chain_name]
    for residue in list(chain):
        het, seqid, _ = residue.id
        if het.strip() or seqid not in wanted:
            chain.detach_child(residue.id)

    ShrakeRupley().compute(chain, level="R")
    out = {}
    for residue in chain:
        name = residue.get_resname()
        if name in MAX_ASA:
            out[residue.id[1]] = residue.sasa / MAX_ASA[name]
    return out


def ligand_footprint(structure: gemmi.Structure, receptor: Receptor,
                     target_chain: str, radius: float = AVOID_RADIUS) -> dict:
    """Auth seq ids on the target chain within `radius` of any native ligand atom."""
    model = structure[0]
    ligand_atoms: list[tuple[str, gemmi.Position]] = []

    for chain_name in receptor.ligand_chains:
        if chain_name not in {c.name for c in model}:
            continue
        for residue in model[chain_name]:
            if residue.name == "HOH":
                continue
            for atom in residue:
                ligand_atoms.append((f"chain_{chain_name}", atom.pos))

    for chain in model:
        for residue in chain:
            if residue.name in receptor.ligand_het:
                for atom in residue:
                    ligand_atoms.append((residue.name, atom.pos))

    if not ligand_atoms:
        return {"by_ligand": {}, "union": []}

    residues = polymer_residues(structure, target_chain)
    coords = [(r.seqid.num, np.array([a.pos.tolist() for a in r])) for r in residues]

    by_ligand: dict[str, set[int]] = {}
    for label, pos in ligand_atoms:
        p = np.array(pos.tolist())
        for seqid, atoms in coords:
            if np.linalg.norm(atoms - p, axis=1).min() <= radius:
                by_ligand.setdefault(label, set()).add(seqid)

    union = sorted(set().union(*by_ligand.values())) if by_ligand else []
    return {"by_ligand": {k: sorted(v) for k, v in by_ligand.items()}, "union": union}


def artefact_residues(residues: list) -> set[int]:
    """Auth ids that look exposed only because the construct ends or breaks there."""
    excluded: set[int] = set()
    if len(residues) > 2 * TERMINAL_EXCLUDE:
        excluded |= {r.seqid.num for r in residues[:TERMINAL_EXCLUDE]}
        excluded |= {r.seqid.num for r in residues[-TERMINAL_EXCLUDE:]}
    for i in range(len(residues) - 1):
        a, b = residues[i].seqid.num, residues[i + 1].seqid.num
        if b != a + 1:  # chain break: an unmodelled loop sits between them
            for j in range(max(0, i - GAP_FLANK_EXCLUDE + 1), i + 1):
                excluded.add(residues[j].seqid.num)
            for j in range(i + 1, min(len(residues), i + 1 + GAP_FLANK_EXCLUDE)):
                excluded.add(residues[j].seqid.num)
    return excluded


def nominate_patches(residues: list, rsa: dict[int, float], avoid: set[int],
                     subdomain: tuple[int, int] | None = None,
                     top_n: int = 5) -> list[dict]:
    """Propose designable surface patches: exposed, clustered, clear of ligand sites.

    Each exposed residue seeds a patch of its exposed neighbours within
    PATCH_RADIUS. Patches are scored on exposure, size, and the hydrophobic
    content a designed interface needs, then greedily de-overlapped. Construct
    termini and chain-break edges are excluded first: they reliably top an
    exposure ranking and reliably are not epitopes.
    """
    artefacts = artefact_residues(residues)
    candidates = []
    for residue in residues:
        seqid = residue.seqid.num
        if seqid in avoid or seqid in artefacts:
            continue
        if rsa.get(seqid, 0.0) < EXPOSURE_CUTOFF:
            continue
        if subdomain and not (subdomain[0] <= seqid <= subdomain[1]):
            continue
        ca = residue.find_atom("CA", "*")
        if ca is not None:
            candidates.append((seqid, residue.name, np.array(ca.pos.tolist())))

    if len(candidates) < MIN_PATCH_SIZE:
        return []

    coords = np.array([c[2] for c in candidates])
    tree = cKDTree(coords)

    patches = []
    for i, (seqid, _, centre) in enumerate(candidates):
        members = tree.query_ball_point(centre, PATCH_RADIUS)
        if len(members) < MIN_PATCH_SIZE:
            continue
        ids = [candidates[j][0] for j in members]
        names = [candidates[j][1] for j in members]
        exposure = float(np.mean([rsa[s] for s in ids]))
        hydrophobic = sum(1 for n in names if n in HYDROPHOBIC) / len(names)
        # Reward exposure and size; reward hydrophobic content up to ~40%, which is
        # roughly what a designable interface core needs, and stop rewarding beyond.
        score = exposure * math.log(len(ids)) * (1.0 + min(hydrophobic, 0.4))
        patches.append({
            "seed_residue_auth": seqid,
            "size": len(ids),
            "residues_auth": sorted(ids),
            "mean_relative_sasa": round(exposure, 3),
            "hydrophobic_fraction": round(hydrophobic, 3),
            "score": round(score, 3),
            "centre_xyz": [round(float(x), 1) for x in centre],
        })

    patches.sort(key=lambda p: -p["score"])
    selected: list[dict] = []
    for patch in patches:
        if all(len(set(patch["residues_auth"]) & set(s["residues_auth"]))
               <= 0.4 * min(patch["size"], s["size"]) for s in selected):
            selected.append(patch)
        if len(selected) >= top_n:
            break
    return selected


def check_published_epitope(residues: list, published: dict) -> dict:
    """Scan numbering offsets for the one that reproduces the published residues.

    A published epitope quoted in the wrong frame points at an unrelated surface,
    so the offset is recovered rather than assumed.
    """
    quoted = published.get("residues_as_published")
    if not quoted:
        return {"status": "no_residues_published"}

    expected = published.get("expected_names")
    by_id = {r.seqid.num: r.name for r in residues}
    results = []
    for offset in range(-60, 61):
        shifted = [q + offset for q in quoted]
        present = sum(1 for s in shifted if s in by_id)
        if present != len(quoted):
            continue
        entry = {
            "offset": offset,
            "mapped_auth": shifted,
            "residue_names": [by_id[s] for s in shifted],
        }
        if expected:
            entry["name_matches"] = sum(
                1 for s, e in zip(shifted, expected) if by_id[s] == e)
        results.append(entry)

    if expected:
        results.sort(key=lambda r: -r.get("name_matches", 0))
    return {"status": "scanned", "candidates": results[:4],
            "quoted_frame": published.get("frame")}


JUNCTION_RADIUS = 12.0  # A; a residue this close to the partner subdomain sits at
                        # the inter-domain junction a bivalent binder would clamp


def clamp_analysis(structure: gemmi.Structure, chain: str,
                   subdomains: dict[str, tuple[int, int]],
                   ligand_footprint_by_subdomain: dict[str, list[int]],
                   rsa: dict[int, float]) -> dict:
    """Geometry of a cross-domain clamp, and whether an orthogonal one can exist.

    When a native ligand bridges two subdomains it is itself the proven bivalent
    binder, so its two footprints -- not the best-exposed patch on each subdomain
    independently -- give the span a bivalent construct has to cover. Picking the
    best patch per subdomain can land on opposite outer faces and overstate the
    span several-fold.

    The junction region and the native-ligand footprint tend to coincide, which is
    the central design tension: the validated geometry is also the competitive one.
    This quantifies the overlap so the trade-off is a number rather than a guess.
    """
    residues = [r for r in polymer_residues(structure, chain)]
    by_id = {r.seqid.num: r for r in residues}

    def atoms_of(ids):
        pts = [a.pos.tolist() for i in ids if i in by_id for a in by_id[i]]
        return np.array(pts) if pts else np.empty((0, 3))

    def span(ids_a, ids_b):
        pa, pb = atoms_of(ids_a), atoms_of(ids_b)
        if len(pa) == 0 or len(pb) == 0:
            return None
        dist = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)
        centroid = float(np.linalg.norm(pa.mean(0) - pb.mean(0)))
        return {
            "centroid_distance_A": round(centroid, 1),
            "closest_approach_A": round(float(dist.min()), 1),
            "farthest_A": round(float(dist.max()), 1),
            "gs_linker_residues_to_span": int(math.ceil(centroid / 3.5)),
        }

    out: dict = {}
    names = list(subdomains)
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            (s1, e1), (s2, e2) = subdomains[first], subdomains[second]
            ids1 = [r.seqid.num for r in residues if s1 <= r.seqid.num <= e1]
            ids2 = [r.seqid.num for r in residues if s2 <= r.seqid.num <= e2]
            if not ids1 or not ids2:
                continue

            # Junction: residues of each subdomain close to the partner subdomain.
            p1, p2 = atoms_of(ids1), atoms_of(ids2)
            t2, t1 = cKDTree(p2), cKDTree(p1)
            junction1 = sorted(i_ for i_ in ids1
                               if t2.query(np.array([a.pos.tolist() for a in by_id[i_]]))[0].min()
                               <= JUNCTION_RADIUS)
            junction2 = sorted(i_ for i_ in ids2
                               if t1.query(np.array([a.pos.tolist() for a in by_id[i_]]))[0].min()
                               <= JUNCTION_RADIUS)

            fp1 = set(ligand_footprint_by_subdomain.get(first, []))
            fp2 = set(ligand_footprint_by_subdomain.get(second, []))

            def orthogonal_exposed(junction, footprint):
                return sorted(i_ for i_ in junction
                              if i_ not in footprint
                              and rsa.get(i_, 0.0) >= EXPOSURE_CUTOFF)

            orth1 = orthogonal_exposed(junction1, fp1)
            orth2 = orthogonal_exposed(junction2, fp2)

            entry = {
                "subdomain_span": span(ids1, ids2),
                "native_ligand_clamp": {
                    "note": "span between the two native-ligand footprints: the "
                            "proven bivalent geometry",
                    f"{first}_footprint": sorted(fp1),
                    f"{second}_footprint": sorted(fp2),
                    "span": span(sorted(fp1), sorted(fp2)) if fp1 and fp2 else None,
                },
                "junction": {
                    "radius_A": JUNCTION_RADIUS,
                    f"{first}_junction_residues": junction1,
                    f"{second}_junction_residues": junction2,
                    "overlap_with_native_ligand": {
                        first: {
                            "n_junction": len(junction1),
                            "n_also_footprint": len(set(junction1) & fp1),
                            "fraction": round(len(set(junction1) & fp1) / len(junction1), 3)
                            if junction1 else None,
                        },
                        second: {
                            "n_junction": len(junction2),
                            "n_also_footprint": len(set(junction2) & fp2),
                            "fraction": round(len(set(junction2) & fp2) / len(junction2), 3)
                            if junction2 else None,
                        },
                    },
                },
                "orthogonal_clamp_option": {
                    "note": "junction residues that are exposed AND outside the "
                            "native-ligand footprint: where a non-competitive "
                            "bivalent binder would have to bind. Whether enough "
                            "surface survives is threshold-dependent, so the "
                            "sweep below is the result -- not a single verdict.",
                    "at_default_thresholds": {
                        "junction_radius_A": JUNCTION_RADIUS,
                        "exposure_cutoff": EXPOSURE_CUTOFF,
                        f"{first}_residues": orth1,
                        f"{second}_residues": orth2,
                        "span": span(orth1, orth2) if orth1 and orth2 else None,
                    },
                    "limiting_subdomain": first if len(orth1) < len(orth2) else second,
                },
            }
            out[f"{first}__{second}"] = entry
    return out


def clamp_threshold_sweep(structure: gemmi.Structure, receptor: Receptor,
                          rsa: dict[int, float]) -> dict:
    """Does an orthogonal cross-domain clamp survive a range of threshold choices?

    How much surface counts as "clear of the native ligand" depends on three
    arbitrary numbers: how far the footprint is dilated, how wide the junction is,
    and how exposed a residue must be. A conclusion that flips across reasonable
    settings is a conclusion about the thresholds, so the sweep is reported and
    the caller decides.
    """
    global JUNCTION_RADIUS, EXPOSURE_CUTOFF
    saved = (JUNCTION_RADIUS, EXPOSURE_CUTOFF)
    rows: dict[str, list[dict]] = {}
    try:
        for fp_radius in (4.5, 6.0, 8.0):
            probe = ligand_footprint(structure, receptor, receptor.chain, radius=fp_radius)
            blocked = set(probe["union"])
            by_sub = {name: sorted(a for a in blocked if rng[0] <= a <= rng[1])
                      for name, rng in receptor.subdomains.items()}
            for junction_radius in (8.0, 12.0, 15.0, 20.0):
                for exposure in (0.15, 0.25):
                    JUNCTION_RADIUS, EXPOSURE_CUTOFF = junction_radius, exposure
                    for pair, entry in clamp_analysis(
                            structure, receptor.chain, receptor.subdomains,
                            by_sub, rsa).items():
                        default = entry["orthogonal_clamp_option"]["at_default_thresholds"]
                        counts = {k.replace("_residues", ""): len(v)
                                  for k, v in default.items()
                                  if k.endswith("_residues")}
                        rows.setdefault(pair, []).append({
                            "footprint_dilation_A": fp_radius,
                            "junction_radius_A": junction_radius,
                            "exposure_cutoff": exposure,
                            "orthogonal_residues": counts,
                            "min_side": min(counts.values()) if counts else 0,
                        })
    finally:
        JUNCTION_RADIUS, EXPOSURE_CUTOFF = saved

    summary = {}
    for pair, entries in rows.items():
        mins = [e["min_side"] for e in entries]
        summary[pair] = {
            "n_settings_tested": len(entries),
            "limiting_side_residue_count": {"min": min(mins), "max": max(mins)},
            "fraction_of_settings_with_at_least_4": round(
                sum(1 for m in mins if m >= 4) / len(mins), 2),
            "rows": entries,
        }
    return summary


def multivalency_geometry(structure: gemmi.Structure, chain_a: str, chain_b: str) -> dict | None:
    """Spacing between two copies of a receptor chain -> multivalent linker span.

    EndoTag's ASGPR gain came from 2- and 3-valent constructs, and its IGF2R
    constructs died outside a narrow linker window, so this span is a design
    input rather than a detail.

    Only meaningful between two copies of the SAME entity. A ligand chain sitting
    in slot B would otherwise be reported as a second receptor copy.
    """
    names = {c.name for c in structure[0]}
    if chain_a not in names or chain_b not in names:
        return None
    seq_a = gemmi.one_letter_code([r.name for r in polymer_residues(structure, chain_a)])
    seq_b = gemmi.one_letter_code([r.name for r in polymer_residues(structure, chain_b)])
    if not seq_a or not seq_b:
        return None
    shorter, longer = sorted((seq_a, seq_b), key=len)
    if len(shorter) < 0.5 * len(longer):
        return None  # different entities, not two receptor copies
    pa = np.array([a.pos.tolist() for r in polymer_residues(structure, chain_a) for a in r])
    pb = np.array([a.pos.tolist() for r in polymer_residues(structure, chain_b) for a in r])
    if len(pa) == 0 or len(pb) == 0:
        return None
    dist = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)
    return {
        "chains": [chain_a, chain_b],
        "centroid_distance_A": round(float(np.linalg.norm(pa.mean(0) - pb.mean(0))), 1),
        "closest_approach_A": round(float(dist.min()), 1),
        "farthest_A": round(float(dist.max()), 1),
    }


def subdomain_geometry(structure: gemmi.Structure, chain: str,
                       subdomains: dict[str, tuple[int, int]],
                       patches: dict[str, list[dict]]) -> dict:
    """Distances between subdomains, and between their best candidate patches.

    The patch-to-patch distance is the number a rigid or linked two-domain
    construct actually has to span -- the quantity EndoTag found empirically by
    screening linker lengths.
    """
    residues = polymer_residues(structure, chain)
    out: dict = {}
    names = list(subdomains)
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            (s1, e1), (s2, e2) = subdomains[first], subdomains[second]
            p1 = np.array([a.pos.tolist() for r in residues if s1 <= r.seqid.num <= e1 for a in r])
            p2 = np.array([a.pos.tolist() for r in residues if s2 <= r.seqid.num <= e2 for a in r])
            if len(p1) == 0 or len(p2) == 0:
                continue
            dist = np.linalg.norm(p1[:, None, :] - p2[None, :, :], axis=2)
            entry = {
                "centroid_distance_A": round(float(np.linalg.norm(p1.mean(0) - p2.mean(0))), 1),
                "closest_approach_A": round(float(dist.min()), 1),
                "farthest_A": round(float(dist.max()), 1),
            }
            best1 = (patches.get(first) or [{}])[0].get("centre_xyz")
            best2 = (patches.get(second) or [{}])[0].get("centre_xyz")
            if best1 and best2:
                span = float(np.linalg.norm(np.array(best1) - np.array(best2)))
                entry["best_patch_centre_distance_A"] = round(span, 1)
                # A Gly-Ser linker contributes roughly 3.5 A of span per residue
                # when extended; this is an order-of-magnitude floor, not a design.
                entry["gs_linker_residues_to_span"] = int(math.ceil(span / 3.5))
                entry["caveat"] = ("independently best-scoring patches can sit on "
                                   "opposite outer faces and overstate the span; "
                                   "prefer the clamp_analysis spans")
            out[f"{first}__{second}"] = entry
    return out


def analyse(receptor: Receptor, struct_dir: Path, out_dir: Path,
            extra_ligand_structs: dict[str, list[str]] | None = None) -> dict:
    cif_path = struct_dir / f"{receptor.pdb_id}.cif"
    structure = gemmi.read_structure(str(cif_path))
    structure.setup_entities()

    clean_path = out_dir.parent / "structures" / f"{receptor.key}_{receptor.pdb_id}_{receptor.chain}.cif"
    residues = write_clean_chain(structure, receptor.chain, clean_path)
    index_of = {r.seqid.num: i for i, r in enumerate(residues)}

    rsa = relative_sasa(cif_path, receptor.chain, residues)
    footprint = ligand_footprint(structure, receptor, receptor.chain)
    avoid = set(footprint["union"])

    # Fold in native-ligand sites that only appear in other crystal forms of the
    # same receptor. A site absent from one structure is still a site.
    cross_structure: dict[str, dict] = {}
    for pdb_id, ligand_spec in (extra_ligand_structs or {}).items():
        extra_path = struct_dir / f"{pdb_id}.cif"
        if not extra_path.exists():
            continue
        extra = gemmi.read_structure(str(extra_path))
        extra.setup_entities()
        chains = [c for c in ligand_spec if c in {ch.name for ch in extra[0]}]
        hets = [h for h in ligand_spec if h not in chains]
        probe = Receptor(key=pdb_id, pdb_id=pdb_id, chain=receptor.chain,
                         description="", ligand_chains=chains, ligand_het=hets)
        extra_fp = ligand_footprint(extra, probe, receptor.chain)
        # Same UniProt entry and same offset across these entries, so auth ids
        # transfer -- but the entries do not cover the same residue range, and the
        # primary structure may have chain breaks the others lack. Keep only ids
        # this structure actually resolves, and record what was dropped: an auth id
        # with no coordinates here cannot be given a 0-based index, and silently
        # carrying it desynchronises the auth and index lists.
        transferable = {a for a in extra_fp["union"] if a in index_of}
        dropped = sorted(set(extra_fp["union"]) - transferable)
        cross_structure[pdb_id] = {
            "by_ligand": extra_fp["by_ligand"],
            "transferred_auth": sorted(transferable),
            "dropped_not_resolved_here": dropped,
        }
        avoid |= transferable

    scopes = receptor.subdomains or {"whole_chain": (
        residues[0].seqid.num, residues[-1].seqid.num)}
    patches = {name: nominate_patches(residues, rsa, avoid, rng)
               for name, rng in scopes.items()}

    # How sensitive is the avoid set to the cutoff? A small domain crowded with
    # ions can lose most of its surface at a generous radius, so report the curve
    # instead of leaving one hard-coded number unexamined.
    exposed = {s for s, v in rsa.items() if v >= EXPOSURE_CUTOFF}
    sensitivity = {}
    for radius in (4.0, 5.0, 6.0, 8.0, 10.0):
        probe = ligand_footprint(structure, receptor, receptor.chain, radius=radius)
        blocked = set(probe["union"])
        sensitivity[f"{radius:g}A"] = {
            "n_avoided": len(blocked),
            "n_exposed_still_free": len(exposed - blocked),
        }

    avoid_auth = sorted(a for a in avoid if a in index_of)
    footprint_by_subdomain = {
        name: [a for a in avoid_auth if rng[0] <= a <= rng[1]]
        for name, rng in receptor.subdomains.items()
    }
    clamp = (clamp_analysis(structure, receptor.chain, receptor.subdomains,
                            footprint_by_subdomain, rsa)
             if len(receptor.subdomains) >= 2 else {})
    if clamp:
        for pair, sweep in clamp_threshold_sweep(
                structure, receptor, rsa).items():
            clamp[pair]["orthogonal_clamp_option"]["threshold_sweep"] = sweep

        # The clamp span above is measured between 8 A-dilated footprints. The
        # undilated contact footprint is the tighter, more literal answer to "how
        # far apart are the two halves of the native bivalent ligand", so report
        # both rather than leaving the quoted span dependent on the dilation.
        raw = ligand_footprint(structure, receptor, receptor.chain, radius=4.5)
        raw_by_sub = {name: [a for a in sorted(raw["union"]) if rng[0] <= a <= rng[1]]
                      for name, rng in receptor.subdomains.items()}
        for pair, entry in clamp_analysis(structure, receptor.chain,
                                          receptor.subdomains, raw_by_sub, rsa).items():
            if pair in clamp:
                clamp[pair]["native_ligand_clamp_raw_contacts"] = {
                    "contact_radius_A": 4.5,
                    **entry["native_ligand_clamp"],
                }

    spec = {
        "receptor": receptor.key,
        "description": receptor.description,
        "structure": {
            "pdb_id": receptor.pdb_id,
            "chain": receptor.chain,
            "clean_file": str(clean_path.relative_to(out_dir.parent)),
            "n_residues": len(residues),
            "auth_range": [residues[0].seqid.num, residues[-1].seqid.num],
        },
        "numbering": {
            "from_mmcif_header": read_uniprot_mapping(cif_path).get(receptor.chain),
            "warning": "0-based indices below are relative to clean_file, which is "
                       "the file a design run must upload. They are NOT auth numbers "
                       "and NOT UniProt numbers.",
        },
        "published_epitope_check": check_published_epitope(residues, receptor.published_epitope),
        "native_ligand_footprint": {
            "radius_A": AVOID_RADIUS,
            "this_structure": footprint["by_ligand"],
            "other_structures": cross_structure,
            # Both lists are derived from one filtered, sorted sequence so that
            # element i of each always refers to the same residue. Building them
            # independently lets a single unresolved id shift one list against the
            # other and silently relocate the whole constraint.
            "union_auth": avoid_auth,
            "union_0based": [index_of[a] for a in avoid_auth],
            "radius_sensitivity": sensitivity,
            "by_subdomain": footprint_by_subdomain,
        },
        "clamp_analysis": clamp,
        "candidate_epitopes": {
            name: [
                {**patch,
                 "residues_0based": [index_of[a] for a in patch["residues_auth"] if a in index_of]}
                for patch in patch_list
            ]
            for name, patch_list in patches.items()
        },
        "geometry": {
            "subdomains": subdomain_geometry(structure, receptor.chain,
                                             receptor.subdomains, patches)
            if receptor.subdomains else {},
            "inter_chain": multivalency_geometry(structure, "A", "B"),
        },
        "surface_summary": {
            "n_exposed": sum(1 for v in rsa.values() if v >= EXPOSURE_CUTOFF),
            "n_avoided": len(avoid_auth),
            "n_excluded_as_artefact": len(artefact_residues(residues)),
            "artefact_residues_auth": sorted(artefact_residues(residues)),
            "exposure_cutoff": EXPOSURE_CUTOFF,
            "terminal_exclude": TERMINAL_EXCLUDE,
            "gap_flank_exclude": GAP_FLANK_EXCLUDE,
        },
        "notes": receptor.notes,
    }
    return spec


def validate_spec(spec: dict, out_root: Path) -> None:
    """Re-read the written coordinate file and assert every 0-based index resolves.

    The 0-based indices are the entire contract with a design run: they decide
    which residues become the epitope and which become non-binding. An off-by-one
    relocates a constraint onto an unrelated surface and the run still succeeds,
    reporting scores for a design against the wrong site. So the indices are
    checked against the file as written, not against the arrays in memory.
    """
    structure = gemmi.read_structure(str(out_root / spec["structure"]["clean_file"]))
    structure.setup_entities()
    residues = polymer_residues(structure, spec["structure"]["chain"])

    n = spec["structure"]["n_residues"]
    if len(residues) != n:
        raise AssertionError(
            f"{spec['receptor']}: clean file has {len(residues)} residues, spec says {n}")
    if [residues[0].seqid.num, residues[-1].seqid.num] != spec["structure"]["auth_range"]:
        raise AssertionError(f"{spec['receptor']}: auth range mismatch")

    def check(auth_list, index_list, label):
        if len(auth_list) != len(index_list):
            raise AssertionError(
                f"{spec['receptor']} {label}: {len(auth_list)} auth ids vs "
                f"{len(index_list)} indices -- lists are desynchronised")
        for auth, idx in zip(auth_list, index_list):
            if not 0 <= idx < len(residues):
                raise AssertionError(f"{spec['receptor']} {label}: index {idx} out of range")
            if residues[idx].seqid.num != auth:
                raise AssertionError(
                    f"{spec['receptor']} {label}: index {idx} is auth "
                    f"{residues[idx].seqid.num}, expected {auth}")
        return len(auth_list)

    checked = check(spec["native_ligand_footprint"]["union_auth"],
                    spec["native_ligand_footprint"]["union_0based"],
                    "non_binding_residues")
    for scope, patches in spec["candidate_epitopes"].items():
        for i, patch in enumerate(patches):
            checked += check(patch["residues_auth"], patch["residues_0based"],
                             f"epitope {scope}[{i}]")
    spec["validation"] = {"index_assertions_passed": checked, "status": "ok"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--struct-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=Path("specs"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Native-ligand sites for sortilin that only show up in other crystal forms:
    # the second neurotensin site (4PO7 chain P) and the progranulin-site compound
    # UMJ (6X3L). All three entries share Q99523 with the same offset.
    extras = {
        "sortilin": {"4PO7": ["N", "P"], "6X3L": ["UMJ"]},
    }

    specs = {}
    for receptor in RECEPTORS:
        print(f"\n{'=' * 72}\n{receptor.key}  ({receptor.pdb_id} chain {receptor.chain})")
        spec = analyse(receptor, args.struct_dir, args.out_dir,
                       extras.get(receptor.key))
        validate_spec(spec, args.out_dir.parent)
        specs[receptor.key] = spec
        print(f"  validated: {spec['validation']['index_assertions_passed']} "
              f"index assertions passed")

        mapping = spec["numbering"]["from_mmcif_header"]
        if mapping:
            m = mapping[0]
            print(f"  UniProt {m['accession']}  auth-minus-uniprot = {m['auth_minus_uniprot']:+d}")
        print(f"  residues {spec['structure']['n_residues']}  "
              f"exposed {spec['surface_summary']['n_exposed']}  "
              f"avoided {spec['surface_summary']['n_avoided']}")

        check = spec["published_epitope_check"]
        for cand in check.get("candidates", [])[:2]:
            match = cand.get("name_matches")
            tail = f"  names {match}/{len(cand['mapped_auth'])}" if match is not None else ""
            print(f"  published epitope offset {cand['offset']:+d} -> "
                  f"{cand['mapped_auth']} {cand['residue_names']}{tail}")

        for scope, patch_list in spec["candidate_epitopes"].items():
            if not patch_list:
                print(f"  {scope}: no patch passed the filters")
            for i, patch in enumerate(patch_list[:3]):
                print(f"  {scope} patch {i}: seed {patch['seed_residue_auth']} "
                      f"size {patch['size']} rSASA {patch['mean_relative_sasa']} "
                      f"phobic {patch['hydrophobic_fraction']} score {patch['score']}")

        for pair, geom in spec["geometry"]["subdomains"].items():
            print(f"  geometry {pair}: centroid {geom['centroid_distance_A']} A, "
                  f"closest {geom['closest_approach_A']} A"
                  + (f"   [naive best-patch span {geom['best_patch_centre_distance_A']} A "
                     f"-- see clamp_analysis]"
                     if "best_patch_centre_distance_A" in geom else ""))

        for pair, clamp in spec.get("clamp_analysis", {}).items():
            native = clamp["native_ligand_clamp"]["span"]
            if native:
                print(f"  CLAMP {pair}: native-ligand footprints "
                      f"{native['centroid_distance_A']} A apart "
                      f"(closest {native['closest_approach_A']} A) "
                      f"-> ~{native['gs_linker_residues_to_span']} GS residues")
            for dom, ov in clamp["junction"]["overlap_with_native_ligand"].items():
                print(f"    {dom} junction {ov['n_junction']} res, "
                      f"{ov['n_also_footprint']} also native-ligand "
                      f"({ov['fraction']})")
            orth = clamp["orthogonal_clamp_option"]
            default = orth["at_default_thresholds"]
            span = default["span"]
            counts = {k.replace("_residues", ""): len(v)
                      for k, v in default.items() if k.endswith("_residues")}
            print(f"    orthogonal clamp surface at default thresholds: {counts}"
                  + (f", span {span['centroid_distance_A']} A "
                     f"-> ~{span['gs_linker_residues_to_span']} GS residues"
                     if span else ""))
            sweep = orth.get("threshold_sweep")
            if sweep:
                lim = sweep["limiting_side_residue_count"]
                print(f"    limiting side ({orth['limiting_subdomain']}) across "
                      f"{sweep['n_settings_tested']} threshold settings: "
                      f"{lim['min']}-{lim['max']} residues; "
                      f"{sweep['fraction_of_settings_with_at_least_4']:.0%} of settings "
                      f"leave >=4 -> THRESHOLD-DEPENDENT, not a clean verdict")
        if spec["geometry"]["inter_chain"]:
            ic = spec["geometry"]["inter_chain"]
            print(f"  inter-chain {ic['chains']}: centroid {ic['centroid_distance_A']} A, "
                  f"closest {ic['closest_approach_A']} A")

        out_file = args.out_dir / f"{receptor.key}.json"
        out_file.write_text(json.dumps(spec, indent=2))
        print(f"  -> {out_file}")

    (args.out_dir / "all_receptors.json").write_text(json.dumps(specs, indent=2))
    print(f"\nwrote {args.out_dir / 'all_receptors.json'}")


if __name__ == "__main__":
    main()
