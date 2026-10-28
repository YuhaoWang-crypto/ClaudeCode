"""Structure layer: fetch a TCR-pMHC complex, type its chains, locate the CDR
loops and measure the interface geometry.

Everything in this module is deterministic measurement of a deposited
structure. No energy model and no prediction is involved, so the numbers it
produces are as good as the crystal structure it was given.
"""

from __future__ import annotations

import math
import re
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

from .constants import (
    BACKBONE_ATOMS,
    CHARGED_SIDECHAIN_ATOMS,
    MAX_SASA,
    POLAR_SIDECHAIN_ATOMS,
)

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# Sequence fingerprints used to type chains. Each is a stretch of constant
# domain or conserved framework that is effectively invariant across alleles.
FINGERPRINTS = {
    "mhc_i": ("GSHSMRY", "SHSMRYF", "GSHSLRY", "SHSLKYF"),
    "b2m": ("IQRTPKIQVYSRHPAEN", "RTPKIQVYSRHPAEN", "SRHPAENGKSNFL"),
    "trac": ("DIQNPDPAVYQL", "IQNPDPAVYQ", "NIQNPDPAVYQL"),
    "trbc": ("VFPPEVAVF", "DLKNVFPPEV", "EDLNKVFPPEV"),
}

# Canonical class I helices that form the walls of the peptide groove. Used as
# the geometric fallback for telling TCR alpha from TCR beta: in the conserved
# diagonal docking mode alpha sits over the alpha2 helix and beta over alpha1.
MHC_A1_HELIX = range(57, 85)
MHC_A2_HELIX = range(138, 181)


def _one_letter(residue) -> str:
    return protein_letters_3to1.get(residue.get_resname(), "X")


def fetch_pdb(pdb_id: str, cache_dir: Path) -> Path:
    """Download a PDB entry, caching it under ``cache_dir``."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb_id.lower()}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path
    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read()
    path.write_bytes(data)
    return path


@dataclass
class Residue:
    """A single amino acid with the measurements the scoring layer needs."""

    chain: str
    number: int
    insertion: str
    aa: str
    index: int  # position within its chain's extracted sequence
    loop: str | None = None  # CDR1/CDR2/CDR3/HV4 where assigned

    sasa_complex: float = 0.0
    sasa_free: float = 0.0
    contacts_peptide: int = 0
    contacts_mhc: int = 0
    contacts_self_chain: int = 0
    nearest_partner_distance: float = math.inf
    void_slots: int = 0
    growth_reaches_partner: bool = False
    partner_charge: float = 0.0
    partner_polar_atoms: int = 0
    backbone_nh_bonded: bool = False
    phi: float | None = None
    in_disulfide: bool = False

    @property
    def label(self) -> str:
        return f"{self.chain}:{self.aa}{self.number}{self.insertion.strip()}"

    @property
    def contacts_partner(self) -> int:
        return self.contacts_peptide + self.contacts_mhc

    @property
    def buried_fraction(self) -> float:
        """Fraction of this residue's accessible surface lost on binding."""
        if self.sasa_free <= 0.0:
            return 0.0
        return max(0.0, (self.sasa_free - self.sasa_complex) / self.sasa_free)

    @property
    def relative_sasa_free(self) -> float:
        return self.sasa_free / MAX_SASA.get(self.aa, 200.0)

    @property
    def mhc_contact_fraction(self) -> float:
        """Share of partner contacts made to MHC rather than to the peptide.

        This is the specificity dial. Affinity bought at MHC-contacting
        positions increases peptide-independent recognition, which is the
        documented route to off-target toxicity.
        """
        total = self.contacts_partner
        if total == 0:
            return 0.0
        return self.contacts_mhc / total


@dataclass
class Complex:
    """A typed TCR-pMHC complex."""

    pdb_id: str
    path: Path
    chains: dict[str, str]  # role -> chain id
    sequences: dict[str, str]  # chain id -> one-letter sequence
    residues: list[Residue]
    loops: dict[str, tuple[str, int, int]] = field(default_factory=dict)
    typing_notes: list[str] = field(default_factory=list)

    @property
    def peptide(self) -> str:
        return self.sequences[self.chains["peptide"]]

    def tcr_chain_ids(self) -> list[str]:
        return [self.chains[r] for r in ("tcr_alpha", "tcr_beta") if r in self.chains]

    def partner_chain_ids(self) -> list[str]:
        roles = ("mhc", "b2m", "peptide")
        return [self.chains[r] for r in roles if r in self.chains]

    def paratope(self, min_contacts: int = 1) -> list[Residue]:
        """TCR residues that actually touch the pMHC surface."""
        tcr = set(self.tcr_chain_ids())
        hits = [
            r for r in self.residues
            if r.chain in tcr and r.contacts_partner >= min_contacts
        ]
        return sorted(hits, key=lambda r: (r.chain, r.number))

    def residue_by_label(self, chain: str, number: int) -> Residue | None:
        for r in self.residues:
            if r.chain == chain and r.number == number:
                return r
        return None


def _extract_chain_sequences(model) -> tuple[dict[str, str], dict[str, list]]:
    sequences: dict[str, str] = {}
    residue_objects: dict[str, list] = {}
    for chain in model:
        keep = []
        for residue in chain:
            # Standard amino acids only: skip waters, ions and ligands.
            if residue.id[0] != " ":
                continue
            if residue.get_resname() not in protein_letters_3to1:
                continue
            if "CA" not in residue:
                continue
            keep.append(residue)
        if not keep:
            continue
        sequences[chain.id] = "".join(_one_letter(r) for r in keep)
        residue_objects[chain.id] = keep
    return sequences, residue_objects


def _match_fingerprint(sequence: str, key: str) -> bool:
    return any(motif in sequence for motif in FINGERPRINTS[key])


def type_chains(
    sequences: dict[str, str],
    residue_objects: dict[str, list],
    override: dict[str, str] | None = None,
) -> tuple[dict[str, str], list[str]]:
    """Assign each chain a role in the complex.

    Sequence fingerprints are tried first. TCR alpha/beta assignment falls back
    to the conserved diagonal docking geometry when the constant domains are
    absent from the model, which happens in single-domain constructs.
    """
    notes: list[str] = []
    roles: dict[str, str] = {}
    override = override or {}

    peptide_candidates = [c for c, s in sequences.items() if 7 <= len(s) <= 20]
    tcr_candidates: list[str] = []

    for chain_id, seq in sequences.items():
        if chain_id in peptide_candidates:
            continue
        if _match_fingerprint(seq, "b2m") or (90 <= len(seq) <= 110 and "b2m" not in roles):
            if _match_fingerprint(seq, "b2m"):
                roles["b2m"] = chain_id
                notes.append(f"chain {chain_id} -> beta-2 microglobulin (sequence fingerprint)")
                continue
        if _match_fingerprint(seq, "mhc_i"):
            roles["mhc"] = chain_id
            notes.append(f"chain {chain_id} -> MHC class I heavy chain (sequence fingerprint)")
            continue
        if _match_fingerprint(seq, "trac"):
            roles["tcr_alpha"] = chain_id
            notes.append(f"chain {chain_id} -> TCR alpha (TRAC constant-domain fingerprint)")
            continue
        if _match_fingerprint(seq, "trbc"):
            roles["tcr_beta"] = chain_id
            notes.append(f"chain {chain_id} -> TCR beta (TRBC constant-domain fingerprint)")
            continue
        tcr_candidates.append(chain_id)

    # Chains that look like an immunoglobulin V domain but had no constant
    # domain to fingerprint.
    leftover_v = [c for c in tcr_candidates if re.search(r"FG.G", sequences[c])]
    for chain_id in leftover_v:
        if "tcr_alpha" not in roles:
            roles["tcr_alpha"] = chain_id
        elif "tcr_beta" not in roles:
            roles["tcr_beta"] = chain_id

    if leftover_v and "mhc" in roles:
        notes.append(
            "TCR alpha/beta resolved by docking geometry: "
            + ", ".join(sorted(leftover_v))
        )
        roles = _orient_tcr_chains_by_geometry(roles, residue_objects, leftover_v, notes)

    if "mhc" not in roles:
        # Fall back to the longest chain, which for class I complexes is the
        # heavy chain.
        longest = max(sequences, key=lambda c: len(sequences[c]))
        roles["mhc"] = longest
        notes.append(f"chain {longest} -> MHC (fallback: longest chain)")

    if peptide_candidates:
        # The peptide is the short chain with the most contacts to the MHC.
        mhc_atoms = _heavy_atoms(residue_objects[roles["mhc"]])
        best, best_count = None, -1
        for chain_id in peptide_candidates:
            count = _count_close_atoms(_heavy_atoms(residue_objects[chain_id]), mhc_atoms, 4.5)
            if count > best_count:
                best, best_count = chain_id, count
        roles["peptide"] = best
        notes.append(
            f"chain {best} -> peptide ({len(sequences[best])} aa, "
            f"{best_count} heavy-atom contacts to MHC)"
        )

    roles.update(override)
    for role, chain_id in override.items():
        notes.append(f"chain {chain_id} -> {role} (user override)")
    return roles, notes


def _orient_tcr_chains_by_geometry(
    roles: dict[str, str],
    residue_objects: dict[str, list],
    candidates: list[str],
    notes: list[str],
) -> dict[str, str]:
    """Resolve TCR alpha vs beta from the canonical docking polarity."""
    if len(candidates) != 2 or "mhc" not in roles:
        return roles
    mhc_residues = residue_objects[roles["mhc"]]
    a1 = _heavy_atoms([r for r in mhc_residues if r.id[1] in MHC_A1_HELIX])
    a2 = _heavy_atoms([r for r in mhc_residues if r.id[1] in MHC_A2_HELIX])
    scores = {}
    for chain_id in candidates:
        atoms = _heavy_atoms(residue_objects[chain_id])
        n_a1 = _count_close_atoms(atoms, a1, 5.0)
        n_a2 = _count_close_atoms(atoms, a2, 5.0)
        scores[chain_id] = n_a2 - n_a1
    # Larger alpha2 preference means TCR alpha.
    ranked = sorted(candidates, key=lambda c: scores[c], reverse=True)
    roles["tcr_alpha"], roles["tcr_beta"] = ranked[0], ranked[1]
    notes.append(
        "docking polarity: "
        + ", ".join(f"{c} alpha2-minus-alpha1 contacts = {scores[c]:+d}" for c in ranked)
    )
    return roles


def _heavy_atoms(residues) -> np.ndarray:
    coords = [
        atom.coord
        for residue in residues
        for atom in residue
        if atom.element != "H"
    ]
    if not coords:
        return np.zeros((0, 3))
    return np.asarray(coords, dtype=float)


def _count_close_atoms(a: np.ndarray, b: np.ndarray, cutoff: float) -> int:
    if len(a) == 0 or len(b) == 0:
        return 0
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return int((distances < cutoff).sum())


def annotate_loops(
    sequence: str, chain_id: str
) -> dict[str, tuple[str, int, int]]:
    """Locate CDR1 and CDR3 from the conserved V-domain anchors.

    CDR3 uses the IMGT junction definition: the stretch bounded by the second
    conserved cysteine and the J-region FGXG motif. CDR1 is placed relative to
    the first conserved cysteine and the framework-2 tryptophan. Both are the
    standard anchors; CDR2 is filled in later from the contact map because
    locating it from sequence alone needs real IMGT numbering.
    """
    loops: dict[str, tuple[str, int, int]] = {}

    junction = [m.start() for m in re.finditer(r"FG.G", sequence)]
    if junction:
        end = junction[0]
        start = sequence.rfind("C", 0, end)
        if start != -1 and 4 <= end - start <= 30:
            loops["CDR3"] = (chain_id, start, end + 1)

    first_cys = sequence.find("C")
    trp = -1
    for match in re.finditer(r"W[YFH][RQKL]Q", sequence):
        if match.start() > first_cys:
            trp = match.start()
            break
    if first_cys != -1 and trp != -1 and trp - first_cys > 8:
        loops["CDR1"] = (chain_id, first_cys + 4, trp - 2)

    return loops


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def _probe_void(
    residue, environment: np.ndarray, partner: np.ndarray
) -> tuple[int, bool]:
    """Measure whether a side chain at this position has room to grow.

    Probe points are placed along the CA->CB direction at 2.0, 3.5 and 5.0 A
    beyond CB. A slot counts as free when no heavy atom sits within 3.0 A of
    it. ``growth_reaches_partner`` records whether that new volume would end up
    against the binding partner rather than pointing into bulk solvent.
    """
    if "CA" not in residue:
        return 0, False
    ca = residue["CA"].coord
    if "CB" in residue:
        direction = residue["CB"].coord - ca
        origin = residue["CB"].coord
    else:
        # Glycine: approximate the missing CB direction from the backbone.
        if "N" not in residue or "C" not in residue:
            return 0, False
        n, c = residue["N"].coord, residue["C"].coord
        bisector = -((n - ca) / np.linalg.norm(n - ca) + (c - ca) / np.linalg.norm(c - ca))
        direction = bisector
        origin = ca + 1.53 * bisector / np.linalg.norm(bisector)
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        return 0, False
    direction = direction / norm

    free_slots = 0
    reaches = False
    for step in (2.0, 3.5, 5.0):
        probe = origin + step * direction
        if len(environment):
            d = np.linalg.norm(environment - probe, axis=1)
            if d.min() < 3.0:
                continue
        free_slots += 1
        if len(partner):
            if np.linalg.norm(partner - probe, axis=1).min() < 6.0:
                reaches = True
    return free_slots, reaches


def load_complex(
    pdb_id: str,
    cache_dir: Path,
    contact_cutoff: float = 5.0,
    chain_override: dict[str, str] | None = None,
) -> Complex:
    """Load a TCR-pMHC complex and measure its interface."""
    path = fetch_pdb(pdb_id, cache_dir)
    structure = PDBParser(QUIET=True).get_structure(pdb_id, str(path))
    model = next(iter(structure))

    sequences, residue_objects = _extract_chain_sequences(model)
    roles, notes = type_chains(sequences, residue_objects, chain_override)

    tcr_ids = [roles[r] for r in ("tcr_alpha", "tcr_beta") if r in roles]
    partner_ids = [roles[r] for r in ("mhc", "b2m", "peptide") if r in roles]
    if not tcr_ids:
        raise ValueError(f"{pdb_id}: no TCR chains identified; pass an explicit override")
    if "peptide" not in roles:
        raise ValueError(f"{pdb_id}: no peptide chain identified; pass an explicit override")

    # Solvent accessibility in the complex, and for the TCR on its own. The
    # difference is the buried surface area, which is a direct measurement.
    sr = ShrakeRupley()
    sr.compute(model, level="R")
    sasa_complex = {
        (chain.id, residue.id[1]): residue.sasa
        for chain in model for residue in chain if hasattr(residue, "sasa")
    }

    tcr_only = _isolate_chains(structure, tcr_ids)
    sr.compute(tcr_only, level="R")
    tcr_only_model = next(iter(tcr_only))
    sasa_free = {
        (chain.id, residue.id[1]): residue.sasa
        for chain in tcr_only_model for residue in chain if hasattr(residue, "sasa")
    }

    peptide_atoms = _heavy_atoms(residue_objects[roles["peptide"]])
    mhc_atoms = _heavy_atoms(
        [r for rid in partner_ids if rid != roles["peptide"] for r in residue_objects[rid]]
    )
    all_partner_atoms = (
        np.vstack([peptide_atoms, mhc_atoms])
        if len(mhc_atoms) and len(peptide_atoms)
        else (peptide_atoms if len(peptide_atoms) else mhc_atoms)
    )

    loops: dict[str, tuple[str, int, int]] = {}
    for role in ("tcr_alpha", "tcr_beta"):
        if role not in roles:
            continue
        chain_id = roles[role]
        prefix = "A" if role == "tcr_alpha" else "B"
        for name, span in annotate_loops(sequences[chain_id], chain_id).items():
            loops[f"{prefix}-{name}"] = span

    residues: list[Residue] = []
    for chain_id, objects in residue_objects.items():
        is_tcr = chain_id in tcr_ids
        self_atoms = _heavy_atoms(objects)
        for index, residue_obj in enumerate(objects):
            aa = _one_letter(residue_obj)
            record = Residue(
                chain=chain_id,
                number=residue_obj.id[1],
                insertion=residue_obj.id[2],
                aa=aa,
                index=index,
                sasa_complex=sasa_complex.get((chain_id, residue_obj.id[1]), 0.0),
                sasa_free=sasa_free.get((chain_id, residue_obj.id[1]), 0.0)
                if is_tcr
                else sasa_complex.get((chain_id, residue_obj.id[1]), 0.0),
            )
            side_chain = np.asarray(
                [a.coord for a in residue_obj if a.get_id() not in BACKBONE_ATOMS and a.element != "H"],
                dtype=float,
            )
            probe_target = side_chain if len(side_chain) else _heavy_atoms([residue_obj])

            if is_tcr and len(probe_target):
                record.contacts_peptide = _count_close_atoms(
                    probe_target, peptide_atoms, contact_cutoff
                )
                record.contacts_mhc = _count_close_atoms(
                    probe_target, mhc_atoms, contact_cutoff
                )
                if len(all_partner_atoms):
                    d = np.linalg.norm(
                        all_partner_atoms[None, :, :] - probe_target[:, None, :], axis=-1
                    )
                    record.nearest_partner_distance = float(d.min())

                neighbours = _neighbour_atoms(residue_objects, chain_id, residue_obj, 12.0)
                record.void_slots, record.growth_reaches_partner = _probe_void(
                    residue_obj, neighbours, all_partner_atoms
                )
                record.partner_charge = _partner_charge_near(
                    probe_target, residue_objects, partner_ids
                )
                record.partner_polar_atoms = _partner_polar_near(
                    probe_target, residue_objects, partner_ids
                )
                record.backbone_nh_bonded = _backbone_nh_bonded(
                    residue_obj, residue_objects
                )
                record.phi = _phi_angle(objects, index)
                record.in_disulfide = aa == "C" and _in_disulfide(
                    residue_obj, residue_objects
                )
            residues.append(record)

    # Fill CDR2/HV4 labels from the contact map: contiguous paratope segments
    # that fall between the anchored CDR1 and CDR3 spans.
    complex_obj = Complex(
        pdb_id=pdb_id.upper(),
        path=path,
        chains=roles,
        sequences=sequences,
        residues=residues,
        loops=loops,
        typing_notes=notes,
    )
    _assign_loop_labels(complex_obj)
    return complex_obj


def _isolate_chains(structure, chain_ids: list[str]):
    import copy

    clone = copy.deepcopy(structure)
    model = next(iter(clone))
    for chain in list(model):
        if chain.id not in chain_ids:
            model.detach_child(chain.id)
    return clone


def _neighbour_atoms(residue_objects, chain_id, residue_obj, radius: float) -> np.ndarray:
    """Heavy atoms near a residue, excluding the residue's own side chain."""
    centre = residue_obj["CA"].coord if "CA" in residue_obj else None
    if centre is None:
        return np.zeros((0, 3))
    coords = []
    for other_chain, objects in residue_objects.items():
        for other in objects:
            if other_chain == chain_id and other.id == residue_obj.id:
                continue
            if "CA" in other and np.linalg.norm(other["CA"].coord - centre) > radius + 8:
                continue
            for atom in other:
                if atom.element == "H":
                    continue
                coords.append(atom.coord)
    if not coords:
        return np.zeros((0, 3))
    arr = np.asarray(coords, dtype=float)
    return arr[np.linalg.norm(arr - centre, axis=1) <= radius]


def _partner_charge_near(probe: np.ndarray, residue_objects, partner_ids) -> float:
    """Net formal charge of partner side chains within 8 A."""
    total = 0.0
    for chain_id in partner_ids:
        for residue in residue_objects[chain_id]:
            atoms = CHARGED_SIDECHAIN_ATOMS.get(residue.get_resname())
            if not atoms:
                continue
            coords = [residue[a].coord for a in atoms if a in residue]
            if not coords:
                continue
            centre = np.mean(coords, axis=0)
            distance = float(np.linalg.norm(probe - centre, axis=1).min())
            if distance < 8.0:
                sign = -1.0 if residue.get_resname() in ("ASP", "GLU") else 1.0
                weight = 1.0 if residue.get_resname() != "HIS" else 0.1
                total += sign * weight * math.exp(-(distance - 3.0) / 3.0)
    return total


def _partner_polar_near(probe: np.ndarray, residue_objects, partner_ids) -> int:
    """Partner nitrogen/oxygen atoms within hydrogen-bonding reach."""
    count = 0
    for chain_id in partner_ids:
        for residue in residue_objects[chain_id]:
            candidates = list(POLAR_SIDECHAIN_ATOMS.get(residue.get_resname(), ()))
            candidates += ["N", "O"]
            for name in candidates:
                if name not in residue:
                    continue
                if float(np.linalg.norm(probe - residue[name].coord, axis=1).min()) < 5.5:
                    count += 1
    return count


def _backbone_nh_bonded(residue_obj, residue_objects) -> bool:
    """Whether the backbone amide is donating a hydrogen bond.

    A position whose backbone NH is engaged cannot tolerate proline, which has
    no amide proton.
    """
    if "N" not in residue_obj:
        return False
    n = residue_obj["N"].coord
    for objects in residue_objects.values():
        for other in objects:
            if other.id == residue_obj.id:
                continue
            for name in ("O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"):
                if name not in other:
                    continue
                if float(np.linalg.norm(other[name].coord - n)) < 3.4:
                    return True
    return False


def _phi_angle(objects, index: int) -> float | None:
    if index == 0:
        return None
    previous, current = objects[index - 1], objects[index]
    try:
        return _dihedral(
            previous["C"].coord, current["N"].coord, current["CA"].coord, current["C"].coord
        )
    except KeyError:
        return None


def _in_disulfide(residue_obj, residue_objects) -> bool:
    if "SG" not in residue_obj:
        return False
    sg = residue_obj["SG"].coord
    for objects in residue_objects.values():
        for other in objects:
            if other.id == residue_obj.id or "SG" not in other:
                continue
            if float(np.linalg.norm(other["SG"].coord - sg)) < 2.5:
                return True
    return False


def _assign_loop_labels(complex_obj: Complex) -> None:
    """Label paratope residues with their loop of origin."""
    spans: dict[str, list[tuple[str, int, int]]] = {}
    for name, (chain_id, start, end) in complex_obj.loops.items():
        spans.setdefault(chain_id, []).append((name, start, end))

    for residue in complex_obj.residues:
        for name, start, end in spans.get(residue.chain, []):
            if start <= residue.index < end:
                residue.loop = name
                break

    # Contiguous paratope stretches with no anchored label sit in CDR2 or HV4.
    for chain_id in complex_obj.tcr_chain_ids():
        prefix = "A" if chain_id == complex_obj.chains.get("tcr_alpha") else "B"
        unlabelled = [
            r for r in complex_obj.residues
            if r.chain == chain_id and r.loop is None and r.contacts_partner > 0
        ]
        for run in _contiguous_runs(unlabelled):
            cdr1 = complex_obj.loops.get(f"{prefix}-CDR1")
            cdr3 = complex_obj.loops.get(f"{prefix}-CDR3")
            midpoint = sum(r.index for r in run) / len(run)
            if cdr1 and cdr3 and cdr1[2] < midpoint < cdr3[1]:
                # Between CDR1 and CDR3: CDR2 comes first, HV4 sits after it.
                name = f"{prefix}-CDR2" if midpoint < (cdr1[2] + cdr3[1]) / 2 else f"{prefix}-HV4"
            else:
                name = f"{prefix}-FR"
            for residue in run:
                residue.loop = name


def _contiguous_runs(residues: list[Residue], max_gap: int = 2) -> list[list[Residue]]:
    runs: list[list[Residue]] = []
    for residue in sorted(residues, key=lambda r: r.index):
        if runs and residue.index - runs[-1][-1].index <= max_gap:
            runs[-1].append(residue)
        else:
            runs.append([residue])
    return runs
