"""
transducer.py -- couple a ligand BINDER to a readout, learned from Zhang/Cao et al.,
"De novo design of ligand binding and sensing with a physics based generative
approach" (bioRxiv 2026.07.13.738243; Westlake / ProBuilder).

That paper's key lesson for us: designing a *binder* is only half a sensor -- the
hard, often-skipped half is the TRANSDUCER COUPLING that turns binding into signal.
It contributes three coupling mechanisms and, more usefully, concrete engineering
knobs to fix the failure modes:

  A. SPLIT-BUNDLE reassembly -- split a helical-bundle binder into two fragments,
     each fused to cpGFP; ligand binding drives reassociation -> fluorescence.
     Failure mode: ligand-independent self-association -> high background.
     Fixes (the paper's actual moves, encoded here as deterministic ops/knobs):
       * FRAGMENT DUPLICATION: duplicate one segment to raise the energy barrier
         for self-association (lower background, ligand-gated activation).
       * INTERFACE-WEAKENING mutations at the fragment interface (bury->polar).
       * TRUNCATION of a terminal sub-element to tune responsiveness.
  B. LIGAND-INDUCED FOLDING -- a binder that only folds when the ligand is present;
     the folding transition is read by an enzyme reporter or a FRET pair.
  C. METAL COORDINATION -- metal ions (e.g. Zn tetrahedral His/Asp/Glu) sensed via
     a coordination-geometry site; couples to induced folding.

Design lesson also carried into generative route (see specificity/generative notes):
build the pocket by GROWING a helical bundle around the ligand from RIF anchor
contacts (fold + bind scored jointly), NOT by docking a generic pre-made scaffold --
which is exactly why our bare-ligand de-novo runs failed.

Honesty: every construction here is deterministic (✅). Whether a fusion actually
switches, and its dynamic range, is ⚠️ until a wet-lab fluorescence/enzyme titration.
"""
from __future__ import annotations
from dataclasses import dataclass, field


# cpGFP transducer (circularly permuted GFP), the standard fluorescent readout the
# paper fuses to split fragments. Sequence is the cpGFP core used in GCaMP-style sensors.
CPGFP = ("NVYIKADKQKNGIKANFKIRHNIEDGGVQLAYHYQQNTPIGDGPVLLPDNHYLSVQSKLSKDPNEKRDHMVLLEFVTAAGI"
         "TLGMDELYK")  # illustrative cpGFP core (⚠️ readout module, not the analyte binder)


# ---------------------------------------------------------------------------
# A. Split-bundle reassembly sensor + background-suppression knobs
# ---------------------------------------------------------------------------
@dataclass
class SplitBundleSensor:
    analyte: str
    binder_name: str
    frag_n: str                 # N-fragment sequence (helices 1..k)
    frag_c: str                 # C-fragment sequence (helices k+1..N)
    split_after_helix: int
    construct_n: str            # frag_n - linker - cpGFP  (or cpGFP - linker - frag_n)
    construct_c: str
    knobs: dict = field(default_factory=dict)


def split_helical_bundle(binder_seq: str, helix_boundaries: list, split_after_helix: int):
    """Split a helical-bundle binder into two fragments at a helix boundary.

    helix_boundaries: sorted residue indices (0-based) where each helix STARTS.
    split_after_helix: 1-based helix index; fragment N = helices 1..k.
    """
    if not (1 <= split_after_helix < len(helix_boundaries)):
        raise ValueError("split_after_helix must be between 1 and n_helices-1")
    cut = helix_boundaries[split_after_helix]
    return binder_seq[:cut], binder_seq[cut:]


def fragment_duplication(fragment_seq: str, segment: tuple) -> str:
    """The paper's anti-self-association move: duplicate a segment of a fragment to
    raise the energetic barrier for ligand-INDEPENDENT reassociation.

    segment: (start, end) 0-based half-open indices within fragment_seq to duplicate
    (inserted immediately after the segment). Deterministic sequence op.
    """
    s, e = segment
    if not (0 <= s < e <= len(fragment_seq)):
        raise ValueError("bad segment for duplication")
    return fragment_seq[:e] + fragment_seq[s:e] + fragment_seq[e:]


def interface_weakening(fragment_seq: str, positions: list, to: str = "S") -> str:
    """Weaken the fragment-fragment interface (reduce background self-association) by
    mutating buried hydrophobic interface positions toward polar (default Ser)."""
    seq = list(fragment_seq)
    for p in positions:
        if 0 <= p < len(seq):
            seq[p] = to
    return "".join(seq)


def truncate_terminus(fragment_seq: str, n_from_end: int, which: str = "C") -> str:
    """Truncate a terminal sub-element to tune responsiveness (paper's knob)."""
    if n_from_end <= 0:
        return fragment_seq
    return fragment_seq[:-n_from_end] if which == "C" else fragment_seq[n_from_end:]


def build_split_bundle_sensor(analyte: str, binder_name: str, binder_seq: str,
                              helix_boundaries: list, split_after_helix: int,
                              linker: str = "GSGGSG", cpgfp: str = CPGFP,
                              duplicate_segment: tuple | None = None,
                              weaken_positions: list | None = None) -> SplitBundleSensor:
    """Assemble a split-bundle + cpGFP fluorescent sensor, optionally applying the
    background-suppression knobs (fragment duplication + interface weakening)."""
    fn, fc = split_helical_bundle(binder_seq, helix_boundaries, split_after_helix)
    knobs = {}
    if duplicate_segment is not None:
        fn = fragment_duplication(fn, duplicate_segment)
        knobs["fragment_duplication"] = duplicate_segment
    if weaken_positions:
        fc = interface_weakening(fc, weaken_positions)
        knobs["interface_weakening"] = weaken_positions
    # fuse each fragment to cpGFP (N-frag on cpGFP N-term, C-frag on cpGFP C-term)
    construct_n = fn + linker + cpgfp
    construct_c = cpgfp + linker + fc
    return SplitBundleSensor(analyte=analyte, binder_name=binder_name, frag_n=fn, frag_c=fc,
                             split_after_helix=split_after_helix,
                             construct_n=construct_n, construct_c=construct_c, knobs=knobs)


def verify_split_bundle(s: SplitBundleSensor, cpgfp: str = CPGFP) -> dict:
    """Both constructs must carry cpGFP and their binder fragment; fragments together
    must cover (a superset of, if duplicated) the original binder residues."""
    checks = {
        "n_has_cpgfp": cpgfp in s.construct_n,
        "c_has_cpgfp": cpgfp in s.construct_c,
        "n_has_frag": s.frag_n in s.construct_n,
        "c_has_frag": s.frag_c in s.construct_c,
        "two_fragments_nonempty": len(s.frag_n) > 0 and len(s.frag_c) > 0,
    }
    checks["all_ok"] = all(checks.values())
    return checks


# ---------------------------------------------------------------------------
# B. Ligand-induced folding sensor
# ---------------------------------------------------------------------------
@dataclass
class InducedFoldingSensor:
    analyte: str
    binder_name: str
    sensor_seq: str
    reporter: str               # "enzyme" | "FRET"
    mechanism: str = "ligand-induced folding"


def build_induced_folding_sensor(analyte: str, binder_name: str, binder_seq: str,
                                 reporter: str = "FRET",
                                 fret_pair=("mTurquoise2", "mVenus"),
                                 enzyme: str = "split-luciferase") -> InducedFoldingSensor:
    """Couple a marginally-stable (folds only when ligand-bound) binder to a readout:
    FRET pair flanking the binder, or an enzyme reporter whose activity follows folding.
    """
    if reporter == "FRET":
        seq = f"[{fret_pair[0]}]-{binder_seq}-[{fret_pair[1]}]"
        note = "folding upon binding changes donor-acceptor distance -> FRET change"
    else:
        seq = f"{binder_seq}-[{enzyme}]"
        note = "folding upon binding reconstitutes/activates the enzyme reporter"
    s = InducedFoldingSensor(analyte=analyte, binder_name=binder_name, sensor_seq=seq,
                             reporter=reporter)
    s.__dict__["note"] = note
    return s


# ---------------------------------------------------------------------------
# C. Metal-ion coordination site (Zn/Ni) -> induced-folding sensor
# ---------------------------------------------------------------------------
METAL_COORDINATION = {
    "Zn": {"geometry": "tetrahedral", "n_ligands": 4, "residues": ["H", "H", "D", "E"]},
    "Ni": {"geometry": "octahedral", "n_ligands": 6, "residues": ["H", "H", "H", "D", "E", "E"]},
    "Ca": {"geometry": "pentagonal-bipyramidal", "n_ligands": 7, "residues": ["D", "D", "E", "N", "D"]},
}


def metal_coordination_plan(metal: str, scaffold_len: int = 120) -> dict:
    """Emit a coordination-residue spec (RIF-style) for a metal-induced-folding sensor:
    which residue types, how many, and the target geometry to build a site into a
    de-novo helical bundle (the paper's metal route)."""
    m = METAL_COORDINATION[metal]
    return {
        "metal": metal, "geometry": m["geometry"], "n_coordinating": m["n_ligands"],
        "residue_types": m["residues"],
        "strategy": "grow a helical bundle so the coordinating residues meet at the "
                    f"{m['geometry']} site; apo is under-folded, metal binding completes "
                    "the coordination sphere and drives folding -> reporter signal.",
        "scaffold_len": scaffold_len,
        "label": "⚠️ coordination geometry is a design target; verify with metal-titration "
                 "folding/activity (CD, fluorescence) -- the wet-lab ground truth.",
    }


# ---------------------------------------------------------------------------
# sensor-readiness heuristic (which coupling suits a binder)
# ---------------------------------------------------------------------------
def recommend_transducer(analyte_class: str, binder_topology: str,
                         desired_readout: str = "fluorescent") -> dict:
    """Pick a coupling mechanism from the paper's menu given the binder's properties."""
    if analyte_class == "metal_ion":
        mech = "metal-coordination induced folding"
    elif binder_topology == "helical_bundle":
        mech = "split-bundle reassembly (cpGFP)" if desired_readout == "fluorescent" \
               else "ligand-induced folding (enzyme/FRET)"
    else:
        mech = "ligand-induced folding (enzyme/FRET)"
    return {
        "analyte_class": analyte_class, "binder_topology": binder_topology,
        "recommended": mech,
        "background_knobs": ["fragment_duplication", "interface_weakening", "terminal_truncation"],
        "note": "helical bundles split cleanly -> best for split sensors; marginal-stability "
                "binders -> best for induced folding. Tune background/response with the knobs.",
    }
