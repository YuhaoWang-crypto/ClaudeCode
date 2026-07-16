"""
systems.py -- concrete, fully public building blocks for the reproduction.

All sequences are REAL and fetched from the RCSB PDB (see provenance below).
The secondary-structure annotations and candidate permutation sites were
computed once with ``structure.py`` (biotite ``annotate_sse``) and are frozen
here so the pipeline runs offline; ``structure.py`` recomputes them on demand
for full reproducibility.

Provenance
----------
Reporter
  TEM-1 beta-lactamase              PDB 3GMW chain A (E.coli class-A beta-lactamase)
      catalytic residues (Ambler)  S70, K73, S130, D131, N132, E166
      author numbering == Ambler   (res_id 28..288, no gaps in this entry)
      paper's insertion loop        Ambler 253  (a surface coil Gly)  -> "BLA-253"
      YES/AND-gate loop             Ambler 196/197

Receptors (ML / de-novo designed minimal ligand-binding domains)
  DIG10.3  digoxigenin binder       PDB 4J9A  (Tinberg et al., Nature 2013)   -> reproduction
  mFAP1    fluorogen (DFHBI) binder PDB 6CZI  (Dou et al., Nature 2018)        -> different-analyte validation

Both receptors are small, designed, and used here as stand-ins for the paper's
17-OHP / cortisol binders (Lee et al., Nat.Commun. 2026), whose exact sequences
are not conveniently machine-fetchable.  The workflow is identical regardless
of which binder is used -- that is precisely the generality the paper claims.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional


# ---------------------------------------------------------------------------
# Reporter enzymes
# ---------------------------------------------------------------------------
@dataclass
class Reporter:
    name: str
    pdb: str
    seq: str
    # 0-indexed catalytic residues (must stay intact for the ON state)
    catalytic: dict            # label -> 0-idx
    # named insertion loops: label -> 0-index cut point (insert BEFORE this idx)
    insertion_sites: dict
    readout: str = "colorimetric"        # colorimetric | electrochemical | luminescent
    mode: str = "insertion"              # insertion | cp_reporter_terminal
    # cp_reporter_terminal only:
    cp_site: int = None                  # 0-idx loop residue kept as the CP point
    terminal_linker: str = "GGSGG"
    reporter_cp_linker: str = "GGGGSGGGGS"


TEM1 = Reporter(
    name="TEM1-BLA",
    pdb="3GMW_A",
    seq=("ETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRR"
         "IHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTR"
         "LDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAG"
         "WFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"),
    # located by unambiguous sequence motifs, cross-checked vs Ambler numbering
    catalytic={"S70": 42, "K73": 45, "S130": 102, "D131": 103, "N132": 104, "E166": 138},
    # cut points = 0-idx of the residue AFTER the named Ambler position
    insertion_sites={"253": 226, "197": 170, "196": 169},
)

# PQQ-glucose dehydrogenase (soluble, A. calcoaceticus), PDB 1CQ1 — ELECTROCHEMICAL
# reporter. Insertion topology at loops 330 (β4/5) and S403-N405 (β5/6).
# Sites + functional residues from Guo et al., JACS 2019 (10.1021/jacs.8b12298);
# integrated from the biosensor-chimera-design skill (verify vs RCSB before wet-lab).
GDH = Reporter(
    name="PQQ-GDH",
    pdb="1CQ1",
    seq=("DVPLTPSQFAKAKSENFDKKVILSNLNKPHALLWGPDNQIWLTERATGKILRVNPESGSVKTVFQVPEIVNDADGQNGLLGFAF"
         "HPDFKNNPYIYISGTFKNPKSTDKELPNQTIIRRYTYNKSTDTLEKPVDLLAGLPSSKDHQSGRLVIGPDQKIYYTIGDQGRNQ"
         "LAYLFLPNQAQHTPTQQELNGKDYHTYMGKVLRLNLDGSIPKDNPSFNGVVSHIYTLGHRNPQGLAFTPNGKLLQSEQGPNSDD"
         "EINLIVKGGNYGWPNVAGYKDDSGYAYANYSAAANKSIKDLAQNGVKVAAGVPVTKESEWTGKNFVPPLKTLYTVQDTYNYNDPT"
         "CGEMTYICWPTVAPSSAYVYKGGKKAITGWENTLLVPSLKRGVIFRIKLDPTYSTTYDDAVPMFKSNNRYRDVIASPDGNVLYVL"
         "TDTAGNVQKDDGSVTNTLENPGSLIKFTYKAK"),
    catalytic={"W346": 345, "T348": 347, "R406": 405, "R408": 407},  # glucose/PQQ, must stay intact
    insertion_sites={"330": 330, "S403-N405": 403},
    readout="electrochemical",
)

# NanoLuc luciferase (O. gracilirostris), PDB 5IBO — LUMINESCENT reporter.
# CP-REPORTER + TERMINAL-FUSION topology: permute NanoLuc at loop 161 and fuse the
# binder at a new terminus. From Guo et al., Nat. Commun. 2022 (10.1038/s41467-022-28425-2).
NANOLUC = Reporter(
    name="NanoLuc",
    pdb="5IBO",
    seq=("SDNMVFTLEDFVGDWRQTAGYNLDQVLEQGGVSSLFQNLGVSVTPIQRIVLSGENGLKIDIHVIIPYEGLSGDQMGQIEKIFKVV"
         "YPVDDHHFKVILHYGTLVIDGVTPNMIDYFGRPYEGIAVFDGKKITVTGTLWNGNKIIDERLINPDGSLLFRVTINGVTGWRLCERILA"),
    catalytic={},                        # whole β-barrel; no single catalytic residue
    insertion_sites={},
    readout="luminescent",
    mode="cp_reporter_terminal",
    cp_site=160,                         # 0-idx == position 161 (paper's CP site)
    terminal_linker="GGSGG",
    reporter_cp_linker="GGGGSGGGGS",
)

REPORTERS = {TEM1.name: TEM1, GDH.name: GDH, NANOLUC.name: NANOLUC}


# ---------------------------------------------------------------------------
# Receptor (ligand-binding) domains
# ---------------------------------------------------------------------------
@dataclass
class Receptor:
    name: str
    pdb: str
    seq: str
    ligand_name: str
    ligand_smiles: str
    ligand_formula: str
    # structure-derived candidate permutation sites (0-idx loop centers),
    # computed by structure.py -> annotate_sse; <10 as in the paper.
    loop_sites: list
    modeled_offset: int = 0    # offset of the modeled structure within seq
    note: str = ""


DIG103 = Receptor(
    name="DIG10.3",
    pdb="4J9A",
    seq=("MNAKEIVVHALRLLENGDARGWSDLFHPEGVLEYPYPPPGYKTRFEGRETIWAHMRLFPEYMTIRFT"
         "DVQFYETADPDLAIGEFHGDGVLTASGGKLAYDYIAVWRTRDGQILLYRLFFNPLRVLEPLGLE"),
    ligand_name="digoxigenin",
    ligand_smiles="CC12CCC(CC1CCC3C2CC(C4(C3(CCC4C5=CC(=O)OC5)O)C)O)O",
    ligand_formula="C23H34O5",
    loop_sites=[36, 58, 74, 92, 110],
    modeled_offset=1,
    note="de-novo digoxigenin binder; His6 tag removed.",
)

MFAP1 = Receptor(
    name="mFAP1",
    pdb="6CZI",
    seq=("MSRAAQLLPGTWQVTMTNEDGQTSQGQWHFQPRSPYTLDIVAQGTISDGRPITGYGKATVKTDDTLH"
         "ANLTYPSLGNIKAQGQITYDSPTQFTWNSTTSDGKKLTGTLQRQE"),
    ligand_name="DFHBI",
    ligand_smiles="CC1=NC(=CC2=CC(=C(C(=C2)F)O)F)C(=O)N1C",
    ligand_formula="C12H10F2N2O2",
    loop_sites=[18, 34, 46, 62, 74, 88, 98],
    modeled_offset=1,
    note="de-novo fluorogen-activating beta-barrel (mini-fluorescence-activating protein).",
)

CDL2_2 = Receptor(
    name="CDL2.2",
    pdb="5IEN",
    # de-novo designed vitamin-D binder (Lee et al. designed-binder family);
    # His6 tag removed. Mined by discover.py from the analyte SMILES.
    seq=("GQSAKEAIEAALADFVKAYNSKDAAGVASKYMDDAAIFPLDMARVDGRQNIQKLWQGLMDMGVSEP"
         "KLTTLDVQESGDFAFESGSISLKGPGKDSKLVDIAGKYVEVWRKGQDGGWKLYRTIANLDPAKLE"),
    ligand_name="vitamin_D3",
    ligand_smiles="CC(CCCC(C)(C)O)C1CCC2C1(CCCC2=CC=C3CC(CCC3=C)O)C",  # VDY (bound vitamin-D analog)
    ligand_formula="C27H44O2",
    loop_sites=[27, 39, 62, 76, 93, 103, 112],
    modeled_offset=1,
    note="de-novo vitamin-D3 binder CDL2.2; discovered via Tier-A PDB mining.",
)

RECEPTORS = {DIG103.name: DIG103, MFAP1.name: MFAP1, CDL2_2.name: CDL2_2}


# ---------------------------------------------------------------------------
# A "system" = one receptor + one reporter + the analyte it detects
# ---------------------------------------------------------------------------
@dataclass
class System:
    key: str
    role: str                  # "reproduction" | "validation"
    receptor: Receptor
    reporter: Reporter
    primary_insertion: str     # key into reporter.insertion_sites
    description: str = ""


SYSTEMS = {
    "dig": System(
        key="dig",
        role="reproduction",
        receptor=DIG103,
        reporter=TEM1,
        primary_insertion="253",
        description="Reproduction: steroid-glycoside sensor (digoxigenin) in the "
                    "paper's TEM-1 / BLA-253 architecture.",
    ),
    "dfhbi": System(
        key="dfhbi",
        role="validation",
        receptor=MFAP1,
        reporter=TEM1,
        primary_insertion="253",
        description="Validation on a DIFFERENT analyte: a small-molecule fluorogen "
                    "(DFHBI) bound by a structurally unrelated beta-barrel receptor, "
                    "same reporter and same design recipe.",
    ),
    "vitd": System(
        key="vitd",
        role="validation",
        receptor=CDL2_2,
        reporter=TEM1,
        primary_insertion="253",
        description="Validation on vitamin D3: a de-novo designed secosteroid binder "
                    "(CDL2.2), mined from the PDB by discover.py, grafted into TEM-1 "
                    "by the same recipe.",
    ),
    # alternative readouts for the same analyte (integrated multi-reporter support)
    "dig-gdh": System(
        key="dig-gdh",
        role="readout-variant",
        receptor=DIG103,
        reporter=GDH,
        primary_insertion="330",
        description="ELECTROCHEMICAL digoxigenin sensor: cpDIG10.3 inserted into "
                    "PQQ-GDH loop 330 (amperometric readout).",
    ),
    "dig-nluc": System(
        key="dig-nluc",
        role="readout-variant",
        receptor=DIG103,
        reporter=NANOLUC,
        primary_insertion="cp161",
        description="LUMINESCENT digoxigenin sensor: NanoLuc circularly permuted at "
                    "loop 161 with cpDIG10.3 fused at a new terminus.",
    ),
}


def get_system(key: str) -> System:
    if key not in SYSTEMS:
        raise KeyError(f"unknown system {key!r}; choose from {list(SYSTEMS)}")
    return SYSTEMS[key]
