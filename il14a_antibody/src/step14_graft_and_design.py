"""Step 14 - Graft TXLNA onto a real bound-peptide backbone, then redesign the
paratope with ProteinMPNN.

WHY A GRAFT AND NOT DE NOVO
---------------------------
ProteinMPNN is an inverse-folding model: it reads a BACKBONE and writes a
sequence. It cannot invent the backbone, and on this target nothing else in
reach can either - the epitope is disordered, so there is no complex structure
to start from and no hydrophobic core for a diffusion model to nucleate on.

Step 13 supplies the missing backbone empirically. Of 1759 Ig-annotated PDB
entries, 96 are antibody-peptide complexes with a measured paratope contact,
and the closest analogues of ERRPEGPGAQAPSSPRVTEAPC are the anti-MUC1
antibody SM3 bound to the tandem repeat - aromatic-free, one third Pro+Gly,
disordered in isolation, and gripped anyway. SM3's bound peptide backbone is
an experimental answer to "what conformation does this kind of chain adopt in
a paratope", and it is exactly the input ProteinMPNN needs.

WHAT IS AND IS NOT ASSUMED
--------------------------
ProteinMPNN reads N, CA, C and O only. So swapping the template peptide's
sequence for a TXLNA window while keeping its backbone is not an approximation
inside the model - the model sees precisely what we intend it to see. The
assumption sits one level up, in whether TXLNA CAN adopt that backbone, and
that is testable rather than hopeful: proline fixes phi near -63 and only
glycine tolerates positive phi, so a threading is admissible only if every
proline in the TXLNA window lands on a template phi that accepts one and no
non-glycine lands on a positive-phi position. Windows that fail are dropped
before any design is run, which is the step the cheap version of this
pipeline skips.
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS
from step13_template_mining import read_atoms, SOLVENT, get_structure

TARGET = "ERRPEGPGAQAPSSPRVTEAPC"        # TXLNA 502-523 (1C6 immunogen)
TARGET_START = 502

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "CIR": "R", "SEP": "S", "TPO": "T", "PTR": "Y",
}
ONE2THREE = {v: k for k, v in list(THREE2ONE.items())[:20]}

# Templates carried forward from step 13. SM3 appears as both a Fab (1SM3,
# VH+VL, 13-residue epitope) and an scFv at higher resolution (5A2J). The Fab
# is the design scaffold because a VH/VL pair is what a therapeutic is built
# from; 2W65 is kept because it is the one verified PTM-specific paratope on a
# Pro/Gly-rich epitope, which is the closest thing to a template for the
# phospho arm of the question.
TEMPLATES = [
    {"pdb": "1SM3", "ab": ["H", "L"], "pep": ["P"], "note": "SM3 Fab / MUC1"},
    {"pdb": "5A2J", "ab": ["H"], "pep": ["P"], "note": "scFv-SM3 / APDTRP"},
    {"pdb": "2W65", "ab": ["A", "B"], "pep": ["E"],
     "note": "ACC4 Fab / citrullinated collagen II"},
]


def residues_of(atoms, chains):
    """Ordered backbone residues: {(chain, resseq): {atom: xyz}, comp}."""
    out = {}
    for ch, rid, comp, xyz in atoms:
        if ch not in chains or comp in SOLVENT:
            continue
        out.setdefault((ch, rid), {"comp": comp, "at": {}})
    return out


def backbone(cif_atoms, chains, want=("N", "CA", "C", "O")):
    """Return ordered [(chain, resseq, comp, {atom:xyz})] for the chains."""
    # read_atoms drops the atom name, so re-read with names here.
    raise NotImplementedError


def read_atoms_named(cif):
    atoms = []
    with open(cif) as fh:
        cols, active = {}, False
        for line in fh:
            if line.startswith("loop_"):
                cols, active = {}, False
                continue
            if line.startswith("_atom_site."):
                cols[line.strip().split(".", 1)[1]] = len(cols)
                active = True
                continue
            if not active or not cols:
                continue
            p = line.split()
            if not p or p[0] not in ("ATOM", "HETATM"):
                if atoms and line.startswith("#"):
                    break
                continue
            if len(p) <= max(cols.values()):
                continue
            alt = p[cols["label_alt_id"]]
            if alt not in (".", "A"):
                continue
            if p[cols.get("pdbx_PDB_model_num", 0)] not in ("1",) and \
               "pdbx_PDB_model_num" in cols:
                continue
            try:
                xyz = (float(p[cols["Cartn_x"]]), float(p[cols["Cartn_y"]]),
                       float(p[cols["Cartn_z"]]))
            except ValueError:
                continue
            # Antibodies are deposited in Kabat/Chothia numbering, which
            # uses insertion codes: 1SM3 carries H52A/B/C in CDR-H2, L27A/B/C
            # in CDR-L1 and H82A/B/C. Keying a residue on auth_seq_id alone
            # merges each insertion into its parent - silently, and precisely
            # in the loops this whole exercise is about.
            icode = p[cols["pdbx_PDB_ins_code"]] if "pdbx_PDB_ins_code" in cols else "."
            atoms.append({
                "chain": p[cols["auth_asym_id"]],
                "resseq": p[cols["auth_seq_id"]] + ("" if icode in (".", "?")
                                                    else icode),
                "comp": p[cols["label_comp_id"]],
                "name": p[cols["auth_atom_id"]] if "auth_atom_id" in cols
                        else p[cols["label_atom_id"]],
                "elem": p[cols["type_symbol"]],
                "xyz": xyz,
            })
    return atoms


def chain_residues(atoms, chain):
    order, res = [], {}
    for a in atoms:
        if a["chain"] != chain or a["comp"] in SOLVENT:
            continue
        key = a["resseq"]
        if key not in res:
            res[key] = {"comp": a["comp"], "at": {}, "resseq": key}
            order.append(key)
        res[key]["at"][a["name"].strip('"')] = a["xyz"]
    return [res[k] for k in order]


def dihedral(p0, p1, p2, p3):
    def sub(a, b):
        return (a[0] - b[0], a[1] - b[1], a[2] - b[2])

    def cross(a, b):
        return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0])

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    b0, b1, b2 = sub(p0, p1), sub(p2, p1), sub(p3, p2)
    n = math.sqrt(dot(b1, b1))
    b1n = tuple(x / n for x in b1)
    v = tuple(b0[i] - dot(b0, b1n) * b1n[i] for i in range(3))
    w = tuple(b2[i] - dot(b2, b1n) * b1n[i] for i in range(3))
    return math.degrees(math.atan2(dot(cross(b1n, v), w), dot(v, w)))


def phi_psi(res):
    out = []
    for i, r in enumerate(res):
        phi = psi = None
        try:
            if i > 0 and "C" in res[i - 1]["at"]:
                phi = dihedral(res[i - 1]["at"]["C"], r["at"]["N"],
                               r["at"]["CA"], r["at"]["C"])
            if i + 1 < len(res) and "N" in res[i + 1]["at"]:
                psi = dihedral(r["at"]["N"], r["at"]["CA"], r["at"]["C"],
                               res[i + 1]["at"]["N"])
        except KeyError:
            pass
        out.append((phi, psi))
    return out


# Backbone admissibility. Proline's ring fixes phi; glycine is the only
# residue that lives comfortably at positive phi. Everything else is a soft
# preference and is deliberately NOT scored here, because a soft preference
# summed over 13 positions will outvote a hard clash.
PRO_PHI = (-105.0, -35.0)


def threading_penalty(window, phipsi):
    """Hard, interpretable reasons a TXLNA window cannot take this backbone."""
    bad = []
    for i, aa in enumerate(window):
        phi, psi = phipsi[i]
        if phi is None:
            continue
        if aa == "P" and not (PRO_PHI[0] <= phi <= PRO_PHI[1]):
            bad.append(f"P@{i+1} needs phi in {PRO_PHI}, template has {phi:.0f}")
        if aa != "G" and phi > 0:
            bad.append(f"{aa}@{i+1} at positive phi {phi:.0f} (Gly only)")
    return bad


def burial(res_pep, ab_atoms, cutoff=4.5):
    """Per-epitope-residue heavy-atom contact count with the antibody."""
    out = []
    c2 = cutoff ** 2
    for r in res_pep:
        n = 0
        for nm, (x, y, z) in r["at"].items():
            for (ax, ay, az) in ab_atoms:
                if (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2 <= c2:
                    n += 1
                    break
        out.append(n)
    return out


def main():
    print("=" * 100)
    print("TEMPLATE BACKBONES AND WHICH TXLNA WINDOWS CAN TAKE THEM")
    print("=" * 100)
    report = []
    for t in TEMPLATES:
        cif = get_structure(t["pdb"])
        atoms = read_atoms_named(cif)
        pep = [r for c in t["pep"] for r in chain_residues(atoms, c)]
        pep = [r for r in pep if {"N", "CA", "C"} <= set(r["at"])]
        ab = [a["xyz"] for a in atoms
              if a["chain"] in t["ab"] and a["comp"] not in SOLVENT
              and a["elem"] != "H"]
        pp = phi_psi(pep)
        seq = "".join(THREE2ONE.get(r["comp"], "X") for r in pep)
        bur = burial(pep, ab)
        L = len(pep)
        print(f"\n{t['pdb']}  {t['note']}   epitope {L} aa  '{seq}'")
        print(f"  {'i':3} {'res':5} {'phi':>8} {'psi':>8}  contacts")
        for i, r in enumerate(pep):
            phi, psi = pp[i]
            f = f"{phi:8.1f}" if phi is not None else "      --"
            s = f"{psi:8.1f}" if psi is not None else "      --"
            print(f"  {i+1:<3} {r['comp']:5} {f} {s}  {'*' * min(bur[i], 20)}")

        ok = []
        for s in range(0, len(TARGET) - L + 1):
            w = TARGET[s:s + L]
            bad = threading_penalty(w, pp)
            # score admissible windows by how well buried positions match the
            # template's own residue class - a graft that puts a charge where
            # the template buries a small neutral is admissible but poor.
            score = sum(b for i, b in enumerate(bur)
                        if w[i] == seq[i]) if len(seq) == L else 0
            ok.append({"offset_in_peptide": s,
                       "txlna_range": [TARGET_START + s, TARGET_START + s + L - 1],
                       "window": w, "admissible": not bad,
                       "violations": bad, "identity_to_template_seq": round(
                           sum(a == b for a, b in zip(w, seq)) / L, 3),
                       "buried_match_score": score})
        n_ok = sum(1 for o in ok if o["admissible"])
        print(f"  -> {n_ok}/{len(ok)} TXLNA {L}-mer windows are backbone-admissible")
        for o in ok:
            flag = "OK " if o["admissible"] else "no "
            print(f"     {flag} {o['txlna_range']} {o['window']}"
                  + ("" if o["admissible"] else f"   {o['violations'][0]}"))
        report.append({"template": t, "epitope_seq": seq, "epitope_len": L,
                       "phi_psi": [[None if x is None else round(x, 1)
                                    for x in p] for p in pp],
                       "contacts_per_residue": bur, "windows": ok})

    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "step14_threading.json").write_text(json.dumps(report, indent=2))
    print("\nwrote results/step14_threading.json")


if __name__ == "__main__":
    main()
