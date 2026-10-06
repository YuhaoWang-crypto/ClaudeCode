"""Audit an AlphaProtein Novo run package before it reaches a GPU.

Every finding below is re-derived when this runs; nothing is quoted. The point
is that an AP Novo manifest can pass its own schema validation and still be
unrunnable, because the schema checks types and the *grammar* of `motif_str` is
checked much later -- on a GPU node, after the weights have loaded.

## The grammar rule that decides everything

`alphaprotein_novo/data/residue_mapping.py`:

    def _is_segment_designable(segment_str: str) -> bool:
      return segment_str[0].isnumeric()

A comma-separated segment starting with a **letter** is a fixed motif residue;
one starting with a **digit** is a designable length. So `A300` fixes chain A
residue 300, while a bare `302` asks for a designable segment of exactly 302
residues. Writing a residue list as `A300,302,305,...` -- the natural way to
abbreviate -- therefore fixes **one** residue and silently requests thousands of
designed ones.

The official grammar interleaves them, each fixed residue carrying its chain:

    "20-100,A56,20-70,A100,20-70,A125,20-100/A2001/A3001"

or lifts the motif out with `|` and `{}` placeholders, or -- the mode meant for
exactly this problem -- puts the residues in `unindexed_motif_residues` and
leaves `motif_str` as a plain length range.

Run:  python3 -m trnaplat.apnovo.audit_apnovo_package <package_dir>
      python3 -m trnaplat.apnovo.audit_apnovo_package <package_dir> --repo <checkout>

Without `--repo` the checks use the documented grammar rule, reimplemented here
in three lines. With `--repo` they additionally call AP Novo's own
`Manifest.from_file` and `sample_designable_lengths`, so the verdict is the
upstream code's, not this module's.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

#: Wild-type *Methanosarcina mazei* PylRS, UniProt Q8PWY1, 454 aa. The same
#: constant `trnaplat.demo1_multisite` uses, so a motif checked here and a
#: transplant checked there are in one numbering.
MM_PYLRS_WT = (
    "MDKKPLNTLISATGLWMSRTGTIHKIKHHEVSRSKIYIEMACGDHLVVNNSRSSRTARAL"
    "RHHKYRKTCKRCRVSDEDLNKFLTKANEDQTSVKVKVVSAPTRTKKAMPKSVARAPKPLE"
    "NTEAAQAQPSGSKFSPAIPVSTQESVSVPASVSTSISSISTGATASALVKGNTNPITSMS"
    "APVQASAPALTKSQTDRLEVLLNPKDEISLNSGKPFRELESELLSRRKKDLQQIYAEERE"
    "NYLGKLEREITRFFVDRGFLEIKSPILIPLEYIERMGIDNDTELSKQIFRVDKNFCLRPM"
    "LAPNLYNYLRKLDRALPDPIKIFEIGPCYRKESDGKEHLEEFTMLNFCQMGSGCTRENLE"
    "SIITDFLNHLGIDFKIVGDSCMVYGDTLDVMHGDLELSSAVVGPIPLDREWGIDKPWIGA"
    "GFGLERLLKVKHDFKNIKRAARSESYYNGISTNL"
)

#: Atom names of the AMP half of CCD `YLY`, pyrrolysyl-adenylate.
#: ✅ CCD YLY is "C22 H35 N8 O9 P" (RCSB chemcomp API) and PDB 2Q7H is titled
#: "Pyrrolysyl-tRNA synthetase bound to adenylated pyrrolysine and
#: pyrophosphate" -- so the ligand in a 2Q7H-derived motif is the **adenylate
#: intermediate**, not the free ncAA. `OAY` is the ester oxygen bridging the
#: pyrrolysyl carboxyl to the phosphate and is counted with the AMP half.
AMP_MOIETY_ATOMS = {
    "OAY", "PBN", "OAI", "OAF", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'",
    "C2'", "O2'", "C1'", "N9", "N7", "C2", "N6",
}


# ------------------------------------------------------------------ grammar

def is_designable(segment: str) -> bool:
    """AP Novo's own rule, from residue_mapping._is_segment_designable."""
    return bool(segment) and segment[0].isnumeric()


#: Placeholder in the `motif|...{}...` form. The motifs left of `|` are
#: substituted into these, so a `{}` is neither a fixed residue of its own nor a
#: designable length -- counting it as either double-counts the motif.
PLACEHOLDER = "{}"


def split_designable_chain(chain: str) -> tuple[list[str], list[tuple[int, int]]]:
    """Return (fixed motif segments, (min, max) of each designable segment)."""
    fixed, ranges = [], []
    for segment in chain.split(","):
        if segment == PLACEHOLDER:
            continue
        if not is_designable(segment):
            fixed.append(segment)
            continue
        if "-" in segment:
            lo, hi = segment.split("-", 1)
            ranges.append((int(lo), int(hi)))
        else:
            ranges.append((int(segment), int(segment)))   # exact length
    return fixed, ranges


def check_motif_str(motif_str: str, seq_length: str | None,
                    declared_residues: int | None = None) -> dict:
    """Does this motif_str mean what its author meant, and can it be sampled?"""
    chains = motif_str.split("|")[-1].split("/")
    designable = [c for c in chains if any(is_designable(s) for s in c.split(","))]
    if len(designable) != 1:
        return {"ok": False,
                "error": f"AP Novo requires exactly 1 designable chain, found "
                         f"{len(designable)}"}

    fixed, ranges = split_designable_chain(designable[0])
    prefix = motif_str.split("|")[0] if "|" in motif_str else ""
    prefix_motifs = [s for s in prefix.split(",") if s] if prefix else []
    min_sum = sum(lo for lo, _ in ranges)
    max_sum = sum(hi for _, hi in ranges)

    result = {
        "ok": True, "error": None,
        "fixed_residues": fixed + prefix_motifs,
        "n_fixed": len(fixed) + len(prefix_motifs),
        "designable_min": min_sum, "designable_max": max_sum,
        "exact_length_segments": [lo for lo, hi in ranges if lo == hi],
    }
    if seq_length:
        lo, hi = ((int(x) for x in seq_length.split("-")) if "-" in seq_length
                  else (int(seq_length), int(seq_length)))
        result["seq_length"] = (lo, hi)
        if min_sum > hi or max_sum < lo:
            result["ok"] = False
            result["error"] = (
                f"Cannot sample values in range [{lo}, {hi}]: the range of "
                f"possible sample values is [{min_sum}, {max_sum}]"
            )
    if declared_residues is not None and result["n_fixed"] != declared_residues:
        result["silently_dropped"] = declared_residues - result["n_fixed"]
    return result


# --------------------------------------------------------------------- CIF

class CifHeaderError(ValueError):
    """Raised when the mmCIF data-block header is malformed."""


def check_cif_header(lines: list[str]) -> str | None:
    """The mmCIF data block must open with `data_<name>`, no space.

    ✅ Found the hard way: the audited package opens with `data pylrs_pyl_motif`,
    which every permissive reader here accepted and AlphaFold 3's parser refused
    with `INVALID_ARGUMENT: The CIF file does not start with the data_ field.`
    That error only surfaces once the real pipeline loads the file, so it is
    checked up front instead.
    """
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_") and len(stripped) > len("data_"):
            return None
        return (
            f"first content line is {stripped[:40]!r}; mmCIF requires a data "
            "block header of the form `data_<name>` with no space. AlphaFold 3 "
            "rejects this with 'The CIF file does not start with the data_ "
            "field.'"
        )
    return "file has no content lines"


#: mmCIF categories the repo's own working example carries
#: (`examples/kemp_eliminase/*.cif`). ✅ A motif file is not an atom list: AF3
#: builds `author_naming_scheme` from the entity and `pdbx_*_scheme` tables, and
#: without them a load fails with `KeyError: ('B', 1)` while renumbering the
#: ligand reference in `motif_str`.
REQUIRED_CIF_CATEGORIES = [
    "_entry", "_chem_comp", "_entity", "_entity_poly", "_entity_poly_seq",
    "_pdbx_poly_seq_scheme", "_struct_asym", "_atom_site",
]


def cif_categories(lines: list[str]) -> list[str]:
    seen = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("_") and "." in stripped:
            category = stripped.split(".", 1)[0]
            if category not in seen:
                seen.append(category)
    return seen


def check_cif_shape(lines: list[str]) -> list[str]:
    """Categories a loadable motif needs and this file lacks."""
    present = set(cif_categories(lines))
    return [c for c in REQUIRED_CIF_CATEGORIES if c not in present]


def read_motif_cif(path: pathlib.Path, strict: bool = True) -> dict:
    """Minimal mmCIF atom_site reader -- enough to audit a motif file."""
    lines = path.read_text().splitlines()
    header_problem = check_cif_header(lines)
    if header_problem and strict:
        raise CifHeaderError(f"{path.name}: {header_problem}")
    header = [l.strip().split(".", 1)[1] for l in lines
              if l.strip().startswith("_atom_site.")]
    if not header:
        raise ValueError(f"{path} has no _atom_site loop")
    col = {name: i for i, name in enumerate(header)}
    rows = [l.split() for l in lines if l.startswith(("ATOM", "HETATM"))]

    def unquote(value: str | None) -> str | None:
        """Strip mmCIF quoting.

        ⚠️ Load-bearing: an atom name containing a prime, like O5', is written
        quoted ("O5'") because the bare token would be ambiguous. Without this,
        9 of the adenylate's 17 nucleotide atoms fail to match
        `AMP_MOIETY_ATOMS` and the moiety split silently shifts by three
        residues -- which it did, until the two CIFs were diffed atom by atom.
        """
        if value is None or len(value) < 2:
            return value
        if value[0] == value[-1] and value[0] in "\"'":
            return value[1:-1]
        return value

    def get(row, *names):
        for name in names:
            if name in col:
                return unquote(row[col[name]])
        return None

    residues, ligand_atoms = {}, []
    for row in rows:
        chain = get(row, "auth_asym_id", "label_asym_id")
        res_id = int(get(row, "auth_seq_id", "label_seq_id"))
        comp = get(row, "label_comp_id", "auth_comp_id")
        if comp in THREE_TO_ONE:
            residues[(chain, res_id)] = comp
        else:
            ligand_atoms.append({
                "chain": chain, "comp": comp,
                "atom": get(row, "label_atom_id"),
                "element": get(row, "type_symbol"),
                "xyz": tuple(float(get(row, f"Cartn_{a}")) for a in "xyz"),
            })
    protein_xyz = {}
    for row in rows:
        comp = get(row, "label_comp_id", "auth_comp_id")
        if comp not in THREE_TO_ONE:
            continue
        key = (get(row, "auth_asym_id", "label_asym_id"),
               int(get(row, "auth_seq_id", "label_seq_id")))
        protein_xyz.setdefault(key, []).append(
            tuple(float(get(row, f"Cartn_{a}")) for a in "xyz"))
    return {"residues": residues, "ligand_atoms": ligand_atoms,
            "protein_xyz": protein_xyz, "n_atoms": len(rows)}


def check_numbering(residues: dict, reference: str = MM_PYLRS_WT) -> list[dict]:
    """Does every motif residue carry the reference sequence's own residue?"""
    out = []
    for (chain, res_id), comp in sorted(residues.items(), key=lambda kv: kv[0][1]):
        one = THREE_TO_ONE.get(comp, "?")
        ref = reference[res_id - 1] if res_id <= len(reference) else "?"
        out.append({"chain": chain, "res_id": res_id, "comp": comp,
                    "one_letter": one, "reference": ref, "match": one == ref})
    return out


#: Below this difference the two subsite distances are a tie, and the residue
#: contacts both halves of the adenylate. ⚠️ Chosen, not derived: it is roughly
#: the coordinate error of a 2.1 A structure, and several motif residues sit
#: inside it, so a clean two-way split of this pocket does not exist.
MOIETY_TIE_ANGSTROM = 0.5


def partition_by_moiety(info: dict, tie: float = MOIETY_TIE_ANGSTROM) -> list[dict]:
    """Split motif residues by which half of the adenylate they contact.

    The question: a motif described as "the ncAA binding pocket", extracted from
    a structure whose ligand is the **adenylate**, is partly an ATP-site motif.
    This measures how much -- and refuses to call the residues that touch both.
    """
    import math

    ncaa = [a for a in info["ligand_atoms"] if a["atom"] not in AMP_MOIETY_ATOMS]
    amp = [a for a in info["ligand_atoms"] if a["atom"] in AMP_MOIETY_ATOMS]
    if not ncaa or not amp:
        return []

    def nearest(coords, atoms):
        return min(
            math.dist(c, a["xyz"]) for c in coords for a in atoms
        )

    rows = []
    for (chain, res_id), coords in sorted(info["protein_xyz"].items(),
                                          key=lambda kv: kv[0][1]):
        d_ncaa, d_amp = nearest(coords, ncaa), nearest(coords, amp)
        if abs(d_ncaa - d_amp) < tie:
            side = "both"
        else:
            side = "ncAA" if d_ncaa < d_amp else "AMP"
        rows.append({
            "residue": f"{info['residues'][(chain, res_id)]}{res_id}",
            "res_id": res_id,
            "d_ncaa": round(d_ncaa, 2), "d_amp": round(d_amp, 2),
            "margin": round(abs(d_ncaa - d_amp), 2),
            "side": side,
        })
    return rows


# ------------------------------------------------------- upstream cross-check

def upstream_check(repo: pathlib.Path, manifest: pathlib.Path) -> dict:
    """Call AP Novo's own validator and length sampler, if a checkout is given.

    Imported lazily and guarded: the repo pulls in jax and alphafold3 at import
    time, so this is best-effort. When it cannot run, the documented-grammar
    checks above still stand on their own.
    """
    out = {"available": False}
    sys.path.insert(0, str(repo / "src"))
    try:
        from alphaprotein_novo.data import design_manifest  # noqa: PLC0415
    except Exception as err:                                 # noqa: BLE001
        out["schema_error"] = f"{type(err).__name__}: {err}"
    else:
        out["available"] = True
        try:
            design_manifest.Manifest.from_file(manifest)
            out["schema"] = "PASSED"
        except Exception as err:                             # noqa: BLE001
            out["schema"] = f"FAILED: {type(err).__name__}: {err}"
    try:
        import numpy as np                                   # noqa: PLC0415

        from alphaprotein_novo.data import residue_mapping   # noqa: PLC0415
    except Exception as err:                                 # noqa: BLE001
        out["sampler_error"] = f"{type(err).__name__}: {err}"
        return out
    spec = json.loads(manifest.read_text())
    defaults = spec.get("defaults", {})
    out["sampler"] = {}
    for job in spec.get("designs", []):
        seq_length = job.get("seq_length", defaults.get("seq_length"))
        try:
            residue_mapping.sample_designable_lengths(
                job["motif_str"], np.random.default_rng(0), seq_length)
            out["sampler"][job["name"]] = "PASSED"
        except Exception as err:                             # noqa: BLE001
            out["sampler"][job["name"]] = f"{type(err).__name__}: {err}"
    return out


# --------------------------------------------------------------------- main

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("package", type=pathlib.Path,
                    help="directory holding the manifest and motif CIF")
    ap.add_argument("--manifest", default=None,
                    help="manifest filename inside the package (default: the "
                         "only *.json, or *manifest*.json)")
    ap.add_argument("--repo", type=pathlib.Path, default=None,
                    help="an alphaprotein-novo checkout, to cross-check with "
                         "its own validator and sampler")
    ap.add_argument("--reference-residues", type=int, default=None,
                    help="how many motif residues the package claims to fix")
    args = ap.parse_args(argv)

    if args.manifest:
        manifest_path = args.package / args.manifest
    else:
        candidates = sorted(args.package.glob("*manifest*.json")) or \
            sorted(args.package.glob("*.json"))
        if not candidates:
            print(f"no manifest JSON found in {args.package}")
            return 2
        manifest_path = candidates[0]

    spec = json.loads(manifest_path.read_text())
    defaults = spec.get("defaults", {})
    print("=" * 78)
    print(f"AP Novo package audit: {args.package}")
    print("=" * 78)
    print(f"\n  manifest: {manifest_path.name}")
    print(f"  jobs    : {[j['name'] for j in spec.get('designs', [])]}")

    # --- 1. the motif CIF ---------------------------------------------------
    cif_name = defaults.get("input_file") or next(
        (j.get("input_file") for j in spec["designs"] if j.get("input_file")), None)
    info = None
    if cif_name and (args.package / cif_name).exists():
        cif_path = args.package / cif_name
        print("\n" + "=" * 78)
        print(f"[1] The motif CIF: {cif_name}")
        print("=" * 78 + "\n")
        cif_lines = cif_path.read_text().splitlines()
        header_problem = check_cif_header(cif_lines)
        missing_categories = check_cif_shape(cif_lines)
        if header_problem:
            print(f"  ❌ malformed mmCIF header -- {header_problem}")
        else:
            print("  ✅ mmCIF data-block header is well formed.")
        if missing_categories:
            print(f"  ❌ not a loadable mmCIF: missing {missing_categories}")
            print(f"     ({len(cif_categories(cif_lines))} categories present; a motif"
                  " file is a complete mmCIF,")
            print("     not an atom list. AF3 builds author_naming_scheme from the")
            print("     entity and pdbx_*_scheme tables, and refuses without them.)")
        else:
            print("  ✅ every mmCIF category the working example carries is present.")
        if header_problem or missing_categories:
            print("     Everything below is read with a permissive parser instead.")
        print()
        info = read_motif_cif(cif_path, strict=False)
        print(f"  {info['n_atoms']} atoms | {len(info['residues'])} protein residues"
              f" | {len(info['ligand_atoms'])} ligand atoms")
        numbering = check_numbering(info["residues"])
        bad = [r for r in numbering if not r["match"]]
        print(f"  residues: "
              f"{', '.join(r['comp'] + str(r['res_id']) for r in numbering)}")
        if bad:
            print(f"\n  ❌ {len(bad)}/{len(numbering)} do NOT match wild-type MmPylRS"
                  " (Q8PWY1):")
            for r in bad:
                print(f"       CIF has {r['comp']}{r['res_id']}"
                      f" ({r['one_letter']}), Q8PWY1 has {r['reference']}")
        else:
            print(f"\n  ✅ {len(numbering)}/{len(numbering)} match wild-type MmPylRS"
                  " (Q8PWY1) at their stated positions.")
            print("     The numbering is full-length Mm numbering, so these"
                  " coordinates are")
            print("     directly comparable with the COM1 transplant positions in"
                  " demo1_multisite.")

        ligand_comps = sorted({a["comp"] for a in info["ligand_atoms"]})
        if ligand_comps:
            print(f"\n  ligand: {ligand_comps} ({len(info['ligand_atoms'])} atoms)")
            moiety = partition_by_moiety(info)
            if moiety:
                groups = {side: [r["residue"] for r in moiety if r["side"] == side]
                          for side in ("ncAA", "AMP", "both")}
                print("\n  ⚠️ CCD YLY is pyrrolysyl-ADENYLATE (C22 H35 N8 O9 P), not the")
                print("     free ncAA -- 2Q7H is 'bound to adenylated pyrrolysine'. So a")
                print("     motif taken from it is partly an ATP-site motif:\n")
                print("       residue   d(ncAA)  d(AMP)  margin  nearer")
                for r in moiety:
                    print(f"       {r['residue']:<9s} {r['d_ncaa']:7.2f}"
                          f"  {r['d_amp']:6.2f}  {r['margin']:6.2f}  {r['side']}")
                for side, label in (("ncAA", "the amino-acid half"),
                                    ("AMP", "the nucleotide half"),
                                    ("both", f"BOTH, within {MOIETY_TIE_ANGSTROM} A")):
                    members = groups[side]
                    print(f"\n     nearer {label} ({len(members)}):"
                          f" {', '.join(members) or '-'}")
                not_ncaa = len(groups["AMP"]) + len(groups["both"])
                print(f"\n  ❌ only {len(groups['ncAA'])}/{len(moiety)} of a motif"
                      " described as 'the ncAA binding pocket' are")
                print(f"     unambiguously in it; {not_ncaa} touch the nucleotide half or"
                      " both. Conditioning")
                print("     on this motif asks the model to rebuild the amino-acid pocket")
                print("     AND the ATP site at once.")
                if groups["both"]:
                    print(f"  ⚠️ {len(groups['both'])} residue(s) sit within"
                          f" {MOIETY_TIE_ANGSTROM} A of both halves, so a clean")
                    print("     two-way split of this pocket does not exist. Any 'ncAA")
                    print("     pocket only' motif is a judgement call, not a measurement.")

    # --- 2. the motif strings ----------------------------------------------
    print("\n" + "=" * 78)
    print("[2] Does each motif_str mean what it looks like?")
    print("=" * 78)
    declared = args.reference_residues
    if declared is None and info is not None:
        declared = len(info["residues"])

    failures = 0
    for job in spec.get("designs", []):
        seq_length = job.get("seq_length", defaults.get("seq_length"))
        unindexed = job.get("unindexed_motif_residues")
        print(f"\n  --- {job['name']} ---")
        print(f"  motif_str : {job['motif_str'][:95]}"
              f"{'...' if len(job['motif_str']) > 95 else ''}")
        print(f"  seq_length: {seq_length}")
        if unindexed:
            n_unindexed = len([s for s in unindexed.split(",") if s])
            print(f"  unindexed_motif_residues: {n_unindexed} residues"
                  " (correct place for a pocket)")
        check = check_motif_str(job["motif_str"], seq_length,
                                None if unindexed else declared)
        print(f"  AP Novo reads {check['n_fixed']} fixed motif residue(s)"
              f": {check['fixed_residues'] or '[]'}")
        if check["exact_length_segments"]:
            segs = check["exact_length_segments"]
            print(f"  ❌ {len(segs)} bare integer(s) read as EXACT designable lengths,"
                  f" summing to {sum(segs)}:")
            print(f"       {segs}")
        print(f"  designable total is forced into [{check['designable_min']},"
              f" {check['designable_max']}]")
        if check.get("silently_dropped"):
            print(f"  ❌ {check['silently_dropped']} of the {declared} residues in the"
                  " CIF are NOT conditioned on.")
        if check["ok"]:
            print("  ✅ satisfiable")
        else:
            failures += 1
            print(f"  ❌ UNSATISFIABLE -> {check['error']}")

    # --- 2b. cross-field consistency ---------------------------------------
    resequence = spec.get("resequence")
    folding_inputs = (spec.get("folding") or {}).get("inputs")
    print("\n" + "=" * 78)
    print("[2b] Cross-field consistency")
    print("=" * 78 + "\n")
    if resequence is None:
        failures += 1
        print("  ❌ no `resequence` block. run_pipeline.py refuses to guess whether")
        print("     designs are resequenced before folding -- set it to true or false.")
    elif folding_inputs and "resequenced" in folding_inputs \
            and not resequence.get("enabled"):
        failures += 1
        print('  ❌ folding.inputs includes "resequenced" while resequence.enabled is')
        print("     false. The pipeline refuses this pair: enable resequencing or fold")
        print('     the generated structure ("generated") instead.')
    else:
        print("  ✅ resequence.enabled and folding.inputs agree.")
    if resequence and resequence.get("enabled"):
        print("  ⚠️ resequence.enabled is true, so --ligandmpnn_dir and")
        print("     --ligandmpnn_python are required at launch or run_pipeline.py")
        print("     raises app.UsageError before any GPU work starts.")

    # --- 3. upstream cross-check -------------------------------------------
    if args.repo:
        print("\n" + "=" * 78)
        print(f"[3] Cross-check against the real code in {args.repo}")
        print("=" * 78 + "\n")
        up = upstream_check(args.repo, manifest_path)
        if "schema" in up:
            print(f"  Manifest.from_file : {up['schema']}")
        if up.get("schema_error"):
            print(f"  Manifest.from_file : unavailable ({up['schema_error'][:90]})")
        for name, verdict in (up.get("sampler") or {}).items():
            print(f"  sample_designable_lengths[{name}]: {verdict}")
        if up.get("sampler_error"):
            print(f"  sampler unavailable: {up['sampler_error'][:90]}")
            print("  (the repo imports jax and alphafold3 at import time; the"
                  " grammar checks")
            print("   above do not need either)")
        if up.get("schema") == "PASSED" and failures:
            print("\n  ⚠️ The schema PASSES while the campaign cannot run. That is the")
            print("     whole reason for this audit: pydantic validates types, and the")
            print("     motif grammar is only checked once the sampler runs, on a GPU")
            print("     node, after the weights have loaded.")

    print("\n" + "=" * 78)
    print("Verdict")
    print("=" * 78)
    total = len(spec.get("designs", []))
    if failures:
        print(f"  ❌ {failures}/{total} job(s) cannot run as written.")
        print("  See manifest_fixed.json in this directory for the corrected form,")
        print("  and README.md for why each change was needed.")
        return 1
    print(f"  ✅ all {total} job(s) are satisfiable as written.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
