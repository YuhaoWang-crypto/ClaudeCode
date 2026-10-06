"""Build a runnable AP Novo package from the audited one.

Three changes, each forced by a finding in `audit_apnovo_package.py`:

1. **Put pocket residues where AP Novo reads them.** A residue list written as
   `A300,302,305,...` fixes one residue and requests ~6,500 designed ones,
   because any segment starting with a digit is a designable length. Both
   corrected forms are emitted: the `unindexed_motif_residues` field (the mode
   the repo documents for enzyme design, where the model places the motif), and
   a correctly interleaved indexed `motif_str` for comparison.

2. **Make the length fields consistent.** `seq_length` bounds the **sum of the
   designable segments**, not the finished protein -- `sample_designable_lengths`
   passes only the designable ranges into `sample_values_with_sum_in_range`. So
   `motif_str: "60-250"` with `seq_length: "260-360"` is a contradiction and
   raises. The generated manifests keep the two in agreement by construction.

3. **Split the motif by what the ligand half actually is.** The 2Q7H ligand is
   pyrrolysyl-adenylate (CCD `YLY`, C22 H35 N8 O9 P), so 10 of the 19 pocket
   residues are nearer the nucleotide than the amino acid. A job aimed at ncAA
   specificity gets the 9-residue motif; the full 19 stays available for a job
   that really means "rebuild the whole active site".

⚠️ The ligand is NOT trimmed. Removing the AMP half would leave an atom set that
no longer matches CCD `YLY`, and `folding.states` builds its ligand from that
code -- so a trimmed ligand needs a custom CCD entry (or a user-supplied
component CIF), which is a decision for whoever runs this, not a default. The
consequence is stated in the README: even the 9-residue job still designs around
the adenylate.

Run:  python3 -m trnaplat.apnovo.build_fixed_package <source_package> [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent.parent))

from trnaplat.apnovo.audit_apnovo_package import (  # noqa: E402
    partition_by_moiety,
    read_motif_cif,
)

#: Designable linker range placed between consecutive fixed residues in the
#: indexed example. Wide enough that the sum can reach the target length, narrow
#: enough that the motif is not scattered across an arbitrary fold.
LINKER = (10, 45)


def ncaa_side_residues(info: dict) -> list[int]:
    """Residue ids unambiguously nearer the adenylate's amino-acid half.

    Residues within `MOIETY_TIE_ANGSTROM` of both halves are excluded: they
    contact the nucleotide too, so including them would smuggle the ATP site
    into a motif advertised as the ncAA pocket.
    """
    return [r["res_id"] for r in partition_by_moiety(info) if r["side"] == "ncAA"]


def normalise_header(lines: list[str]) -> list[str]:
    """Repair a `data <name>` header into the mmCIF-legal `data_<name>`.

    ✅ The audited package ships `data pylrs_pyl_motif`. AlphaFold 3's parser
    refuses it with "The CIF file does not start with the data_ field", so every
    job fails at structure load — after the manifest has validated and the
    motif grammar has been checked. Repaired here rather than passed through.
    """
    out = list(lines)
    for i, line in enumerate(out):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("data_") and len(stripped) > len("data_"):
            return out
        parts = stripped.split(None, 1)
        if parts and parts[0] == "data" and len(parts) == 2:
            out[i] = f"data_{parts[1].strip()}"
        else:
            out[i] = "data_motif"
            out.insert(i + 1 if stripped else i, line)
        return out
    return ["data_motif"]


#: The `_atom_site` columns AlphaFold 3's parser reads, plus the author and
#: bookkeeping columns a well-formed mmCIF carries. ✅ Enumerated from
#: alphafold3/structure/parsing.py rather than guessed: the audited package
#: supplies only 12 columns and is missing `label_entity_id`, `occupancy`,
#: `pdbx_PDB_ins_code` and `pdbx_PDB_model_num`, which fails with
#: `KeyError: '_atom_site.pdbx_PDB_model_num'` after the data_ header is fixed.
ATOM_SITE_COLUMNS = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
    "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
    "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
    "auth_atom_id", "pdbx_PDB_model_num",
]

AMINO_ACIDS = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP "
    "TYR VAL".split()
)


def write_motif_cif(source: pathlib.Path, target: pathlib.Path,
                    keep_res_ids: set[int] | None, entry_id: str) -> int:
    """Rewrite a motif as a complete mmCIF, optionally keeping only some residues.

    Rebuilt rather than patched, because the source is missing four columns the
    parser requires and `label_seq_id` has to be renumbered: AF3 wants
    `label_seq_id` contiguous per entity while the motif's identity lives in
    `auth_seq_id` (the manifest sets `is_author_naming: true` and names A300,
    A302, …). Copying the author numbering into both fields is what a naive
    patch would do, and it is wrong.
    """
    lines = source.read_text().splitlines()
    header = [l.strip().split(".", 1)[1] for l in lines
              if l.strip().startswith("_atom_site.")]
    col = {name: i for i, name in enumerate(header)}

    def get(row, *names, default=None):
        for name in names:
            if name in col:
                return row[col[name]]
        return default

    rows = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        row = line.split()
        comp = get(row, "label_comp_id", "auth_comp_id")
        auth_seq = int(get(row, "auth_seq_id", "label_seq_id"))
        if comp in AMINO_ACIDS and keep_res_ids is not None \
                and auth_seq not in keep_res_ids:
            continue
        rows.append({
            "group_PDB": get(row, "group_PDB"),
            "type_symbol": get(row, "type_symbol"),
            "label_atom_id": get(row, "label_atom_id"),
            "label_comp_id": comp,
            "label_asym_id": get(row, "label_asym_id", "auth_asym_id"),
            "auth_asym_id": get(row, "auth_asym_id", "label_asym_id"),
            "auth_seq_id": auth_seq,
            "Cartn_x": get(row, "Cartn_x"), "Cartn_y": get(row, "Cartn_y"),
            "Cartn_z": get(row, "Cartn_z"),
            "B_iso_or_equiv": get(row, "B_iso_or_equiv", default="1.00"),
            "is_polymer": comp in AMINO_ACIDS,
        })

    # label_entity_id per chain, label_seq_id contiguous within each polymer
    # chain, "." for non-polymer atoms.
    chains = sorted({r["label_asym_id"] for r in rows})
    entity_of = {chain: str(i + 1) for i, chain in enumerate(chains)}
    seq_of: dict[tuple[str, int], int] = {}
    for chain in chains:
        residues = sorted({r["auth_seq_id"] for r in rows
                           if r["label_asym_id"] == chain and r["is_polymer"]})
        for index, auth in enumerate(residues, start=1):
            seq_of[(chain, auth)] = index

    out = [f"data_{entry_id}", f"_entry.id {entry_id}", "#", "loop_"]
    out += [f"_atom_site.{name}" for name in ATOM_SITE_COLUMNS]
    for serial, r in enumerate(rows, start=1):
        label_seq = (str(seq_of[(r["label_asym_id"], r["auth_seq_id"])])
                     if r["is_polymer"] else ".")
        values = {
            "group_PDB": r["group_PDB"], "id": str(serial),
            "type_symbol": r["type_symbol"],
            "label_atom_id": r["label_atom_id"], "label_alt_id": ".",
            "label_comp_id": r["label_comp_id"],
            "label_asym_id": r["label_asym_id"],
            "label_entity_id": entity_of[r["label_asym_id"]],
            "label_seq_id": label_seq, "pdbx_PDB_ins_code": "?",
            "Cartn_x": r["Cartn_x"], "Cartn_y": r["Cartn_y"],
            "Cartn_z": r["Cartn_z"], "occupancy": "1.00",
            "B_iso_or_equiv": r["B_iso_or_equiv"],
            "auth_seq_id": str(r["auth_seq_id"]),
            "auth_comp_id": r["label_comp_id"],
            "auth_asym_id": r["auth_asym_id"],
            "auth_atom_id": r["label_atom_id"], "pdbx_PDB_model_num": "1",
        }
        out.append(" ".join(values[name] for name in ATOM_SITE_COLUMNS))
    out.append("#")
    target.write_text("\n".join(out) + "\n")
    return len(rows)


def indexed_motif_str(res_ids: list[int], chain: str = "A",
                      ligand_chain: str = "B") -> tuple[str, tuple[int, int]]:
    """A correctly interleaved indexed motif_str, plus its designable (min, max).

    Alternates designable linker / fixed residue / ... / designable linker, so
    every fixed residue carries its chain letter and no bare integer is left to
    be read as a length.
    """
    lo, hi = LINKER
    parts = [f"{lo}-{hi}"]
    for res_id in res_ids:
        parts.append(f"{chain}{res_id}")
        parts.append(f"{lo}-{hi}")
    n_linkers = len(res_ids) + 1
    return "/".join([",".join(parts), f"{ligand_chain}1"]), (lo * n_linkers,
                                                             hi * n_linkers)


def clamp(target: tuple[int, int], feasible: tuple[int, int]) -> str:
    """A seq_length inside both the requested window and what the ranges allow."""
    lo = max(target[0], feasible[0])
    hi = min(target[1], feasible[1])
    if lo > hi:
        lo, hi = feasible
    return f"{lo}-{hi}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", type=pathlib.Path,
                    help="the audited package directory")
    ap.add_argument("--manifest", default="pylrs_pyl_manifest.json")
    ap.add_argument("--cif", default="pylrs_pyl_motif.cif")
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "fixed")
    ap.add_argument("--target-length", default="260-360")
    ap.add_argument("--pdb", type=pathlib.Path, default=None,
                    help="deposited mmCIF to re-extract the motif from "
                         "(e.g. 2q7h.cif). Strongly recommended: without it the "
                         "output carries only an _atom_site loop, which AF3 "
                         "cannot load.")
    ap.add_argument("--ligand", default="YLY")
    args = ap.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    src_cif = args.source / args.cif
    # strict=False: this builder repairs a malformed data-block header rather
    # than refusing it, so it must be able to read the broken original.
    info = read_motif_cif(src_cif, strict=False)
    all_ids = sorted(r for _, r in info["residues"])
    ncaa_ids = ncaa_side_residues(info)
    target = tuple(int(x) for x in args.target_length.split("-"))

    # --- the motif files ---------------------------------------------------
    # Rebuilt from the deposited entry when one is given, because a loadable
    # motif is a complete mmCIF (entity + pdbx_*_scheme categories), not an
    # atom list -- see extract_motif for the three failures that established it.
    if args.pdb:
        from trnaplat.apnovo import extract_motif as em
        plans = [("motif_active_site_full.cif", all_ids, "pylrs_active_site"),
                 ("motif_ncaa_pocket.cif", ncaa_ids, "pylrs_ncaa_pocket")]
        print(f"motif files in {args.out}, extracted from {args.pdb.name} with gemmi:")
        for name, ids, entry in plans:
            structure, meta = em.extract(args.pdb, ids, args.ligand, entry)
            em.write(structure, args.out / name)
            missing = [c for c in em.REQUIRED_CATEGORIES
                       if c not in em.categories(args.out / name)]
            status = "✅" if not missing else f"❌ missing {missing}"
            print(f"  {name:28s} {len(ids)} residues,"
                  f" {meta['protein_atoms']} protein + {meta['ligand_atoms']}"
                  f" ligand atoms  {status}")
    else:
        n_full = write_motif_cif(src_cif, args.out / "motif_active_site_full.cif",
                                 None, "pylrs_active_site")
        n_kept = write_motif_cif(src_cif, args.out / "motif_ncaa_pocket.cif",
                                 set(ncaa_ids), "pylrs_ncaa_pocket")
        print(f"⚠️ no --pdb given, so the motifs are rewritten from"
              f" {src_cif.name} with an")
        print(f"   {len(ATOM_SITE_COLUMNS)}-column _atom_site loop only. AF3 needs the"
              " entity and")
        print("   pdbx_*_scheme categories too, so pass --pdb for a loadable file.")
        print(f"  motif_active_site_full.cif  {len(all_ids)} residues, {n_full} atoms")
        print(f"  motif_ncaa_pocket.cif       {len(ncaa_ids)} residues, {n_kept} atoms")
    print(f"  ncAA-side residues (unambiguous): {ncaa_ids}")

    # --- the fixed manifest ------------------------------------------------
    unindexed_full = ",".join(f"A{r}" for r in all_ids)
    unindexed_ncaa = ",".join(f"A{r}" for r in ncaa_ids)
    idx_motif, idx_feasible = indexed_motif_str(ncaa_ids)
    idx_seq_length = clamp(target, idx_feasible)
    un_seq_length = f"{target[0]}-{target[1]}"

    # ⚠️ No top-level comment key: `Manifest` rejects unknown top-level keys
    # ("Unknown top level manifest keys: ['_comment']"), which the upstream
    # cross-check in audit_apnovo_package.py caught. Rationale lives in the
    # README and in each job's `description`, which is a real schema field.
    manifest = {
        "settings": {"model_dir": "./models/apnovo_generator"},
        "defaults": {
            "is_author_naming": True,
            "num_sampling_steps": 1000,
            "num_designs": 50,
        },
        "designs": [
            {
                "name": "ncaa_pocket_unindexed",
                "description": (
                    "9 residues nearer the amino-acid half of the adenylate. "
                    "Unindexed: the model places the motif."
                ),
                "input_file": "motif_ncaa_pocket.cif",
                "motif_str": f"{un_seq_length}/B1",
                "unindexed_motif_residues": unindexed_ncaa,
                "seq_length": un_seq_length,
            },
            {
                "name": "active_site_full_unindexed",
                "description": (
                    "All 19 residues, i.e. amino-acid pocket AND adenylate site."
                ),
                "input_file": "motif_active_site_full.cif",
                "motif_str": f"{un_seq_length}/B1",
                "unindexed_motif_residues": unindexed_full,
                "seq_length": un_seq_length,
            },
            {
                "name": "ncaa_pocket_indexed",
                "description": (
                    "Same 9 residues, indexed: fixed residues interleaved with "
                    "designable linkers. Shows the grammar the broken manifest "
                    "was reaching for."
                ),
                "input_file": "motif_ncaa_pocket.cif",
                "motif_str": idx_motif,
                "seq_length": idx_seq_length,
            },
        ],
        "resequence": {"enabled": True, "temperature": 0.1},
        "folding": {
            "inputs": ["resequenced"],
            "seeds": [230],
            "states": [
                {"name": "apo_monomer", "ligands": []},
                {"name": "adenylate_complex", "ligands": [
                    {"id": "B", "ccd_code": "YLY"}]},
            ],
        },
    }
    (args.out / "manifest_fixed.json").write_text(
        json.dumps(manifest, indent=2) + "\n")

    partial = {
        "settings": {"model_dir": "./models/apnovo_generator"},
        "defaults": {"is_author_naming": True, "num_designs": 50},
        "designs": [
            {
                "name": "pylrs_partial_diffusion",
                "description": (
                    "Partial diffusion from the real structure: preserves the "
                    "measured pocket geometry instead of asking a de novo "
                    "scaffold to rediscover it. Point "
                    "partial_diffusion_input_file at a full 2Q7H chain A, not "
                    "at the motif fragment."
                ),
                "input_file": "motif_ncaa_pocket.cif",
                "motif_str": f"{un_seq_length}/B1",
                "unindexed_motif_residues": unindexed_ncaa,
                "seq_length": un_seq_length,
                "partial_diffusion_input_file": "REPLACE_WITH_2q7h_chainA.cif",
                "partial_diffusion_num_steps": 200,
            }
        ],
        "resequence": {"enabled": True, "temperature": 0.1},
        "folding": {
            "inputs": ["resequenced"],
            "seeds": [230],
            "states": [
                {"name": "apo_monomer", "ligands": []},
                {"name": "adenylate_complex", "ligands": [
                    {"id": "B", "ccd_code": "YLY"}]},
            ],
        },
    }
    (args.out / "manifest_partial_diffusion.json").write_text(
        json.dumps(partial, indent=2) + "\n")

    commands = f"""#!/bin/bash
# Corrected AP Novo launch. Differences from the audited script, each one a
# failure it would have hit:
#   * LigandMPNN is installed and its two flags are passed, because the manifest
#     sets resequence.enabled -- run_pipeline.py raises app.UsageError otherwise.
#   * the "quick" run really is quick: --only_stage=generation on a cut-down
#     manifest, instead of 1000 sampling steps x 50 designs.
set -euo pipefail

# --- 1. AP Novo environment (per the repo's own README) ---------------------
uv venv --python 3.12 .venv
source .venv/bin/activate
uv pip install -U "jax[cuda12]>=0.4.30"      # CPU-only: uv pip install -U jax
uv pip install git+https://github.com/google-deepmind/alphafold3.git
uv pip install -e .
build_data                                    # compiles ccd.pickle

# --- 2. LigandMPNN environment (separate, PyTorch) -------------------------
# Required because manifest.resequence.enabled is true. Drop this section only
# if you also set {{"resequence": {{"enabled": false}}}}.
git clone https://github.com/dauparas/LigandMPNN.git
( cd LigandMPNN \\
  && uv venv --python 3.12 .venv_ligandmpnn \\
  && uv pip install --python .venv_ligandmpnn/bin/python -r requirements.txt \\
  && uv pip install --python .venv_ligandmpnn/bin/python "setuptools<82" \\
  && bash get_model_params.sh ./model_params )
export LIGANDMPNN_DIR="$PWD/LigandMPNN"
export LIGANDMPNN_PYTHON="$PWD/LigandMPNN/.venv_ligandmpnn/bin/python"

# --- 3. Weights ------------------------------------------------------------
# ✅ verified reachable: 547,140,096 B and 1,020,518,661 B, HTTP 200.
mkdir -p models/apnovo_generator models/af3_la
wget -nc -P models/apnovo_generator \\
  https://storage.googleapis.com/alphaprotein_novo/generator.bin.zst
wget -nc -P models/af3_la \\
  https://storage.googleapis.com/alphafold3/af3_leaving_atom.bin.zst

# --- 4. Smoke test: generation only, few steps, one design ----------------
python - <<'PY'
import json, pathlib
m = json.loads(pathlib.Path('manifest_fixed.json').read_text())
m['defaults'].update(num_sampling_steps=50, num_designs=1)
m['designs'] = m['designs'][:1]
pathlib.Path('manifest_smoke.json').write_text(json.dumps(m, indent=2))
PY
python run_pipeline.py --manifest=manifest_smoke.json --output_dir=./out_smoke \\
  --only_stage=generation --apn_model_dir=./models/apnovo_generator

# --- 5. Full campaign ------------------------------------------------------
python run_pipeline.py --manifest=manifest_fixed.json --output_dir=./out_full \\
  --apn_model_dir=./models/apnovo_generator \\
  --af3_model_dir=./models/af3_la \\
  --ligandmpnn_dir="$LIGANDMPNN_DIR" \\
  --ligandmpnn_python="$LIGANDMPNN_PYTHON"
"""
    (args.out / "run_commands_fixed.sh").write_text(commands)
    (args.out / "run_commands_fixed.sh").chmod(0o755)

    print(f"\nmanifests in {args.out}:")
    for job in manifest["designs"]:
        print(f"  {job['name']:28s} motif_str={job['motif_str'][:58]}"
              f"{'...' if len(job['motif_str']) > 58 else ''}")
        print(f"  {'':28s} seq_length={job['seq_length']}")
    print(f"  manifest_partial_diffusion.json  (needs a full 2Q7H chain A CIF)")
    print(f"  run_commands_fixed.sh")
    print("\nVerify with:")
    print(f"  python3 -m trnaplat.apnovo.audit_apnovo_package {args.out}"
          " --manifest manifest_fixed.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
