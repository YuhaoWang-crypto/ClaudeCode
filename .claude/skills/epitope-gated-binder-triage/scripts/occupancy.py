#!/usr/bin/env python3
"""Epitope occupancy from predicted binder:target complexes.

Every agent that analysed this campaign independently rewrote this contact
loop, usually correctly and occasionally not. It is bundled so the geometry and
— more importantly — the verification assertions are identical every time.

The assertions are the point. A residue-numbering shift produces occupancy
numbers that are well-formatted, plausible, and wrong, with no error anywhere.
`--target-seq` makes that failure loud.

Usage
-----
    python occupancy.py --cif-dir ./cif \\
        --target-chain A --binder-chain NANO1 \\
        --target-start 47 \\
        --target-seq LPQLCKFCDVRLSTC... \\
        --epitope 50,53,55,72,73,74,75,76,77,78,100,142 \\
        --hotspots 76,142 \\
        --out occupancy.tsv

Add `--paratope` to also emit the binder-side contact residues (needed when
checking that a framework repair falls outside the interface).

Outputs one TSV row per structure. Files failing verification are not scored;
they are counted and listed on stderr so a silent partial analysis is impossible.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Field indices in whitespace-split mmCIF ATOM records as emitted by Boltz-2.
# Verify on the first file of any new run rather than trusting these blindly --
# a different producer may lay the loop out differently.
F_RESNAME, F_RESNUM, F_CHAIN, F_X, F_Y, F_Z = 5, 6, 9, 10, 11, 12


def parse_cif(path: Path):
    """-> {chain: {resnum: (one_letter, [(x,y,z), ...])}}"""
    chains: dict[str, dict[int, list]] = defaultdict(dict)
    with path.open() as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            f = line.split()
            resname = f[F_RESNAME]
            if resname not in THREE_TO_ONE:
                continue
            try:
                resnum = int(f[F_RESNUM])
                xyz = (float(f[F_X]), float(f[F_Y]), float(f[F_Z]))
            except (ValueError, IndexError):
                continue
            entry = chains[f[F_CHAIN]].setdefault(resnum, [THREE_TO_ONE[resname], []])
            entry[1].append(xyz)
    return chains


def verify(chains, target_chain, binder_chain, target_start, target_seq):
    """Return None if the structure is usable, else a reason string."""
    if target_chain not in chains:
        return f"missing target chain {target_chain!r}"
    if binder_chain not in chains:
        return f"missing binder chain {binder_chain!r}"
    ids = sorted(chains[target_chain])
    if target_start is not None and ids[0] != target_start:
        return f"target starts at {ids[0]}, expected {target_start}"
    if target_seq:
        observed = "".join(chains[target_chain][i][0] for i in ids)
        if observed != target_seq:
            return (
                f"target sequence mismatch (len {len(observed)} vs "
                f"{len(target_seq)}) -- numbering or crop differs"
            )
    return None


def contacts(chains, target_chain, binder_chain, cutoff, want_paratope=False):
    """-> (target_residue -> n_contacting_atoms, binder_residues_in_contact)"""
    c2 = cutoff * cutoff
    binder_atoms = [xyz for _, atoms in chains[binder_chain].values() for xyz in atoms]
    tgt_counts: dict[int, int] = {}
    paratope: set[int] = set()
    for resnum, (_aa, atoms) in chains[target_chain].items():
        n = 0
        for a in atoms:
            for b in binder_atoms:
                dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
                if dx * dx + dy * dy + dz * dz < c2:
                    n += 1
                    break
        if n:
            tgt_counts[resnum] = n
    if want_paratope:
        tgt_atoms = [xyz for _, atoms in chains[target_chain].values() for xyz in atoms]
        for resnum, (_aa, atoms) in chains[binder_chain].items():
            hit = False
            for a in atoms:
                for b in tgt_atoms:
                    dx, dy, dz = a[0] - b[0], a[1] - b[1], a[2] - b[2]
                    if dx * dx + dy * dy + dz * dz < c2:
                        hit = True
                        break
                if hit:
                    break
            if hit:
                paratope.add(resnum)
    return tgt_counts, paratope


def classify_off_epitope(off, epitope, chains, target_chain, rim_seq=2, rim_ang=8.0):
    """Split off-epitope contacts into rim (adjacent to the epitope) and remote.

    A rim contact is the normal edge of a real interface. A remote contact is a
    second binding patch elsewhere on the surface -- a different binding mode.
    """
    rim, remote = [], []
    ep_atoms = [
        xyz
        for e in epitope
        if e in chains[target_chain]
        for xyz in chains[target_chain][e][1]
    ]
    for r in sorted(off):
        if any(abs(r - e) <= rim_seq for e in epitope):
            rim.append(r)
            continue
        near = False
        for a in chains[target_chain][r][1]:
            for b in ep_atoms:
                if math.dist(a, b) <= rim_ang:
                    near = True
                    break
            if near:
                break
        (rim if near else remote).append(r)
    return rim, remote


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cif-dir", required=True, type=Path)
    p.add_argument("--glob", default="*.cif")
    p.add_argument("--target-chain", default="A")
    p.add_argument("--binder-chain", default="NANO1")
    p.add_argument("--target-start", type=int, default=None,
                   help="expected first residue number of the target chain")
    p.add_argument("--target-seq", default="",
                   help="expected target chain sequence; mismatch = not scored")
    p.add_argument("--epitope", required=True,
                   help="comma-separated epitope residue numbers (target numbering)")
    p.add_argument("--hotspots", default="",
                   help="comma-separated subset of the epitope treated as hotspots")
    p.add_argument("--cutoff", type=float, default=4.5)
    p.add_argument("--paratope", action="store_true",
                   help="also emit binder-side contact residues")
    p.add_argument("--out", required=True, type=Path)
    a = p.parse_args()

    epitope = {int(x) for x in a.epitope.split(",") if x.strip()}
    hotspots = {int(x) for x in a.hotspots.split(",") if x.strip()}
    if not hotspots <= epitope:
        print("warning: hotspots not a subset of the epitope", file=sys.stderr)

    cols = [
        "structure", "n_contact", "n_epi", "recall", "purity", "atom_purity",
        "n_pairs", "hotspots_contacted", "n_hotspots", "n_rim", "n_remote",
        "contacted", "off_epitope", "remote_residues",
    ]
    if a.paratope:
        cols.append("paratope")

    files = sorted(a.cif_dir.glob(a.glob))
    if not files:
        print(f"no files matching {a.glob} in {a.cif_dir}", file=sys.stderr)
        return 1

    failures: list[tuple[str, str]] = []
    rows: list[list[str]] = []

    for path in files:
        chains = parse_cif(path)
        why = verify(chains, a.target_chain, a.binder_chain, a.target_start, a.target_seq)
        if why:
            failures.append((path.stem, why))
            continue

        tgt_counts, paratope = contacts(
            chains, a.target_chain, a.binder_chain, a.cutoff, want_paratope=a.paratope
        )
        contacted = set(tgt_counts)
        hit = contacted & epitope
        off = contacted - epitope
        rim, remote = classify_off_epitope(off, epitope, chains, a.target_chain)
        epi_atoms = sum(tgt_counts[r] for r in hit)
        all_atoms = sum(tgt_counts.values()) or 1
        hs = sorted(hit & hotspots)

        row = [
            path.stem,
            str(len(contacted)),
            str(len(hit)),
            f"{len(hit)/len(epitope):.4f}",
            f"{len(hit)/max(1, len(contacted)):.4f}",
            f"{epi_atoms/all_atoms:.4f}",
            str(all_atoms),
            ";".join(map(str, hs)),
            str(len(hs)),
            str(len(rim)),
            str(len(remote)),
            ";".join(map(str, sorted(contacted))),
            ";".join(map(str, sorted(off))),
            ";".join(map(str, remote)),
        ]
        if a.paratope:
            row.append(";".join(map(str, sorted(paratope))))
        rows.append(row)

    with a.out.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"scored {len(rows)}/{len(files)} structures -> {a.out}", file=sys.stderr)
    if failures:
        print(f"VERIFICATION FAILURES: {len(failures)} (not scored)", file=sys.stderr)
        for name, why in failures[:20]:
            print(f"  {name}: {why}", file=sys.stderr)
        if len(failures) > 20:
            print(f"  ... and {len(failures)-20} more", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
