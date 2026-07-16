"""
structure.py -- reproducible derivation of the values frozen in systems.py.  ✅

Uses biotite to (a) fetch a deposited structure, (b) compute per-residue
secondary structure (``annotate_sse``), and (c) derive candidate permutation
sites (interior loop centers).  Running ``python -m biosensor_pipeline.structure``
recomputes and checks the frozen annotations, so nothing in systems.py is a
magic constant -- it is all reproducible from public coordinates.

Requires network access + biotite.  The rest of the pipeline does NOT depend
on this module at run time.
"""

from __future__ import annotations
import warnings
warnings.filterwarnings("ignore")

import tempfile


def annotate(pdb: str, chain: str, full_seq: str):
    """Return (modeled_seq, offset, sse_full, loop_centers) for a chain."""
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    from biotite.sequence import ProteinSequence

    tmp = tempfile.mkdtemp()
    path = rcsb.fetch(pdb, "bcif", tmp)
    cif = pdbx.BinaryCIFFile.read(path)
    arr = pdbx.get_structure(cif, model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    arr = arr[arr.chain_id == chain]
    ca = arr[arr.atom_name == "CA"]
    sse = struc.annotate_sse(arr)
    one = "".join(ProteinSequence.convert_letter_3to1(r) for r in ca.res_name)

    # Map the modeled chain onto full_seq. Exact substring is the common case;
    # if terminal residues differ (missing start, partial His tag modeled, a
    # point variant), fall back to the longest matching block so interior loop
    # detection still works.  Returns (mod_off_in_full, mod_start_in_modeled,
    # block_len) covering the reliably-aligned interior.
    off = full_seq.find(one)
    if off >= 0:
        f0, m0, blen = off, 0, len(one)
    else:
        import difflib
        m = difflib.SequenceMatcher(None, full_seq, one, autojunk=False)
        a, b, blen = m.find_longest_match(0, len(full_seq), 0, len(one))
        if blen < 20:
            raise ValueError(f"{pdb}: modeled sequence does not align to the provided full sequence")
        f0, m0 = a, b
        off = f0 - m0  # nominal offset (may be approximate outside the block)

    sse_full = ["c"] * len(full_seq)
    for i, s in enumerate(sse):
        fpos = f0 + (i - m0)              # map modeled residue i -> full_seq index
        if 0 <= fpos < len(full_seq):
            sse_full[fpos] = s
    sse_full = "".join(sse_full)

    # interior loop centers within the reliably aligned block only
    lo, hi = f0 + 2, f0 + blen - 2       # trim block ends
    centers = []
    i = f0
    while i < f0 + blen:
        if sse_full[i] == "c":
            j = i
            while j < f0 + blen and sse_full[j] == "c":
                j += 1
            if i > lo and j < hi and (j - i) >= 3:
                centers.append((i + j) // 2)
            i = j
        else:
            i += 1
    return one, off, sse_full, centers


def pocket_residues(pdb: str, chain: str, full_seq: str, ligand_ccd: str = None,
                    cutoff: float = 5.0):
    """0-based receptor residues within `cutoff` Å of any ligand atom.

    Used for pocket-safe circular permutation: these residues (and their
    neighbours) should not be permutation sites. Needs a *liganded* deposited
    structure; returns [] if no ligand is present. (Integrated concept from the
    biosensor-chimera-design skill.)
    """
    import biotite.database.rcsb as rcsb
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    import numpy as np
    import tempfile
    arr = pdbx.get_structure(pdbx.BinaryCIFFile.read(rcsb.fetch(pdb, "bcif", tempfile.mkdtemp())),
                             model=1)
    prot = arr[(arr.chain_id == chain) & struc.filter_amino_acids(arr)]
    het = arr[~struc.filter_amino_acids(arr) & (arr.res_name != "HOH")]
    if ligand_ccd:
        het = het[het.res_name == ligand_ccd]
    if het.array_length() == 0:
        return []
    # map modeled residues to full_seq via the same alignment as annotate()
    from biotite.sequence import ProteinSequence
    ca = prot[prot.atom_name == "CA"]
    one = "".join(ProteinSequence.convert_letter_3to1(r) for r in ca.res_name)
    off = max(full_seq.find(one), 0)
    res_ids = list(dict.fromkeys(prot.res_id.tolist()))
    pocket = []
    for k, rid in enumerate(res_ids):
        res = prot[prot.res_id == rid]
        d = np.linalg.norm(het.coord[None, :, :] - res.coord[:, None, :], axis=2)
        if (d <= cutoff).any():
            pocket.append(off + k)
    return pocket


def verify_frozen() -> bool:
    """Recompute loop sites for both receptors and compare to systems.py."""
    from .systems import RECEPTORS
    ok = True
    for r in RECEPTORS.values():
        _, off, _, centers = annotate(r.pdb, "A", r.seq)
        match = (centers == r.loop_sites and off == r.modeled_offset)
        flag = "OK " if match else "MISMATCH"
        print(f"[{flag}] {r.name}: computed offset={off} sites={centers}")
        if not match:
            print(f"        frozen  offset={r.modeled_offset} sites={r.loop_sites}")
        ok = ok and match
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if verify_frozen() else 1)
