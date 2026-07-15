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

    off = full_seq.find(one)
    if off < 0:
        raise ValueError(f"{pdb}: modeled sequence is not a substring of the provided full sequence")

    sse_full = ["c"] * len(full_seq)
    for i, s in enumerate(sse):
        sse_full[off + i] = s
    sse_full = "".join(sse_full)

    # interior loop centers (flanked by SSE both sides, length>=3, within modeled span)
    centers = []
    i = 0
    while i < len(sse_full):
        if sse_full[i] == "c":
            j = i
            while j < len(sse_full) and sse_full[j] == "c":
                j += 1
            if i > off + 2 and j < off + len(one) - 2 and (j - i) >= 3:
                centers.append((i + j) // 2)
            i = j
        else:
            i += 1
    return one, off, sse_full, centers


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
