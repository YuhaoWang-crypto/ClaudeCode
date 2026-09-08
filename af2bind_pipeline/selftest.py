"""Offline self-test: everything except the GPU pass.

    python -m af2bind_pipeline.selftest

Checks the head loads, the scoring math matches the upstream reference formula
bit-for-bit, the feature layout is the one the head expects, the PDB parser and
ligand-contact logic behave, and the ranking metrics are correct on constructed
cases. Needs numpy and one download of the weights archive; no GPU, no jax.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from . import core, pockets, structure, validate, weights

MINI_PDB = """\
ATOM      1  N   ALA A   1      -1.000   0.000   0.000  1.00 90.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C
ATOM      3  C   ALA A   1       1.000   0.000   0.000  1.00 90.00           C
ATOM      4  O   ALA A   1       1.500   1.000   0.000  1.00 90.00           O
ATOM      5  CB  ALA A   1       0.000   1.500   0.000  1.00 90.00           C
ATOM      6  N   GLY A   2       3.000   0.000   0.000  1.00 40.00           N
ATOM      7  CA  GLY A   2       4.000   0.000   0.000  1.00 40.00           C
ATOM      8  C   GLY A   2       5.000   0.000   0.000  1.00 40.00           C
ATOM      9  N   LEU A   3      20.000   0.000   0.000  1.00 95.00           N
ATOM     10  CA  LEU A   3      21.000   0.000   0.000  1.00 95.00           C
ATOM     11  C   LEU A   3      22.000   0.000   0.000  1.00 95.00           C
ATOM     12  CB  LEU A   3      21.000   1.500   0.000  1.00 95.00           C
HETATM   13  C1  LIG A 101       0.500   2.500   0.000  1.00 20.00           C
HETATM   14  C2  LIG A 101       1.500   2.500   0.000  1.00 20.00           C
HETATM   15  C3  LIG A 101       2.500   2.500   0.000  1.00 20.00           C
HETATM   16  C4  LIG A 101       0.500   3.500   0.000  1.00 20.00           C
HETATM   17  C5  LIG A 101       1.500   3.500   0.000  1.00 20.00           C
HETATM   18  C6  LIG A 101       2.500   3.500   0.000  1.00 20.00           C
HETATM   19  O   HOH A 201      50.000  50.000  50.000  1.00 30.00           O
END
"""

_checks: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    _checks.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))


def test_head_and_math() -> None:
    print("head + scoring math")
    head = weights.load_head(mask_sidechains=True, seed=0)
    check(
        "head shapes",
        head["w"].shape == (core.N_FEATURES, 1)
        and head["mean"].shape == (core.N_FEATURES,)
        and head["std"].shape == (core.N_FEATURES,)
        and head["b"].shape == (1,),
        f"w{head['w'].shape}",
    )
    check(
        "nosc and plain heads differ",
        not np.allclose(head["w"], weights.load_head(False, 0)["w"]),
    )
    check(
        "seeds are distinct folds",
        not np.allclose(head["w"], weights.load_head(True, 1)["w"]),
    )

    rng = np.random.default_rng(0)
    L = 23
    x = rng.normal(size=(L, core.N_FEATURES)).astype(np.float32)
    got = core.predict(x, [head])

    # Reference formula transcribed from the upstream notebook's af2bind().
    z = (x - head["mean"]) / head["std"]
    z = (z * head["w"][:, 0]) + (head["b"] / x.shape[-1])
    ref_aa = z.reshape(z.shape[0], 2, 20, -1).sum((1, 3))
    ref = 1.0 / (1.0 + np.exp(-ref_aa.sum(-1).astype(np.float64)))
    check(
        "p_bind matches upstream formula",
        np.allclose(got["p_bind"], ref, atol=1e-12),
        f"max diff {np.abs(got['p_bind'] - ref).max():.2e}",
    )
    check("per-bait logits match", np.allclose(got["p_bind_aa"], ref_aa, atol=1e-5))
    check(
        "logits sum to the total",
        np.allclose(got["p_bind_aa"].sum(-1), got["logit"], atol=1e-5),
    )
    check(
        "sigmoid is stable at extremes",
        np.allclose(core.sigmoid(np.array([-1e4, 0.0, 1e4])), [0.0, 0.5, 1.0]),
    )
    ens = core.predict(x, weights.load_heads(True, range(10)))
    check(
        "ensemble stays in range and differs from seed 0",
        (ens["p_bind"] >= 0).all()
        and (ens["p_bind"] <= 1).all()
        and not np.allclose(ens["p_bind"], got["p_bind"]),
    )


def test_features() -> None:
    print("feature extraction")
    rng = np.random.default_rng(1)
    L = 15
    pair = rng.normal(size=(L + 20, L + 20, 128)).astype(np.float32)
    f = core.features_from_pair(pair)
    check("shape is (L, 5120)", f.shape == (L, core.N_FEATURES), str(f.shape))
    check(
        "first half is target->bait",
        np.allclose(f[:, : 20 * 128].reshape(L, 20, 128), pair[:-20, -20:]),
    )
    check(
        "second half is bait->target, transposed",
        np.allclose(
            f[:, 20 * 128 :].reshape(L, 20, 128), pair[-20:, :-20].swapaxes(0, 1)
        ),
    )
    bad = False
    try:
        core.features_from_pair(rng.normal(size=(10, 10, 64)))
    except ValueError:
        bad = True
    check("rejects a wrong-width pair tensor", bad)


def test_structure() -> None:
    print("structure parsing + ligand contacts")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "mini.pdb"
        p.write_text(MINI_PDB)
        res = structure.parse_pdb(p)
        prot = structure.protein_residues(res)
        check("3 protein residues", len(prot) == 3, str(len(prot)))
        check("sequence read", "".join(r.one_letter for r in prot) == "AGL")
        lig = structure.ligands(res)
        check(
            "water excluded, 6-atom ligand kept",
            len(lig) == 1 and lig[0].resn == "LIG",
            str([r.resn for r in lig]),
        )
        contacts = structure.contact_residues(prot, lig, cutoff=5.0)
        check(
            "near residues contact, far one does not",
            contacts == {("A", 1), ("A", 2)},
            str(sorted(contacts)),
        )
        check(
            "additives excluded by default",
            structure.ligands(structure.parse_pdb(p), exclude={"LIG"}) == [],
        )

        out = Path(d) / "trim.pdb"
        kept = structure.strip_low_plddt(p, out, min_plddt=70.0)
        check("pLDDT trim drops the low-confidence residue", kept == 2, str(kept))

        bf = Path(d) / "bf.pdb"
        structure.write_bfactor_pdb(p, bf, {("A", 1): 0.75}, scale=100.0)
        line = [
            l for l in bf.read_text().splitlines()
            if l.startswith("ATOM") and l[22:26].strip() == "1"
        ][0]
        check("b-factor column rewritten", abs(float(line[60:66]) - 75.0) < 1e-6,
              line[60:66])


def _atom(serial, name, resn, chain, resi, xyz, hetero=False, bfac=90.0):
    rec = "HETATM" if hetero else "ATOM  "
    return (
        f"{rec}{serial:5d} {name:<4s} {resn:>3s} {chain}{resi:4d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1.00:6.2f}{bfac:6.2f}"
    )


def _two_chain_pdb() -> str:
    """Chains A and B, plus three HETATM groups probing the selection rules.

    INT sits midway between the chains and contacts three residues of each, the
    way a ligand bound at a homodimer interface does. BRS leans on chain B and
    only brushes one residue of A, the way a neighbouring copy's ligand does
    through crystal packing. NAG is a glycan hanging off chain A.
    """
    lines, serial = [], 1
    for chain, x0 in (("A", 0.0), ("B", 8.0)):
        for i, y in enumerate((0.0, 4.0, 8.0), start=1):
            for name in ("N", "CA", "C"):
                dx = {"N": -0.5, "CA": 0.0, "C": 0.5}[name]
                lines.append(
                    _atom(serial, name, "ALA", chain, i, (x0 + dx, y, 0.0))
                )
                serial += 1
    for name, xyz in zip(
        ("C1", "C2", "C3", "C4", "C5", "C6"),
        ((4.0, 0.0, 0.0), (4.0, 2.0, 0.0), (4.0, 4.0, 0.0),
         (4.0, 6.0, 0.0), (4.0, 8.0, 0.0), (4.0, 1.0, 0.0)),
    ):
        lines.append(_atom(serial, name, "INT", "A", 301, xyz, hetero=True))
        serial += 1
    for name, xyz in zip(
        ("C1", "C2", "C3", "C4", "C5", "C6"),
        ((4.5, 0.0, 0.0), (8.0, 2.0, 0.0), (8.0, 6.0, 0.0),
         (9.0, 0.0, 0.0), (9.0, 4.0, 0.0), (9.0, 8.0, 0.0)),
    ):
        lines.append(_atom(serial, name, "BRS", "B", 302, xyz, hetero=True))
        serial += 1
    for name, xyz in zip(
        ("C1", "C2", "C3", "C4", "C5", "C6", "C7"),
        ((0.0, -3.0, 0.0), (0.5, -3.5, 0.0), (1.0, -4.0, 0.0),
         (1.5, -4.5, 0.0), (2.0, -5.0, 0.0), (2.5, -5.5, 0.0),
         (3.0, -6.0, 0.0)),
    ):
        lines.append(_atom(serial, name, "NAG", "A", 401, xyz, hetero=True))
        serial += 1
    return "\n".join(lines) + "\nEND\n"


def test_ligand_selection() -> None:
    print("ligand selection + ground truth")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "two.pdb"
        p.write_text(_two_chain_pdb())
        res = structure.parse_pdb(p)
        prot_a = structure.protein_residues(res, "A")
        prot_b = structure.protein_residues(res, "B")
        check("both chains parsed", len(prot_a) == 3 and len(prot_b) == 3,
              f"A={len(prot_a)} B={len(prot_b)}")

        auto = {r.resn for r in structure.ligands(res)}
        check("glycan excluded by default", "NAG" not in auto, str(sorted(auto)))
        check("real ligands kept", {"INT", "BRS"} <= auto, str(sorted(auto)))

        forced = {r.resn for r in structure.ligands(res, only={"NAG"})}
        check("explicit `only` overrides the exclusion list", forced == {"NAG"},
              str(sorted(forced)))

        lig = structure.ligands(res)
        keep_a = {r.resn for r in structure.ligands_near_chain(res, lig, "A")}
        keep_b = {r.resn for r in structure.ligands_near_chain(res, lig, "B")}
        check("interface ligand counts for both chains",
              "INT" in keep_a and "INT" in keep_b, f"A={keep_a} B={keep_b}")
        check("packing-neighbour ligand dropped for the far chain",
              "BRS" not in keep_a and "BRS" in keep_b, f"A={keep_a} B={keep_b}")

        pos = structure.contact_residues(prot_a, [r for r in lig if r.resn == "INT"], 5.0)
        check("interface ligand yields chain-A ground truth", len(pos) == 3,
              str(sorted(pos)))

        check("a few known additive classes are covered",
              {"CLR", "OLC", "BOG", "NAG", "BU1", "HOH"}
              <= structure.CRYSTALLIZATION_ADDITIVES)
        check("genuine cofactors are NOT excluded",
              not ({"HEM", "FAD", "NAP", "ATP"} & structure.CRYSTALLIZATION_ADDITIVES))


def test_pockets() -> None:
    print("pocket clustering")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "mini.pdb"
        p.write_text(MINI_PDB)
        prot = structure.protein_residues(structure.parse_pdb(p))
    scores = {("A", 1): 0.9, ("A", 2): 0.8, ("A", 3): 0.7}
    one = pockets.cluster_pockets(prot, scores, top_n=3, link_cutoff=10.0, min_size=2)
    check(
        "residue 3 is 20 A away, so two clusters, one too small",
        len(one) == 1 and one[0].size == 2,
        f"{[c.size for c in one]}",
    )
    two = pockets.cluster_pockets(prot, scores, top_n=3, link_cutoff=30.0, min_size=2)
    check(
        "a generous cutoff merges them",
        len(two) == 1 and two[0].size == 3,
        f"{[c.size for c in two]}",
    )
    check(
        "pymol selection is well formed",
        two[0].pymol_selection().startswith("select pocket_1, (chain A and resi "),
        two[0].pymol_selection(),
    )


def test_metrics() -> None:
    print("ranking metrics")
    labels = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    perfect = np.array([1.0, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    check("perfect ranking gives AUC 1", validate.roc_auc(labels, perfect) == 1.0)
    check(
        "reversed ranking gives AUC 0",
        validate.roc_auc(labels, -perfect) == 0.0,
    )
    check(
        "all-tied scores give AUC 0.5",
        validate.roc_auc(labels, np.ones(10)) == 0.5,
    )
    check("perfect AP is 1", validate.average_precision(labels, perfect) == 1.0)
    # AP for [hit, miss, hit] = (1/1 + 2/3) / 2
    lab2 = np.array([1, 0, 1])
    sc2 = np.array([0.9, 0.8, 0.7])
    check(
        "AP matches hand calculation",
        abs(validate.average_precision(lab2, sc2) - (1.0 + 2 / 3) / 2) < 1e-12,
    )
    keys = [("A", i) for i in range(10)]
    ev = validate.evaluate(keys, perfect, {("A", 0), ("A", 1)}, top_ns=(2,))
    check(
        "evaluate reports enrichment over the base rate",
        ev["top_n"]["2"]["enrichment_over_random"] == 5.0,
        str(ev["top_n"]["2"]),
    )


def main() -> int:
    for fn in (test_head_and_math, test_features, test_structure,
               test_ligand_selection, test_pockets, test_metrics):
        fn()
    failed = [n for n, ok, _ in _checks if not ok]
    print(f"\n{len(_checks) - len(failed)}/{len(_checks)} checks passed")
    if failed:
        print("failed: " + ", ".join(failed))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
