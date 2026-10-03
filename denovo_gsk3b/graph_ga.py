"""A graph genetic algorithm over molecular graphs (Jensen-style).

Crossover is a single-bond cut-and-join on acyclic single bonds: both parents
are fragmented at a random non-ring single bond and one fragment from each is
re-joined at the cut points. Mutation applies one of several local graph edits
(atom append / replace / delete, bond-order change, atom insertion, ring
attachment from a medicinal-chemistry fragment library).

This reproduces the operator family of Jensen (Chem. Sci. 2019) and the
"graph GA" baseline in GuacaMol (Brown et al., J. Chem. Inf. Model. 2019),
implemented directly on RDKit rather than taken from either codebase.
"""

from __future__ import annotations

import random

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# Fragments appended by the ring-attachment mutation. Common, benign,
# synthetically accessible medicinal-chemistry rings and small groups.
FRAGMENTS = [
    "c1ccccc1", "c1ccncc1", "c1ccc(F)cc1", "c1ccc(Cl)cc1", "c1cc[nH]c1",
    "c1cn[nH]c1", "c1cnc[nH]1", "c1ccsc1", "c1ccoc1", "c1ncccn1",
    "C1CCNCC1", "C1CCOCC1", "C1CCNC1", "C1COCCN1", "C1CC1", "C1CCCC1",
    "C1CCCCC1", "c1ccc2[nH]ccc2c1", "c1ccc2ncccc2c1",
]
APPEND_GROUPS = ["C", "N", "O", "F", "Cl", "C(=O)N", "C(=O)O", "S(=O)(=O)N",
                 "C#N", "OC", "NC", "C(F)(F)F", "O", "N"]
ELEMENTS = ["C", "N", "O", "S"]

_JOIN_RXN = AllChem.ReactionFromSmarts("[*:1]-[1*].[1*]-[*:2]>>[*:1]-[*:2]")


def _sanitise(mol: Chem.Mol) -> str | None:
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    smi = Chem.MolToSmiles(mol)
    if not smi or "." in smi:
        return None
    # round-trip to be certain the string parses back
    return smi if Chem.MolFromSmiles(smi) is not None else None


def _cut(mol: Chem.Mol, rng: random.Random):
    """Fragment at a random acyclic single bond; return the two pieces."""
    bonds = [
        b.GetIdx()
        for b in mol.GetBonds()
        if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE
        and b.GetBeginAtom().GetDegree() > 1 and b.GetEndAtom().GetDegree() > 1
    ]
    if not bonds:
        return None
    frag = Chem.FragmentOnBonds(mol, [rng.choice(bonds)], addDummies=True,
                                dummyLabels=[(1, 1)])
    try:
        pieces = Chem.GetMolFrags(frag, asMols=True, sanitizeFrags=True)
    except Exception:
        return None
    return pieces if len(pieces) == 2 else None


def crossover(smi_a: str, smi_b: str, rng: random.Random, tries: int = 12):
    """Cut both parents, join a piece of A to a piece of B."""
    mol_a, mol_b = Chem.MolFromSmiles(smi_a), Chem.MolFromSmiles(smi_b)
    if mol_a is None or mol_b is None:
        return None
    for _ in range(tries):
        pa, pb = _cut(mol_a, rng), _cut(mol_b, rng)
        if pa is None or pb is None:
            continue
        try:
            products = _JOIN_RXN.RunReactants((rng.choice(pa), rng.choice(pb)))
        except Exception:
            continue
        if not products:
            continue
        child = _sanitise(rng.choice(products)[0])
        if child:
            return child
    return None


def _mut_append(mol, rng):
    """Attach a group or ring to an atom that has a free hydrogen."""
    cand = [a.GetIdx() for a in mol.GetAtoms() if a.GetTotalNumHs() > 0]
    if not cand:
        return None
    group = rng.choice(FRAGMENTS + APPEND_GROUPS)
    combo = Chem.RWMol(Chem.CombineMols(mol, Chem.MolFromSmiles(group)))
    combo.AddBond(rng.choice(cand), mol.GetNumAtoms(), Chem.BondType.SINGLE)
    return combo.GetMol()


def _mut_replace_atom(mol, rng):
    """Change one heavy atom's element."""
    cand = [a.GetIdx() for a in mol.GetAtoms() if not a.GetIsAromatic()]
    if not cand:
        return None
    rw = Chem.RWMol(mol)
    atom = rw.GetAtomWithIdx(rng.choice(cand))
    new = rng.choice([e for e in ELEMENTS if e != atom.GetSymbol()])
    atom.SetAtomicNum(Chem.GetPeriodicTable().GetAtomicNumber(new))
    atom.SetNoImplicit(False)
    atom.SetNumExplicitHs(0)
    return rw.GetMol()


def _mut_delete_atom(mol, rng):
    """Remove a terminal heavy atom."""
    cand = [a.GetIdx() for a in mol.GetAtoms()
            if a.GetDegree() == 1 and not a.IsInRing()]
    if not cand:
        return None
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(rng.choice(cand))
    return rw.GetMol()


def _mut_bond_order(mol, rng):
    """Flip an acyclic bond between single and double."""
    cand = [b.GetIdx() for b in mol.GetBonds()
            if not b.IsInRing() and not b.GetIsAromatic()
            and b.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE)]
    if not cand:
        return None
    rw = Chem.RWMol(mol)
    bond = rw.GetBondWithIdx(rng.choice(cand))
    bond.SetBondType(Chem.BondType.DOUBLE
                     if bond.GetBondType() == Chem.BondType.SINGLE
                     else Chem.BondType.SINGLE)
    return rw.GetMol()


def _mut_insert_atom(mol, rng):
    """Insert a heavy atom into an acyclic bond (A-B becomes A-X-B)."""
    cand = [b.GetIdx() for b in mol.GetBonds()
            if not b.IsInRing() and b.GetBondType() == Chem.BondType.SINGLE]
    if not cand:
        return None
    rw = Chem.RWMol(mol)
    bond = rw.GetBondWithIdx(rng.choice(cand))
    i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    rw.RemoveBond(i, j)
    new_idx = rw.AddAtom(Chem.Atom(rng.choice(ELEMENTS)))
    rw.AddBond(i, new_idx, Chem.BondType.SINGLE)
    rw.AddBond(new_idx, j, Chem.BondType.SINGLE)
    return rw.GetMol()


MUTATIONS = [
    (_mut_append, 0.35),
    (_mut_replace_atom, 0.20),
    (_mut_delete_atom, 0.15),
    (_mut_bond_order, 0.10),
    (_mut_insert_atom, 0.20),
]


def mutate(smi: str, rng: random.Random, tries: int = 12):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    ops, weights = zip(*MUTATIONS)
    for _ in range(tries):
        op = rng.choices(ops, weights=weights, k=1)[0]
        try:
            out = op(mol, rng)
        except Exception:
            continue
        child = _sanitise(out)
        if child and child != smi:
            return child
    return None


def reproduce(parents: list[str], rng: random.Random, mutation_rate: float = 0.5):
    """Produce one child: crossover two parents, then mutate with some rate."""
    a, b = rng.sample(parents, 2) if len(parents) >= 2 else (parents[0], parents[0])
    child = crossover(a, b, rng)
    if child is None:
        child = mutate(a, rng)
    elif rng.random() < mutation_rate:
        child = mutate(child, rng) or child
    return child
