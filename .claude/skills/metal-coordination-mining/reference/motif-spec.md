# Motif spec reference

A motif spec is one JSON file. It is the whole interface: a new mechanism is a
new spec, never new code. Specs live in `assets/motifs/`.

## Full schema

```jsonc
{
  "name": "my_motif",                    // used in output filenames
  "description": "Mechanism in prose — state the anchor AND the discriminator.",

  "anchor": {
    "residues": ["HIS", "HIS"],          // exactly two; may differ (["HIS","CYS"])
    "anchor_atoms": {                    // side-chain atoms per residue type;
      "HIS": ["NE2"]                     // the closest valid pair defines the probe
    },
    "anchor_distance_max": 4.0,          // Å, side-chain preorganization
    "backbone_hbond_max": 4.0,           // Å, BOTH O(i)-N(j) and N(i)-O(j); null disables
    "same_chain_only": false,            // true to forbid inter-chain sites
    "min_sequence_separation": 0         // e.g. 3 to reject i,i+1/i+2 neighbours
  },

  "atom_set": "proximity",               // "proximity" | "coordinating" | inline {RES:[atoms]}
  "ligand_shell_cutoff": 4.5,            // Å, first-shell membership for shell_residues/typing

  "rules": {
    "require_absent":  [ {"residues": ["ASP","GLU"], "min_distance": 5.5} ],
    "require_present": [ {"residues": ["ALA","GLY"], "max_distance": 7.0} ],
    "sequence_context": {"from": "anchorA", "offset": 2, "allowed": ["ALA","GLY"]},
    "min_protein_length": 200
  },

  "contrast_rules": {                    // optional: the sibling class, counted alongside
    "name": "facial_triad_2His_1AspGlu",
    "require_present": [ {"residues": ["ASP","GLU"], "max_distance": 5.0} ],
    "min_protein_length": 0
  }
}
```

### Rule semantics

- `require_absent` — **every** listed residue type must be farther than
  `min_distance` from the probe. Multiple blocks allow different cutoffs for
  different residue groups (carboxylates at 5.5 Å, other ligands at 4.0 Å).
- `require_present` — **at least one** listed type within `max_distance` (an OR
  within a block; blocks AND together).
- `sequence_context` — residue at a fixed sequence offset from an anchor
  (`anchorA` = the first anchor residue, `anchorB` = the second) must be in
  `allowed`. Cheap orthogonal evidence; drop it if the family's numbering is
  irregular. `"NA"` (offset past the chain end) never matches.
- `min_protein_length` — residue count of the model; a fragment filter.
- A residue type absent from a structure is treated as **distance = ∞**, so
  `require_absent` passes and `require_present` fails. This is the intended
  reading of "nothing of that type is nearby".

### Atom sets

`coordinating` for what IS there (typing, first-shell membership);
`proximity` for what must NOT be there (exclusion rules — it includes side-chain
carbons, so a rotatable carboxylate cannot hide). See `methodology.md` §2.

## Sweep addressing

`benchmark --sweep "<path>:<lo>:<hi>:<step>"` overrides one numeric field, addressed
by dotted path into `rules`, with list indices as integers:

```
require_absent.0.min_distance:4.0:7.0:0.25     # carboxylate exclusion radius
require_absent.1.min_distance:3.0:5.0:0.25     # other-ligand exclusion radius
require_present.0.max_distance:4.0:9.0:0.5     # Ala/Gly presence radius
min_protein_length:100:400:50
```

Anchor geometry is not sweepable this way (it changes which sites exist, not how
they are classified) — edit the spec and re-run `mine` to explore it.

## Worked example 1 — reciprocal motif (find the siblings)

Invert the halogenase spec to enumerate canonical facial-triad hydroxylases:

```jsonc
"rules": {
  "require_present": [ {"residues": ["ASP","GLU"], "max_distance": 5.0} ],
  "min_protein_length": 200
}
```

Useful as an explicit negative set, and to ask which subfamilies contain both
classes (the paper's observation that halogenases and non-halogenases interleave
within the same SSN clusters).

## Worked example 2 — a different metal site (2His-1Cys)

Anchor unchanged, discriminator becomes a required thiolate:

```jsonc
"atom_set": "coordinating",
"rules": {
  "require_present": [ {"residues": ["CYS"], "max_distance": 3.0} ],
  "require_absent":  [ {"residues": ["ASP","GLU"], "min_distance": 5.0} ],
  "min_protein_length": 150
}
```

Discover which environments are worth a spec by running `mcmine.py type` first
and reading the coordination-type census.

## Worked example 3 — a non-metal anchor (His/Cys catalytic dyad)

```jsonc
"anchor": {
  "residues": ["HIS", "CYS"],
  "anchor_atoms": {"HIS": ["NE2", "ND1"], "CYS": ["SG"]},
  "anchor_distance_max": 4.5,
  "backbone_hbond_max": null,          // dyad partners are not strand-paired
  "min_sequence_separation": 3
}
```

With `backbone_hbond_max: null` the anchor is purely side-chain geometry, so
expect more accidental pairs — compensate with a stricter distance, a
`require_present` third member (e.g. Asp/Glu of a triad), and a benchmark set.

## Checklist before running at scale

1. Anchor present in **100%** of known members? (`benchmark` reports positives
   with no anchor site separately.)
2. Discriminator mechanistically **required**, not merely correlated?
3. Negative class named, and members of it in the benchmark set?
4. Recall 1.0 ✅ and specificity reported on the benchmark?
5. Sweep run — plateau or knife edge?
6. Post-mining annotation filter planned for non-enzyme members of the fold?
