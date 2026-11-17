"""Stage 1 CLI: fetch -> classify -> dedup -> split -> write manifest + stats.

Dry/demo run on the seed list (real RCSB data, ~1 network call):
    python -m abproto.data.build --out stage1

Full discovery run (RCSB search for Ab-Ag complexes, may be large/slow):
    python -m abproto.data.build --search --max-rows 500 --out stage1
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .classify import to_complexes
from .oas import default_cdr_length_prior, save_prior
from .rcsb import fetch_entries, search_antibody_antigen
from .seeds import SEED_COMPLEXES
from .splits import make_splits


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="Stage 1 antibody-antigen dataset builder")
    p.add_argument("--search", action="store_true",
                   help="discover ids via RCSB search instead of the seed list")
    p.add_argument("--max-rows", type=int, default=300)
    p.add_argument("--date-cutoff", default="2021-09-30")
    p.add_argument("--antigen-homology", type=float, default=0.4)
    p.add_argument("--out", type=Path, default=Path("stage1"))
    args = p.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)

    if args.search:
        ids = search_antibody_antigen(max_rows=args.max_rows)
        print(f"RCSB search returned {len(ids)} candidate ids")
        if not ids:
            print("search returned nothing (API/network); falling back to seed list")
            ids = SEED_COMPLEXES
    else:
        ids = SEED_COMPLEXES

    print(f"fetching {len(ids)} entries from RCSB GraphQL ...")
    entries = fetch_entries(ids)
    complexes, skipped = to_complexes(entries)
    print(f"  valid Ab-Ag complexes: {len(complexes)}   skipped (not Ab-Ag): {len(skipped)}")
    if skipped:
        print(f"  skipped: {', '.join(skipped)}")

    split = make_splits(complexes, date_cutoff=args.date_cutoff,
                        antigen_homology=args.antigen_homology)

    print("\n=== split stats ===")
    for k, v in split.stats.items():
        print(f"  {k}: {v}")

    # Write manifest + splits + CDR-length prior.
    manifest = {
        "stats": split.stats,
        "train": [c.summary() for c in split.train],
        "test": [c.summary() for c in split.test],
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2))
    _write_fasta(args.out / "train_antibodies.fasta", split.train)
    _write_fasta(args.out / "test_antibodies.fasta", split.test)
    save_prior(default_cdr_length_prior(), str(args.out / "cdr_length_prior.json"))

    print(f"\nwrote {args.out}/manifest.json, train/test FASTA, cdr_length_prior.json")
    print("(CDR-length prior is the approximate human default; replace by running "
          "abproto.data.oas.extract_from_unit on an OAS paired unit.)")
    return 0


def _write_fasta(path: Path, complexes) -> None:
    lines = []
    for c in complexes:
        lines.append(f">{c.pdb_id}|ab heavy={len(c.heavy)} light={len(c.light)} gpcr={c.is_gpcr}")
        lines.append(c.antibody_seq)
        lines.append(f">{c.pdb_id}|antigen {c.antigen[0].description if c.antigen else ''}")
        lines.append(c.antigen_seq)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    sys.exit(main())
