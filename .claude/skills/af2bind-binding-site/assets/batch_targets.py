"""Template: score many targets in one Modal session.

Copy, edit TARGETS, run from the repository root:

    python .claude/skills/af2bind-binding-site/assets/batch_targets.py

Everything stays inside a single `app.run()` block so the container, the loaded
AlphaFold parameters and the jax compilation cache are reused across targets.
Opening a session per target instead is the main way to waste money here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from af2bind_pipeline import structure  # noqa: E402
from af2bind_pipeline.modal_app import app, pair_features  # noqa: E402
from af2bind_pipeline.run import run_from_features  # noqa: E402

# (target, chain) — target is a PDB file path, a 4-char PDB ID, or a UniProt id.
TARGETS = [
    ("6w70", "A"),
    ("1stp", "A"),
]

OUT_ROOT = Path("results/batch")
SEEDS = (0,)
MIN_PLDDT = None  # set e.g. 70.0 when the inputs are AlphaFold models


def main() -> int:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    summary = []

    with app.run():
        for target, chain in TARGETS:
            name = Path(target).stem if Path(target).exists() else target
            out_dir = OUT_ROOT / f"{name}_{chain}"
            out_dir.mkdir(parents=True, exist_ok=True)
            try:
                pdb_path = structure.fetch_structure(target, out_dir)
                if MIN_PLDDT is not None:
                    trimmed = out_dir / "trimmed.pdb"
                    structure.strip_low_plddt(pdb_path, trimmed, MIN_PLDDT)
                    pdb_path = trimmed

                result = pair_features.remote(
                    pdb_text=Path(pdb_path).read_text(),
                    chain=chain,
                    mask_sidechains=True,
                )
                report = run_from_features(
                    result, pdb_path, out_dir, seeds=SEEDS
                )
            except Exception as exc:  # noqa: BLE001
                print(f"[af2bind] {name}/{chain} FAILED: {exc}")
                summary.append({"target": name, "chain": chain, "error": str(exc)})
                continue

            row = {
                "target": name,
                "chain": chain,
                "n_residues": report["n_residues"],
                "max_p_bind": report["p_bind"]["max"],
                "n_pockets": len(report["pockets"]),
                "top_pocket": (
                    [f"{r['resn']}{r['resi']}" for r in report["pockets"][0]["residues"][:8]]
                    if report["pockets"] else []
                ),
            }
            if "ligand_validation" in report:
                row["roc_auc_vs_ligand"] = report["ligand_validation"]["roc_auc"]
            summary.append(row)
            print(f"[af2bind] {name}/{chain}: {json.dumps(row)}")

    (OUT_ROOT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[af2bind] wrote {OUT_ROOT}/summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
