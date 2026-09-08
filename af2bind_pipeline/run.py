"""CLI and orchestration for an AF2BIND run.

    # GPU pass on Modal, then scoring + analysis locally
    python -m af2bind_pipeline.run --target 6w70 --chain A --out results/6w70

    # re-score cached features with the 10-fold ensemble, no GPU needed
    python -m af2bind_pipeline.run --features results/6w70/features.npz \
        --seeds 0,1,2,3,4,5,6,7,8,9 --out results/6w70_ens
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from . import core, pockets, structure, validate, weights


def score_features(
    features: np.ndarray,
    mask_sidechains: bool = True,
    seeds=(0,),
    cache_dir=None,
) -> dict:
    heads = weights.load_heads(mask_sidechains, seeds, cache_dir)
    return core.predict(features, heads)


def _write_csv(path: Path, header: list[str], rows: list[list]) -> None:
    lines = [",".join(header)]
    for r in rows:
        lines.append(",".join("" if v is None else str(v) for v in r))
    path.write_text("\n".join(lines) + "\n")


def run_from_features(
    result: dict,
    pdb_path: str | Path,
    out_dir: str | Path,
    seeds=(0,),
    top_n: int = 15,
    pocket_top_n: int = 25,
    ligand_cutoff: float = 5.0,
    cache_dir=None,
) -> dict:
    """Score a GPU result, write every output, and self-validate if possible.

    `result` is what `modal_app.pair_features` returns (or a reload of it).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = Path(pdb_path)

    meta = dict(result.get("meta", {}))
    features = np.asarray(result["features"], dtype=np.float32)
    chains, resis, resns = result["chain"], result["resi"], result["resn"]
    keys = [(c, int(r)) for c, r in zip(chains, resis)]

    pred = score_features(
        features, meta.get("mask_sidechains", True), seeds, cache_dir
    )
    p_bind = pred["p_bind"]
    score_map = {k: float(p) for k, p in zip(keys, p_bind)}

    np.savez_compressed(
        out_dir / "features.npz",
        features=features,
        chain=np.array(chains),
        resi=np.array(resis),
        resn=np.array(resns),
        plddt=np.asarray(result.get("plddt", [])),
        meta=np.array(json.dumps(meta)),
    )

    order = np.argsort(-p_bind)
    _write_csv(
        out_dir / "results.csv",
        ["rank", "chain", "resi", "resn", "p_bind", "logit"],
        [
            [i + 1, chains[j], resis[j], resns[j],
             round(float(p_bind[j]), 5), round(float(pred["logit"][j]), 4)]
            for i, j in enumerate(order)
        ],
    )

    aa_matrix, aa_labels = core.blosum_reorder(pred["p_bind_aa"])
    _write_csv(
        out_dir / "bait_activations.csv",
        ["chain", "resi", "resn", *aa_labels],
        [
            [chains[j], resis[j], resns[j],
             *[round(float(v), 5) for v in aa_matrix[j]]]
            for j in order[: max(top_n, 25)]
        ],
    )

    structure.write_bfactor_pdb(
        pdb_path, out_dir / "p_bind.pdb", score_map, scale=100.0
    )

    residues = structure.parse_pdb(pdb_path)
    target_res = structure.protein_residues(residues)
    found = pockets.cluster_pockets(target_res, score_map, top_n=pocket_top_n)

    top_keys = [keys[j] for j in order[:top_n]]
    by_chain: dict[str, list[int]] = {}
    for c, r in top_keys:
        by_chain.setdefault(c, []).append(r)
    pymol_top = "select af2bind_top, " + " or ".join(
        f"(chain {c} and resi {'+'.join(str(r) for r in sorted(v))})"
        for c, v in sorted(by_chain.items())
    )

    report = {
        "target": str(pdb_path),
        "meta": meta,
        "seeds": list(seeds),
        "n_residues": len(keys),
        "p_bind": {
            "max": round(float(p_bind.max()), 4),
            "mean": round(float(p_bind.mean()), 4),
            "n_above_0.5": int((p_bind >= 0.5).sum()),
            "n_above_0.9": int((p_bind >= 0.9).sum()),
        },
        "top_residues": [
            {"rank": i + 1, "chain": chains[j], "resi": resis[j],
             "resn": resns[j], "p_bind": round(float(p_bind[j]), 4)}
            for i, j in enumerate(order[:top_n])
        ],
        "pymol_top": pymol_top,
        "pockets": [p.as_dict() for p in found],
    }

    lig = structure.ligands(residues)
    target_chains = {c for c, _ in keys}
    if lig and len(target_chains) == 1:
        lig = structure.assign_ligands_to_chain(
            residues, lig, next(iter(target_chains))
        )
    if lig:
        positives = structure.contact_residues(
            target_res, lig, cutoff=ligand_cutoff
        )
        positives &= set(keys)
        report["ligand_validation"] = {
            "ligands": sorted({r.resn for r in lig}),
            "contact_cutoff_A": ligand_cutoff,
            **validate.evaluate(keys, p_bind, positives),
        }

    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        prog="af2bind", description="AF2BIND small-molecule binding-site prediction"
    )
    ap.add_argument("--target", help="PDB file path, 4-char PDB ID, or UniProt accession")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--out", default="af2bind_out")
    ap.add_argument(
        "--features",
        help="re-score a cached features.npz instead of running the GPU pass",
    )
    ap.add_argument("--seeds", default="0", help="comma-separated fold seeds, 0..9")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--pocket-top-n", type=int, default=25)
    ap.add_argument("--no-mask-sidechains", action="store_true",
                    help="use the head trained on unmasked target side chains")
    ap.add_argument("--af2-params", default="2021-07-14")
    ap.add_argument("--min-plddt", type=float, default=None,
                    help="trim residues below this B-factor/pLDDT before running")
    args = ap.parse_args(argv)

    seeds = tuple(int(s) for s in args.seeds.split(",") if s.strip())
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.features:
        data = np.load(args.features, allow_pickle=False)
        result = {
            "features": data["features"],
            "chain": [str(c) for c in data["chain"]],
            "resi": [int(r) for r in data["resi"]],
            "resn": [str(r) for r in data["resn"]],
            "plddt": data["plddt"],
            "meta": json.loads(str(data["meta"])),
        }
        pdb_path = args.target or result["meta"].get("pdb_path")
        if not pdb_path:
            ap.error("--features also needs --target pointing at the same PDB")
    else:
        if not args.target:
            ap.error("one of --target or --features is required")
        pdb_path = structure.fetch_structure(args.target, out_dir)
        if args.min_plddt is not None:
            trimmed = out_dir / "trimmed.pdb"
            kept = structure.strip_low_plddt(pdb_path, trimmed, args.min_plddt)
            print(f"[af2bind] kept {kept} residues above pLDDT {args.min_plddt}")
            pdb_path = trimmed
        from .modal_app import app as modal_app_handle
        from .modal_app import pair_features

        with modal_app_handle.run():
            result = pair_features.remote(
                pdb_text=Path(pdb_path).read_text(),
                chain=args.chain,
                mask_sidechains=not args.no_mask_sidechains,
                af2_params=args.af2_params,
            )
        result["meta"]["pdb_path"] = str(pdb_path)

    report = run_from_features(
        result, pdb_path, out_dir,
        seeds=seeds, top_n=args.top_n, pocket_top_n=args.pocket_top_n,
    )

    print(f"[af2bind] {report['n_residues']} residues scored, "
          f"max p(bind)={report['p_bind']['max']}")
    print(f"[af2bind] {len(report['pockets'])} candidate pocket(s)")
    for p in report["pockets"]:
        top = ", ".join(
            f"{r['resn']}{r['resi']}" for r in p["residues"][:6]
        )
        print(f"    pocket {p['rank']}: {p['size']} residues, "
              f"score_sum={p['score_sum']}, {top}")
    if "ligand_validation" in report:
        v = report["ligand_validation"]
        print(f"[af2bind] ligand check ({'+'.join(v['ligands'])}): "
              f"ROC-AUC={v['roc_auc']}, AP={v['average_precision']}, "
              f"precision@15={v['top_n']['15']['precision_at_15']}")
    print(f"[af2bind] outputs in {out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
