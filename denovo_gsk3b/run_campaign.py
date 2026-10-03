"""Run the GSK3-beta de novo graph-GA campaign.

Usage:
    python -m denovo_gsk3b.run_campaign --actives gsk3b_actives.json --out results/
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import time

from rdkit import Chem, RDLogger, SimDivFilters
from rdkit.Chem import rdFingerprintGenerator

from . import filters, scoring

RDLogger.DisableLog("rdApp.*")
_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def pick_diverse(smiles: list[str], n: int, seed: int = 42) -> list[str]:
    """MaxMin diversity pick, so the starting population spans chemotypes."""
    fps, keep = [], []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            fps.append(_MORGAN.GetFingerprint(mol))
            keep.append(s)
    if len(keep) <= n:
        return keep
    picker = SimDivFilters.MaxMinPicker()
    idx = picker.LazyBitVectorPick(fps, len(fps), n, seed=seed)
    return [keep[i] for i in idx]


def rank_weights(n: int) -> list[float]:
    """Linear rank weights (best gets the largest share)."""
    return [float(n - i) for i in range(n)]


def run(actives_path: str, out_dir: str, pop_size: int = 200,
        n_generations: int = 40, offspring: int | None = None,
        seed: int = 42, mutation_rate: float = 0.5) -> dict:
    from . import graph_ga

    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    offspring = offspring or pop_size

    actives = json.load(open(actives_path))
    ref_smiles = [v["smiles"] for v in actives.values()]
    ref_smiles = [s for s in (scoring.canonical(s) for s in ref_smiles) if s]
    print(f"reference actives: {len(ref_smiles)}")

    seeds = pick_diverse(ref_smiles, pop_size, seed=seed)
    print(f"seed population (MaxMin diverse): {len(seeds)}")

    t0 = time.time()
    scored = {r["smiles"]: r for r in scoring.score_batch(seeds)}
    population = sorted(scored.values(), key=lambda r: r["score"], reverse=True)
    population = [r for r in population if not r["gated"]][:pop_size]
    if len(population) < 2:
        raise RuntimeError("seed population collapsed under the gate")

    all_seen: dict[str, dict] = dict(scored)
    history: list[dict] = []
    ref_fps = filters.reference_fingerprints(ref_smiles)

    for gen in range(1, n_generations + 1):
        parents = [r["smiles"] for r in population]
        weights = rank_weights(len(parents))
        children: list[str] = []
        attempts = 0
        while len(children) < offspring and attempts < offspring * 8:
            attempts += 1
            pool = rng.choices(parents, weights=weights, k=2)
            child = graph_ga.reproduce(pool, rng, mutation_rate=mutation_rate)
            if child and child not in all_seen:
                children.append(child)

        new_records = scoring.score_batch(children) if children else []
        for r in new_records:
            all_seen[r["smiles"]] = r

        merged = {r["smiles"]: r for r in population}
        for r in new_records:
            if not r["gated"]:
                merged[r["smiles"]] = r
        population = sorted(merged.values(), key=lambda r: r["score"],
                            reverse=True)[:pop_size]

        scores = [r["score"] for r in population]
        novel = [r for r in population
                 if filters.max_similarity(r["smiles"], ref_fps)
                 < filters.NOVELTY_TANIMOTO_MAX]
        row = {
            "generation": gen,
            "best": max(scores),
            "mean": sum(scores) / len(scores),
            "median": sorted(scores)[len(scores) // 2],
            "n_new_scored": len(new_records),
            "n_unique_total": len(all_seen),
            "n_novel_in_pop": len(novel),
            "best_novel": max((r["score"] for r in novel), default=0.0),
            "elapsed_s": round(time.time() - t0, 1),
        }
        history.append(row)
        print(f"gen {gen:3d}  best {row['best']:.4f}  mean {row['mean']:.4f}  "
              f"novel_in_pop {row['n_novel_in_pop']:3d}  "
              f"best_novel {row['best_novel']:.4f}  "
              f"unique {row['n_unique_total']:6d}  {row['elapsed_s']:.0f}s")

    json.dump(history, open(out / "ga_history.json", "w"), indent=1)
    json.dump(list(all_seen.values()), open(out / "ga_all_scored.json", "w"))
    json.dump(seeds, open(out / "ga_seeds.json", "w"), indent=1)
    json.dump({"pop_size": pop_size, "n_generations": n_generations,
               "offspring_per_gen": offspring, "seed": seed,
               "mutation_rate": mutation_rate,
               "n_reference_actives": len(ref_smiles),
               "runtime_s": round(time.time() - t0, 1)},
              open(out / "ga_params.json", "w"), indent=1)
    print(f"\ndone in {time.time()-t0:.0f}s; {len(all_seen)} unique molecules scored")
    return {"history": history, "all_scored": list(all_seen.values()),
            "seeds": seeds, "ref_smiles": ref_smiles}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--actives", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pop-size", type=int, default=200)
    ap.add_argument("--generations", type=int, default=40)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    run(args.actives, args.out, pop_size=args.pop_size,
        n_generations=args.generations, seed=args.seed)


if __name__ == "__main__":
    main()
