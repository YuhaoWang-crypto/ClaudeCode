"""Retrosynthesis on the selected designs via AiZynthFinder (USPTO templates).

Runs in its own virtualenv because AiZynthFinder's dependency pins conflict with
the scikit-learn version required to unpickle the TDC oracles.

Routes are template-based proposals from a model trained on USPTO reaction data.
A "solved" route means every leaf matched the ZINC stock under the search budget
-- it is NOT a validated synthesis, and says nothing about yield or selectivity.
"""

from __future__ import annotations

import argparse
import json
import time


def build_finder(configfile: str, time_limit: int):
    from aizynthfinder.aizynthfinder import AiZynthFinder

    finder = AiZynthFinder(configfile=configfile)
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto")
    try:
        finder.filter_policy.select("uspto")
    except Exception as exc:          # filter policy is optional
        print(f"  [warn] filter policy unavailable: {exc}")
    # time limit lives in different places across 3.x / 4.x
    for obj, attr in ((getattr(finder.config, "search", None), "time_limit"),
                      (finder.config, "time_limit")):
        if obj is not None and hasattr(obj, attr):
            setattr(obj, attr, time_limit)
            break
    return finder


def run(designs_path: str, configfile: str, out_path: str,
        time_limit: int = 120) -> list[dict]:
    designs = json.load(open(designs_path))
    finder = build_finder(configfile, time_limit)
    results = []

    for d in designs:
        did, smi = d["design_id"], d["smiles"]
        t0 = time.time()
        rec = {"design_id": did, "smiles": smi}
        try:
            finder.target_smiles = smi
            finder.tree_search()
            finder.build_routes()
            stats = finder.extract_statistics()
            routes = finder.routes
            scores = []
            for i in range(len(routes)):
                try:
                    scores.append(float(routes[i]["all_score"]["state score"]))
                except Exception:
                    try:
                        scores.append(float(routes[i]["score"]))
                    except Exception:
                        pass
            rec.update({
                "solved": bool(stats.get("is_solved")),
                "n_routes": len(routes),
                "best_state_score": max(scores) if scores else None,
                "n_steps_best": stats.get("number_of_steps"),
                "n_precursors": stats.get("number_of_precursors"),
                "precursors_in_stock": stats.get("number_of_precursors_in_stock"),
                "top_route_smiles": stats.get("precursors_in_stock"),
                "search_time_s": round(time.time() - t0, 1),
                "status": "ok",
            })
            print(f"  {did}: solved={rec['solved']} routes={rec['n_routes']} "
                  f"steps={rec['n_steps_best']} "
                  f"score={rec['best_state_score']} ({rec['search_time_s']}s)")
        except Exception as exc:
            rec.update({"status": f"error: {type(exc).__name__}: {exc}",
                        "solved": None,
                        "search_time_s": round(time.time() - t0, 1)})
            print(f"  {did}: FAILED {exc}")
        results.append(rec)

    json.dump(results, open(out_path, "w"), indent=1)
    n_ok = sum(1 for r in results if r.get("solved"))
    print(f"\nsolved {n_ok}/{len(results)} designs")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--designs", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--time-limit", type=int, default=120)
    args = ap.parse_args()
    run(args.designs, args.config, args.out, time_limit=args.time_limit)


if __name__ == "__main__":
    main()
