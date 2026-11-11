"""Post-training analysis.

    python -m alphago_zero.analyze --run alphago_zero/runs/7x7

1. Round-robin tournament between saved generations -> Elo by maximum
   likelihood (the chained gating Elo is biased upward by the selection).
2. Final strength vs fixed opponents (random player, untrained network).
3. Training-curve figure and a few SGF game records of the final network.
"""
import argparse
import glob
import itertools
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import torch

from .go import GoState
from .mcts import Node, NetEvaluator, search
from .model import load
from .selfplay import arena_worker, init_worker


def spec_of(path, sims):
    ck = torch.load(path, map_location="cpu", weights_only=False)
    return dict(kind="mcts", net_cfg=ck["net_cfg"], state_dict=ck["state_dict"],
                sims=sims), ck


def fit_elo(n_players, results, anchor=0):
    """results: list of (i, j, score_i, games). Bradley-Terry MLE, +0.5 pseudo-wins."""
    W = np.zeros((n_players, n_players))
    G = np.zeros((n_players, n_players))
    for i, j, s, g in results:
        W[i, j] += s + 0.5
        W[j, i] += g - s + 0.5
        G[i, j] += g + 1
        G[j, i] += g + 1
    r = np.zeros(n_players)
    k = np.log(10) / 400
    for _ in range(20000):
        p = 1 / (1 + np.exp(-k * (r[:, None] - r[None, :])))
        grad = (W - G * p).sum(1)
        r += 2.0 * grad / max(G.sum(1).max(), 1)
        r -= r[anchor]
        if np.abs(grad).max() < 1e-6:
            break
    return r


def to_sgf(moves, n, komi, result, names):
    def c(m):
        if m == n * n:
            return ""
        r, col = divmod(m, n)
        return "abcdefghijklmnopqrs"[col] + "abcdefghijklmnopqrs"[r]
    body = "".join(f";{'B' if i % 2 == 0 else 'W'}[{c(m)}]" for i, m in enumerate(moves))
    return (f"(;GM[1]FF[4]SZ[{n}]KM[{komi}]RU[Tromp-Taylor]PB[{names[0]}]"
            f"PW[{names[1]}]RE[{result}]{body})\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="alphago_zero/runs/7x7")
    ap.add_argument("--max-gens", type=int, default=10)
    ap.add_argument("--games-per-pair", type=int, default=20)
    ap.add_argument("--sims", type=int, default=50)
    ap.add_argument("--final-games", type=int, default=100)
    ap.add_argument("--final-sims", type=int, default=200)
    ap.add_argument("--sample-sims", type=int, default=400)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()

    cfg = json.load(open(os.path.join(args.run, "config.json")))
    cfg["eval_open_moves"] = cfg.get("eval_open_moves", 2)
    paths = sorted(glob.glob(os.path.join(args.run, "gens", "gen_*.pt")))
    pick = np.unique(np.linspace(0, len(paths) - 1, min(args.max_gens, len(paths)))
                     .round().astype(int))
    paths = [paths[i] for i in pick]
    specs, cks = zip(*[spec_of(p, args.sims) for p in paths])
    gens = [ck["generation"] for ck in cks]
    games = [ck.get("games", 0) or 0 for ck in cks]
    print("ladder generations:", gens, flush=True)

    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(args.workers, mp_context=ctx, initializer=init_worker) as pool:
        # ---- round robin
        jobs = {}
        for (i, j) in itertools.combinations(range(len(paths)), 2):
            half = args.games_per_pair // 2
            for side in (True, False):
                jobs[pool.submit(arena_worker, specs[i], specs[j], cfg, half,
                                 side, 1000 * i + j + side)] = (i, j)
        res = {}
        for f, (i, j) in jobs.items():
            out = f.result()
            s = sum(g["winner"] == g["a_color"] for g in out)
            a, b = res.get((i, j), (0, 0))
            res[(i, j)] = (a + s, b + len(out))
        results = [(i, j, s, g) for (i, j), (s, g) in res.items()]
        elo = fit_elo(len(paths), results)
        matrix = [[None] * len(paths) for _ in paths]
        for i, j, s, g in results:
            matrix[i][j] = round(s / g, 3)
            matrix[j][i] = round(1 - s / g, 3)
        print("elo:", dict(zip(gens, elo.round(0))), flush=True)

        # ---- final strength vs fixed opponents (more search)
        best = spec_of(os.path.join(args.run, "best.pt"), args.final_sims)[0]
        final = {}
        for name, opp in (("random", dict(kind="random")), ("gen0", spec_of(paths[0], args.final_sims)[0])):
            futs = [pool.submit(arena_worker, best, opp, cfg, args.final_games // args.workers,
                                k % 2 == 0, 777 + k) for k in range(args.workers)]
            out = [g for f in futs for g in f.result()]
            final[name] = dict(games=len(out),
                               win_rate=round(float(np.mean([g["winner"] == g["a_color"] for g in out])), 3))
        print("final:", final, flush=True)

    lad = dict(generations=gens, games=games, elo=[round(float(e), 1) for e in elo],
               games_per_pair=args.games_per_pair, sims=args.sims, matrix=matrix,
               final=final)
    json.dump(lad, open(os.path.join(args.run, "ladder.json"), "w"), indent=1)

    # ---- sample games of the final network against itself
    torch.set_num_threads(4)
    net, ck = load(os.path.join(args.run, "best.pt"))
    rng = np.random.default_rng(7)
    ev = NetEvaluator(net, cfg["history"], True, rng)
    sgf = ""
    for g in range(3):
        s = GoState.new(cfg["board_size"], cfg["komi"])
        moves = []
        while not s.is_terminal():
            root = Node(s)
            search([root], ev, args.sample_sims, cfg["c_puct"], None, rng)
            a = (int(rng.choice(len(root.N), p=root.visit_policy())) if s.move_num < 2
                 else int(np.argmax(root.N)))
            moves.append(a)
            s = s.play(a)
        sc = s.score()
        res = f"{'B' if sc > 0 else 'W'}+{abs(sc)}"
        name = f"zero-gen{ck['generation']}"
        sgf += to_sgf(moves, s.n, s.komi, res, (name, name))
        print(f"sample game {g}: {len(moves)} moves, {res}")
        print(s)
    open(os.path.join(args.run, "sample_games.sgf"), "w").write(sgf)

    plot(args.run, lad)


def plot(run, lad):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = [json.loads(l) for l in open(os.path.join(run, "metrics.jsonl"))]
    g = [r["games"] for r in rows]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
    ax[0, 0].plot(lad["games"], lad["elo"], "o-", color="#23705f", label="round-robin MLE")
    ax[0, 0].plot(g, [r["elo"] for r in rows], "--", color="#999", label="gating chain")
    ax[0, 0].set(title="Elo of the best player (gen 0 = 0)", xlabel="self-play games")
    ax[0, 0].legend()
    pr = [r for r in rows if r.get("best_vs_random") is not None]
    ax[0, 1].plot([r["games"] for r in pr], [100 * r["best_vs_random"] for r in pr], "o-",
                  color="#23705f", label="vs uniform random")
    ax[0, 1].plot([r["games"] for r in pr], [100 * r["best_vs_gen0"] for r in pr], "s--",
                  color="#a2402d", label="vs untrained net + MCTS")
    ax[0, 1].set(title="Win rate of the best player (%)", xlabel="self-play games", ylim=(0, 102))
    ax[0, 1].legend()
    ax[1, 0].plot(g, [r["loss_p"] for r in rows], color="#23705f", label="policy  -pi.log p")
    ax[1, 0].plot(g, [r["loss_v"] for r in rows], color="#a2402d", label="value  (z-v)^2")
    ax[1, 0].set(title="Training loss", xlabel="self-play games")
    ax[1, 0].legend()
    ax[1, 1].plot(g, [r["mean_len"] for r in rows], color="#555", label="mean game length")
    ax[1, 1].plot(g, [100 * r["resigned"] for r in rows], color="#c58b2a", label="% games resigned")
    ax[1, 1].set(title="Self-play games", xlabel="self-play games")
    ax[1, 1].legend()
    for a in ax.flat:
        a.grid(alpha=.3)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
    os.makedirs(out, exist_ok=True)
    name = os.path.basename(os.path.normpath(run))
    fig.savefig(os.path.join(out, f"training_{name}.png"), dpi=130)
    print("saved figure")


if __name__ == "__main__":
    main()
