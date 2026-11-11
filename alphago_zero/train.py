"""AlphaGo Zero training pipeline: self-play -> optimise -> evaluate (gate).

    python -m alphago_zero.train                      # 7x7, ~2.5 h on 4 CPU cores
    python -m alphago_zero.train --time-budget-min 10 --sims 32   # smoke test

The paper runs three components asynchronously.  Here they run in turn,
which is the same algorithm with a fixed schedule:

  1. Self-play  : the current *best* player plays games against itself with
                  MCTS; every position yields (s_t, pi_t, z_t).
  2. Optimise   : the *training* network is updated by SGD on positions sampled
                  from the most recent games, loss = (z-v)^2 - pi.log p + c|theta|^2.
  3. Evaluate   : the new checkpoint plays the best player; if it wins >= 55%
                  it becomes the new best player (and generates future data).
"""
import argparse
import collections
import json
import math
import multiprocessing as mp
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import torch
import torch.nn.functional as F

from .config import Config
from .go import sym_planes, sym_policy
from .model import build
from .selfplay import arena_worker, init_worker, selfplay_worker


def cpu_state(net):
    return {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}


def save_ckpt(path, net_cfg, state_dict, cfg, **extra):
    torch.save(dict(net_cfg=net_cfg, state_dict=state_dict, cfg=cfg.to_dict(),
                    **extra), path)


def split(total, parts):
    base = [total // parts] * parts
    for i in range(total % parts):
        base[i] += 1
    return [b for b in base if b > 0]


def run_selfplay(pool, net_cfg, sd, cfg, v_resign, seed):
    d = cfg.to_dict()
    futs = [pool.submit(selfplay_worker, net_cfg, sd, d, k, seed + i, v_resign)
            for i, k in enumerate(split(cfg.games_per_iter, cfg.workers))]
    recs, evals = [], 0
    for f in futs:
        r, e = f.result()
        recs.extend(r)
        evals += e
    return recs, evals


def run_arena(pool, spec_a, spec_b, cfg, games, seed):
    """Returns the score of A (wins / games)."""
    d = cfg.to_dict()
    futs = [pool.submit(arena_worker, spec_a, spec_b, d, k, i % 2 == 0, seed + i)
            for i, k in enumerate(split(games, cfg.workers))]
    res = [g for f in futs for g in f.result()]
    wins = sum(g["winner"] == g["a_color"] for g in res)
    a_black = [g for g in res if g["a_color"] == 1]
    return wins / len(res), dict(
        games=len(res), wins=wins,
        wins_as_black=sum(g["winner"] == 1 for g in a_black),
        games_as_black=len(a_black),
        mean_len=float(np.mean([g["length"] for g in res])))


def elo_diff(score, games):
    s = min(max(score, 0.5 / games), 1 - 0.5 / games)
    return 400 * math.log10(s / (1 - s))


class Replay:
    """Positions from the most recent `window_games` self-play games."""

    def __init__(self, window_games):
        self.games = collections.deque(maxlen=window_games)
        self._cache = None

    def add(self, records):
        self.games.extend(records)
        self._cache = None

    # Each iteration's games are also written to disk so that a restarted
    # run (--resume) continues with the same training window.
    @staticmethod
    def save_chunk(out_dir, it, records):
        d = os.path.join(out_dir, "buffer")
        os.makedirs(d, exist_ok=True)
        np.savez_compressed(
            os.path.join(d, f"iter_{it:05d}.npz"),
            feats=np.concatenate([r["feats"] for r in records]),
            pis=np.concatenate([r["pis"] for r in records]).astype(np.float16),
            z=np.concatenate([r["z"] for r in records]).astype(np.int8),
            lengths=np.array([r["length"] for r in records]))

    def load_chunks(self, out_dir, upto_iter, games_per_chunk):
        d = os.path.join(out_dir, "buffer")
        if not os.path.isdir(d):
            return 0
        need = -(-self.games.maxlen // games_per_chunk) + 1
        files = sorted(f for f in os.listdir(d)
                       if int(f[5:10]) <= upto_iter)[-need:]
        for f in files:
            c = np.load(os.path.join(d, f))
            ends = np.cumsum(c["lengths"])
            for a, b in zip(np.r_[0, ends[:-1]], ends):
                self.games.append(dict(feats=c["feats"][a:b],
                                       pis=c["pis"][a:b].astype(np.float32),
                                       z=c["z"][a:b].astype(np.float32)))
        self._cache = None
        return len(files)

    @staticmethod
    def prune(out_dir, it, keep):
        d = os.path.join(out_dir, "buffer")
        for f in os.listdir(d):
            if int(f[5:10]) <= it - keep:
                os.remove(os.path.join(d, f))

    def arrays(self):
        if self._cache is None:
            self._cache = (np.concatenate([g["feats"] for g in self.games]),
                           np.concatenate([g["pis"] for g in self.games]),
                           np.concatenate([g["z"] for g in self.games]))
        return self._cache

    def sample(self, batch, rng, n):
        X, P, Z = self.arrays()
        idx = rng.integers(0, len(Z), batch)
        x = X[idx].astype(np.float32)
        p = P[idx].copy()
        z = Z[idx]
        ks = rng.integers(0, 8, batch)        # data augmentation: 8 symmetries
        for k in range(1, 8):
            m = ks == k
            if m.any():
                x[m] = sym_planes(x[m], k)
                p[m] = sym_policy(p[m], n, k)
        return (torch.from_numpy(np.ascontiguousarray(x)),
                torch.from_numpy(np.ascontiguousarray(p)),
                torch.from_numpy(z))


def train_steps(net, opt, replay, cfg, steps, rng):
    net.train()
    tot = collections.Counter()
    for _ in range(steps):
        x, pi, z = replay.sample(cfg.batch_size, rng, cfg.board_size)
        logits, v = net(x)
        loss_v = F.mse_loss(v, z)                                   # (z - v)^2
        loss_p = -(pi * F.log_softmax(logits, dim=1)).sum(1).mean()  # -pi.log p
        loss = loss_v + loss_p                  # + c|theta|^2 via weight_decay
        opt.zero_grad()
        loss.backward()
        opt.step()
        tot["v"] += loss_v.item()
        tot["p"] += loss_p.item()
    net.eval()
    return {k: v / max(steps, 1) for k, v in tot.items()}


def calibrate_resign(min_vals, target):
    """Pick v_resign so that < `target` of eventual winners would have resigned."""
    if len(min_vals) < 50:
        return None
    v = float(np.quantile(np.array(min_vals), target))
    return float(np.clip(v - 0.02, -0.99, -0.8))


def parse():
    ap = argparse.ArgumentParser()
    c = Config()
    for k, v in c.to_dict().items():
        flag = "--" + k.replace("_", "-")
        if isinstance(v, bool):
            ap.add_argument(flag, type=lambda s: s.lower() in ("1", "true", "yes"),
                            default=v)
        elif isinstance(v, tuple):
            ap.add_argument(flag, type=float, nargs="+", default=v)
        else:
            ap.add_argument(flag, type=type(v), default=v)
    ap.add_argument("--resume", action="store_true")
    a = vars(ap.parse_args())
    resume = a.pop("resume")
    a = {k: tuple(v) if isinstance(v, list) else v for k, v in a.items()}
    return Config(**a), resume


def main():
    cfg, resume = parse()
    os.makedirs(os.path.join(cfg.out_dir, "gens"), exist_ok=True)
    json.dump(cfg.to_dict(), open(os.path.join(cfg.out_dir, "config.json"), "w"),
              indent=1)
    torch.manual_seed(cfg.seed)
    rng = np.random.default_rng(cfg.seed)
    torch.set_num_threads(cfg.workers)

    net = build(cfg)
    net_cfg = net.cfg
    opt = torch.optim.SGD(net.parameters(), lr=cfg.lr, momentum=cfg.momentum,
                          weight_decay=cfg.weight_decay)
    state = dict(iteration=0, generation=0, elo=0.0, games=0, positions=0,
                 steps=0, elapsed=0.0, v_resign=None)
    best_sd = cpu_state(net)
    gen0_sd = best_sd
    latest = os.path.join(cfg.out_dir, "latest.pt")
    if resume and os.path.exists(latest):
        ck = torch.load(latest, weights_only=False)
        net.load_state_dict(ck["state_dict"])
        opt.load_state_dict(ck["optim"])
        state = ck["state"]
        best_sd = torch.load(os.path.join(cfg.out_dir, "best.pt"),
                             weights_only=False)["state_dict"]
        gen0_sd = torch.load(os.path.join(cfg.out_dir, "gens", "gen_000.pt"),
                             weights_only=False)["state_dict"]
        print(f"resumed at iteration {state['iteration']}", flush=True)
    else:
        save_ckpt(os.path.join(cfg.out_dir, "gens", "gen_000.pt"), net_cfg,
                  best_sd, cfg, generation=0, elo=0.0)
        save_ckpt(os.path.join(cfg.out_dir, "best.pt"), net_cfg, best_sd, cfg,
                  generation=0, elo=0.0)
    replay = Replay(cfg.window_games)
    keep_chunks = -(-cfg.window_games // cfg.games_per_iter) + 1
    if resume:
        n = replay.load_chunks(cfg.out_dir, state["iteration"], cfg.games_per_iter)
        print(f"reloaded {len(replay.games)} games from {n} buffer files", flush=True)
    resign_vals = collections.deque(maxlen=600)
    log = open(os.path.join(cfg.out_dir, "metrics.jsonl"), "a")
    budget = cfg.time_budget_min * 60
    t_start = time.time() - state["elapsed"]

    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(cfg.workers, mp_context=ctx,
                             initializer=init_worker) as pool:
        while state["iteration"] < cfg.iterations:
            it = state["iteration"] = state["iteration"] + 1
            row = dict(iteration=it)
            seed = cfg.seed * 1_000_003 + it * 1000

            # ---- 1. self-play with the best player ------------------------
            t0 = time.time()
            v_res = state["v_resign"] if it > cfg.resign_start_iter else None
            recs, evals = run_selfplay(pool, net_cfg, best_sd, cfg, v_res, seed)
            replay.add(recs)
            Replay.save_chunk(cfg.out_dir, it, recs)
            Replay.prune(cfg.out_dir, it, keep_chunks)
            for r in recs:
                if r["winner_min_v"] is not None:
                    resign_vals.append(r["winner_min_v"])
            new_v = calibrate_resign(resign_vals, cfg.resign_fp_target)
            if new_v is not None:
                state["v_resign"] = new_v
            calib = [r for r in recs if r["no_resign"] and v_res is not None]
            new_pos = int(sum(r["length"] for r in recs))
            state["games"] += len(recs)
            state["positions"] += new_pos
            row.update(
                selfplay_s=round(time.time() - t0, 1), games=state["games"],
                positions=state["positions"], evals=evals,
                mean_len=round(float(np.mean([r["length"] for r in recs])), 1),
                black_win=round(float(np.mean([r["winner"] == 1 for r in recs])), 3),
                resigned=round(float(np.mean([r["resigned"] for r in recs])), 3),
                v_resign=v_res,
                resign_fp=(round(float(np.mean([r["winner_min_v"] < v_res
                                                for r in calib])), 3)
                           if calib else None))

            # ---- 2. optimise the training network -------------------------
            t0 = time.time()
            frac = (time.time() - t_start) / budget
            lr = cfg.lr * 0.1 ** sum(frac >= m for m in cfg.lr_milestones)
            for g in opt.param_groups:
                g["lr"] = lr
            steps = max(1, round(cfg.sample_reuse * new_pos / cfg.batch_size))
            losses = train_steps(net, opt, replay, cfg, steps, rng)
            state["steps"] += steps
            row.update(train_s=round(time.time() - t0, 1), steps=state["steps"],
                       lr=lr, loss_v=round(losses["v"], 4),
                       loss_p=round(losses["p"], 4),
                       window_positions=len(replay.arrays()[2]))

            # ---- 3. evaluate: candidate vs best (gatekeeper) --------------
            t0 = time.time()
            cand_sd = cpu_state(net)
            spec = lambda sd: dict(kind="mcts", net_cfg=net_cfg, state_dict=sd,
                                   sims=cfg.eval_sims)
            if cfg.gating:
                score, info = run_arena(pool, spec(cand_sd), spec(best_sd), cfg,
                                        cfg.eval_games, seed + 500)
                promote = score >= cfg.gate_threshold
                row.update(gate_score=round(score, 3), gate=info)
            else:
                score, promote = None, True
            if promote:
                if score is not None:
                    state["elo"] += elo_diff(score, cfg.eval_games)
                state["generation"] += 1
                best_sd = cand_sd
                gpath = os.path.join(cfg.out_dir, "gens",
                                     f"gen_{state['generation']:03d}.pt")
                for p in (gpath, os.path.join(cfg.out_dir, "best.pt")):
                    save_ckpt(p, net_cfg, best_sd, cfg, iteration=it,
                              generation=state["generation"], elo=state["elo"],
                              games=state["games"])
            row.update(promoted=promote, generation=state["generation"],
                       elo=round(state["elo"], 1))

            # ---- fixed yardsticks: random player and the untrained net ----
            if it % cfg.random_eval_every == 0:
                s_rand, _ = run_arena(pool, spec(best_sd), dict(kind="random"),
                                      cfg, cfg.random_eval_games, seed + 700)
                s_gen0, _ = run_arena(pool, spec(best_sd), spec(gen0_sd), cfg,
                                      cfg.random_eval_games, seed + 800)
                row.update(best_vs_random=round(s_rand, 3),
                           best_vs_gen0=round(s_gen0, 3))
            row["eval_s"] = round(time.time() - t0, 1)

            state["elapsed"] = time.time() - t_start
            row["elapsed_min"] = round(state["elapsed"] / 60, 1)
            torch.save(dict(state_dict=cpu_state(net), optim=opt.state_dict(),
                            state=state, net_cfg=net_cfg, cfg=cfg.to_dict()),
                       latest)
            log.write(json.dumps(row) + "\n")
            log.flush()
            print(json.dumps(row), flush=True)
            if state["elapsed"] >= budget:
                print("time budget reached", flush=True)
                break


if __name__ == "__main__":
    main()
