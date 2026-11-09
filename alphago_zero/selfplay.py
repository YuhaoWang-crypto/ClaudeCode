"""Self-play (data generation) and arena matches (evaluation / gating).

Both run inside worker processes; each worker plays many games at once so
that one network call evaluates one leaf from every game (see mcts.py).
"""
import numpy as np
import torch

from .go import GoState
from .mcts import Node, NetEvaluator, search
from .model import AlphaZeroNet


def init_worker():
    torch.set_num_threads(1)


def _net(net_cfg, state_dict):
    net = AlphaZeroNet(**net_cfg)
    net.load_state_dict(state_dict)
    return net.eval()


# ------------------------------------------------------------------ self-play
def selfplay_worker(net_cfg, state_dict, cfg, num_games, seed, v_resign):
    """Play `num_games` games of best-player vs itself.

    Returns one record per game with the training triples (s_t, pi_t, z_t).
    v_resign=None disables resignation.
    """
    rng = np.random.default_rng(seed)
    ev = NetEvaluator(_net(net_cfg, state_dict), cfg["history"], True, rng)
    n2 = cfg["board_size"] ** 2
    games = []
    for _ in range(num_games):
        games.append(dict(
            root=Node(GoState.new(cfg["board_size"], cfg["komi"])),
            feats=[], pis=[], players=[],
            no_resign=(v_resign is None) or rng.random() < cfg["no_resign_frac"],
            min_v={1: 1.0, -1: 1.0},
            winner=None, resigned=False, score=None))
    active = list(games)
    while active:
        search([g["root"] for g in active], ev, cfg["sims"], cfg["c_puct"],
               noise=(cfg["dir_alpha"], cfg["dir_eps"]), rng=rng)
        still = []
        for g in active:
            root = g["root"]
            s = root.state
            pi = root.visit_policy()
            g["feats"].append(s.features(cfg["history"]).astype(np.uint8))
            g["pis"].append(pi.astype(np.float32))
            g["players"].append(s.to_play)

            # resignation: root value AND best-child value below v_resign
            best = int(np.argmax(root.N))
            v_root = root.value()
            v_best = root.W[best] / max(root.N[best], 1)
            v_hi = max(v_root, v_best)
            g["min_v"][s.to_play] = min(g["min_v"][s.to_play], v_hi)
            if v_resign is not None and v_hi < v_resign and not g["no_resign"]:
                g["winner"], g["resigned"] = -s.to_play, True
                continue

            # move choice: tau = 1 (sample) early, tau -> 0 (argmax) later
            if s.move_num < cfg["temp_moves"]:
                a = int(rng.choice(n2 + 1, p=pi))
            else:
                top = np.flatnonzero(root.N == root.N.max())
                a = int(rng.choice(top))
            g["root"] = root.child(a)          # tree reuse: subtree survives
            if g["root"].state.is_terminal():
                g["winner"] = g["root"].state.winner()
                g["score"] = g["root"].state.score()
                continue
            still.append(g)
        active = still

    records = []
    for g in games:
        players = np.array(g["players"], dtype=np.int8)
        z = np.where(players == g["winner"], 1.0, -1.0).astype(np.float32)
        records.append(dict(
            feats=np.stack(g["feats"]), pis=np.stack(g["pis"]), z=z,
            winner=g["winner"], length=len(players), resigned=g["resigned"],
            no_resign=g["no_resign"], score=g["score"],
            # lowest value the eventual winner ever saw; in games played to
            # the end this calibrates v_resign (false positive = winner dipped
            # below the threshold at some point)
            winner_min_v=g["min_v"][g["winner"]] if g["no_resign"] else None))
    return records, ev.positions


# --------------------------------------------------------------------- arena
def random_move(state, rng):
    """Uniform random legal move that does not fill its own eye."""
    n2 = state.n * state.n
    cand = [m for m in state.legal_moves()
            if m != n2 and not state.is_eye(m, state.to_play)]
    return int(rng.choice(cand)) if cand else n2


def arena_worker(spec_a, spec_b, cfg, num_games, a_black_first, seed):
    """Play A vs B.  spec = {'kind': 'mcts'|'random', 'net_cfg', 'state_dict', 'sims'}.

    Colours alternate between games.  No Dirichlet noise; moves are the most
    visited child (tau -> 0), except the first `eval_open_moves` moves which
    are sampled (tau = 1) so that the small number of games are not copies.
    Returns list of dicts {a_color, winner, length}.
    """
    rng = np.random.default_rng(seed)
    evs = {}
    for name, spec in (("a", spec_a), ("b", spec_b)):
        if spec["kind"] == "mcts":
            evs[name] = NetEvaluator(_net(spec["net_cfg"], spec["state_dict"]),
                                     cfg["history"], True, rng)
    games = []
    for i in range(num_games):
        a_color = 1 if (i % 2 == 0) == a_black_first else -1
        games.append(dict(state=GoState.new(cfg["board_size"], cfg["komi"]),
                          a_color=a_color, winner=None))
    open_moves = cfg.get("eval_open_moves", 2)
    active = list(games)
    while active:
        for name, spec in (("a", spec_a), ("b", spec_b)):
            mine = [g for g in active if not g["state"].is_terminal()
                    and (g["state"].to_play == g["a_color"]) == (name == "a")]
            if not mine:
                continue
            if spec["kind"] == "random":
                for g in mine:
                    g["state"] = g["state"].play(random_move(g["state"], rng))
                continue
            roots = [Node(g["state"]) for g in mine]
            search(roots, evs[name], spec["sims"], cfg["c_puct"], None, rng)
            for g, r in zip(mine, roots):
                if g["state"].move_num < open_moves:
                    a = int(rng.choice(len(r.N), p=r.visit_policy()))
                else:
                    a = int(rng.choice(np.flatnonzero(r.N == r.N.max())))
                g["state"] = g["state"].play(a)
        still = []
        for g in active:
            if g["state"].is_terminal():
                g["winner"] = g["state"].winner()
            else:
                still.append(g)
        active = still
    return [dict(a_color=g["a_color"], winner=g["winner"],
                 length=g["state"].move_num) for g in games]
