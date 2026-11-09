"""PUCT Monte-Carlo tree search, exactly the AlphaGo Zero variant.

Each edge (s, a) stores N(s,a), W(s,a), Q(s,a)=W/N and prior P(s,a).
One simulation:
  1. Select   a = argmax Q(s,a) + U(s,a),
              U(s,a) = c_puct * P(s,a) * sqrt(sum_b N(s,b)) / (1 + N(s,a))
  2. Expand   the leaf s_L with ONE network call (p, v) = f(d_i(s_L)),
              d_i a random one of the 8 board symmetries.
  3. Backup   N += 1, W += v along the path (sign flips every ply).
No rollouts: the value head replaces the random playouts of the 2016 AlphaGo.

Engineering detail: instead of the paper's virtual-loss threads we run one
simulation in each of many independent games at the same time and send all
their leaves to the network as one batch.  Same search, better CPU usage.
"""
import math
import numpy as np
import torch

from .go import sym_planes, inv_sym_policy


class NetEvaluator:
    """Wraps the network: list[GoState] -> (priors (B, n*n+1), values (B,))."""

    def __init__(self, net, history=8, use_symmetry=True, rng=None):
        self.net = net.eval()
        self.history = history
        self.use_symmetry = use_symmetry
        self.rng = rng or np.random.default_rng()
        self.calls = 0
        self.positions = 0

    def __call__(self, states):
        n = states[0].n
        x = np.stack([s.features(self.history) for s in states])
        B = len(states)
        ks = (self.rng.integers(0, 8, B) if self.use_symmetry
              else np.zeros(B, dtype=int))
        if self.use_symmetry:
            x = np.stack([sym_planes(x[i], ks[i]) for i in range(B)])
        with torch.inference_mode():
            logits, v = self.net(torch.from_numpy(np.ascontiguousarray(x)))
        logits = logits.numpy().astype(np.float64)
        if self.use_symmetry:
            logits = np.stack([inv_sym_policy(logits[i], n, ks[i])
                               for i in range(B)])
        legal = np.stack([s.legal_mask() for s in states])
        logits = np.where(legal, logits, -np.inf)
        logits -= logits.max(axis=1, keepdims=True)
        p = np.exp(logits)
        p /= p.sum(axis=1, keepdims=True)
        self.calls += 1
        self.positions += B
        return p, v.numpy().astype(np.float64)


class Node:
    __slots__ = ("state", "P", "N", "W", "Nsum", "children", "pen", "expanded")

    def __init__(self, state):
        self.state = state
        self.children = {}
        self.expanded = False

    def expand(self, priors):
        legal = self.state.legal_mask()
        self.P = priors
        self.N = np.zeros(len(priors))
        self.W = np.zeros(len(priors))
        self.Nsum = 0
        self.pen = np.where(legal, 0.0, -np.inf)   # forbid illegal moves
        self.expanded = True

    def select(self, c_puct):
        N = self.N
        q = self.W / np.maximum(N, 1.0)            # Q = 0 for unvisited (paper)
        u = c_puct * self.P * (math.sqrt(max(self.Nsum, 1)) / (1.0 + N))
        return int(np.argmax(q + u + self.pen))

    def add_dirichlet(self, alpha, eps, rng):
        idx = np.flatnonzero(self.pen == 0.0)
        if len(idx) > 1:
            eta = rng.dirichlet([alpha] * len(idx))
            self.P = self.P.copy()
            self.P[idx] = (1 - eps) * self.P[idx] + eps * eta

    # statistics used for move choice / training targets / resignation
    def visit_policy(self):
        return self.N / max(self.Nsum, 1)

    def value(self):
        """Mean value of the root from the perspective of the player to move."""
        return float(self.W.sum() / self.Nsum) if self.Nsum else 0.0

    def child(self, a):
        c = self.children.get(a)
        if c is None:
            c = Node(self.state.play(a))
        return c


def terminal_value(state):
    """+1 if the player to move at this finished position has won."""
    return 1.0 if state.winner() == state.to_play else -1.0


def _backup(path, v):
    # v is from the leaf player's view; each edge stores its parent's view.
    for node, a in reversed(path):
        v = -v
        node.N[a] += 1
        node.W[a] += v
        node.Nsum += 1


def search(roots, evaluator, num_sims, c_puct=1.5, noise=None, rng=None):
    """Run `num_sims` simulations from every root (batched across roots)."""
    todo = [r for r in roots if not r.expanded and not r.state.is_terminal()]
    if todo:
        P, _ = evaluator([r.state for r in todo])
        for r, p in zip(todo, P):
            r.expand(p)
    if noise is not None:
        alpha, eps = noise
        for r in roots:
            if r.expanded:
                r.add_dirichlet(alpha, eps, rng)

    live = [r for r in roots if r.expanded]
    for _ in range(num_sims):
        leaves, paths = [], []
        for root in live:
            node, path = root, []
            while node.expanded:                       # 1. select
                a = node.select(c_puct)
                path.append((node, a))
                nxt = node.children.get(a)
                if nxt is None:
                    nxt = Node(node.state.play(a))
                    node.children[a] = nxt
                node = nxt
            if node.state.is_terminal():               # game over: exact value
                _backup(path, terminal_value(node.state))
            else:
                leaves.append(node)
                paths.append(path)
        if leaves:                                     # 2. evaluate (batched)
            P, V = evaluator([l.state for l in leaves])
            for leaf, path, p, v in zip(leaves, paths, P, V):
                leaf.expand(p)                         #    expand
                _backup(path, v)                       # 3. backup
    return roots
