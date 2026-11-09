"""Go rules engine (Tromp-Taylor area scoring, positional superko).

This is the *only* domain knowledge AlphaGo Zero is given: the rules.
No patterns, no openings, no handcrafted features.

Board points are indexed p = row * n + col, row 0 at the top.
The pass move is index n*n.  Colours: BLACK = 1, WHITE = -1, EMPTY = 0.
States are immutable: play() returns a new state, so MCTS can share them.
"""
import numpy as np

BLACK, WHITE, EMPTY = 1, -1, 0
DEFAULT_KOMI = {5: 0.5, 7: 8.5, 9: 7.5, 13: 7.5, 19: 7.5}
COLS = "ABCDEFGHJKLMNOPQRST"  # Go convention skips "I"

_GEOM = {}


def geometry(n):
    """Neighbour tables + Zobrist keys for an n x n board (cached)."""
    if n in _GEOM:
        return _GEOM[n]
    nbrs, diags = [], []
    for p in range(n * n):
        r, c = divmod(p, n)
        nbrs.append(tuple(
            (r + dr) * n + (c + dc)
            for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1))
            if 0 <= r + dr < n and 0 <= c + dc < n))
        diags.append(tuple(
            (r + dr) * n + (c + dc)
            for dr, dc in ((-1, -1), (-1, 1), (1, -1), (1, 1))
            if 0 <= r + dr < n and 0 <= c + dc < n))
    rng = np.random.default_rng(20171019)  # AGZ publication date
    keys = rng.integers(1, 2 ** 62, size=(2, n * n), dtype=np.int64)
    zob = {BLACK: [int(k) for k in keys[0]], WHITE: [int(k) for k in keys[1]]}
    _GEOM[n] = (tuple(nbrs), tuple(diags), zob)
    return _GEOM[n]


class GoState:
    __slots__ = ("n", "komi", "board", "to_play", "passes", "hash", "seen",
                 "history", "move_num", "max_moves", "last_move", "_groups",
                 "_legal")

    @classmethod
    def new(cls, n=7, komi=None, max_moves=None):
        s = cls()
        s.n = n
        s.komi = DEFAULT_KOMI.get(n, 7.5) if komi is None else komi
        s.board = (EMPTY,) * (n * n)
        s.to_play = BLACK
        s.passes = 0
        s.hash = 0
        s.seen = frozenset((0,))
        s.history = (s.board,)          # most recent first, length <= T
        s.move_num = 0
        s.max_moves = max_moves or 2 * n * n   # AGZ: 19*19*2
        s.last_move = None
        s._groups = None
        s._legal = None
        return s

    # ------------------------------------------------------------------ groups
    def groups(self):
        """gid[p] -> group index; libs[g] = set of liberties; stones[g]."""
        if self._groups is None:
            n2 = self.n * self.n
            nbrs = geometry(self.n)[0]
            b = self.board
            gid = [-1] * n2
            libs, stones = [], []
            for p in range(n2):
                if b[p] == EMPTY or gid[p] >= 0:
                    continue
                col, g = b[p], len(libs)
                st, lb, stack = [p], set(), [p]
                gid[p] = g
                while stack:
                    q = stack.pop()
                    for r in nbrs[q]:
                        v = b[r]
                        if v == EMPTY:
                            lb.add(r)
                        elif v == col and gid[r] < 0:
                            gid[r] = g
                            st.append(r)
                            stack.append(r)
                libs.append(lb)
                stones.append(st)
            self._groups = (gid, libs, stones)
        return self._groups

    def _analyse(self, p):
        """For a stone move at p: (is_legal, captured_points, new_hash)."""
        if self.board[p] != EMPTY:
            return False, None, None
        nbrs, _, zob = geometry(self.n)
        gid, libs, stones = self.groups()
        c, b = self.to_play, self.board
        has_lib, cap_groups = False, set()
        for q in nbrs[p]:
            v = b[q]
            if v == EMPTY:
                has_lib = True
            elif v == c:
                if len(libs[gid[q]]) > 1:
                    has_lib = True
            elif len(libs[gid[q]]) == 1:
                cap_groups.add(gid[q])
        if not has_lib and not cap_groups:
            return False, None, None            # suicide
        h = self.hash ^ zob[c][p]
        captured = []
        for g in cap_groups:
            for s in stones[g]:
                h ^= zob[-c][s]
                captured.append(s)
        if h in self.seen:
            return False, None, None            # positional superko
        return True, captured, h

    # ------------------------------------------------------------------ moves
    def legal_mask(self):
        """Boolean array of length n*n+1 (last entry = pass, always legal)."""
        if self._legal is None:
            n2 = self.n * self.n
            m = np.zeros(n2 + 1, dtype=bool)
            m[n2] = True
            for p in range(n2):
                if self.board[p] == EMPTY and self._analyse(p)[0]:
                    m[p] = True
            self._legal = m
        return self._legal

    def legal_moves(self):
        return np.flatnonzero(self.legal_mask()).tolist()

    def play(self, move, history_len=8):
        n2 = self.n * self.n
        s = GoState.__new__(GoState)
        s.n, s.komi, s.max_moves = self.n, self.komi, self.max_moves
        s.to_play = -self.to_play
        s.move_num = self.move_num + 1
        s.last_move = move
        s._groups = None
        s._legal = None
        if move == n2:                          # pass
            s.board, s.hash, s.seen = self.board, self.hash, self.seen
            s.passes = self.passes + 1
        else:
            ok, captured, h = self._analyse(move)
            if not ok:
                raise ValueError(f"illegal move {self.move_to_str(move)}")
            b = list(self.board)
            b[move] = self.to_play
            for q in captured:
                b[q] = EMPTY
            s.board = tuple(b)
            s.hash = h
            s.seen = self.seen | {h}
            s.passes = 0
        s.history = (s.board,) + self.history[:history_len - 1]
        return s

    def is_terminal(self):
        return self.passes >= 2 or self.move_num >= self.max_moves

    # ---------------------------------------------------------------- scoring
    def ownership(self):
        """Tromp-Taylor: each point owned by BLACK/WHITE/EMPTY(dame)."""
        nbrs = geometry(self.n)[0]
        b = self.board
        own = list(b)
        seen = set()
        for p in range(self.n * self.n):
            if b[p] != EMPTY or p in seen:
                continue
            region, border, stack = [p], set(), [p]
            seen.add(p)
            while stack:
                q = stack.pop()
                for r in nbrs[q]:
                    if b[r] == EMPTY:
                        if r not in seen:
                            seen.add(r)
                            region.append(r)
                            stack.append(r)
                    else:
                        border.add(b[r])
            owner = border.pop() if len(border) == 1 else EMPTY
            for q in region:
                own[q] = owner
        return own

    def score(self):
        """Black area - white area - komi (positive = black ahead)."""
        own = self.ownership()
        return own.count(BLACK) - own.count(WHITE) - self.komi

    def winner(self):
        return BLACK if self.score() > 0 else WHITE

    # --------------------------------------------------------------- features
    def features(self, history_len=8):
        """AGZ input: [X_t, Y_t, X_{t-1}, Y_{t-1}, ..., C]  -> (2T+1, n, n).

        X = stones of the player to move, Y = opponent stones,
        C = all ones if black to move else zeros.
        """
        n, T = self.n, history_len
        planes = np.zeros((2 * T + 1, n * n), dtype=np.float32)
        h = np.asarray(self.history[:T], dtype=np.int8)
        k = h.shape[0]
        planes[0:2 * k:2] = h == self.to_play
        planes[1:2 * k:2] = h == -self.to_play
        if self.to_play == BLACK:
            planes[2 * T] = 1.0
        return planes.reshape(2 * T + 1, n, n)

    # ------------------------------------------------------------------ utils
    def move_to_str(self, move):
        n = self.n
        if move == n * n:
            return "pass"
        r, c = divmod(move, n)
        return f"{COLS[c]}{n - r}"

    def str_to_move(self, txt):
        txt = txt.strip().upper()
        if txt in ("PASS", "P"):
            return self.n * self.n
        c = COLS.index(txt[0])
        r = self.n - int(txt[1:])
        if not (0 <= c < self.n and 0 <= r < self.n):
            raise ValueError(txt)
        return r * self.n + c

    def is_eye(self, p, color):
        """Cheap 'true eye' test used only by the random baseline player."""
        nbrs, diags, _ = geometry(self.n)
        b = self.board
        if b[p] != EMPTY or any(b[q] != color for q in nbrs[p]):
            return False
        bad = sum(1 for q in diags[p] if b[q] == -color)
        edge = len(diags[p]) < 4
        return bad == 0 if edge else bad <= 1

    def __str__(self):
        n = self.n
        lines = []
        for r in range(n):
            row = []
            for c in range(n):
                v = self.board[r * n + c]
                ch = "X" if v == BLACK else "O" if v == WHITE else "."
                if self.last_move == r * n + c:
                    ch = ch.lower() if ch != "." else ch
                row.append(ch)
            lines.append(f"{n - r:2d} " + " ".join(row))
        lines.append("   " + " ".join(COLS[:n]))
        return "\n".join(lines)


# ------------------------------------------------------------------ symmetry
def sym_planes(x, k):
    """Apply dihedral symmetry k (0..7) to (..., n, n) arrays."""
    x = np.rot90(x, k % 4, axes=(-2, -1))
    if k >= 4:
        x = np.flip(x, axis=-1)
    return x


def sym_policy(pi, n, k):
    """Apply the same symmetry to a policy vector of length n*n+1."""
    board = sym_planes(pi[..., : n * n].reshape(pi.shape[:-1] + (n, n)), k)
    return np.concatenate(
        [board.reshape(pi.shape[:-1] + (n * n,)), pi[..., n * n:]], axis=-1)


def inv_sym_policy(pi, n, k):
    """Undo symmetry k on a policy vector (flip first, then rotate back)."""
    board = pi[..., : n * n].reshape(pi.shape[:-1] + (n, n))
    if k >= 4:
        board = np.flip(board, axis=-1)
    board = np.rot90(board, -(k % 4), axes=(-2, -1))
    return np.concatenate(
        [board.reshape(pi.shape[:-1] + (n * n,)), pi[..., n * n:]], axis=-1)
