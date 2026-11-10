"""
Reference (numpy-only) port of the chess-alpha-zero encoding and network.

Mirrors zeta36/chess-alpha-zero exactly:
  * 1968 move labels (config.create_uci_labels) and their rank-mirrored flip
  * 18x8x8 input planes, always from the side-to-move's point of view
    (env.chess_env.canon_input_planes, including its quirks)
  * the residual policy/value network (agent.model_chess.ChessModel.build)

It is used to convert the weights and to produce test vectors that the
browser engine is checked against.
"""
import numpy as np

PIECES_ORDER = "KQRBNPkqrbnp"
BN_EPS = 1e-3  # Keras BatchNormalization default used by the original model


# ---------------------------------------------------------------- labels
def create_uci_labels():
    labels = []
    letters = "abcdefgh"
    numbers = "12345678"
    for l1 in range(8):
        for n1 in range(8):
            dest = [(t, n1) for t in range(8)] + \
                   [(l1, t) for t in range(8)] + \
                   [(l1 + t, n1 + t) for t in range(-7, 8)] + \
                   [(l1 + t, n1 - t) for t in range(-7, 8)] + \
                   [(l1 + a, n1 + b) for (a, b) in
                    [(-2, -1), (-1, -2), (-2, 1), (1, -2), (2, -1), (-1, 2), (2, 1), (1, 2)]]
            for (l2, n2) in dest:
                if (l1, n1) != (l2, n2) and 0 <= l2 < 8 and 0 <= n2 < 8:
                    labels.append(letters[l1] + numbers[n1] + letters[l2] + numbers[n2])
    for l1 in range(8):
        l = letters[l1]
        for p in "qrbn":
            labels.append(l + "2" + l + "1" + p)
            labels.append(l + "7" + l + "8" + p)
            if l1 > 0:
                ll = letters[l1 - 1]
                labels.append(l + "2" + ll + "1" + p)
                labels.append(l + "7" + ll + "8" + p)
            if l1 < 7:
                lr = letters[l1 + 1]
                labels.append(l + "2" + lr + "1" + p)
                labels.append(l + "7" + lr + "8" + p)
    return labels


LABELS = create_uci_labels()
assert len(LABELS) == 1968


def mirror_uci(m):
    """Mirror a UCI move top<->bottom (a2a4 -> a7a5)."""
    return "".join(str(9 - int(c)) if c.isdigit() else c for c in m)


# ---------------------------------------------------------------- planes
def _swapcase(s):
    return "".join(c.lower() if c.isupper() else c.upper() if c.islower() else c for c in s)


def maybe_flip_fen(fen, flip):
    if not flip:
        return fen
    f = fen.split(" ")
    rows = f[0].split("/")
    return "/".join(_swapcase(r) for r in reversed(rows)) + " " + \
        ("w" if f[1] == "b" else "b") + " " + "".join(sorted(_swapcase(f[2]))) + \
        " " + f[3] + " " + f[4] + " " + f[5]   # NB: ep square is NOT mirrored (original quirk)


def canon_input_planes(fen):
    """(18, 8, 8) planes; row 0 = rank 8 of the (possibly flipped) board."""
    fen = maybe_flip_fen(fen, fen.split(" ")[1] == "b")
    f = fen.split(" ")
    planes = np.zeros((18, 8, 8), dtype=np.float32)
    board = f[0]
    for d in "2345678":
        board = board.replace(d, "1" * int(d))
    board = board.replace("/", "")
    for i, c in enumerate(board):
        if c.isalpha():
            planes[PIECES_ORDER.index(c), i // 8, i % 8] = 1
    for k, c in enumerate("KQkq"):
        planes[12 + k] = float(c in f[2])
    planes[16] = float(int(f[4]))           # raw half-move clock
    if f[3] != "-":
        planes[17, 8 - int(f[3][1]), ord(f[3][0]) - 97] = 1
    return planes


# ---------------------------------------------------------------- network
def load_h5(path):
    import h5py
    w = {}
    with h5py.File(path, "r") as f:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                layer = name.split("/")[0]
                w.setdefault(layer, {})[name.split("/")[-1].split(":")[0]] = np.array(obj)
        f.visititems(visit)
    return w


def _conv_same(x, k):
    """x: (C,8,8), k: (kh,kw,Cin,Cout) Keras layout -> (Cout,8,8)."""
    kh, kw, cin, cout = k.shape
    ph, pw = kh // 2, kw // 2
    xp = np.pad(x, ((0, 0), (ph, ph), (pw, pw)))
    out = np.zeros((cout, 8, 8), dtype=np.float64)
    for i in range(kh):
        for j in range(kw):
            patch = xp[:, i:i + 8, j:j + 8]              # (Cin,8,8)
            out += np.einsum("chw,co->ohw", patch, k[i, j])
    return out


def _bn(x, p):
    s = p["gamma"] / np.sqrt(p["moving_variance"] + BN_EPS)
    return x * s[:, None, None] + (p["beta"] - p["moving_mean"] * s)[:, None, None]


def find(w, prefix):
    return w[next(k for k in w if k.startswith(prefix))]


def forward_h5(w, planes):
    """Exact float forward of the original Keras graph. Returns (policy[1968], value)."""
    relu = lambda a: np.maximum(a, 0)
    x = relu(_bn(_conv_same(planes, find(w, "input_conv")["kernel"]), w["input_batchnorm"]))
    for i in range(1, 8):
        y = relu(_bn(_conv_same(x, find(w, f"res{i}_conv1")["kernel"]), w[f"res{i}_batchnorm1"]))
        y = _bn(_conv_same(y, find(w, f"res{i}_conv2")["kernel"]), w[f"res{i}_batchnorm2"])
        x = relu(x + y)
    p = relu(_bn(_conv_same(x, find(w, "policy_conv")["kernel"]), w["policy_batchnorm"])).reshape(-1)
    logits = p @ w["policy_out"]["kernel"] + w["policy_out"]["bias"]
    e = np.exp(logits - logits.max())
    policy = e / e.sum()
    v = relu(_bn(_conv_same(x, find(w, "value_conv")["kernel"]), w["value_batchnorm"])).reshape(-1)
    v = relu(v @ w["value_dense"]["kernel"] + w["value_dense"]["bias"])
    value = np.tanh(v @ w["value_out"]["kernel"] + w["value_out"]["bias"])[0]
    return policy, float(value)
