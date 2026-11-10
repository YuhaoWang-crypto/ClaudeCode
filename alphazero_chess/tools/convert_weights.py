"""
Convert the pretrained chess-alpha-zero Keras weights (model_best_weight.h5)
into a compact file the browser engine can load.

  * BatchNorm is folded into the preceding conv (kernel * s, bias = beta - mean * s)
  * weights are stored as float16 (2x smaller than the original, ~lossless)
    and the file is split into base64 text parts of ~8 MB for easy hosting
  * the policy dense rows are reordered from channels-first to channels-last
    flatten order so the browser can use NHWC tensors directly

Usage:
  python convert_weights.py path/to/model_best_weight.h5 ../web/model/weights \
         --vectors ../tests/vectors.json
"""
import argparse
import base64
import json
import struct
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import az_reference as ref  # noqa: E402

PART_BYTES = 6_000_000  # 8 MB per part once base64-encoded

TEST_FENS = [
    "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
    "rnbqkbnr/pppppppp/8/8/4P3/8/PPPP1PPP/RNBQKBNR b KQkq e3 0 1",
    "r1bqkbnr/pppp1ppp/2n5/4p3/4P3/5N2/PPPP1PPP/RNBQKB1R w KQkq - 2 3",
    "r1bqkb1r/pppp1ppp/2n2n2/4p2Q/2B1P3/8/PPPP1PPP/RNB1K1NR w KQkq - 4 4",
    "rnb1kbnr/pppp1ppp/8/4p3/6Pq/5P2/PPPPP2P/RNBQKBNR w KQkq - 1 3",
    "r3k2r/ppp2ppp/2n1bn2/3qp3/3P4/2N1BN2/PPP2PPP/R2QK2R b KQkq - 3 9",
    "8/5k2/8/3K4/8/8/4Q3/8 w - - 0 60",
    "6k1/5ppp/8/8/8/8/5PPP/3R2K1 w - - 0 30",
]


def fold(w, conv_prefix, bn_name):
    k = ref.find(w, conv_prefix)["kernel"].astype(np.float64)
    p = w[bn_name]
    s = p["gamma"] / np.sqrt(p["moving_variance"] + ref.BN_EPS)
    return k * s, p["beta"] - p["moving_mean"] * s


def build_layers(w):
    layers = [("input", *fold(w, "input_conv", "input_batchnorm"))]
    for i in range(1, 8):
        layers.append((f"res{i}a", *fold(w, f"res{i}_conv1", f"res{i}_batchnorm1")))
        layers.append((f"res{i}b", *fold(w, f"res{i}_conv2", f"res{i}_batchnorm2")))
    layers.append(("policy_conv", *fold(w, "policy_conv", "policy_batchnorm")))
    layers.append(("value_conv", *fold(w, "value_conv", "value_batchnorm")))
    # channels-first flatten index c*64+hw  ->  channels-last index hw*2+c
    pk = w["policy_out"]["kernel"]
    perm = np.array([c * 64 + hw for hw in range(64) for c in range(2)])
    layers.append(("policy_fc", pk[perm], w["policy_out"]["bias"]))
    layers.append(("value_fc1", w["value_dense"]["kernel"], w["value_dense"]["bias"]))
    layers.append(("value_fc2", w["value_out"]["kernel"], w["value_out"]["bias"]))
    return layers


def write_bin(layers, out_path):
    header, blobs, offset = [], [], 0

    def add(arr):
        nonlocal offset
        b = arr.tobytes()
        pad = (-len(b)) % 4
        start = offset
        blobs.append(b + b"\0" * pad)
        offset += len(b) + pad
        return start

    stored = []
    for name, k, b in layers:
        k16 = k.astype(np.float16)
        stored.append((name, k16.astype(np.float64), b))
        header.append({
            "name": name, "shape": list(k.shape),
            "w": add(k16), "bias": add(b.astype(np.float32)),
        })
    hj = json.dumps({"format": "AZC1-f16", "layers": header}).encode()
    hj += b" " * ((-len(hj)) % 4)
    data = b"AZC1" + struct.pack("<I", len(hj)) + hj + b"".join(blobs)
    # Parts are base64 text so that any static host / CDN serves them as-is.
    parts = []
    for i in range(0, len(data), PART_BYTES):
        name = f"{Path(out_path).name}.b64.{len(parts)}.txt"
        Path(out_path).with_name(name).write_text(base64.b64encode(data[i:i + PART_BYTES]).decode())
        parts.append(name)
    Path(f"{out_path}.json").write_text(json.dumps(
        {"parts": parts, "bytes": len(data), "encoding": "base64"}))
    return stored


def forward_folded(layers, planes):
    """Forward pass with folded (and float16-rounded) weights, NHWC flatten order."""
    L = {n: (k, b) for n, k, b in layers}
    relu = lambda a: np.maximum(a, 0)

    def conv(x, n):
        k, b = L[n]
        return ref._conv_same(x, k) + b[:, None, None]

    x = relu(conv(planes, "input"))
    for i in range(1, 8):
        x = relu(x + conv(relu(conv(x, f"res{i}a")), f"res{i}b"))
    p = relu(conv(x, "policy_conv")).transpose(1, 2, 0).reshape(-1)
    logits = p @ L["policy_fc"][0] + L["policy_fc"][1]
    e = np.exp(logits - logits.max())
    v = relu(conv(x, "value_conv")).transpose(1, 2, 0).reshape(-1)
    v = relu(v @ L["value_fc1"][0] + L["value_fc1"][1])
    return e / e.sum(), float(np.tanh(v @ L["value_fc2"][0] + L["value_fc2"][1])[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h5")
    ap.add_argument("out")
    ap.add_argument("--vectors")
    a = ap.parse_args()

    w = ref.load_h5(a.h5)
    layers = build_layers(w)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    deq = write_bin(layers, a.out)
    print(f"wrote {a.out}.json + parts")

    vectors = []
    for fen in TEST_FENS:
        planes = ref.canon_input_planes(fen)
        p_ref, v_ref = ref.forward_h5(w, planes)
        p_q, v_q = forward_folded(deq, planes)
        top = np.argsort(-p_ref)[:8]
        # policy is in the canonical (side-to-move = white) frame; map back to real moves
        black = fen.split(" ")[1] == "b"
        name = lambda i: ref.mirror_uci(ref.LABELS[i]) if black else ref.LABELS[i]
        print(f"{fen[:40]:40s} v={v_ref:+.3f} (f16 {v_q:+.3f})  "
              f"top={[name(i) for i in top[:4]]}  "
              f"max|dp|={np.abs(p_ref - p_q).max():.4f}")
        vectors.append({"fen": fen, "value": v_ref,
                        "top": [{"move": name(i), "p": float(p_ref[i])} for i in top]})
    if a.vectors:
        Path(a.vectors).parent.mkdir(parents=True, exist_ok=True)
        Path(a.vectors).write_text(json.dumps(vectors, indent=1))
        print(f"wrote {a.vectors}")


if __name__ == "__main__":
    main()
