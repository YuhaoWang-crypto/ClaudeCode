"""Export a checkpoint to the self-contained browser game.

    python -m alphago_zero.export_web --ckpt alphago_zero/runs/7x7/best.pt

BatchNorm is folded into the preceding convolution (inference only), weights
are packed as base64 float32, and engine.js + model are inlined into
web/template.html -> web/go_zero.html (open it in any browser, no server).
"""
import argparse
import base64
import json
import os

import numpy as np
import torch

from .model import load

HERE = os.path.dirname(os.path.abspath(__file__))


def fold(conv, bn):
    w = conv.weight.detach().numpy().astype(np.float64)
    g = bn.weight.detach().numpy()
    b = bn.bias.detach().numpy()
    m = bn.running_mean.numpy()
    v = bn.running_var.numpy()
    s = g / np.sqrt(v + bn.eps)
    return (w * s[:, None, None, None]).astype(np.float32), (b - m * s).astype(np.float32)


def export_model(net, meta):
    layers, chunks, off = [], [], 0

    def add(name, arr):
        nonlocal off
        a = np.ascontiguousarray(arr, dtype=np.float32).reshape(-1)
        layers.append(dict(name=name, offset=off, size=int(a.size)))
        chunks.append(a)
        off += a.size

    for nm, (w, b) in [("stem", fold(net.stem, net.stem_bn))]:
        add(nm + ".w", w); add(nm + ".b", b)
    for i, blk in enumerate(net.tower):
        for j, (c, bn) in enumerate(((blk.c1, blk.b1), (blk.c2, blk.b2)), 1):
            w, b = fold(c, bn)
            add(f"b{i}.c{j}.w", w); add(f"b{i}.c{j}.b", b)
    w, b = fold(net.p_conv, net.p_bn); add("p_conv.w", w); add("p_conv.b", b)
    add("p_fc.w", net.p_fc.weight.detach()); add("p_fc.b", net.p_fc.bias.detach())
    w, b = fold(net.v_conv, net.v_bn); add("v_conv.w", w); add("v_conv.b", b)
    add("v_fc1.w", net.v_fc1.weight.detach()); add("v_fc1.b", net.v_fc1.bias.detach())
    add("v_fc2.w", net.v_fc2.weight.detach()); add("v_fc2.b", net.v_fc2.bias.detach())
    data = np.concatenate(chunks).astype("<f4").tobytes()
    return dict(meta=meta, layers=layers,
                weights=base64.b64encode(data).decode("ascii"))


def load_metrics(run_dir):
    path = os.path.join(run_dir, "metrics.jsonl")
    if not os.path.exists(path):
        return []
    keep = ("iteration", "games", "elo", "generation", "best_vs_random",
            "best_vs_gen0", "loss_p", "loss_v", "mean_len", "elapsed_min")
    rows = [json.loads(l) for l in open(path)]
    return [{k: r.get(k) for k in keep} for r in rows]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default="alphago_zero/runs/7x7/best.pt")
    ap.add_argument("--out", default=os.path.join(HERE, "web", "go_zero.html"))
    ap.add_argument("--json", default=None, help="also write the model JSON here")
    ap.add_argument("--ladder", default=None, help="ladder.json from analyze.py")
    args = ap.parse_args()

    net, ck = load(args.ckpt)
    cfg = ck["cfg"]
    meta = dict(net.cfg, history=cfg["history"], komi=cfg["komi"],
                c_puct=cfg["c_puct"], generation=ck.get("generation"),
                elo=ck.get("elo"), games=ck.get("games"))
    model = export_model(net, meta)
    run_dir = os.path.dirname(args.ckpt)
    extra = dict(metrics=load_metrics(run_dir), config=cfg)
    lad = args.ladder or os.path.join(run_dir, "ladder.json")
    if os.path.exists(lad):
        extra["ladder"] = json.load(open(lad))
    if args.json:
        json.dump(model, open(args.json, "w"))
    tpl = open(os.path.join(HERE, "web", "template.html")).read()
    engine = open(os.path.join(HERE, "web", "engine.js")).read()
    html = (tpl.replace("/*__ENGINE__*/", engine)
               .replace("/*__MODEL__*/", "const MODEL = " + json.dumps(model) + ";")
               .replace("/*__RUN__*/", "const RUN = " + json.dumps(extra) + ";"))
    open(args.out, "w").write(html)
    print(f"wrote {args.out} ({len(html) / 1e6:.2f} MB), generation {meta['generation']}")


if __name__ == "__main__":
    main()
