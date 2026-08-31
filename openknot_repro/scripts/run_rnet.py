"""Run RibonanzaNet on OpenKnotBench designs -> predicted 2A3 -> in-silico OpenKnot score.

Reproduces the in-silico screening step of Townley et al. Science 2026
(gRNAde src/evaluator.py::openknot_score_ribonanzanet).

Predictions are checkpointed to disk during inference, and scoring runs as a
separate resumable pass, so a crash in scoring never loses the (expensive) NN pass.
"""
import argparse, os, pickle, sys, time
from collections import defaultdict

import numpy as np
import torch

sys.path.insert(0, "/home/user/work/shim")
sys.path.append("/home/user/work/gRNAde")
os.environ.setdefault("PROJECT_PATH", "/home/user/work/gRNAde/")
os.chdir("/home/user/work/gRNAde")

from tools.ribonanzanet.network import RibonanzaNet
from src.openknot_score import get_openknot_score

ap = argparse.ArgumentParser()
ap.add_argument("--rounds", type=int, nargs="+", default=[1, 3])
ap.add_argument("--batch", type=int, default=16)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--out", default="/home/user/work/rnet_preds.pkl")
ap.add_argument("--ckpt", default=None, help="raw-prediction checkpoint (default: <out>.raw)")
args = ap.parse_args()
ckpt_path = args.ckpt or args.out + ".raw"

torch.set_num_threads(args.threads)

rows = pickle.load(open("/home/user/work/okb_designs.pkl", "rb"))
work = [r for r in rows if r["round"] in args.rounds]
print(f"{len(work)} designs, rounds {args.rounds}", flush=True)

# ---------------- pass 1: RibonanzaNet inference (checkpointed) ----------------
pred = [None] * len(work)
if os.path.exists(ckpt_path):
    saved = pickle.load(open(ckpt_path, "rb"))
    if saved["ids"] == [r["id"] for r in work]:
        pred = saved["pred"]
        print(f"resumed {sum(p is not None for p in pred)}/{len(work)} from {ckpt_path}", flush=True)

todo = [i for i, p in enumerate(pred) if p is None]
if todo:
    model = RibonanzaNet("tools/ribonanzanet/config.yaml",
                         "checkpoints/ribonanzanet/ribonanzanet.pt", "cpu").eval()
    by_len = defaultdict(list)
    for i in todo:
        by_len[len(work[i]["seq"])].append(i)

    t0, done, last_ckpt = time.time(), 0, 0
    with torch.no_grad():
        for L, idxs in sorted(by_len.items()):
            for b in range(0, len(idxs), args.batch):
                chunk = idxs[b:b + args.batch]
                out = model.predict([work[i]["seq"] for i in chunk]).numpy()[:, :, 0]  # ch0 = 2A3
                for k, i in enumerate(chunk):
                    pred[i] = out[k].astype(np.float32)
                done += len(chunk)
                if done - last_ckpt >= 1000:
                    pickle.dump({"ids": [r["id"] for r in work], "pred": pred},
                                open(ckpt_path, "wb"))
                    last_ckpt = done
                    el = time.time() - t0
                    print(f"\r{done}/{len(todo)}  {el:.0f}s  "
                          f"eta {el / done * (len(todo) - done) / 60:.1f}min", end="", flush=True)
    pickle.dump({"ids": [r["id"] for r in work], "pred": pred}, open(ckpt_path, "wb"))
    print(f"\ninference done, checkpoint saved -> {ckpt_path}", flush=True)

# ---------------- pass 2: OpenKnot scoring (cheap, fault-tolerant) ----------------
results, skipped = [], []
for r, p in zip(work, pred):
    try:
        _, _, _, ok_pred = get_openknot_score(r["struct"], p[None, :])
        exp = np.nan_to_num(r["exp_react"], nan=0.0)[None, :]
        _, _, _, ok_exp = get_openknot_score(r["struct"], exp)
    except Exception as e:  # e.g. malformed target dot-bracket in the published data
        skipped.append((r["id"], r["puzzle"], repr(e)[:90]))
        continue
    results.append(dict(
        id=r["id"], round=r["round"], puzzle=r["puzzle"], method=r["method"],
        sn=r["sn"], sn_filter=r["sn_filter"], reads=r["reads"],
        rnet_f1=r["rnet_f1"], rnet_f1_cp=r["rnet_f1_cp"],
        oks_exp_published=r["oks"], oks_exp_recomputed=float(ok_exp[0]),
        oks_insilico=float(ok_pred[0]),
        pred_react=p, exp_react=r["exp_react"], seq=r["seq"], struct=r["struct"],
    ))

pickle.dump(results, open(args.out, "wb"))
print(f"saved {len(results)} -> {args.out}")
if skipped:
    print(f"skipped {len(skipped)} designs with unscoreable target structures:")
    for s in skipped[:5]:
        print("   ", s)
