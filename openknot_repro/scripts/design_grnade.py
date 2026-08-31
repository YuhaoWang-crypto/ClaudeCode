"""gRNAde 2D-mode inverse folding for an OpenKnot target + in-silico OpenKnot screening.

Mirrors projects/openknot_benchmark/design.ipynb (MODE='2d'), which needs no PDB input.
"""
import argparse, os, random, sys, time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

sys.path.insert(0, "/home/user/work/shim"); sys.path.append("/home/user/work/gRNAde")
os.environ.setdefault("PROJECT_PATH", "/home/user/work/gRNAde/")
os.environ.setdefault("DATA_PATH", "/home/user/work/gRNAde/data/")
os.chdir("/home/user/work/gRNAde")

from src.data.featurizer import RNAGraphFeaturizer
from src.models import gRNAde
from src.evaluator import openknot_score_ribonanzanet, self_consistency_score_ribonanzanet_sec_struct
from src.constants import NUM_TO_LETTER, FILL_VALUE
from tools.ribonanzanet.network import RibonanzaNet
from tools.ribonanzanet_sec_struct.network import RibonanzaNetSS

ap = argparse.ArgumentParser()
ap.add_argument("--puzzle", default="P20")
ap.add_argument("--set", default="a", choices=["a", "b"])
ap.add_argument("--total", type=int, default=2048)
ap.add_argument("--batch", type=int, default=64)
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--seed", type=int, default=42)
ap.add_argument("--out", default="/home/user/work/grnade_designs.csv")
args = ap.parse_args()

torch.set_num_threads(args.threads)
random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)

meta = pd.read_csv(f"projects/openknot_benchmark/metadata_7{args.set}.csv")
row = meta[meta["Puzzle"] == args.puzzle].iloc[0]
native_seq, target_ss = row["Sequence"], row["Dot-bracket"]
print(f"{args.puzzle} ({row['Title']}) L={len(target_ss)}")
print(f"target: {target_ss}")

featurizer = RNAGraphFeaturizer(split="test_2d", radius=0.0, top_k=32, num_rbf=32,
                                num_posenc=32, max_num_conformers=1, noise_scale=0.1,
                                drop_prob_3d=0.5)
model = gRNAde(node_in_dim=(15, 4), node_h_dim=(128, 16), edge_in_dim=(132, 3),
               edge_h_dim=(64, 4), num_layers=4, drop_rate=0.5, out_dim=4)
model.load_state_dict(torch.load("checkpoints/gRNAde_drop3d@0.75_maxlen@500.h5",
                                 map_location="cpu"))
model.eval()

rnet = RibonanzaNet("tools/ribonanzanet/config.yaml",
                    "checkpoints/ribonanzanet/ribonanzanet.pt", "cpu").eval()
rnet_ss = RibonanzaNetSS("tools/ribonanzanet_sec_struct/config.yaml",
                         "checkpoints/ribonanzanet_sec_struct/ribonanzanet_ss.pt", "cpu").eval()

data = featurizer({"sequence": native_seq,
                   "coords_list": [torch.ones(len(target_ss), 3, 3) * FILL_VALUE],
                   "sec_struct_list": [target_ss]})

recs, t0 = [], time.time()
with torch.no_grad():
    while len(recs) < args.total:
        n = min(args.batch, args.total - len(recs))
        temperature = float(np.random.uniform(0.1, 1.0))
        samples, logits = model.sample(data, n, temperature, return_logits=True)
        nn_ = logits.shape[1]
        ppl = torch.exp(F.cross_entropy(logits.view(n * nn_, model.out_dim),
                                        samples.view(n * nn_).long(), reduction="none")
                        .view(n, nn_).mean(dim=1)).cpu().numpy()
        s = samples.cpu().numpy()
        mask = data.mask_seq.cpu().numpy()
        oks = openknot_score_ribonanzanet(s, target_ss, mask, rnet)
        sc = self_consistency_score_ribonanzanet_sec_struct(s, target_ss, mask, rnet_ss)
        for i in range(n):
            recs.append(dict(sequence="".join(NUM_TO_LETTER[c] for c in s[i][mask]),
                             temperature=round(temperature, 3),
                             perplexity=float(ppl[i]),
                             openknot_score=float(oks[i]),
                             sc_score_rnet_ss=float(sc[i])))
        el = time.time() - t0
        print(f"\r{len(recs)}/{args.total}  {el:.0f}s  eta {el/len(recs)*(args.total-len(recs))/60:.1f}min",
              end="", flush=True)
print()

df = pd.DataFrame(recs).sort_values("openknot_score", ascending=False)
df.insert(0, "puzzle", args.puzzle)
df.to_csv(args.out, index=False)
print(df.head(10).to_string(index=False))
print(f"\nin-silico OpenKnot score: mean {df.openknot_score.mean():.1f}  "
      f"max {df.openknot_score.max():.1f}  frac>90 {(df.openknot_score > 90).mean():.3f}")
print(f"saved -> {args.out}")
