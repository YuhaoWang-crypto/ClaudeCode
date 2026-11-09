"""Check that web/engine.js reproduces the Python rules and network exactly.

    python -m alphago_zero.web.check_parity   (needs node)
"""
import json, os, subprocess, sys, tempfile
import numpy as np, torch
from ..go import GoState
from ..model import load
from ..export_web import export_model
from ..selfplay import random_move

HERE = os.path.dirname(os.path.abspath(__file__))

def main(ckpt="alphago_zero/runs/7x7/best.pt", n_games=30):
    torch.set_num_threads(1)
    net, ck = load(ckpt)
    cfg = ck["cfg"]
    rng = np.random.default_rng(1)
    cases = []
    for g in range(n_games):
        s = GoState.new(cfg["board_size"], cfg["komi"])
        moves = []
        while not s.is_terminal():
            m = random_move(s, rng) if rng.random() > 0.03 else s.n * s.n
            moves.append(m)
            s = s.play(m)
            if rng.random() < 0.15:
                x = torch.from_numpy(s.features(cfg["history"])[None])
                with torch.no_grad():
                    lg, v = net(x)
                cases.append(dict(moves=list(moves), legal=s.legal_mask().astype(int).tolist(),
                                  score=s.score(), logits=lg[0].tolist(), value=float(v[0])))
    with tempfile.TemporaryDirectory() as d:
        mpath, cpath = os.path.join(d, "m.json"), os.path.join(d, "c.json")
        json.dump(export_model(net, dict(net.cfg, history=cfg["history"])), open(mpath, "w"))
        json.dump(dict(komi=cfg["komi"], cases=cases), open(cpath, "w"))
        js = f"""
const G = require({json.dumps(os.path.join(HERE, 'engine.js'))});
const fs = require('fs');
const model = JSON.parse(fs.readFileSync({json.dumps(mpath)}));
const data = JSON.parse(fs.readFileSync({json.dumps(cpath)}));
const net = new G.Net(model);
let maxLogit = 0, maxV = 0, legalBad = 0, scoreBad = 0, symMax = 0;
const t0 = Date.now(); let evals = 0;
for (const c of data.cases) {{
  let s = G.State.initial(net.n, data.komi);
  for (const m of c.moves) s = s.play(m);
  const lm = s.legalMask();
  for (let i = 0; i < lm.length; i++) if (lm[i] !== c.legal[i]) legalBad++;
  if (Math.abs(s.score() - c.score) > 1e-9) scoreBad++;
  const out = net.forward(s.features(net.history)); evals++;
  for (let i = 0; i < out.logits.length; i++) maxLogit = Math.max(maxLogit, Math.abs(out.logits[i] - c.logits[i]));
  maxV = Math.max(maxV, Math.abs(out.value - c.value));
  // symmetry consistency: value should be (nearly) invariant only for a trained net, so just check it runs
  for (let k = 0; k < 8; k++) {{ net.evaluate(s, k); evals++; }}
}}
console.log(JSON.stringify({{cases: data.cases.length, maxLogitErr: maxLogit, maxValueErr: maxV, legalMismatch: legalBad, scoreMismatch: scoreBad,
  msPerEval: (Date.now() - t0) / evals}}));
"""
        out = subprocess.run(["node", "-e", js], capture_output=True, text=True)
        print(out.stdout, out.stderr)
        r = json.loads(out.stdout)
        ok = r["maxLogitErr"] < 1e-3 and r["maxValueErr"] < 1e-4 and not r["legalMismatch"] and not r["scoreMismatch"]
        print("PARITY OK" if ok else "PARITY FAILED")
        return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))
