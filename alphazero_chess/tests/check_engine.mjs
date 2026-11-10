// Node check: the JS network must reproduce the Python reference outputs
// (tests/vectors.json, produced from the original Keras weights), and the
// search must find simple mates.
//
//   npm install chess.js@1.4.0 @tensorflow/tfjs@4.22.0
//   node alphazero_chess/tests/check_engine.mjs
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import * as tf from "@tensorflow/tfjs";
// (swap in "@tensorflow/tfjs-node" for a much faster run if it is installed)
import { Chess } from "chess.js";
import { AlphaZeroNet, LABEL_INDEX, mirrorUci, search } from "../web/engine.js";

const here = path.dirname(fileURLToPath(import.meta.url));
const modelDir = path.join(here, "../web/model");
const manifest = JSON.parse(fs.readFileSync(path.join(modelDir, "weights.json")));
const raw = manifest.parts.map((p) => fs.readFileSync(path.join(modelDir, p)));
const buf = manifest.encoding === "base64"
  ? Buffer.concat(raw.map((b) => Buffer.from(b.toString().trim(), "base64")))
  : Buffer.concat(raw);
const net = AlphaZeroNet.fromBuffer(tf, buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length));

let failures = 0;
const check = (ok, msg) => { console.log((ok ? "PASS " : "FAIL ") + msg); if (!ok) failures++; };

// 1. network vs reference
const vectors = JSON.parse(fs.readFileSync(path.join(here, "vectors.json")));
const out = await net.predict(vectors.map((v) => v.fen));
vectors.forEach((v, i) => {
  const black = v.fen.split(" ")[1] === "b";
  const prob = (m) => out[i].policy[LABEL_INDEX.get(black ? mirrorUci(m) : m)];
  const dv = Math.abs(out[i].value - v.value);
  const dp = Math.max(...v.top.map((t) => Math.abs(prob(t.move) - t.p)));
  check(dv < 0.01 && dp < 0.01, `${v.fen.padEnd(62)} |dv|=${dv.toFixed(4)} max|dp|=${dp.toFixed(4)}`);
});

// 2. search finds mate in one
for (const [fen, mate] of [
  ["r1bqkb1r/pppp1ppp/2n2n2/4p2Q/2B1P3/8/PPPP1PPP/RNB1K1NR w KQkq - 4 4", "h5f7"],
  ["6k1/5ppp/8/8/8/8/5PPP/3R2K1 w - - 0 30", "d1d8"],
  ["3r2k1/8/8/8/8/8/5PPP/6K1 b - - 0 30", "d8d1"],
]) {
  const t = performance.now();
  const r = await search({ net, Chess, startFen: fen, sims: 64, batch: 16 });
  check(r.best.uci === mate,
    `mate search ${fen.slice(0, 30)}… -> ${r.best.uci} (n=${r.best.n}, q=${r.best.q.toFixed(2)}, ${Math.round(performance.now() - t)} ms)`);
}

// 3. search from a history (repetition handling) runs and returns a legal move
const r = await search({ net, Chess, startFen: new Chess().fen(), moves: ["e2e4", "e7e5", "g1f3"], sims: 32 });
check(r.stats.length === new Chess("rnbqkbnr/pppp1ppp/8/4p3/4P3/5N2/PPPP1PPP/RNBQKB1R b KQkq - 1 2").moves().length,
  `search after 1.e4 e5 2.Nf3 -> ${r.best.san}, ${r.stats.length} legal moves`);

console.log(failures ? `\n${failures} failure(s)` : "\nall checks passed");
process.exit(failures ? 1 : 0);
