// Engine-vs-engine match to measure whether a search setting makes play stronger.
// Each opening is played twice with colours swapped, so neither side gets an
// easier set of positions.
//
//   npm install chess.js@1.4.0 @tensorflow/tfjs-node@4.22.0
//   node alphazero_chess/tests/match.mjs '{"solver":true}' '{}' 10 100
//   args: engine A options, engine B options, number of openings, simulations per move
import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { Chess } from "chess.js";
import { AlphaZeroNet, search } from "../web/engine.js";

let tf;
try { tf = await import("@tensorflow/tfjs-node"); } catch { tf = await import("@tensorflow/tfjs"); }

const here = path.dirname(fileURLToPath(import.meta.url));
const modelDir = path.join(here, "../web/model");
const manifest = JSON.parse(fs.readFileSync(path.join(modelDir, "weights.json")));
const raw = manifest.parts.map((p) => fs.readFileSync(path.join(modelDir, p)));
const buf = manifest.encoding === "base64"
  ? Buffer.concat(raw.map((b) => Buffer.from(b.toString().trim(), "base64")))
  : Buffer.concat(raw);
const net = AlphaZeroNet.fromBuffer(tf, buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length));

const A = JSON.parse(process.argv[2] || "{}");
const B = JSON.parse(process.argv[3] || "{}");
const OPENINGS = Number(process.argv[4] || 10);
const SIMS = Number(process.argv[5] || 100);
const START = new Chess().fen();
const VAL = { p: 1, n: 3, b: 3, r: 5, q: 9, k: 0 };

// deterministic random openings: 6 plies sampled from the network's own policy
let seed = 12345;
const rand = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
async function makeOpening() {
  const hist = [];
  for (let ply = 0; ply < 6; ply++) {
    const r = await search({ net, Chess, startFen: START, moves: hist, sims: 1 });
    const total = r.stats.reduce((s, x) => s + x.p, 0);
    let x = rand() * total, pick = r.stats[0];
    for (const m of r.stats) { x -= m.p; if (x <= 0) { pick = m; break; } }
    hist.push(pick.uci);
  }
  return hist;
}

function material(game) {
  let s = 0;
  for (const row of game.board()) for (const p of row) if (p) s += (p.color === "w" ? 1 : -1) * VAL[p.type];
  return s;
}

async function playGame(opening, white, black) {
  const game = new Chess();
  const hist = [];
  for (const u of opening) { game.move({ from: u.slice(0, 2), to: u.slice(2, 4), promotion: u[4] }); hist.push(u); }
  const state = { w: { cfg: white, tree: null, cache: new Map() }, b: { cfg: black, tree: null, cache: new Map() } };
  while (!game.isGameOver() && hist.length < 200) {
    const s = state[game.turn()];
    const { reuse, ...cfg } = s.cfg;
    const r = await search({ net, Chess, startFen: START, moves: hist, sims: SIMS, ...cfg,
                             tree: reuse ? s.tree : null, cache: reuse ? s.cache : null });
    s.tree = r.tree;
    const u = r.best.uci;
    game.move({ from: u.slice(0, 2), to: u.slice(2, 4), promotion: u[4] });
    hist.push(u);
  }
  if (game.isCheckmate()) return game.turn() === "w" ? -1 : 1;
  if (game.isGameOver()) return 0;
  const m = material(game);            // adjudicate long games on material
  return m >= 3 ? 1 : m <= -3 ? -1 : 0;
}

let score = 0, games = 0, w = 0, d = 0, l = 0;
const t0 = Date.now();
for (let i = 0; i < OPENINGS; i++) {
  const op = await makeOpening();
  for (const aWhite of [true, false]) {
    const res = await playGame(op, aWhite ? A : B, aWhite ? B : A);
    const forA = aWhite ? res : -res;
    games++; score += (forA + 1) / 2;
    if (forA > 0) w++; else if (forA < 0) l++; else d++;
  }
  const s = score / games;
  const elo = s <= 0 ? -999 : s >= 1 ? 999 : Math.round(-400 * Math.log10(1 / s - 1));
  console.log(`after ${games} games: A +${w} =${d} -${l}  score ${(100 * s).toFixed(0)}%  ~${elo} Elo  (${Math.round((Date.now() - t0) / 1000)} s)`);
}
