// AlphaZero-style chess engine for the browser.
//
// A JavaScript port of the method used in zeta36/chess-alpha-zero:
//   1. Board encoding  - 18 planes of 8x8, always seen from the side to move
//   2. Neural network  - 7 residual blocks (256 filters) with a policy head
//                        (1968 move probabilities) and a value head (-1..+1)
//   3. Search          - PUCT Monte-Carlo tree search guided by the network,
//                        leaves evaluated in batches using "virtual loss"
// The network weights are the pretrained ones shipped with that repository.
//
// The TensorFlow.js namespace and the chess.js Chess class are passed in,
// so this file runs both in the browser and in Node tests.

// ------------------------------------------------------------ move labels
function createLabels() {
  const L = "abcdefgh", N = "12345678", out = [];
  for (let l1 = 0; l1 < 8; l1++) {
    for (let n1 = 0; n1 < 8; n1++) {
      const dest = [];
      for (let t = 0; t < 8; t++) dest.push([t, n1]);
      for (let t = 0; t < 8; t++) dest.push([l1, t]);
      for (let t = -7; t < 8; t++) dest.push([l1 + t, n1 + t]);
      for (let t = -7; t < 8; t++) dest.push([l1 + t, n1 - t]);
      for (const [a, b] of [[-2, -1], [-1, -2], [-2, 1], [1, -2], [2, -1], [-1, 2], [2, 1], [1, 2]])
        dest.push([l1 + a, n1 + b]);
      for (const [l2, n2] of dest) {
        if ((l1 !== l2 || n1 !== n2) && l2 >= 0 && l2 < 8 && n2 >= 0 && n2 < 8)
          out.push(L[l1] + N[n1] + L[l2] + N[n2]);
      }
    }
  }
  for (let l1 = 0; l1 < 8; l1++) {
    const l = L[l1];
    for (const p of "qrbn") {
      out.push(l + "2" + l + "1" + p, l + "7" + l + "8" + p);
      if (l1 > 0) out.push(l + "2" + L[l1 - 1] + "1" + p, l + "7" + L[l1 - 1] + "8" + p);
      if (l1 < 7) out.push(l + "2" + L[l1 + 1] + "1" + p, l + "7" + L[l1 + 1] + "8" + p);
    }
  }
  return out;
}

export const LABELS = createLabels();
export const LABEL_INDEX = new Map(LABELS.map((m, i) => [m, i]));
export const mirrorUci = (m) => m.replace(/[1-8]/g, (d) => String(9 - Number(d)));

// ------------------------------------------------------------ input planes
const PIECES = "KQRBNPkqrbnp";
const swapCase = (s) => s.replace(/[a-zA-Z]/g, (c) => (c === c.toUpperCase() ? c.toLowerCase() : c.toUpperCase()));

// Same transform as the original maybe_flip_fen (note: en-passant square is kept as-is).
export function canonicalFen(fen) {
  const f = fen.split(" ");
  if (f[1] !== "b") return fen;
  const rows = f[0].split("/").reverse().map(swapCase).join("/");
  const castling = swapCase(f[2]).split("").sort().join("");
  return [rows, "w", castling, f[3], f[4], f[5]].join(" ");
}

// Float32Array in NHWC order: index = (row*8 + file)*18 + plane, row 0 = rank 8.
export function inputPlanes(fen, out = new Float32Array(8 * 8 * 18), offset = 0) {
  const f = canonicalFen(fen).split(" ");
  out.fill(0, offset, offset + 8 * 8 * 18);
  let sq = 0;
  for (const c of f[0]) {
    if (c === "/") continue;
    if (c >= "1" && c <= "8") { sq += Number(c); continue; }
    out[offset + sq * 18 + PIECES.indexOf(c)] = 1;
    sq++;
  }
  const clock = Number(f[4]);
  let ep = -1;
  if (f[3] !== "-") ep = (8 - Number(f[3][1])) * 8 + (f[3].charCodeAt(0) - 97);
  for (let s = 0; s < 64; s++) {
    const o = offset + s * 18;
    out[o + 12] = f[2].includes("K") ? 1 : 0;
    out[o + 13] = f[2].includes("Q") ? 1 : 0;
    out[o + 14] = f[2].includes("k") ? 1 : 0;
    out[o + 15] = f[2].includes("q") ? 1 : 0;
    out[o + 16] = clock;
    out[o + 17] = s === ep ? 1 : 0;
  }
  return out;
}

// ------------------------------------------------------------ network
function halfToFloat(h) {
  const s = (h & 0x8000) ? -1 : 1, e = (h >> 10) & 0x1f, m = h & 0x3ff;
  if (e === 0) return s * m * 2 ** -24;
  if (e === 31) return m ? NaN : s * Infinity;
  return s * (1 + m / 1024) * 2 ** (e - 15);
}

export class AlphaZeroNet {
  constructor(tf, layers) {
    this.tf = tf;
    this.L = layers;
    this.padBatches = false;
  }

  // buffer: the concatenated weight file (ArrayBuffer)
  static fromBuffer(tf, buffer) {
    const dv = new DataView(buffer);
    const magic = String.fromCharCode(...new Uint8Array(buffer, 0, 4));
    if (magic !== "AZC1") throw new Error("bad weight file");
    const hlen = dv.getUint32(4, true);
    const header = JSON.parse(new TextDecoder().decode(new Uint8Array(buffer, 8, hlen)));
    const base = 8 + hlen;
    const table = new Float32Array(65536);
    for (let i = 0; i < 65536; i++) table[i] = halfToFloat(i);
    const layers = {};
    for (const l of header.layers) {
      const size = l.shape.reduce((a, b) => a * b, 1);
      const raw = new Uint16Array(buffer, base + l.w, size);
      const w = new Float32Array(size);
      for (let i = 0; i < size; i++) w[i] = table[raw[i]];
      const nout = l.shape[l.shape.length - 1];
      const bias = new Float32Array(buffer.slice(base + l.bias, base + l.bias + nout * 4));
      layers[l.name] = { w: tf.tensor(w, l.shape), b: tf.tensor1d(bias) };
    }
    return new AlphaZeroNet(tf, layers);
  }

  static async load(tf, baseUrl, onProgress) {
    const manifest = await (await fetch(new URL("weights.json", baseUrl))).json();
    const buf = new Uint8Array(manifest.bytes);
    let off = 0;
    for (const part of manifest.parts) {
      const res = await fetch(new URL(part, baseUrl));
      if (!res.ok) throw new Error(`failed to load ${part}: ${res.status}`);
      if (manifest.encoding === "base64") {
        // for hosts that only serve text files
        const bin = atob((await res.text()).trim());
        for (let i = 0; i < bin.length; i++) buf[off + i] = bin.charCodeAt(i);
        off += bin.length;
        if (onProgress) onProgress(off / manifest.bytes);
        continue;
      }
      const reader = res.body.getReader();
      for (;;) {
        const { done, value } = await reader.read();
        if (done) break;
        buf.set(value, off);
        off += value.length;
        if (onProgress) onProgress(off / manifest.bytes);
      }
    }
    return AlphaZeroNet.fromBuffer(tf, buf.buffer);
  }

  // fens: array of FEN strings. Returns [{policy: Float32Array(1968), value}] in the
  // canonical frame (side to move plays "up the board" as White).
  // With padBatches on, batches are rounded up to a power of two so WebGL only
  // compiles shaders for a handful of shapes.
  async predict(fens) {
    const tf = this.tf, L = this.L;
    const n = this.padBatches ? 2 ** Math.ceil(Math.log2(fens.length)) : fens.length;
    const data = new Float32Array(n * 8 * 8 * 18);
    fens.forEach((fen, i) => inputPlanes(fen, data, i * 8 * 8 * 18));
    const [pT, vT] = tf.tidy(() => {
      const conv = (x, name, activation) =>
        tf.fused.conv2d({ x, filter: L[name].w, strides: 1, pad: "same", bias: L[name].b, activation });
      let x = conv(tf.tensor4d(data, [n, 8, 8, 18]), "input", "relu");
      for (let i = 1; i <= 7; i++) {
        const y = conv(conv(x, `res${i}a`, "relu"), `res${i}b`, "linear");
        x = tf.relu(tf.add(x, y));
      }
      const dense = (a, name, activation) =>
        tf.fused.matMul({ a, b: L[name].w, bias: L[name].b, activation });
      const p = tf.softmax(dense(conv(x, "policy_conv", "relu").reshape([n, 128]), "policy_fc", "linear"));
      const h = dense(conv(x, "value_conv", "relu").reshape([n, 64]), "value_fc1", "relu");
      const v = tf.tanh(dense(h, "value_fc2", "linear"));
      return [p, v];
    });
    const [pd, vd] = await Promise.all([pT.data(), vT.data()]);
    pT.dispose(); vT.dispose();
    return fens.map((_, i) => ({ policy: pd.subarray(i * 1968, (i + 1) * 1968), value: vd[i] }));
  }
}

// ------------------------------------------------------------ MCTS
const uciOf = (m) => m.from + m.to + (m.promotion || "");

class Node {
  constructor(move, prior) {
    this.move = move;     // chess.js move object leading here (null at root)
    this.P = prior;       // prior probability from the policy head
    this.N = 0;           // visit count
    this.W = 0;           // total value, from the view of the player who made `move`
    this.children = null; // null = not expanded yet
    this.pending = false; // waiting for a network evaluation
  }
  get Q() { return this.N ? this.W / this.N : 0; }
}

// Prior over the legal moves of `game`, renormalised (as in the original
// select_action_q_and_u). Policy is canonical: mirror moves when Black is to move.
function priorsFor(policy, legal, blackToMove) {
  const pri = legal.map((m) => {
    const u = uciOf(m);
    const idx = LABEL_INDEX.get(blackToMove ? mirrorUci(u) : u);
    return idx === undefined ? 0 : policy[idx];
  });
  const tot = pri.reduce((a, b) => a + b, 1e-8);
  return pri.map((p) => p / tot);
}

function terminalValue(game) {
  // value from the view of the side to move
  if (game.isCheckmate()) return -1;
  if (game.isDraw() || game.isStalemate()) return 0;
  return null;
}

/**
 * Run a search and return the chosen move plus statistics.
 * opts: { net, Chess, startFen, moves (array of uci), sims, batch, cpuct,
 *         onProgress(info), shouldStop() }
 */
export async function search(opts) {
  const { net, Chess, startFen, moves = [], sims = 200, batch = 8, cpuct = 1.0,
          virtualLoss = 1, onProgress, shouldStop } = opts;
  const game = new Chess(startFen);
  for (const u of moves) game.move({ from: u.slice(0, 2), to: u.slice(2, 4), promotion: u[4] });
  if (terminalValue(game) !== null) return null;

  const root = new Node(null, 1);
  const expand = (node, policy) => {
    const legal = game.moves({ verbose: true });
    const pri = priorsFor(policy, legal, game.turn() === "b");
    node.children = legal.map((m, i) => new Node(m, pri[i]));
  };

  const [rootEval] = await net.predict([game.fen()]);
  expand(root, rootEval.policy);
  root.N = 1;
  const rawValue = rootEval.value;

  const backup = (path, v) => {
    // v is from the view of the side to move at the leaf
    for (let i = path.length - 1; i >= 0; i--) {
      v = -v;
      path[i].N += 1;
      path[i].W += v;
    }
    root.N += 1;
  };

  let done = 1;
  const t0 = performance.now();
  while (done < sims) {
    if (shouldStop && shouldStop()) break;
    const leaves = [];
    const want = Math.min(batch, sims - done);
    let attempts = 0;
    while (leaves.length < want && attempts < want * 2) {
      attempts++;
      // ---- select
      const path = [];
      let node = root;
      while (node.children && node.children.length) {
        const sq = Math.sqrt(node.N + 1);
        let best = null, bestS = -Infinity;
        for (const c of node.children) {
          const s = c.Q + cpuct * c.P * sq / (1 + c.N);
          if (s > bestS) { bestS = s; best = c; }
        }
        node = best;
        path.push(node);
        game.move(node.move);
        if (node.pending) break;
      }
      const tv = terminalValue(game);
      if (node.pending) {
        // already queued in this batch; stop collecting to avoid duplicate work
        for (let i = 0; i < path.length; i++) game.undo();
        break;
      }
      if (tv !== null) {
        for (let i = 0; i < path.length; i++) game.undo();
        backup(path, tv);
        done++;
        continue;
      }
      // ---- virtual loss: make other threads in this batch look elsewhere
      for (const p of path) { p.N += virtualLoss; p.W -= virtualLoss; }
      node.pending = true;
      leaves.push({ path, node, fen: game.fen(), legal: game.moves({ verbose: true }), black: game.turn() === "b" });
      for (let i = 0; i < path.length; i++) game.undo();
    }
    if (!leaves.length) continue; // only terminal leaves this round (already counted)

    // ---- evaluate + expand + backup
    const evals = await net.predict(leaves.map((l) => l.fen));
    leaves.forEach((leaf, i) => {
      const pri = priorsFor(evals[i].policy, leaf.legal, leaf.black);
      leaf.node.children = leaf.legal.map((m, k) => new Node(m, pri[k]));
      leaf.node.pending = false;
      for (const p of leaf.path) { p.N -= virtualLoss; p.W += virtualLoss; }
      backup(leaf.path, evals[i].value);
      done++;
    });
    if (onProgress) onProgress(summary());
  }

  function summary() {
    const stats = root.children
      .map((c) => ({ uci: uciOf(c.move), san: c.move.san, n: c.N, q: c.Q, p: c.P }))
      .sort((a, b) => b.n - a.n || b.p - a.p);
    const best = stats[0];
    return {
      stats, best, sims: done, total: sims,
      value: best && best.n ? best.q : rawValue, rawValue,
      ms: performance.now() - t0, turn: game.turn(),
    };
  }
  return summary();
}
