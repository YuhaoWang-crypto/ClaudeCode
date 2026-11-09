// AlphaGo Zero mini — browser engine (no DOM).
// A line-by-line port of go.py (rules), model.py (forward pass, BN folded
// into the convolutions at export time) and mcts.py (PUCT search).
const GoZero = (() => {
  "use strict";
  const BLACK = 1, WHITE = -1, EMPTY = 0;
  const COLS = "ABCDEFGHJKLMNOPQRST";

  // ------------------------------------------------------------ geometry
  const geomCache = {};
  function geometry(n) {
    if (geomCache[n]) return geomCache[n];
    const nbrs = [], diags = [];
    for (let p = 0; p < n * n; p++) {
      const r = (p / n) | 0, c = p % n, a = [], d = [];
      for (const [dr, dc] of [[-1, 0], [1, 0], [0, -1], [0, 1]])
        if (r + dr >= 0 && r + dr < n && c + dc >= 0 && c + dc < n) a.push((r + dr) * n + c + dc);
      for (const [dr, dc] of [[-1, -1], [-1, 1], [1, -1], [1, 1]])
        if (r + dr >= 0 && r + dr < n && c + dc >= 0 && c + dc < n) d.push((r + dr) * n + c + dc);
      nbrs.push(a); diags.push(d);
    }
    // Zobrist keys: two 32-bit halves combined into one 53-bit number
    let seed = 20171019 >>> 0;
    const rnd = () => { seed ^= seed << 13; seed >>>= 0; seed ^= seed >>> 17; seed ^= seed << 5; return seed >>> 0; };
    const z = {};
    for (const col of [BLACK, WHITE]) {
      z[col] = { lo: new Uint32Array(n * n), hi: new Uint32Array(n * n) };
      for (let p = 0; p < n * n; p++) { z[col].lo[p] = rnd(); z[col].hi[p] = rnd() & 0x1fffff; }
    }
    // 8 board symmetries as permutations: t[p] = src[perm[p]]
    const perms = [];
    for (let k = 0; k < 8; k++) {
      const perm = new Int32Array(n * n);
      for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) {
        let rr = r, cc = c;
        if (k >= 4) cc = n - 1 - cc;
        for (let i = 0; i < k % 4; i++) { const t = rr; rr = cc; cc = n - 1 - t; }
        perm[r * n + c] = rr * n + cc;
      }
      perms.push(perm);
    }
    return (geomCache[n] = { nbrs, diags, z, perms });
  }

  // --------------------------------------------------------------- state
  class State {
    constructor(n, komi, board, toPlay, passes, lo, hi, prev, moveNum, lastMove) {
      Object.assign(this, { n, komi, board, toPlay, passes, lo, hi, prev, moveNum, lastMove });
      this.maxMoves = 2 * n * n;
      this._groups = null; this._legal = null; this._seen = null;
    }
    static initial(n, komi) {
      return new State(n, komi, new Int8Array(n * n), BLACK, 0, 0, 0, null, 0, null);
    }
    get key() { return this.lo + this.hi * 4294967296; }
    seen() {
      if (!this._seen) {
        this._seen = new Set();
        for (let s = this; s; s = s.prev) this._seen.add(s.key);
      }
      return this._seen;
    }
    groups() {
      if (this._groups) return this._groups;
      const n2 = this.n * this.n, { nbrs } = geometry(this.n), b = this.board;
      const gid = new Int16Array(n2).fill(-1), libs = [], stones = [];
      const mark = new Int32Array(n2).fill(-1);
      for (let p = 0; p < n2; p++) {
        if (b[p] === EMPTY || gid[p] >= 0) continue;
        const col = b[p], g = libs.length, st = [p], stack = [p];
        let nl = 0; gid[p] = g;
        while (stack.length) {
          const q = stack.pop();
          for (const r of nbrs[q]) {
            const v = b[r];
            if (v === EMPTY) { if (mark[r] !== g) { mark[r] = g; nl++; } }
            else if (v === col && gid[r] < 0) { gid[r] = g; st.push(r); stack.push(r); }
          }
        }
        libs.push(nl); stones.push(st);
      }
      return (this._groups = { gid, libs, stones });
    }
    analyse(p) {
      if (this.board[p] !== EMPTY) return null;
      const { nbrs, z } = geometry(this.n), { gid, libs, stones } = this.groups();
      const c = this.toPlay, b = this.board;
      let hasLib = false; const cap = new Set();
      for (const q of nbrs[p]) {
        const v = b[q];
        if (v === EMPTY) hasLib = true;
        else if (v === c) { if (libs[gid[q]] > 1) hasLib = true; }
        else if (libs[gid[q]] === 1) cap.add(gid[q]);
      }
      if (!hasLib && cap.size === 0) return null;            // suicide
      let lo = this.lo ^ z[c].lo[p], hi = this.hi ^ z[c].hi[p];
      const captured = [];
      for (const g of cap) for (const s of stones[g]) { lo ^= z[-c].lo[s]; hi ^= z[-c].hi[s]; captured.push(s); }
      lo >>>= 0; hi >>>= 0;
      if (this.seen().has(lo + hi * 4294967296)) return null; // superko
      return { captured, lo, hi };
    }
    legalMask() {
      if (!this._legal) {
        const n2 = this.n * this.n, m = new Uint8Array(n2 + 1);
        m[n2] = 1;
        for (let p = 0; p < n2; p++) if (this.board[p] === EMPTY && this.analyse(p)) m[p] = 1;
        this._legal = m;
      }
      return this._legal;
    }
    play(move) {
      const n2 = this.n * this.n;
      if (move === n2)
        return new State(this.n, this.komi, this.board, -this.toPlay, this.passes + 1,
                         this.lo, this.hi, this, this.moveNum + 1, move);
      const a = this.analyse(move);
      if (!a) throw new Error("illegal move");
      const b = this.board.slice();
      b[move] = this.toPlay;
      for (const q of a.captured) b[q] = EMPTY;
      const s = new State(this.n, this.komi, b, -this.toPlay, 0, a.lo, a.hi, this, this.moveNum + 1, move);
      s.captures = a.captured.length;
      return s;
    }
    isTerminal() { return this.passes >= 2 || this.moveNum >= this.maxMoves; }
    ownership() {
      const { nbrs } = geometry(this.n), b = this.board, n2 = this.n * this.n;
      const own = Int8Array.from(b), seen = new Uint8Array(n2);
      for (let p = 0; p < n2; p++) {
        if (b[p] !== EMPTY || seen[p]) continue;
        const region = [p], stack = [p], border = new Set(); seen[p] = 1;
        while (stack.length) {
          const q = stack.pop();
          for (const r of nbrs[q]) {
            if (b[r] === EMPTY) { if (!seen[r]) { seen[r] = 1; region.push(r); stack.push(r); } }
            else border.add(b[r]);
          }
        }
        const owner = border.size === 1 ? [...border][0] : EMPTY;
        for (const q of region) own[q] = owner;
      }
      return own;
    }
    score() {
      let s = -this.komi;
      for (const v of this.ownership()) s += v;
      return s;
    }
    winner() { return this.score() > 0 ? BLACK : WHITE; }
    features(T) {
      const n2 = this.n * this.n, x = new Float32Array((2 * T + 1) * n2);
      let s = this;
      for (let i = 0; i < T && s; i++, s = s.prev) {
        const b = s.board, o1 = 2 * i * n2, o2 = o1 + n2;
        for (let p = 0; p < n2; p++) {
          if (b[p] === this.toPlay) x[o1 + p] = 1;
          else if (b[p] === -this.toPlay) x[o2 + p] = 1;
        }
      }
      if (this.toPlay === BLACK) x.fill(1, 2 * T * n2);
      return x;
    }
    moveName(m) {
      if (m === this.n * this.n) return "pass";
      return COLS[m % this.n] + (this.n - ((m / this.n) | 0));
    }
  }

  // ------------------------------------------------------------- network
  function b64ToFloat32(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return new Float32Array(u8.buffer);
  }

  // 3x3 convolution as im2col + dot products (contiguous inner loop, fast in JS)
  const nb9Cache = {};
  function nb9(n) {
    if (nb9Cache[n]) return nb9Cache[n];
    const t = new Int32Array(n * n * 9);
    for (let p = 0; p < n * n; p++) {
      const y = (p / n) | 0, x = p % n;
      for (let k = 0; k < 9; k++) {
        const yy = y + ((k / 3) | 0) - 1, xx = x + (k % 3) - 1;
        t[p * 9 + k] = yy >= 0 && yy < n && xx >= 0 && xx < n ? yy * n + xx : -1;
      }
    }
    return (nb9Cache[n] = t);
  }
  let colBuf = new Float32Array(0);
  function conv3x3(inp, cin, W, B, cout, n, out) {
    const n2 = n * n, K = cin * 9, nb = nb9(n);
    if (colBuf.length < n2 * K) colBuf = new Float32Array(n2 * K);
    const col = colBuf;
    for (let p = 0; p < n2; p++) {
      const cb = p * K, pb = p * 9;
      for (let ci = 0; ci < cin; ci++) {
        const ib = ci * n2, o = cb + ci * 9;
        for (let k = 0; k < 9; k++) { const q = nb[pb + k]; col[o + k] = q >= 0 ? inp[ib + q] : 0; }
      }
    }
    for (let co = 0; co < cout; co++) {
      const wb = co * K, bias = B[co], ob = co * n2;
      for (let p = 0; p < n2; p++) {
        let s = bias;
        const cb = p * K;
        for (let j = 0; j < K; j++) s += W[wb + j] * col[cb + j];
        out[ob + p] = s;
      }
    }
  }
  function conv1x1(inp, cin, W, B, cout, n2, out) {
    for (let co = 0; co < cout; co++) {
      const ob = co * n2;
      out.fill(B[co], ob, ob + n2);
      for (let ci = 0; ci < cin; ci++) {
        const w = W[co * cin + ci], ib = ci * n2;
        for (let p = 0; p < n2; p++) out[ob + p] += w * inp[ib + p];
      }
    }
  }
  function fc(inp, W, B, nout, nin, out) {
    for (let j = 0; j < nout; j++) {
      let s = B[j];
      const wb = j * nin;
      for (let i = 0; i < nin; i++) s += W[wb + i] * inp[i];
      out[j] = s;
    }
  }
  const relu = (a) => { for (let i = 0; i < a.length; i++) if (a[i] < 0) a[i] = 0; return a; };

  class Net {
    constructor(model) {
      Object.assign(this, model.meta);
      const data = b64ToFloat32(model.weights);
      this.L = {};
      for (const l of model.layers) this.L[l.name] = data.subarray(l.offset, l.offset + l.size);
      const n2 = this.n * this.n, F = this.filters;
      this.a = new Float32Array(F * n2); this.b = new Float32Array(F * n2); this.c = new Float32Array(F * n2);
    }
    forward(x) {
      const n = this.n, n2 = n * n, F = this.filters, L = this.L;
      let { a, b, c } = this;
      conv3x3(x, this.in_planes, L["stem.w"], L["stem.b"], F, n, a); relu(a);
      for (let i = 0; i < this.blocks; i++) {
        conv3x3(a, F, L[`b${i}.c1.w`], L[`b${i}.c1.b`], F, n, b); relu(b);
        conv3x3(b, F, L[`b${i}.c2.w`], L[`b${i}.c2.b`], F, n, c);
        for (let j = 0; j < a.length; j++) { const v = a[j] + c[j]; a[j] = v > 0 ? v : 0; }
      }
      const ph = new Float32Array(2 * n2), logits = new Float32Array(n2 + 1);
      conv1x1(a, F, L["p_conv.w"], L["p_conv.b"], 2, n2, ph); relu(ph);
      fc(ph, L["p_fc.w"], L["p_fc.b"], n2 + 1, 2 * n2, logits);
      const vh = new Float32Array(n2), h = new Float32Array(this.value_hidden), v = new Float32Array(1);
      conv1x1(a, F, L["v_conv.w"], L["v_conv.b"], 1, n2, vh); relu(vh);
      fc(vh, L["v_fc1.w"], L["v_fc1.b"], this.value_hidden, n2, h); relu(h);
      fc(h, L["v_fc2.w"], L["v_fc2.b"], 1, this.value_hidden, v);
      return { logits, value: Math.tanh(v[0]) };
    }
    // (priors over legal moves, value) for a state, under symmetry k
    evaluate(state, k = 0) {
      const n2 = this.n * this.n, perm = geometry(this.n).perms[k];
      const x = state.features(this.history), xt = new Float32Array(x.length);
      for (let ch = 0; ch < this.in_planes; ch++)
        for (let p = 0; p < n2; p++) xt[ch * n2 + p] = x[ch * n2 + perm[p]];
      const { logits, value } = this.forward(xt);
      const lg = new Float64Array(n2 + 1);
      for (let p = 0; p < n2; p++) lg[perm[p]] = logits[p];
      lg[n2] = logits[n2];
      const legal = state.legalMask();
      let mx = -Infinity;
      for (let i = 0; i <= n2; i++) if (legal[i] && lg[i] > mx) mx = lg[i];
      let sum = 0;
      const P = new Float64Array(n2 + 1);
      for (let i = 0; i <= n2; i++) if (legal[i]) { P[i] = Math.exp(lg[i] - mx); sum += P[i]; }
      for (let i = 0; i <= n2; i++) P[i] /= sum;
      return { P, value };
    }
  }

  // ---------------------------------------------------------------- MCTS
  class Node {
    constructor(state) { this.state = state; this.children = new Map(); this.expanded = false; }
    expand(P) {
      const A = P.length;
      this.P = P; this.N = new Float64Array(A); this.W = new Float64Array(A); this.Nsum = 0;
      this.legal = this.state.legalMask(); this.expanded = true;
    }
    select(cPuct) {
      const s = Math.sqrt(Math.max(this.Nsum, 1));
      let best = -1, bv = -Infinity;
      for (let a = 0; a < this.P.length; a++) {
        if (!this.legal[a]) continue;
        const n = this.N[a];
        const v = (n > 0 ? this.W[a] / n : 0) + cPuct * this.P[a] * s / (1 + n);
        if (v > bv) { bv = v; best = a; }
      }
      return best;
    }
    value() { let w = 0; for (const x of this.W) w += x; return this.Nsum ? w / this.Nsum : 0; }
  }

  class Search {
    constructor(net, cPuct = 1.5) { this.net = net; this.cPuct = cPuct; }
    evalNode(node) {
      const k = (Math.random() * 8) | 0;      // random symmetry, as in the paper
      const r = this.net.evaluate(node.state, k);
      node.expand(r.P);
      node.netValue = r.value;
      return r.value;
    }
    simulate(root) {
      if (!root.expanded) this.evalNode(root);
      let node = root; const path = [];
      while (node.expanded) {
        const a = node.select(this.cPuct);
        path.push([node, a]);
        let nxt = node.children.get(a);
        if (!nxt) { nxt = new Node(node.state.play(a)); node.children.set(a, nxt); }
        node = nxt;
      }
      let v;
      if (node.state.isTerminal()) v = node.state.winner() === node.state.toPlay ? 1 : -1;
      else v = this.evalNode(node);
      for (let i = path.length - 1; i >= 0; i--) {
        const [nd, a] = path[i];
        v = -v; nd.N[a] += 1; nd.W[a] += v; nd.Nsum += 1;
      }
    }
  }

  return { BLACK, WHITE, EMPTY, COLS, State, Net, Node, Search, geometry };
})();
if (typeof module !== "undefined") module.exports = GoZero;
