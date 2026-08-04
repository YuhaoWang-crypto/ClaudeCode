// ---------------------------------------------------------------------------
//  engine.js - CPU Particle Life engine.
//
//  Same update order as the platform's WebGPU pipeline:
//      1. bin particles into a uniform grid of cell size max(beta)
//      2. force pass:  velocitySum = sum_j f_{ab}(r) * rhat
//                      v += velocitySum * (forceFactor * dt * 60)
//      3. advance:     v *= (1 - frictionFactor)^(dt*60);   x += v * dt
//      4. wrap / reflect at the walls
//
//  Everything lives in flat typed arrays; the per-ordered-pair parameters are
//  packed into S*S Float32Arrays so the inner loop is a single index lookup.
// ---------------------------------------------------------------------------

import { DISTRIBUTIONS } from "./model.js";

export class CpuEngine {
  constructor(opts) {
    this.W = opts.W; this.H = opts.H;
    this.S = opts.S;
    this.repel = opts.repel ?? 1;
    this.forceFactor = opts.forceFactor ?? 1;
    this.frictionFactor = opts.frictionFactor ?? 0.3;
    this.wallMode = opts.wallMode ?? "wrap";   // "wrap" | "repel" | "none"
    this.setMatrices(opts.A, opts.alpha, opts.beta);
    this.setCount(opts.N ?? 4000, opts.distribution ?? "random");
    this.steps = 0;
  }

  setMatrices(A, alpha, beta) {
    const S = this.S;
    this.A = A; this.alpha = alpha; this.beta = beta;
    this.fA = new Float32Array(S * S);
    this.fAl = new Float32Array(S * S);
    this.fBe = new Float32Array(S * S);
    let rmax = 1;
    for (let a = 0; a < S; a++) for (let b = 0; b < S; b++) {
      const k = a * S + b;
      this.fA[k] = A[a][b];
      this.fAl[k] = Math.max(0.5, alpha[a][b]);
      this.fBe[k] = Math.max(this.fAl[k] * 1.02, beta[a][b]);
      rmax = Math.max(rmax, this.fBe[k]);
    }
    this.rmax = rmax;
    this._buildGrid();
  }

  _buildGrid() {
    this.cell = Math.max(2, this.rmax);
    this.nx = Math.max(1, Math.ceil(this.W / this.cell));
    this.ny = Math.max(1, Math.ceil(this.H / this.cell));
    this.nCells = this.nx * this.ny;
    this.cellStart = new Int32Array(this.nCells + 1);
    this.cellCount = new Int32Array(this.nCells);
    this.order = new Int32Array(this.n || 0);
  }

  setCount(N, distribution = "random") {
    this.n = N;
    this.px = new Float32Array(N); this.py = new Float32Array(N);
    this.vx = new Float32Array(N); this.vy = new Float32Array(N);
    this.type = new Int32Array(N);
    this.order = new Int32Array(N);
    this.cellOf = new Int32Array(N);
    this.respawn(distribution);
  }

  respawn(distribution = "random") {
    const gen = DISTRIBUTIONS[distribution] || DISTRIBUTIONS.random;
    const pts = gen(this.n, this.S, this.W, this.H);
    for (let i = 0; i < this.n; i++) {
      this.px[i] = pts[i][0]; this.py[i] = pts[i][1];
      this.vx[i] = 0; this.vy[i] = 0;
      this.type[i] = pts[i][2] % this.S;
    }
    this.steps = 0;
  }

  resize(W, H) {
    const sx = W / this.W, sy = H / this.H;
    for (let i = 0; i < this.n; i++) { this.px[i] *= sx; this.py[i] *= sy; }
    this.W = W; this.H = H;
    this._buildGrid();
  }

  // ------------------------------------------------------------------ binning
  _bin() {
    const { nx, ny, cell, nCells } = this;
    this.cellCount.fill(0);
    for (let i = 0; i < this.n; i++) {
      let cx = Math.floor(this.px[i] / cell), cy = Math.floor(this.py[i] / cell);
      cx = cx < 0 ? 0 : (cx >= nx ? nx - 1 : cx);
      cy = cy < 0 ? 0 : (cy >= ny ? ny - 1 : cy);
      const c = cy * nx + cx;
      this.cellOf[i] = c;
      this.cellCount[c]++;
    }
    let acc = 0;
    for (let c = 0; c < nCells; c++) { this.cellStart[c] = acc; acc += this.cellCount[c]; }
    this.cellStart[nCells] = acc;
    const cursor = this.cellStart.slice(0, nCells);
    for (let i = 0; i < this.n; i++) this.order[cursor[this.cellOf[i]]++] = i;
  }

  // -------------------------------------------------------------------- step
  step(dt) {
    this._bin();
    const { nx, ny, S, W, H } = this;
    const wrap = this.wallMode === "wrap";
    const px = this.px, py = this.py, vx = this.vx, vy = this.vy, type = this.type;
    const fA = this.fA, fAl = this.fAl, fBe = this.fBe, repel = this.repel;
    const ff = this.forceFactor * dt * 60;
    const halfW = W * 0.5, halfH = H * 0.5;

    for (let idx = 0; idx < this.n; idx++) {
      const i = this.order[idx];
      const xi = px[i], yi = py[i], ti = type[i] * S;
      let sx = 0, sy = 0;
      const cx = this.cellOf[i] % nx, cy = (this.cellOf[i] / nx) | 0;

      for (let oy = -1; oy <= 1; oy++) {
        let gy = cy + oy;
        if (wrap) gy = (gy + ny) % ny; else if (gy < 0 || gy >= ny) continue;
        for (let ox = -1; ox <= 1; ox++) {
          let gx = cx + ox;
          if (wrap) gx = (gx + nx) % nx; else if (gx < 0 || gx >= nx) continue;
          const c = gy * nx + gx;
          const s0 = this.cellStart[c], s1 = this.cellStart[c + 1];
          for (let q = s0; q < s1; q++) {
            const j = this.order[q];
            if (j === i) continue;
            let dx = px[j] - xi, dy = py[j] - yi;
            if (wrap) {
              if (dx > halfW) dx -= W; else if (dx < -halfW) dx += W;
              if (dy > halfH) dy -= H; else if (dy < -halfH) dy += H;
            }
            const d2 = dx * dx + dy * dy;
            const k = ti + type[j];
            const be = fBe[k];
            if (d2 <= 1e-12 || d2 >= be * be) continue;
            const d = Math.sqrt(d2);
            const al = fAl[k];
            let f;
            if (d < al) {
              f = (d / al - 1) * repel;
            } else {
              const m = 0.5 * (al + be);
              f = fA[k] * (1 - Math.abs(d - m) / (m - al));
            }
            if (f !== 0) { const g = f / d; sx += dx * g; sy += dy * g; }
          }
        }
      }
      vx[i] += sx * ff;
      vy[i] += sy * ff;
    }

    // ---- advance
    const drag = Math.pow(1 - this.frictionFactor, dt * 60);
    for (let i = 0; i < this.n; i++) {
      vx[i] *= drag; vy[i] *= drag;
      px[i] += vx[i] * dt; py[i] += vy[i] * dt;
      if (this.wallMode === "wrap") {
        if (px[i] < 0) px[i] += W; else if (px[i] >= W) px[i] -= W;
        if (py[i] < 0) py[i] += H; else if (py[i] >= H) py[i] -= H;
      } else if (this.wallMode === "repel") {
        if (px[i] < 0) { px[i] = 0; vx[i] = -vx[i]; }
        else if (px[i] > W) { px[i] = W; vx[i] = -vx[i]; }
        if (py[i] < 0) { py[i] = 0; vy[i] = -vy[i]; }
        else if (py[i] > H) { py[i] = H; vy[i] = -vy[i]; }
      }
    }
    this.steps++;
  }

  // --------------------------------------------------------------- diagnostics
  meanSpeed() {
    let s = 0;
    for (let i = 0; i < this.n; i++) s += Math.hypot(this.vx[i], this.vy[i]);
    return s / Math.max(1, this.n);
  }

  /** Brush: push or pull particles near (x, y). */
  brush(x, y, radius, strength) {
    const r2 = radius * radius;
    for (let i = 0; i < this.n; i++) {
      let dx = this.px[i] - x, dy = this.py[i] - y;
      if (this.wallMode === "wrap") {
        if (dx > this.W / 2) dx -= this.W; else if (dx < -this.W / 2) dx += this.W;
        if (dy > this.H / 2) dy -= this.H; else if (dy < -this.H / 2) dy += this.H;
      }
      const d2 = dx * dx + dy * dy;
      if (d2 > r2 || d2 < 1e-6) continue;
      const d = Math.sqrt(d2), w = (1 - d / radius) * strength;
      this.vx[i] += (dx / d) * w;
      this.vy[i] += (dy / d) * w;
    }
  }
}
