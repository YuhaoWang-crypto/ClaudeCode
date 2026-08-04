// ---------------------------------------------------------------------------
//  theory.js - live continuum theory for the running simulation.
//
//  Implements, in the browser, the same linear-stability calculation as
//  theory/plife/stability.py:
//
//      Uhat_ab(k) = 2 pi \int_0^beta r J0(kr) U_ab(r) dr
//      sigma_n(k) = eig[ diag(rho) Uhat(k) ]
//      lambda_+   = ( -gamma + sqrt(gamma^2 - 240 kappa k^2 sigma_n) ) / 2
//
//  plus the measured counterparts (structure factor peak, motility) so the UI
//  can show prediction and measurement side by side.
// ---------------------------------------------------------------------------

import { pairPotential } from "./model.js";

// --------------------------------------------------------------------------- //
//  Bessel J0  (Abramowitz & Stegun 9.4.1 / 9.4.3, |err| < 1e-8)
// --------------------------------------------------------------------------- //
export function besselJ0(x) {
  const ax = Math.abs(x);
  if (ax < 8) {
    const y = x * x;
    const p = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 +
      y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
    const q = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 +
      y * (59272.64853 + y * (267.8532712 + y))));
    return p / q;
  }
  const z = 8 / ax, y = z * z, xx = ax - 0.785398164;
  const p = 1 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 +
    y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
  const q = -0.1562499995e-1 + y * (0.1430488765e-3 +
    y * (-0.6911147651e-5 + y * (0.7621095161e-6 + y * (-0.934935152e-7))));
  return Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * p - z * Math.sin(xx) * q);
}

// --------------------------------------------------------------------------- //
//  Hankel transform of the pair potential.
//  Split at the two kinks (alpha and the tent apex) so composite Simpson is
//  applied only to smooth pieces -> near machine accuracy with few nodes.
// --------------------------------------------------------------------------- //
function simpson(f, a, b, n) {
  if (b <= a) return 0;
  if (n % 2) n++;
  const h = (b - a) / n;
  let s = f(a) + f(b);
  for (let i = 1; i < n; i++) s += f(a + i * h) * (i % 2 ? 4 : 2);
  return (s * h) / 3;
}

export function hankelU(k, A, alpha, beta, repel, n = 60) {
  const m = 0.5 * (alpha + beta);
  const g = (r) => r * besselJ0(k * r) * pairPotential(r, A, alpha, beta, repel);
  return 2 * Math.PI * (simpson(g, 0, alpha, n) + simpson(g, alpha, m, n) +
    simpson(g, m, beta, n));
}

/** Uhat(k) for every ordered species pair -> S x S array. */
export function hankelMatrix(k, A, alpha, beta, repel) {
  const S = A.length;
  const out = [];
  for (let a = 0; a < S; a++) {
    const row = [];
    for (let b = 0; b < S; b++) row.push(hankelU(k, A[a][b], alpha[a][b], beta[a][b], repel));
    out.push(row);
  }
  return out;
}

// --------------------------------------------------------------------------- //
//  eigenvalues of a small real matrix: Hessenberg reduction + shifted QR.
//  Returns [{re, im}, ...].  Adequate for S <= 16.
// --------------------------------------------------------------------------- //
export function eigenvalues(Min) {
  const n = Min.length;
  const a = Min.map((r) => r.slice());

  // --- Householder reduction to upper Hessenberg form
  for (let k = 1; k < n - 1; k++) {
    let scale = 0;
    for (let i = k; i < n; i++) scale += Math.abs(a[i][k - 1]);
    if (scale === 0) continue;
    let h = 0;
    const v = new Array(n).fill(0);
    for (let i = k; i < n; i++) { v[i] = a[i][k - 1] / scale; h += v[i] * v[i]; }
    let g = v[k] >= 0 ? -Math.sqrt(h) : Math.sqrt(h);
    h -= v[k] * g;
    v[k] -= g;
    for (let j = 0; j < n; j++) {          // A <- A (I - vv^T/h)
      let f = 0;
      for (let i = k; i < n; i++) f += v[i] * a[j][i];
      f /= h;
      for (let i = k; i < n; i++) a[j][i] -= f * v[i];
    }
    for (let i = 0; i < n; i++) {          // A <- (I - vv^T/h) A
      let f = 0;
      for (let j = k; j < n; j++) f += v[j] * a[j][i];
      f /= h;
      for (let j = k; j < n; j++) a[j][i] -= f * v[j];
    }
    a[k][k - 1] = scale * g;
    for (let i = k + 1; i < n; i++) a[i][k - 1] = 0;
  }

  // --- Francis-style QR with Wilkinson shift on the trailing submatrix
  const ev = [];
  let hi = n - 1, iter = 0;
  const eps = 1e-13;
  while (hi >= 0) {
    // find a small subdiagonal element
    let lo = hi;
    while (lo > 0) {
      const s = Math.abs(a[lo - 1][lo - 1]) + Math.abs(a[lo][lo]);
      if (Math.abs(a[lo][lo - 1]) <= eps * (s || 1)) { a[lo][lo - 1] = 0; break; }
      lo--;
    }
    if (lo === hi) { ev.push({ re: a[hi][hi], im: 0 }); hi--; iter = 0; continue; }
    if (lo === hi - 1) {                    // 2x2 block -> closed form
      const p = 0.5 * (a[hi - 1][hi - 1] - a[hi][hi]);
      const q = p * p + a[hi][hi - 1] * a[hi - 1][hi];
      if (q >= 0) {
        const z = Math.sqrt(q) * (p >= 0 ? 1 : -1);
        ev.push({ re: a[hi][hi] + p + z, im: 0 });
        ev.push({ re: a[hi][hi] + p - (z !== 0 ? (p * p + q - z * z) / z || 0 : 0), im: 0 });
        ev[ev.length - 1] = { re: a[hi - 1][hi - 1] + a[hi][hi] - ev[ev.length - 2].re, im: 0 };
      } else {
        const c = a[hi][hi] + p, d = Math.sqrt(-q);
        ev.push({ re: c, im: d });
        ev.push({ re: c, im: -d });
      }
      hi -= 2; iter = 0; continue;
    }
    if (++iter > 90) {                       // give up on this block
      ev.push({ re: a[hi][hi], im: 0 }); hi--; iter = 0; continue;
    }

    // Wilkinson shift
    const d = 0.5 * (a[hi - 1][hi - 1] - a[hi][hi]);
    const bb = a[hi][hi - 1] * a[hi - 1][hi];
    let mu = a[hi][hi];
    const disc = d * d + bb;
    if (disc >= 0) {
      const sq = Math.sqrt(disc);
      mu = a[hi][hi] + d - (d >= 0 ? sq : -sq);
    }
    for (let i = lo; i <= hi; i++) a[i][i] -= mu;

    // explicit QR of the active block via Givens rotations
    const cs = [], sn = [];
    for (let i = lo; i < hi; i++) {
      const x = a[i][i], y = a[i + 1][i];
      const r = Math.hypot(x, y) || 1;
      const c = x / r, s = y / r;
      cs.push(c); sn.push(s);
      for (let j = i; j <= hi; j++) {
        const t1 = a[i][j], t2 = a[i + 1][j];
        a[i][j] = c * t1 + s * t2;
        a[i + 1][j] = -s * t1 + c * t2;
      }
    }
    for (let i = lo; i < hi; i++) {
      const c = cs[i - lo], s = sn[i - lo];
      for (let j = lo; j <= Math.min(hi, i + 2); j++) {
        const t1 = a[j][i], t2 = a[j][i + 1];
        a[j][i] = c * t1 + s * t2;
        a[j][i + 1] = -s * t1 + c * t2;
      }
    }
    for (let i = lo; i <= hi; i++) a[i][i] += mu;
  }
  return ev;
}

// --------------------------------------------------------------------------- //
//  the dispersion relation
// --------------------------------------------------------------------------- //
function csqrt(re, im) {
  const r = Math.hypot(re, im);
  const a = Math.sqrt(Math.max(0, (r + re) / 2));
  let b = Math.sqrt(Math.max(0, (r - re) / 2));
  if (im < 0) b = -b;
  return [a, b];
}

/**
 * Full linear analysis of the uniform state.
 *
 *  @returns {kStar, growth, omega, length, oscGrowth, oscK, oscSpeed, regime, curve}
 */
export function analyseStability(A, alpha, beta, repel, rho, gamma, kappa, opts = {}) {
  const S = A.length;
  const nk = opts.nk ?? 90;
  let alMin = Infinity;
  for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) alMin = Math.min(alMin, alpha[i][j]);
  const kMax = opts.kMax ?? (5 * Math.PI / alMin);

  const curve = [];
  let best = { growth: -Infinity }, bestOsc = { growth: -Infinity };

  for (let ik = 1; ik <= nk; ik++) {
    const k = (kMax * ik) / nk;
    const U = hankelMatrix(k, A, alpha, beta, repel);
    const M = U.map((row, a) => row.map((v) => rho[a] * v));
    const sig = eigenvalues(M);
    let maxRe = -Infinity, maxIm = 0;
    for (const s of sig) {
      // lambda^2 + gamma lambda + 240/4 * ... => lambda = (-g + sqrt(g^2 - 240 kappa k^2 sigma))/2
      const zr = gamma * gamma - 240 * kappa * k * k * s.re;
      const zi = -240 * kappa * k * k * s.im;
      const [sr, si] = csqrt(zr, zi);
      for (const sgn of [1, -1]) {
        const lr = 0.5 * (-gamma + sgn * sr), li = 0.5 * (sgn * si);
        if (lr > maxRe) { maxRe = lr; maxIm = Math.abs(li); }
        if (Math.abs(li) > 1e-9 && lr > bestOsc.growth) {
          bestOsc = { growth: lr, k, omega: Math.abs(li) };
        }
      }
    }
    curve.push({ k, re: maxRe, im: maxIm });
    if (maxRe > best.growth) best = { growth: maxRe, k, omega: maxIm };
  }

  const oscGrowth = bestOsc.growth === -Infinity ? -Infinity : bestOsc.growth;
  const travelling = best.growth > 0 && best.omega > 1e-9;
  let regime = "uniform";
  if (best.growth > 0) regime = travelling ? "travelling" : (oscGrowth > 0 ? "condense+move" : "static");

  return {
    kStar: best.k, growth: best.growth, omega: best.omega,
    length: best.k > 0 ? (2 * Math.PI) / best.k : Infinity,
    oscGrowth, oscK: bestOsc.k ?? NaN,
    oscSpeed: bestOsc.k ? bestOsc.omega / bestOsc.k : 0,
    travelling, regime, curve,
  };
}

// --------------------------------------------------------------------------- //
//  measurements on the live particle cloud
// --------------------------------------------------------------------------- //
function fft(re, im, inverse = false) {
  const n = re.length;
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (inverse ? 2 : -2) * Math.PI / len;
    const wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      for (let j = 0; j < len / 2; j++) {
        const ur = re[i + j], ui = im[i + j];
        const vr = re[i + j + len / 2] * cr - im[i + j + len / 2] * ci;
        const vi = re[i + j + len / 2] * ci + im[i + j + len / 2] * cr;
        re[i + j] = ur + vr; im[i + j] = ui + vi;
        re[i + j + len / 2] = ur - vr; im[i + j + len / 2] = ui - vi;
        const nr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = nr;
      }
    }
  }
}

/** Radially averaged structure factor of the particle cloud (grid FFT). */
export function structureFactor(px, py, n, W, H, G = 64) {
  const dens = new Float64Array(G * G);
  for (let i = 0; i < n; i++) {
    let gx = Math.floor((px[i] / W) * G) % G, gy = Math.floor((py[i] / H) * G) % G;
    if (gx < 0) gx += G; if (gy < 0) gy += G;
    dens[gy * G + gx] += 1;
  }
  const mean = n / (G * G);
  for (let i = 0; i < G * G; i++) dens[i] -= mean;

  const re = new Float64Array(G * G), im = new Float64Array(G * G);
  re.set(dens);
  const rowR = new Float64Array(G), rowI = new Float64Array(G);
  for (let y = 0; y < G; y++) {
    for (let x = 0; x < G; x++) { rowR[x] = re[y * G + x]; rowI[x] = im[y * G + x]; }
    fft(rowR, rowI);
    for (let x = 0; x < G; x++) { re[y * G + x] = rowR[x]; im[y * G + x] = rowI[x]; }
  }
  for (let x = 0; x < G; x++) {
    for (let y = 0; y < G; y++) { rowR[y] = re[y * G + x]; rowI[y] = im[y * G + x]; }
    fft(rowR, rowI);
    for (let y = 0; y < G; y++) { re[y * G + x] = rowR[y]; im[y * G + x] = rowI[y]; }
  }

  const nb = 40, kmaxBin = Math.PI * G / Math.max(W, H);
  const sum = new Float64Array(nb), cnt = new Float64Array(nb);
  for (let y = 0; y < G; y++) for (let x = 0; x < G; x++) {
    const fx = x <= G / 2 ? x : x - G, fy = y <= G / 2 ? y : y - G;
    const kx = (2 * Math.PI * fx) / W, ky = (2 * Math.PI * fy) / H;
    const k = Math.hypot(kx, ky);
    const b = Math.floor((k / kmaxBin) * nb);
    if (b < 1 || b >= nb) continue;
    sum[b] += (re[y * G + x] ** 2 + im[y * G + x] ** 2) / n;
    cnt[b]++;
  }
  const ks = [], Ss = [];
  for (let b = 1; b < nb; b++) if (cnt[b] > 0) {
    ks.push(((b + 0.5) / nb) * kmaxBin);
    Ss.push(sum[b] / cnt[b]);
  }
  let peak = NaN, bestv = -Infinity;
  for (let i = 0; i < ks.length; i++) if (Ss[i] > bestv) { bestv = Ss[i]; peak = ks[i]; }
  return { k: ks, S: Ss, peak };
}

/** Cluster count and mean cluster centre-of-mass speed, via grid union-find. */
export function clusterStats(px, py, vx, vy, n, W, H, cutoff) {
  const cell = cutoff;
  const nx = Math.max(1, Math.ceil(W / cell)), ny = Math.max(1, Math.ceil(H / cell));
  const heads = new Int32Array(nx * ny).fill(-1);
  const next = new Int32Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    let cx = Math.floor(px[i] / cell) % nx, cy = Math.floor(py[i] / cell) % ny;
    if (cx < 0) cx += nx; if (cy < 0) cy += ny;
    const c = cy * nx + cx;
    next[i] = heads[c]; heads[c] = i;
  }
  const parent = new Int32Array(n);
  for (let i = 0; i < n; i++) parent[i] = i;
  const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  const uni = (x, y) => { const a = find(x), b = find(y); if (a !== b) parent[a] = b; };

  const c2 = cutoff * cutoff;
  for (let cy = 0; cy < ny; cy++) for (let cx = 0; cx < nx; cx++) {
    for (let i = heads[cy * nx + cx]; i !== -1; i = next[i]) {
      for (let oy = -1; oy <= 1; oy++) for (let ox = -1; ox <= 1; ox++) {
        const gx = (cx + ox + nx) % nx, gy = (cy + oy + ny) % ny;
        for (let j = heads[gy * nx + gx]; j !== -1; j = next[j]) {
          if (j <= i) continue;
          let dx = px[j] - px[i], dy = py[j] - py[i];
          if (dx > W / 2) dx -= W; else if (dx < -W / 2) dx += W;
          if (dy > H / 2) dy -= H; else if (dy < -H / 2) dy += H;
          if (dx * dx + dy * dy < c2) uni(i, j);
        }
      }
    }
  }
  const size = new Map(), vsx = new Map(), vsy = new Map();
  for (let i = 0; i < n; i++) {
    const r = find(i);
    size.set(r, (size.get(r) || 0) + 1);
    vsx.set(r, (vsx.get(r) || 0) + vx[i]);
    vsy.set(r, (vsy.get(r) || 0) + vy[i]);
  }
  let nClusters = 0, wsum = 0, vsum = 0, biggest = 0;
  for (const [r, m] of size) {
    if (m < 4) continue;
    nClusters++;
    biggest = Math.max(biggest, m);
    vsum += m * Math.hypot(vsx.get(r) / m, vsy.get(r) / m);
    wsum += m;
  }
  return {
    nClusters, biggest,
    meanCluster: nClusters ? wsum / nClusters : 0,
    clusterSpeed: wsum ? vsum / wsum : 0,
  };
}
