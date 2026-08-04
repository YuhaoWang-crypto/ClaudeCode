// ---------------------------------------------------------------------------
//  model.js - the Particle Life rule set.
//
//  Force law (identical to sandbox-science.com/particle-life's WebGPU kernel):
//
//      m   = (alpha + beta)/2
//      f(r)= repel * (r/alpha - 1)                 0     <= r < alpha
//            A * (1 - |r - m| / (m - alpha))       alpha <= r < beta
//            0                                     r     >= beta
//
//  f > 0 pulls a towards b.  A, alpha and beta are all *ordered-pair* matrices,
//  so f_ab != f_ba in general: Newton's third law is optional here, and that is
//  precisely what makes the structures move.
// ---------------------------------------------------------------------------

export function pairForce(r, A, alpha, beta, repel) {
  if (r >= beta) return 0;
  if (r < alpha) return (r / alpha - 1) * repel;
  const m = 0.5 * (alpha + beta);
  return A * (1 - Math.abs(r - m) / (m - alpha));
}

export function pairPotential(r, A, alpha, beta, repel) {
  const w = 0.5 * (beta - alpha), m = 0.5 * (alpha + beta);
  if (r >= beta) return 0;
  if (r < alpha) return (repel / (2 * alpha)) * (r - alpha) ** 2 - A * w;
  if (r < m) return (-A / (2 * w)) * (2 * w * w - (r - alpha) ** 2);
  return (-A / (2 * w)) * (beta - r) ** 2;
}

// --------------------------------------------------------------------------- //
//  interaction-matrix generators (ports of the platform's preset list)
// --------------------------------------------------------------------------- //
const zeros = (S) => Array.from({ length: S }, () => new Array(S).fill(0));
const full = (S, v) => Array.from({ length: S }, () => new Array(S).fill(v));

export const MATRIX_GENERATORS = {
  random: (S) => Array.from({ length: S }, () =>
    Array.from({ length: S }, () => Math.random() * 2 - 1)),

  symmetric: (S) => {
    const A = MATRIX_GENERATORS.random(S);
    for (let i = 0; i < S; i++) for (let j = i + 1; j < S; j++) A[j][i] = A[i][j];
    return A;
  },

  snake: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) { A[i][i] = 1; A[i][(i + 1) % S] = 0.2; }
    return A;
  },

  chains1: (S) => {
    const A = full(S, -1);
    for (let i = 0; i < S; i++) { A[i][i] = 1; A[i][(i + 1) % S] = 1; A[i][(i + S - 1) % S] = 1; }
    return A;
  },

  chains2: (S) => {
    const A = full(S, -1);
    for (let i = 0; i < S; i++) { A[i][i] = 1; A[i][(i + 1) % S] = 0.2; A[i][(i + S - 1) % S] = 0.2; }
    return A;
  },

  chains3: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) { A[i][i] = 1; A[i][(i + 1) % S] = 0.2; A[i][(i + S - 1) % S] = 0.2; }
    return A;
  },

  rps: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) { A[i][i] = -0.1; A[i][(i + 1) % S] = 0.9; A[i][(i + S - 1) % S] = -0.7; }
    return A;
  },

  bipartite: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++)
      A[i][j] = i === j ? 0.2 : ((i % 2) === (j % 2) ? 0.8 : -0.8);
    return A;
  },

  hub: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++)
      A[i][j] = i === j ? 0 : (i === 0 ? 1 : (j === 0 ? 0.6 : 0));
    return A;
  },

  shells: (S) => {
    const A = full(S, -0.6);
    for (let i = 0; i < S; i++) { A[i][i] = 0.9; A[i][(i + 1) % S] = 0.3; }
    return A;
  },

  checker: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
      if (i === j) A[i][j] = 0.2;
      else if (j === (i + 1) % S || j === (i + S - 1) % S) A[i][j] = -0.8;
      else if (j === (i + 2) % S || j === (i + S - 2) % S) A[i][j] = 0.8;
    }
    return A;
  },

  dimers: (S) => {
    const A = full(S, -0.9);
    for (let i = 0; i < S; i++) { A[i][i] = 0; const p = i ^ 1; if (p < S) A[i][p] = 1; }
    return A;
  },

  triads: (S) => {
    const A = zeros(S), g = (i) => Math.floor(i / 3);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++)
      A[i][j] = i === j ? 0.1 : (g(i) === g(j) ? 0.9 : -0.7);
    return A;
  },

  spiralConveyor: (S) => {
    const A = full(S, -0.6);
    for (let i = 0; i < S; i++) { A[i][i] = -0.1; A[i][(i + 1) % S] = 0.7; A[i][(i + 2) % S] = 0.3; }
    return A;
  },

  patchwork: (S) => {
    const A = zeros(S);
    const h = (i, j) => Math.abs(Math.sin((i * 73856093) ^ (j * 19349663)) * 43758.5453 % 1);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
      if (i === j) continue;
      const v = h(i, j);
      A[i][j] = v < 0.35 ? 0.9 : (v < 0.7 ? -0.9 : 0);
    }
    return A;
  },

  wavefield: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++)
      A[i][j] = i === j ? -0.05 : 0.9 * Math.sin(2 * Math.PI * (j - i) / S);
    return A;
  },

  chiralBandpass: (S) => {
    const A = zeros(S);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
      if (i === j) { A[i][j] = -0.05; continue; }
      const d = (j - i + S) % S, c = 2 * Math.PI * d / S;
      A[i][j] = 0.7 * Math.sin(c) + 0.35 * Math.sin(2 * c) + 0.15 * (d > 0 ? 1 : -1);
    }
    return A;
  },

  rotatingConveyor: (S) => {
    const A = zeros(S), off = 2 * Math.PI / (S * 1.5);
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++)
      A[i][j] = i === j ? -0.05 : 0.8 * Math.sin(2 * Math.PI * (j - i) / S + off * i);
    return A;
  },

  parityVortex: (S) => {
    const A = zeros(S), ev = (i) => i % 2 === 0;
    for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
      if (i === j) A[i][j] = -0.05;
      else if (ev(i) && !ev(j)) A[i][j] = 0.8;
      else if (!ev(i) && ev(j)) A[i][j] = 0.6;
      else A[i][j] = ev(i) ? -0.4 : -0.6;
    }
    return A;
  },

  condensate: (S) => {
    const A = full(S, -0.3);
    for (let i = 0; i < S; i++) A[i][i] = 0.9;
    return A;
  },

  gas: (S) => full(S, -0.6),
};

export const MATRIX_NAMES = {
  random: "Random", symmetric: "Symmetric (reciprocal)", snake: "Snake",
  chains1: "Chains 1", chains2: "Chains 2", chains3: "Chains 3",
  rps: "Rock–Paper–Scissors", bipartite: "Bipartite Alliances",
  hub: "Hub and Spokes", shells: "Concentric Shells", checker: "Checker Offsets",
  dimers: "Dimers & Chains", triads: "Triad Flocks",
  spiralConveyor: "Spiral Conveyor", patchwork: "Patchwork", wavefield: "Wavefield",
  chiralBandpass: "Chiral Bandpass", rotatingConveyor: "Rotating Conveyor",
  parityVortex: "Parity Vortex", condensate: "Condensate (control)",
  gas: "All-repulsive gas (control)",
};

export function generateMatrix(name, S) {
  const gen = MATRIX_GENERATORS[name] || MATRIX_GENERATORS.random;
  return gen(S).map((row) => row.map((v) => Math.max(-1, Math.min(1, Math.round(v * 100) / 100))));
}

// --------------------------------------------------------------------------- //
//  radii
// --------------------------------------------------------------------------- //
export function randomRadii(S, minRange, maxRange, symmetric = false) {
  const rnd = (a, b) => a + Math.random() * (b - a);
  const alpha = Array.from({ length: S }, () => Array.from({ length: S }, () => rnd(...minRange)));
  const beta = Array.from({ length: S }, () => Array.from({ length: S }, () => rnd(...maxRange)));
  for (let i = 0; i < S; i++) for (let j = 0; j < S; j++) {
    if (symmetric && j > i) { alpha[j][i] = alpha[i][j]; beta[j][i] = beta[i][j]; }
    beta[i][j] = Math.max(beta[i][j], alpha[i][j] * 1.05 + 1);
  }
  return { alpha, beta };
}

// --------------------------------------------------------------------------- //
//  reciprocity diagnostics - the single most informative number about a matrix
// --------------------------------------------------------------------------- //
export function nonreciprocity(A) {
  // NOTE: this only sees the force matrix.  alpha and beta are ordered-pair
  // matrices too, so use nonreciprocityForces() for the quantity the
  // reciprocity theorem actually depends on.
  let s = 0, n = 0;
  for (let i = 0; i < A.length; i++) for (let j = 0; j < A.length; j++) {
    s += (A[i][j] + A[j][i]) ** 2;
    n += (A[i][j] - A[j][i]) ** 2;
  }
  s = Math.sqrt(s); n = Math.sqrt(n);
  return s + n > 0 ? n / (s + n) : 0;
}

/**
 * Non-reciprocity of the actual pair forces, over all three matrices:
 *   nu_f = ||f_ab - f_ba|| / (||f_ab + f_ba|| + ||f_ab - f_ba||)
 * Zero iff f_ab === f_ba for every pair - exactly the hypothesis of the
 * reciprocity theorem.  The platform's randomiser samples min/max radius per
 * ordered pair, so a symmetric force matrix alone does NOT make nu_f zero.
 */
export function nonreciprocityForces(A, alpha, beta, repel, nR = 200) {
  const S = A.length;
  let rmax = 0;
  for (let a = 0; a < S; a++) for (let b = 0; b < S; b++) rmax = Math.max(rmax, beta[a][b]);
  let sym2 = 0, asym2 = 0;
  for (let a = 0; a < S; a++) for (let b = a; b < S; b++) {
    for (let i = 0; i < nR; i++) {
      const r = ((i + 0.5) / nR) * rmax;
      const fab = pairForce(r, A[a][b], alpha[a][b], beta[a][b], repel);
      const fba = pairForce(r, A[b][a], alpha[b][a], beta[b][a], repel);
      sym2 += (fab + fba) ** 2;
      asym2 += (fab - fba) ** 2;
    }
  }
  const s = Math.sqrt(sym2), n = Math.sqrt(asym2);
  return s + n > 0 ? n / (s + n) : 0;
}

export function isReciprocal(A, alpha, beta, tol = 1e-9) {
  for (let i = 0; i < A.length; i++) for (let j = 0; j < A.length; j++) {
    if (Math.abs(A[i][j] - A[j][i]) > tol) return false;
    if (Math.abs(alpha[i][j] - alpha[j][i]) > tol) return false;
    if (Math.abs(beta[i][j] - beta[j][i]) > tol) return false;
  }
  return true;
}

// --------------------------------------------------------------------------- //
//  colour palettes
// --------------------------------------------------------------------------- //
function hsl(h, s, l) { return `hsl(${((h % 360) + 360) % 360} ${s}% ${l}%)`; }

export const PALETTES = {
  rainbow: (S) => Array.from({ length: S }, (_, i) => hsl(360 * i / S, 90, 58)),
  neonWarm: (S) => Array.from({ length: S }, (_, i) => hsl(-20 + 120 * i / Math.max(1, S - 1), 95, 60)),
  coldBlue: (S) => Array.from({ length: S }, (_, i) => hsl(180 + 80 * i / Math.max(1, S - 1), 85, 45 + 25 * i / Math.max(1, S - 1))),
  pastel: (S) => Array.from({ length: S }, (_, i) => hsl(360 * i / S, 65, 76)),
  fire: (S) => Array.from({ length: S }, (_, i) => hsl(0 + 55 * i / Math.max(1, S - 1), 100, 40 + 30 * i / Math.max(1, S - 1))),
  vaporwave: (S) => Array.from({ length: S }, (_, i) => hsl(280 + 120 * i / Math.max(1, S - 1), 90, 65)),
  grayscale: (S) => Array.from({ length: S }, (_, i) => hsl(0, 0, 35 + 55 * i / Math.max(1, S - 1))),
  gameboy: (S) => Array.from({ length: S }, (_, i) => hsl(95, 45, 25 + 50 * i / Math.max(1, S - 1))),
  greenMagenta: (S) => Array.from({ length: S }, (_, i) => (i % 2 ? hsl(315, 100, 60) : hsl(130, 90, 50))),
  random: (S) => Array.from({ length: S }, () => hsl(Math.random() * 360, 70 + Math.random() * 30, 45 + Math.random() * 25)),
};

export const PALETTE_NAMES = {
  rainbow: "Rainbow", neonWarm: "Neon Warm", coldBlue: "Cold Blue", pastel: "Pastel",
  fire: "Fire", vaporwave: "Vaporwave", grayscale: "Grayscale", gameboy: "Game Boy",
  greenMagenta: "Green / Magenta", random: "Random",
};

// --------------------------------------------------------------------------- //
//  initial particle distributions
// --------------------------------------------------------------------------- //
export const DISTRIBUTIONS = {
  random: (n, S, W, H) => {
    const p = [];
    for (let i = 0; i < n; i++) p.push([Math.random() * W, Math.random() * H, i % S]);
    return p;
  },
  disk: (n, S, W, H) => {
    const p = [], R = 0.46 * Math.min(W, H);
    for (let i = 0; i < n; i++) {
      const t = Math.random() * 2 * Math.PI, r = R * Math.sqrt(Math.random());
      p.push([W / 2 + r * Math.cos(t), H / 2 + r * Math.sin(t), i % S]);
    }
    return p;
  },
  ring: (n, S, W, H) => {
    const p = [], R = 0.44 * Math.min(W, H);
    for (let i = 0; i < n; i++) {
      const t = Math.random() * 2 * Math.PI, r = R * (1 - 0.08 * Math.random());
      p.push([W / 2 + r * Math.cos(t), H / 2 + r * Math.sin(t), i % S]);
    }
    return p;
  },
  rainbowSpiral: (n, S, W, H) => {
    const p = [], R = 0.45 * Math.min(W, H);
    for (let i = 0; i < n; i++) {
      const s = i % S, f = i / n;
      const t = 6 * Math.PI * f + (2 * Math.PI * s) / S, r = R * Math.sqrt(f);
      p.push([W / 2 + r * Math.cos(t), H / 2 + r * Math.sin(t), s]);
    }
    return p;
  },
  stripes: (n, S, W, H) => {
    const p = [];
    for (let i = 0; i < n; i++) {
      const s = i % S;
      p.push([W * (s + Math.random()) / S, Math.random() * H, s]);
    }
    return p;
  },
  grid: (n, S, W, H) => {
    const p = [], m = Math.ceil(Math.sqrt(n));
    for (let i = 0; i < n; i++) {
      const gx = i % m, gy = Math.floor(i / m);
      p.push([(gx + 0.5) * W / m, (gy + 0.5) * H / m, i % S]);
    }
    return p;
  },
  clusters: (n, S, W, H) => {
    const p = [], nc = Math.max(3, S * 2), cx = [], cy = [];
    for (let c = 0; c < nc; c++) { cx.push(0.1 * W + 0.8 * W * Math.random()); cy.push(0.1 * H + 0.8 * H * Math.random()); }
    const R = 0.06 * Math.min(W, H);
    for (let i = 0; i < n; i++) {
      const c = i % nc, t = Math.random() * 2 * Math.PI, r = R * Math.sqrt(Math.random());
      p.push([cx[c] + r * Math.cos(t), cy[c] + r * Math.sin(t), i % S]);
    }
    return p;
  },
};

export const DISTRIBUTION_NAMES = {
  random: "Random", disk: "Disk", ring: "Ring", rainbowSpiral: "Rainbow Spiral",
  stripes: "Stripes", grid: "Grid", clusters: "Soft Clusters",
};
