// ---------------------------------------------------------------------------
//  main.js - wiring: engine + renderer + UI + live physics panel.
// ---------------------------------------------------------------------------

import {
  MATRIX_NAMES, PALETTES, PALETTE_NAMES, DISTRIBUTIONS, DISTRIBUTION_NAMES,
  generateMatrix, randomRadii, nonreciprocity, pairForce,
} from "./model.js";
import { CpuEngine } from "./engine.js";
import { analyseStability, structureFactor, clusterStats } from "./theory.js";

// --------------------------------------------------------------------------- //
const state = {
  S: 6,
  N: 6000,
  forceFactor: 1.0,
  frictionFactor: 0.12,
  repel: 1.0,
  dt: 1 / 60,
  substeps: 1,
  wallMode: "wrap",
  minRange: [12, 26],
  maxRange: [34, 68],
  symRadii: false,
  matrixPreset: "random",
  palettePreset: "greenMagenta",
  distPreset: "random",
  size: 1.75,
  trail: 0.55,
  glow: true,
  colourBySpeed: false,
  paused: false,
  tab: "force",
  sel: null,
  hover: null,
};

let A, alpha, beta, colors, engine;
const $ = (id) => document.getElementById(id);

// --------------------------------------------------------------------------- //
//  set-up
// --------------------------------------------------------------------------- //
function buildUniverse({ matrix = true, radii = true, palette = true, spawn = true } = {}) {
  if (matrix) A = generateMatrix(state.matrixPreset, state.S);
  if (radii) ({ alpha, beta } = randomRadii(state.S, state.minRange, state.maxRange, state.symRadii));
  if (palette) colors = (PALETTES[state.palettePreset] || PALETTES.rainbow)(state.S);

  const canvas = $("sim");
  const W = canvas.width, H = canvas.height;
  if (!engine) {
    engine = new CpuEngine({
      W, H, S: state.S, N: state.N, A, alpha, beta,
      repel: state.repel, forceFactor: state.forceFactor,
      frictionFactor: state.frictionFactor, wallMode: state.wallMode,
      distribution: state.distPreset,
    });
  } else {
    engine.S = state.S;
    engine.setMatrices(A, alpha, beta);
    engine.repel = state.repel;
    engine.forceFactor = state.forceFactor;
    engine.frictionFactor = state.frictionFactor;
    engine.wallMode = state.wallMode;
    if (spawn || engine.n !== state.N) engine.setCount(state.N, state.distPreset);
  }
  renderMatrix();
  refreshTheory();
}

// --------------------------------------------------------------------------- //
//  matrix grid UI
// --------------------------------------------------------------------------- //
function activeMatrix() {
  return state.tab === "force" ? A : (state.tab === "min" ? alpha : beta);
}

function cellColour(v) {
  if (state.tab === "force") {
    const t = Math.max(-1, Math.min(1, v));
    return t >= 0
      ? `rgb(${Math.round(20 + 20 * t)},${Math.round(60 + 150 * t)},${Math.round(90 + 60 * t)})`
      : `rgb(${Math.round(60 + 150 * -t)},${Math.round(24 + 10 * -t)},${Math.round(38 + 20 * -t)})`;
  }
  const lo = state.tab === "min" ? 2 : 10, hi = state.tab === "min" ? 60 : 140;
  const t = Math.max(0, Math.min(1, (v - lo) / (hi - lo)));
  return `rgb(${Math.round(20 + 40 * t)},${Math.round(50 + 90 * t)},${Math.round(80 + 130 * t)})`;
}

function renderMatrix() {
  const g = $("matrix-grid");
  const S = state.S;
  g.style.gridTemplateColumns = `22px repeat(${S}, 1fr)`;
  g.innerHTML = "";
  const corner = document.createElement("div");
  corner.className = "corner";
  g.appendChild(corner);
  for (let b = 0; b < S; b++) {
    const d = document.createElement("div");
    d.className = "hdr";
    d.style.background = colors[b];
    d.title = `exerted by colour ${b}`;
    g.appendChild(d);
  }
  const Mx = activeMatrix();
  for (let a = 0; a < S; a++) {
    const d = document.createElement("div");
    d.className = "hdr";
    d.style.background = colors[a];
    d.title = `felt by colour ${a}`;
    g.appendChild(d);
    for (let b = 0; b < S; b++) {
      const c = document.createElement("div");
      c.className = "cell" + (state.sel && state.sel[0] === a && state.sel[1] === b ? " sel" : "");
      c.style.background = cellColour(Mx[a][b]);
      c.textContent = state.tab === "force" ? Mx[a][b].toFixed(1) : Math.round(Mx[a][b]);
      c.title = `${a} → ${b}`;
      c.onclick = () => { state.sel = [a, b]; renderMatrix(); showCellEditor(); };
      c.onmouseenter = () => { state.hover = [a, b]; drawForceCurve(); };
      c.onmouseleave = () => { state.hover = null; drawForceCurve(); };
      g.appendChild(c);
    }
  }
  drawForceCurve();
}

function showCellEditor() {
  const box = $("cell-editor");
  if (!state.sel) { box.hidden = true; return; }
  const [a, b] = state.sel;
  const Mx = activeMatrix();
  const isForce = state.tab === "force";
  box.hidden = false;
  $("cell-label").textContent =
    `${isForce ? "force" : state.tab === "min" ? "min radius" : "max radius"}  ${a} → ${b}`;
  const sl = $("cell-slider"), nm = $("cell-number");
  sl.min = isForce ? -1 : 1; sl.max = isForce ? 1 : (state.tab === "min" ? 90 : 220);
  sl.step = isForce ? 0.01 : 1;
  sl.value = Mx[a][b]; nm.value = isForce ? Mx[a][b].toFixed(2) : Math.round(Mx[a][b]);
  const apply = (v) => {
    v = parseFloat(v);
    if (!isFinite(v)) return;
    Mx[a][b] = v;
    if (state.tab !== "force") beta[a][b] = Math.max(beta[a][b], alpha[a][b] * 1.05 + 1);
    engine.setMatrices(A, alpha, beta);
    sl.value = v; nm.value = isForce ? v.toFixed(2) : Math.round(v);
    renderMatrix(); refreshTheory();
  };
  sl.oninput = (e) => apply(e.target.value);
  nm.onchange = (e) => apply(e.target.value);
}

// --------------------------------------------------------------------------- //
//  physics X-ray
// --------------------------------------------------------------------------- //
function ctxOf(id) {
  const c = $(id);
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const hCss = parseFloat(c.dataset.h) || 120;
  const wWant = Math.max(60, Math.round(c.clientWidth * dpr));
  if (c.width !== wWant) {
    c.width = wWant;
    c.height = Math.max(40, Math.round(hCss * dpr));
    c.style.height = hCss + "px";
  }
  const g = c.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, c.width, c.height);
  return { g, w: c.width / dpr, h: c.height / dpr };
}

function drawForceCurve() {
  const { g, w, h } = ctxOf("c-force");
  const pad = 4;
  const pairs = state.hover ? [state.hover] : (state.sel ? [state.sel] : [[0, 0]]);
  const [a, b] = pairs[0];
  const rmax = Math.max(beta[a][b], beta[b][a]) * 1.08;
  const fmax = 1.15 * Math.max(0.2, state.repel, Math.abs(A[a][b]), Math.abs(A[b][a]));

  const X = (r) => pad + (r / rmax) * (w - 2 * pad);
  const Y = (f) => h / 2 - (f / fmax) * (h / 2 - pad);

  g.strokeStyle = "#1b232e"; g.lineWidth = 1;
  g.beginPath(); g.moveTo(pad, Y(0)); g.lineTo(w - pad, Y(0)); g.stroke();

  for (const [p, q, col, dash] of [[a, b, "#39d353", []], [b, a, "#ff5ec4", [3, 3]]]) {
    g.strokeStyle = col; g.lineWidth = 1.6; g.setLineDash(dash);
    g.beginPath();
    for (let i = 0; i <= 220; i++) {
      const r = (i / 220) * rmax;
      const f = pairForce(r, A[p][q], alpha[p][q], beta[p][q], state.repel);
      i ? g.lineTo(X(r), Y(f)) : g.moveTo(X(r), Y(f));
    }
    g.stroke(); g.setLineDash([]);
  }
  g.strokeStyle = "#42506380"; g.setLineDash([2, 3]);
  g.beginPath(); g.moveTo(X(alpha[a][b]), pad); g.lineTo(X(alpha[a][b]), h - pad); g.stroke();
  g.setLineDash([]);
  g.fillStyle = "#8b949e"; g.font = "9px ui-sans-serif";
  g.fillText(`${a}→${b} (solid)   ${b}→${a} (dashed)`, pad + 2, 10);
  g.fillText(`α=${alpha[a][b].toFixed(0)}`, X(alpha[a][b]) + 3, h - 4);
}

let lastTheory = null;

function refreshTheory() {
  const S = state.S;
  const nu = nonreciprocity(A);
  $("st-nu").textContent = nu.toFixed(3);
  $("bar-nu").style.width = `${Math.round(nu * 100)}%`;
  $("note-nu").innerHTML = nu < 1e-6
    ? "ν = 0 — Newton's third law holds. E = K + V is a Lyapunov function, so the system <b>must</b> come to rest. Only static structures are reachable."
    : `ν = ${nu.toFixed(2)} — the third law is broken. Bound clusters feel a non-zero net internal force (translation) and torque (rotation): the structures can stay alive indefinitely.`;

  const area = engine.W * engine.H;
  const rho = new Array(S).fill(engine.n / S / area);
  const gamma = -60 * Math.log(Math.max(1e-9, 1 - state.frictionFactor));
  const res = analyseStability(A, alpha, beta, state.repel, rho, gamma, state.forceFactor, { nk: 70 });
  lastTheory = res;

  $("v-gamma").textContent = gamma.toFixed(1);
  $("st-kstar").textContent = res.growth > 0 ? res.kStar.toFixed(4) : "—";
  $("st-length").textContent = res.growth > 0 ? `${res.length.toFixed(0)} px` : "—";
  $("st-growth").textContent = `${res.growth >= 0 ? "+" : ""}${res.growth.toFixed(3)} 1/s`;
  $("st-osc").textContent = res.oscGrowth > 0
    ? `present (Re λ = +${res.oscGrowth.toFixed(3)})` : "suppressed";
  $("st-speed").textContent = res.oscGrowth > 0 ? `${res.oscSpeed.toFixed(1)} px/s` : "—";

  const pill = $("hud-regime");
  pill.textContent = res.regime;
  pill.style.background = { uniform: "#8b949e", static: "#58a6ff",
    "condense+move": "#f0883e", travelling: "#ff5ec4" }[res.regime] || "#39d353";

  drawDispersion(res);
}

function drawDispersion(res) {
  const { g, w, h } = ctxOf("c-disp");
  const pad = 4;
  const ks = res.curve.map((c) => c.k);
  const re = res.curve.map((c) => c.re);
  const kmax = ks[ks.length - 1];
  const lo = Math.min(-0.2, Math.min(...re)), hi = Math.max(0.2, Math.max(...re));
  const X = (k) => pad + (k / kmax) * (w - 2 * pad);
  const Y = (v) => h - pad - ((v - lo) / (hi - lo)) * (h - 2 * pad);

  g.strokeStyle = "#1b232e"; g.beginPath();
  g.moveTo(pad, Y(0)); g.lineTo(w - pad, Y(0)); g.stroke();

  g.strokeStyle = "#58a6ff"; g.lineWidth = 1.7; g.beginPath();
  res.curve.forEach((c, i) => (i ? g.lineTo(X(c.k), Y(c.re)) : g.moveTo(X(c.k), Y(c.re))));
  g.stroke();

  res.curve.forEach((c) => {
    if (c.im > 1e-9 && c.re > 0) {
      g.fillStyle = "#ff5ec4";
      g.fillRect(X(c.k) - 1, Y(c.re) - 1, 2.5, 2.5);
    }
  });

  if (res.growth > 0) {
    g.strokeStyle = "#39d353"; g.setLineDash([3, 3]); g.beginPath();
    g.moveTo(X(res.kStar), pad); g.lineTo(X(res.kStar), h - pad); g.stroke();
    g.setLineDash([]);
  }
  g.fillStyle = "#8b949e"; g.font = "9px ui-sans-serif";
  g.fillText("Re λ(k)   (pink = oscillatory ⇒ travels)", pad + 2, 10);
}

function drawSk(sk) {
  const { g, w, h } = ctxOf("c-sk");
  const pad = 4;
  if (!sk || !sk.k.length) return;
  const kmax = sk.k[sk.k.length - 1];
  const smax = Math.max(...sk.S) || 1;
  const X = (k) => pad + (k / kmax) * (w - 2 * pad);
  const Y = (v) => h - pad - (v / smax) * (h - 2 * pad);
  g.strokeStyle = "#39d353"; g.lineWidth = 1.5; g.beginPath();
  sk.k.forEach((k, i) => (i ? g.lineTo(X(k), Y(sk.S[i])) : g.moveTo(X(k), Y(sk.S[i]))));
  g.stroke();
  if (lastTheory && lastTheory.growth > 0 && lastTheory.kStar < kmax) {
    g.strokeStyle = "#58a6ff"; g.setLineDash([3, 3]); g.beginPath();
    g.moveTo(X(lastTheory.kStar), pad); g.lineTo(X(lastTheory.kStar), h - pad); g.stroke();
    g.setLineDash([]);
  }
  g.fillStyle = "#8b949e"; g.font = "9px ui-sans-serif";
  g.fillText("measured S(k)   (blue = predicted k*)", pad + 2, 10);
}

function measure() {
  const sk = structureFactor(engine.px, engine.py, engine.n, engine.W, engine.H, 64);
  drawSk(sk);
  $("st-skpeak").textContent = isFinite(sk.peak) && sk.peak > 0
    ? `${(2 * Math.PI / sk.peak).toFixed(0)} px` : "—";

  let alTyp = 0;
  for (let a = 0; a < state.S; a++) for (let b = 0; b < state.S; b++) alTyp += alpha[a][b];
  alTyp /= state.S * state.S;

  const cs = clusterStats(engine.px, engine.py, engine.vx, engine.vy,
    engine.n, engine.W, engine.H, 1.35 * alTyp);
  $("st-speedm").textContent = engine.meanSpeed().toFixed(2);
  $("st-nclust").textContent = cs.nClusters;
  $("st-clsize").textContent = cs.meanCluster.toFixed(0);
  $("st-motility").textContent = (cs.clusterSpeed / alTyp).toFixed(4);
}

// --------------------------------------------------------------------------- //
//  rendering
// --------------------------------------------------------------------------- //
function draw() {
  const c = $("sim"), g = c.getContext("2d");
  g.globalCompositeOperation = "source-over";
  g.fillStyle = `rgba(5,7,12,${1 - state.trail})`;
  g.fillRect(0, 0, c.width, c.height);
  g.globalCompositeOperation = state.glow ? "lighter" : "source-over";

  const r = state.size;
  const n = engine.n, px = engine.px, py = engine.py, ty = engine.type;
  if (state.colourBySpeed) {
    const vx = engine.vx, vy = engine.vy;
    for (let i = 0; i < n; i++) {
      const s = Math.min(1, Math.hypot(vx[i], vy[i]) / 60);
      g.fillStyle = `hsl(${Math.round(210 - 210 * s)} 95% ${Math.round(38 + 32 * s)}%)`;
      g.beginPath(); g.arc(px[i], py[i], r, 0, 6.2832); g.fill();
    }
  } else {
    const buckets = Array.from({ length: state.S }, () => []);
    for (let i = 0; i < n; i++) buckets[ty[i]].push(i);
    for (let s = 0; s < state.S; s++) {
      g.fillStyle = colors[s];
      g.beginPath();
      for (const i of buckets[s]) {
        g.moveTo(px[i] + r, py[i]);
        g.arc(px[i], py[i], r, 0, 6.2832);
      }
      g.fill();
    }
  }
  g.globalCompositeOperation = "source-over";
}

// --------------------------------------------------------------------------- //
//  loop
// --------------------------------------------------------------------------- //
let last = performance.now(), fpsAvg = 60, frame = 0;

function loop() {
  const now = performance.now();
  const ms = now - last; last = now;
  fpsAvg = 0.9 * fpsAvg + 0.1 * (1000 / Math.max(1, ms));

  if (!state.paused) for (let s = 0; s < state.substeps; s++) engine.step(state.dt);
  draw();

  if (++frame % 12 === 0) {
    $("fps").textContent = `${fpsAvg.toFixed(0)} fps`;
    $("hud-n").textContent = `${engine.n.toLocaleString()} particles · ${state.S} species`;
  }
  if (frame % 45 === 0) measure();
  requestAnimationFrame(loop);
}

// --------------------------------------------------------------------------- //
//  controls
// --------------------------------------------------------------------------- //
function fillSelect(id, names, value) {
  const s = $(id);
  s.innerHTML = "";
  for (const [k, label] of Object.entries(names)) {
    const o = document.createElement("option");
    o.value = k; o.textContent = label;
    s.appendChild(o);
  }
  s.value = value;
}

function bindRange(id, key, fmt = (v) => v, after = () => {}) {
  const el = $(id), out = $(id.replace("s-", "v-"));
  el.value = state[key];
  if (out) out.textContent = fmt(state[key]);
  el.oninput = () => {
    state[key] = parseFloat(el.value);
    if (out) out.textContent = fmt(state[key]);
    after();
  };
}

function resize() {
  const c = $("sim");
  const w = window.innerWidth, h = window.innerHeight;
  if (c.width === w && c.height === h) return;
  c.width = w; c.height = h;
  if (engine) engine.resize(w, h);
}

function init() {
  resize();
  fillSelect("sel-matrix", MATRIX_NAMES, state.matrixPreset);
  fillSelect("sel-palette", PALETTE_NAMES, state.palettePreset);
  fillSelect("sel-dist", DISTRIBUTION_NAMES, state.distPreset);
  $("n-min-lo").value = state.minRange[0]; $("n-min-hi").value = state.minRange[1];
  $("n-max-lo").value = state.maxRange[0]; $("n-max-hi").value = state.maxRange[1];
  $("v-minr").textContent = state.minRange.join("–");
  $("v-maxr").textContent = state.maxRange.join("–");

  buildUniverse();

  $("sel-matrix").onchange = (e) => { state.matrixPreset = e.target.value; buildUniverse({ radii: false, palette: false, spawn: false }); };
  $("btn-matrix").onclick = () => buildUniverse({ radii: false, palette: false, spawn: false });
  $("sel-palette").onchange = (e) => { state.palettePreset = e.target.value; buildUniverse({ matrix: false, radii: false, spawn: false }); };
  $("btn-palette").onclick = () => buildUniverse({ matrix: false, radii: false, spawn: false });
  $("sel-dist").onchange = (e) => { state.distPreset = e.target.value; engine.respawn(state.distPreset); };
  $("btn-dist").onclick = () => engine.respawn(state.distPreset);
  $("btn-randomise").onclick = () => buildUniverse();

  document.querySelectorAll(".tab").forEach((t) => {
    t.onclick = () => {
      document.querySelectorAll(".tab").forEach((x) => x.classList.remove("active"));
      t.classList.add("active");
      state.tab = t.dataset.tab;
      renderMatrix(); showCellEditor();
    };
  });

  $("btn-symmetrise").onclick = () => {
    for (let a = 0; a < state.S; a++) for (let b = a + 1; b < state.S; b++) {
      const m = 0.5 * (A[a][b] + A[b][a]); A[a][b] = A[b][a] = Math.round(m * 100) / 100;
      const p = 0.5 * (alpha[a][b] + alpha[b][a]); alpha[a][b] = alpha[b][a] = p;
      const q = 0.5 * (beta[a][b] + beta[b][a]); beta[a][b] = beta[b][a] = q;
    }
    engine.setMatrices(A, alpha, beta);
    renderMatrix(); refreshTheory();
  };

  bindRange("s-species", "S", (v) => v, () => { state.S = Math.round(state.S); buildUniverse(); });
  bindRange("s-count", "N", (v) => Math.round(v).toLocaleString(), () => engine.setCount(Math.round(state.N), state.distPreset));
  bindRange("s-force", "forceFactor", (v) => v.toFixed(2), () => { engine.forceFactor = state.forceFactor; refreshTheory(); });
  bindRange("s-friction", "frictionFactor", (v) => v.toFixed(3), () => { engine.frictionFactor = state.frictionFactor; refreshTheory(); });
  bindRange("s-repel", "repel", (v) => v.toFixed(1), () => { engine.repel = state.repel; refreshTheory(); drawForceCurve(); });
  bindRange("s-dt", "dt", (v) => v.toFixed(3));
  bindRange("s-sub", "substeps", (v) => Math.round(v));
  bindRange("s-size", "size", (v) => v.toFixed(2));
  bindRange("s-trail", "trail", (v) => v.toFixed(2));

  $("sel-wall").value = state.wallMode;
  $("sel-wall").onchange = (e) => { state.wallMode = e.target.value; engine.wallMode = state.wallMode; };
  $("c-glow").onchange = (e) => (state.glow = e.target.checked);
  $("c-vel").onchange = (e) => (state.colourBySpeed = e.target.checked);
  $("c-symradii").onchange = (e) => (state.symRadii = e.target.checked);
  for (const [id, idx, key] of [["n-min-lo", 0, "minRange"], ["n-min-hi", 1, "minRange"],
                                ["n-max-lo", 0, "maxRange"], ["n-max-hi", 1, "maxRange"]]) {
    $(id).onchange = (e) => {
      state[key][idx] = parseFloat(e.target.value);
      $("v-minr").textContent = state.minRange.join("–");
      $("v-maxr").textContent = state.maxRange.join("–");
    };
  }

  $("btn-pause").onclick = () => {
    state.paused = !state.paused;
    $("btn-pause").textContent = state.paused ? "▶ Resume" : "⏸ Pause";
  };
  $("btn-reset").onclick = () => engine.respawn(state.distPreset);

  $("btn-export").onclick = () => {
    const blob = new Blob([JSON.stringify({
      version: 1, species: state.S, matrices: { forces: A, minRadius: alpha, maxRadius: beta },
      settings: {
        numParticles: state.N, forceFactor: state.forceFactor,
        frictionFactor: state.frictionFactor, repel: state.repel,
        dt: state.dt, wallMode: state.wallMode,
      },
      colors,
    }, null, 2)], { type: "application/json" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "particle-life-preset.json";
    a.click();
  };
  $("btn-import").onclick = () => $("file-import").click();
  $("file-import").onchange = async (e) => {
    const f = e.target.files[0];
    if (!f) return;
    const p = JSON.parse(await f.text());
    state.S = p.species ?? p.matrices.forces.length;
    A = p.matrices.forces; alpha = p.matrices.minRadius; beta = p.matrices.maxRadius;
    colors = p.colors || PALETTES.rainbow(state.S);
    Object.assign(state, p.settings || {});
    $("s-species").value = state.S;
    engine.S = state.S;
    engine.setMatrices(A, alpha, beta);
    engine.setCount(state.N, state.distPreset);
    renderMatrix(); refreshTheory();
  };

  $("toggle-left").onclick = () => $("left").classList.toggle("hidden");
  $("toggle-right").onclick = () => $("right").classList.toggle("hidden");

  // pointer brush
  const c = $("sim");
  let down = 0;
  c.onpointerdown = (e) => { down = e.button === 2 ? -1 : 1; };
  window.onpointerup = () => (down = 0);
  c.oncontextmenu = (e) => e.preventDefault();
  c.onpointermove = (e) => {
    if (!down) return;
    engine.brush(e.clientX, e.clientY, 140, down * 30);
  };

  window.addEventListener("resize", resize);
  window.addEventListener("keydown", (e) => {
    if (e.key === " ") { e.preventDefault(); $("btn-pause").click(); }
    if (e.key === "r") engine.respawn(state.distPreset);
    if (e.key === "n") buildUniverse();
  });

  // expose for automated checks
  window.__pl = { state, get engine() { return engine; }, get A() { return A; },
                  get alpha() { return alpha; }, get beta() { return beta; },
                  get theory() { return lastTheory; }, measure, refreshTheory };

  requestAnimationFrame(loop);
}

init();
