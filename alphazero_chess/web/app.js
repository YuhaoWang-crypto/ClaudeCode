// Game interface: board, input, game flow and the live "what the AI is thinking" panel.
import { Chess } from "https://cdn.jsdelivr.net/npm/chess.js@1.4.0/dist/esm/chess.js";
import { AlphaZeroNet, search } from "./engine.js";

const tf = window.tf;
const $ = (id) => document.getElementById(id);
const START = new Chess().fen();
const GLYPH = { k: "♚", q: "♛", r: "♜", b: "♝", n: "♞", p: "♟" };
const FILES = "abcdefgh";
const STORE_KEY = "az-chess-v1";

const S = {
  net: null,
  game: new Chess(),
  history: [],          // uci moves from START
  human: "w",
  sims: 200,
  flipped: false,
  selected: null,
  thinking: false,
  token: 0,             // bumped to cancel a running search
  analysis: null,       // last search summary
  analysisFor: null,    // "ai" | "hint"
  drag: null,
  pendingPromo: null,
};

// ------------------------------------------------------------ persistence
function save() {
  try { localStorage.setItem(STORE_KEY, JSON.stringify({ history: S.history, human: S.human, sims: S.sims })); } catch {}
}
function restore(data) {
  let d = data && data.history ? data : null;
  if (!d) { try { d = JSON.parse(localStorage.getItem(STORE_KEY) || "null"); } catch {} }
  if (!d) return;
  if (d.human === "w" || d.human === "b") S.human = d.human;
  if ([1, 64, 200, 800].includes(d.sims)) { S.sims = d.sims; S.simsSaved = true; }
  for (const u of d.history || []) {
    try { S.game.move(toMove(u)); S.history.push(u); } catch { break; }
  }
}
try { window.claude?.hot?.snapshot?.(() => ({ history: S.history, human: S.human, sims: S.sims })); } catch {}

const toMove = (u) => ({ from: u.slice(0, 2), to: u.slice(2, 4), promotion: u[4] });

// ------------------------------------------------------------ board rendering
const whiteBottom = () => (S.human === "w") !== S.flipped;
function squareAt(row, col) {
  return whiteBottom() ? FILES[col] + (8 - row) : FILES[7 - col] + (row + 1);
}
function squareCenter(sq) {
  const f = FILES.indexOf(sq[0]), r = Number(sq[1]);
  return whiteBottom() ? [f + 0.5, 8 - r + 0.5] : [7 - f + 0.5, r - 1 + 0.5];
}

function renderBoard() {
  const board = $("board");
  const g = S.game;
  const last = S.history.length ? S.history[S.history.length - 1] : null;
  const targets = new Map();
  if (S.selected) for (const m of g.moves({ square: S.selected, verbose: true })) targets.set(m.to, !!m.captured);
  let checkSq = null;
  if (g.inCheck()) {
    for (const row of g.board()) for (const p of row) if (p && p.type === "k" && p.color === g.turn()) checkSq = p.square;
  }
  const frag = document.createDocumentFragment();
  for (let row = 0; row < 8; row++) {
    for (let col = 0; col < 8; col++) {
      const sq = squareAt(row, col);
      const el = document.createElement("div");
      const light = (FILES.indexOf(sq[0]) + Number(sq[1])) % 2 === 1;
      el.className = "sq " + (light ? "l" : "d");
      el.dataset.sq = sq;
      if (last && (sq === last.slice(0, 2) || sq === last.slice(2, 4))) el.classList.add("last");
      if (sq === S.selected) el.classList.add("sel");
      if (sq === checkSq) el.classList.add("check");
      if (targets.has(sq)) { el.classList.add("tgt"); if (targets.get(sq)) el.classList.add("cap"); }
      if (col === 0) el.insertAdjacentHTML("beforeend", `<span class="coord r">${sq[1]}</span>`);
      if (row === 7) el.insertAdjacentHTML("beforeend", `<span class="coord f">${sq[0]}</span>`);
      const p = g.get(sq);
      if (p) {
        const pc = document.createElement("span");
        pc.className = "pc " + p.color + (S.drag && S.drag.from === sq ? " dragging" : "");
        pc.textContent = GLYPH[p.type] + "︎";
        el.appendChild(pc);
      }
      frag.appendChild(el);
    }
  }
  board.replaceChildren(frag);
  $("evalbar").classList.toggle("flipped", !whiteBottom());
  renderArrows();
}

function renderArrows() {
  const svg = $("arrows");
  const a = S.analysis;
  if (!a || !a.stats.length) { svg.innerHTML = ""; return; }
  const total = Math.max(1, a.stats.reduce((s, x) => s + x.n, 0));
  const shown = a.stats.slice(0, 3).filter((x, i) => i === 0 || x.n / total > 0.08 || (a.total <= 1 && x.p > 0.12));
  let out = `<defs><marker id="ah" viewBox="0 0 10 10" refX="5" refY="5" markerWidth="2.6" markerHeight="2.6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="var(--accent)"/></marker></defs>`;
  shown.forEach((m, i) => {
    const [x1, y1] = squareCenter(m.uci.slice(0, 2));
    const [x2, y2] = squareCenter(m.uci.slice(2, 4));
    const len = Math.hypot(x2 - x1, y2 - y1), k = (len - 0.32) / len;
    const share = a.total <= 1 ? m.p : m.n / total;
    const w = 0.07 + 0.13 * share;
    out += `<line x1="${x1}" y1="${y1}" x2="${x1 + (x2 - x1) * k}" y2="${y1 + (y2 - y1) * k}" stroke="var(--accent)" stroke-width="${w}" stroke-linecap="round" opacity="${i === 0 ? 0.85 : 0.45}" marker-end="url(#ah)"/>`;
  });
  svg.innerHTML = out;
}

// ------------------------------------------------------------ side panel
const pct = (v) => Math.round(((v + 1) / 2) * 100);

function renderAnalysis() {
  const a = S.analysis;
  const body = $("candbody");
  if (!a) {
    body.innerHTML = "";
    setWin(0);
    $("an-title").textContent = "AI 的思考";
    $("an-note").hidden = false;
    return;
  }
  $("an-note").hidden = true;
  const who = a.turn === "w" ? "白方" : "黑方";
  $("an-title").textContent = S.analysisFor === "hint"
    ? `给你的提示（${who}）· ${a.sims}/${a.total} 次模拟`
    : `AI（${who}）的思考 · ${a.sims}/${a.total} 次模拟`;
  const whiteV = a.turn === "w" ? a.value : -a.value;
  setWin(whiteV);
  const total = Math.max(1, a.stats.reduce((s, x) => s + x.n, 0));
  body.innerHTML = a.stats.slice(0, 6).map((m, i) => {
    const share = a.total <= 1 ? m.p : m.n / total;
    const q = m.n ? `${pct(m.q)}%` : "–";
    return `<tr class="${i === 0 ? "top" : ""}"><td class="mv">${m.san}</td>` +
      `<td class="vis"><div class="vbar"><i style="width:${(share * 100).toFixed(1)}%"></i></div></td>` +
      `<td>${(m.p * 100).toFixed(1)}%</td><td>${q}</td></tr>`;
  }).join("");
}

function setWin(whiteV) {
  const w = pct(whiteV);
  $("win-w").textContent = `白 ${w}%`;
  $("win-b").textContent = `黑 ${100 - w}%`;
  $("winfill").style.width = w + "%";
  $("evalfill").style.height = w + "%";
}

function renderMoves() {
  const sans = S.game.history();
  const ol = $("moves");
  let html = "";
  for (let i = 0; i < sans.length; i += 2) {
    const cw = i === sans.length - 1 ? "cur" : "", cb = i + 1 === sans.length - 1 ? "cur" : "";
    html += `<li><span class="n">${i / 2 + 1}.</span><span class="${cw}">${sans[i]}</span><span class="${cb}">${sans[i + 1] || ""}</span></li>`;
  }
  ol.innerHTML = html || `<li><span class="n"></span><span style="color:var(--muted)">还没有走棋</span><span></span></li>`;
  ol.scrollTop = ol.scrollHeight;
}

function gameOverText() {
  const g = S.game;
  if (g.isCheckmate()) return g.turn() === S.human ? "将死 · AI 获胜" : "将死 · 你赢了！";
  if (g.isStalemate()) return "逼和 · 和棋";
  if (g.isInsufficientMaterial()) return "子力不足 · 和棋";
  if (g.isThreefoldRepetition()) return "三次重复局面 · 和棋";
  if (g.isDrawByFiftyMoves()) return "50 步规则 · 和棋";
  return null;
}

function renderStatus(extra) {
  const g = S.game, turn = g.turn();
  const over = gameOverText();
  let text;
  if (over) text = over;
  else if (S.thinking) text = "AI 思考中";
  else text = turn === S.human ? "轮到你走" : "轮到 AI";
  if (!over && g.inCheck()) text += " · 将军！";
  $("status").innerHTML = `<span class="dot ${turn}"></span><span>${text}</span>` +
    (extra ? `<span class="think">${extra}</span>` : "");
}

function renderControls() {
  for (const b of document.querySelectorAll("#side button")) b.setAttribute("aria-pressed", String(b.dataset.side === S.human));
  for (const b of document.querySelectorAll("#level button")) b.setAttribute("aria-pressed", String(Number(b.dataset.sims) === S.sims));
  const ready = !!S.net;
  $("undo").disabled = !ready || S.history.length === 0;
  $("hint").disabled = !ready || S.thinking || S.game.turn() !== S.human || !!gameOverText();
}

function renderAll(extra) { renderBoard(); renderAnalysis(); renderMoves(); renderStatus(extra); renderControls(); }

// ------------------------------------------------------------ game flow
function canHumanMove() { return S.net && !S.thinking && !S.pendingPromo && S.game.turn() === S.human && !gameOverText(); }

function play(uci) {
  S.game.move(toMove(uci));
  S.history.push(uci);
  S.selected = null;
  save();
}

function humanMove(from, to) {
  const opts = S.game.moves({ square: from, verbose: true }).filter((m) => m.to === to);
  if (!opts.length) return false;
  if (opts.some((m) => m.promotion)) { askPromotion(from, to); return true; }
  play(from + to);
  S.analysis = null;
  renderAll();
  maybeAiMove();
  return true;
}

function askPromotion(from, to) {
  S.pendingPromo = { from, to };
  const box = $("promobox");
  box.innerHTML = "";
  for (const t of ["q", "r", "b", "n"]) {
    const b = document.createElement("button");
    b.type = "button";
    b.id = "promo-" + t;
    b.setAttribute("aria-label", { q: "后", r: "车", b: "象", n: "马" }[t]);
    b.innerHTML = `<span class="pc ${S.human}">${GLYPH[t]}︎</span>`;
    b.onclick = () => {
      $("promo").hidden = true;
      S.pendingPromo = null;
      play(from + to + t);
      S.analysis = null;
      renderAll();
      maybeAiMove();
    };
    box.appendChild(b);
  }
  $("promo").hidden = false;
  box.querySelector("button").focus();
}

const nextFrame = () => new Promise((r) => setTimeout(r, 30));

async function runSearch(kind) {
  const token = ++S.token;
  S.thinking = kind === "ai";
  S.analysisFor = kind;
  renderAll(kind === "ai" ? "" : "");
  await nextFrame();
  const res = await search({
    net: S.net, Chess, startFen: START, moves: S.history.slice(), sims: S.sims, batch: 8, cpuct: 1.0,
    shouldStop: () => token !== S.token,
    onProgress: (info) => {
      if (token !== S.token) return;
      S.analysis = info;
      renderAnalysis(); renderArrows();
      if (kind === "ai") renderStatus(`${info.sims}/${info.total}`);
    },
  });
  if (token !== S.token) return null;
  S.analysis = res;
  return res;
}

async function maybeAiMove() {
  if (!S.net || gameOverText() || S.game.turn() === S.human) return;
  const res = await runSearch("ai");
  if (!res) return;
  S.thinking = false;
  play(res.best.uci);
  renderAll();
}

async function hint() {
  if (!canHumanMove()) return;
  const res = await runSearch("hint");
  if (!res) return;
  renderAll();
  $("an-note").hidden = false;
  $("an-note").textContent = `建议：${res.best.san}（箭头所示）`;
}

function cancelSearch() { S.token++; S.thinking = false; }

function newGame(human) {
  cancelSearch();
  if (human) S.human = human;
  S.game.reset();
  S.history = [];
  S.selected = null;
  S.analysis = null;
  S.flipped = false;
  $("promo").hidden = true; S.pendingPromo = null;
  save();
  renderAll();
  maybeAiMove();
}

function undo() {
  if (!S.history.length) return;
  cancelSearch();
  do { S.game.undo(); S.history.pop(); } while (S.history.length && S.game.turn() !== S.human);
  S.selected = null;
  S.analysis = null;
  $("promo").hidden = true; S.pendingPromo = null;
  save();
  renderAll();
  maybeAiMove();
}

// ------------------------------------------------------------ input (click + drag)
function squareFromPoint(x, y) {
  const r = $("board").getBoundingClientRect();
  const col = Math.floor(((x - r.left) / r.width) * 8), row = Math.floor(((y - r.top) / r.height) * 8);
  if (col < 0 || col > 7 || row < 0 || row > 7) return null;
  return squareAt(row, col);
}

function setupInput() {
  const board = $("board");
  let ghost = null;
  board.addEventListener("pointerdown", (e) => {
    if (!canHumanMove()) return;
    const sq = squareFromPoint(e.clientX, e.clientY);
    if (!sq) return;
    if (S.selected && S.selected !== sq && humanMove(S.selected, sq)) return;
    const p = S.game.get(sq);
    if (p && p.color === S.human) {
      S.selected = sq;
      S.drag = { from: sq, x0: e.clientX, y0: e.clientY, moved: false };
      board.setPointerCapture(e.pointerId);
    } else {
      S.selected = null;
    }
    renderBoard();
  });
  board.addEventListener("pointermove", (e) => {
    const d = S.drag;
    if (!d) return;
    if (!d.moved && Math.hypot(e.clientX - d.x0, e.clientY - d.y0) > 6) {
      d.moved = true;
      const p = S.game.get(d.from);
      ghost = document.createElement("span");
      ghost.className = "pc ghost " + p.color;
      ghost.textContent = GLYPH[p.type] + "︎";
      const size = board.getBoundingClientRect().width / 8;
      ghost.style.fontSize = size * 0.82 + "px";
      document.body.appendChild(ghost);
      renderBoard();
    }
    if (ghost) { ghost.style.left = e.clientX + "px"; ghost.style.top = e.clientY + "px"; }
  });
  const end = (e) => {
    const d = S.drag;
    if (!d) return;
    S.drag = null;
    if (ghost) { ghost.remove(); ghost = null; }
    if (d.moved && e.type === "pointerup") {
      const to = squareFromPoint(e.clientX, e.clientY);
      if (to && to !== d.from && humanMove(d.from, to)) return;
    }
    renderBoard();
  };
  board.addEventListener("pointerup", end);
  board.addEventListener("pointercancel", end);
}

function setupControls() {
  for (const b of document.querySelectorAll("#side button")) b.onclick = () => newGame(b.dataset.side);
  for (const b of document.querySelectorAll("#level button")) {
    b.onclick = () => { S.sims = Number(b.dataset.sims); save(); renderControls(); };
  }
  $("new").onclick = () => newGame();
  $("undo").onclick = undo;
  $("hint").onclick = hint;
  $("flip").onclick = () => { S.flipped = !S.flipped; renderBoard(); };
  $("pgn").onclick = async () => {
    const btn = $("pgn");
    try { await navigator.clipboard.writeText(S.game.pgn() || "*"); btn.textContent = "已复制"; }
    catch { btn.textContent = "复制失败"; }
    setTimeout(() => { btn.textContent = "复制 PGN"; }, 1500);
  };
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && S.selected) { S.selected = null; renderBoard(); }
  });
}

// ------------------------------------------------------------ boot
async function pickBackend() {
  for (const b of ["webgl", "cpu"]) {
    try { if (await tf.setBackend(b)) { await tf.ready(); return b; } } catch {}
  }
  return tf.getBackend();
}

async function start(data) {
  restore(data);
  setupInput();
  setupControls();
  renderAll();
  try {
    const backend = await pickBackend();
    const base = new URL("model/", location.href);
    S.net = await AlphaZeroNet.load(tf, base, (f) => {
      $("loadbar").style.width = (f * 100).toFixed(0) + "%";
      $("loadsub").textContent = `${(f * 17.2).toFixed(1)} / 17.2 MB`;
    });
    $("loadmsg").textContent = "正在预热网络…";
    $("loadsub").textContent = backend === "webgl" ? "使用 GPU（WebGL）计算" : "使用 CPU 计算（会慢一些）";
    S.net.padBatches = backend === "webgl";
    for (const n of [1, 2, 4, 8]) await S.net.predict(Array(n).fill(START));   // compile shaders
    const t0 = performance.now();
    await S.net.predict(Array(8).fill(START));
    const perSec = 8000 / Math.max(1, performance.now() - t0);
    if (!S.simsSaved && 200 / perSec > 12) S.sims = 64 / perSec > 12 ? 1 : 64;
    $("speed").textContent = `本设备约每秒 ${perSec >= 10 ? Math.round(perSec) : perSec.toFixed(1)} 次模拟（${backend === "webgl" ? "GPU" : "CPU"}）。模拟越多越强，也越慢。`;
    $("loader").hidden = true;
    renderAll();
    maybeAiMove();
  } catch (err) {
    console.error(err);
    $("loadmsg").textContent = "神经网络没能加载";
    $("loadsub").textContent = "请检查网络连接后刷新页面。如果是本地打开，请用 python3 -m http.server 启动后访问。";
  }
}

const hot = window.claude?.hot;
if (hot?.ready) hot.ready(start); else start(hot?.data ?? {});
