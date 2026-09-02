/* VC-FF · 虚拟细胞力场 —— 浏览器端引擎 (与 vcff/ Python 包同一套组合律) */
(function () {
"use strict";

const D = JSON.parse(document.getElementById("vcff-data").textContent);
const $ = (id) => document.getElementById(id);
const NORM = (s) => String(s).toLowerCase().replace(/[^a-z0-9]/g, "");

/* ---------- 索引 ---------- */
const G2I = Object.create(null);
D.tox.genes.forEach((g, i) => { G2I[NORM(g)] = i; });

const SUPP = {
  paclitaxel: ["TUBB", false], docetaxel: ["TUBB", false],
  vincristine: ["TUBB", false], vinblastine: ["TUBB", false],
  doxorubicin: ["TOP2A", false], daunorubicin: ["TOP2A", false],
  etoposide: ["TOP2A", false], methotrexate: ["DHFR", false],
  gemcitabine: ["RRM1", false], hydroxyurea: ["RRM1", false],
  sorafenib: ["BRAF", true], sunitinibmalate: ["KDR", true],
  dasatinib: ["ABL1", true], imatinib: ["ABL1", true],
  nilotinib: ["ABL1", true], midostaurin: ["FLT3", true],
  cabozantinib: ["MET", true], lenvatinib: ["KDR", true]
};

const CMBI = Object.create(null);
D.cmb.pairs.forEach((p, i) => {
  const [a, b] = p.split("+");
  CMBI[[NORM(a), NORM(b)].sort().join("|")] = i;
});
const EPI_MED = (function () {
  const v = D.cmb.epi.slice().sort((x, y) => x - y);
  const n = v.length;
  return n % 2 ? v[(n - 1) / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
})();
const EPI_IQR = (function () {
  const v = D.cmb.epi.slice().sort((x, y) => x - y);
  return [v[Math.floor(v.length / 4)], v[Math.floor(3 * v.length / 4)]];
})();
const FRAC_NOISE = D.cmb.non.filter((v) => v !== null && v > 1).length / D.cmb.non.length;
const NON_MED = (function () {
  const v = D.cmb.non.filter((x) => x !== null).sort((a, b) => a - b);
  const n = v.length;
  return n % 2 ? v[(n - 1) / 2] : (v[n / 2 - 1] + v[n / 2]) / 2;
})();

const CONTEXTS = {
  hepg2: "HepG2 · 肝细胞癌", rpe1: "RPE1 · 正常视网膜色素上皮",
  k562: "K562 · 慢性髓系白血病", jurkat: "Jurkat · T 细胞急性淋巴细胞白血病"
};
const ESS_LINES = ["hepg2", "rpe1"];
const LINE_LABEL = { hepg2: "HepG2（肝癌）", rpe1: "RPE1（正常上皮）" };

/* ---------- 组合律 (physics) ---------- */
const LAM = Math.LN2, TREF = 72;
const occupancy = (c, ic, h) => {
  if (!(c > 0)) return 0;
  if (!(ic > 0)) return 1;
  const x = Math.pow(c / ic, h);
  return x / (1 + x);
};
const killPressure = (dep, th) => Math.max(0, -dep) * th;
const logViab = (ps, hrs, g) =>
  -LAM * ps.reduce((a, b) => a + b, 0) * (hrs / TREF) * g;
const viab = (lv) => Math.exp(lv);
const epiGamma = (epi) => {
  if (!(EPI_MED > 0)) return 1;
  return Math.min(2, Math.max(0.5, 1 + 0.5 * (epi / EPI_MED - 1)));
};

function crossing(ms, vs, level) {
  for (let i = 1; i < ms.length; i++) {
    if (vs[i] <= level && level <= vs[i - 1]) {
      const v0 = vs[i - 1], v1 = vs[i];
      if (v0 === v1) return ms[i];
      const f = (v0 - level) / (v0 - v1);
      return Math.exp(Math.log(ms[i - 1]) + f * (Math.log(ms[i]) - Math.log(ms[i - 1])));
    }
  }
  return null;
}

/* ---------- kernel 查询 ---------- */
function resolveEss(name, target) {
  let via = null, supp = false, multi = false;
  const key = target ? NORM(target) : NORM(name);
  let idx = G2I[key];
  if (idx === undefined && !target) {
    const q = NORM(name);
    let t = D.tox.drug2target[q];
    if (t === undefined && SUPP[q]) { t = SUPP[q][0]; multi = SUPP[q][1]; supp = true; }
    if (t !== undefined) { idx = G2I[NORM(t)]; via = t; }
  }
  if (idx === undefined) return null;
  const h = D.tox.hepg2[idx], r = D.tox.rpe1[idx];
  let klass;
  if (h < -1 && r > -0.5) klass = "HepG2 选择性必需（肝特异脆弱点）";
  else if (h < -1 && r < -1) klass = "两系共同必需（广谱毒性风险）";
  else if (r < -1 && h > -0.5) klass = "RPE1 选择性必需（正常组织脆弱点）";
  else klass = "非必需（两系依赖性均弱）";
  return { gene: D.tox.genes[idx], dep: { hepg2: h, rpe1: r }, klass,
           via, supp, multi };
}
const getPhe = (drug) => {
  const k = D.phe.drugs.find((x) => NORM(x) === NORM(drug));
  return k ? Object.assign({ drug: k }, D.phe.per[k]) : null;
};
const getPrt = (gene, line) => D.prt.combos[gene + "|" + line] || null;
const getEpi = (a, b) => {
  const i = CMBI[[NORM(a), NORM(b)].sort().join("|")];
  if (i === undefined) return null;
  return { pair: D.cmb.pairs[i], non: D.cmb.non[i], epi: D.cmb.epi[i],
           rank: D.cmb.rank[i], ncells: D.cmb.ncells[i] };
};
const getTwn = (line) => {
  const r = D.twn.per[line];
  if (!r) return null;
  const verdict = r.frac >= 0.8 ? "GO · 高置信"
                : r.frac >= 0.7 ? "CAUTION · 需验证" : "GAP · 数据不足";
  return Object.assign({}, r, { verdict,
    note: String(D.twn.rel[line]).replace(/<[^>]+>/g, "") });
};

/* ---------- 引擎 ---------- */
function evaluate(spec) {
  const led = { claims: [], notes: [] };
  const claim = (o) => { led.claims.push(o); };
  const traces = spec.compounds.map((c) => {
    const ic = (c.ic50 === null || c.ic50 === undefined || c.ic50 === "")
      ? null : Number(c.ic50);
    const effIc = ic === null ? 1.0 : ic;
    const th = occupancy(Number(c.conc), effIc, Number(c.hill) || 1);
    const t = { name: c.name, target: c.target || null, via: "", theta: th,
                conc: Number(c.conc), ic50: effIc, placeholder: ic === null,
                hill: Number(c.hill) || 1, dep: null, klass: "",
                press: null, pwSrc: "none", kernels: [], warn: [] };
    if (ic === null) {
      t.warn.push("IC50 未提供，使用 1.0 µM 占位。效价是你的输入参数，引擎不预测效价；" +
                  "占据率与由它推出的一切数值都随该占位值漂移。");
    }
    const hit = resolveEss(c.name, c.target);
    if (!hit) {
      t.via = "未命中";
      t.warn.push("「" + c.name + "」既不在 17,787 基因表中，也无法映射到已知化合物靶点；" +
                  "请直接填靶基因符号。该化合物不进入毒性主轴。");
      return t;
    }
    t.target = hit.gene;
    t.via = hit.via && hit.supp ? "补充映射（人工整理）→ " + hit.via
          : hit.via ? "Tahoe 化合物→靶点映射 → " + hit.via
          : c.target ? "显式指定靶点" : "基因符号直接命中";
    if (hit.multi) {
      t.warn.push(c.name + " 是多靶点化合物，这里只按主靶点 " + hit.gene +
                  " 建模。单靶点近似是本次装配最弱的一环。");
    }
    t.dep = hit.dep; t.klass = hit.klass; t.kernels.push("VC-TOX");
    t.press = {};
    ESS_LINES.forEach((ln) => { t.press[ln] = killPressure(hit.dep[ln], th); });
    if (getPhe(c.name)) { t.pwSrc = "VC-PHE"; t.kernels.push("VC-PHE"); }
    else if (getPrt(hit.gene, spec.context)) { t.pwSrc = "VC-PRT"; t.kernels.push("VC-PRT"); }
    else {
      t.warn.push("通路层未覆盖：不在 PhenoMap 的 8 个化合物中，靶点 " + hit.gene +
                  " 也不在 PerturbLens 的 6 个基因中。本化合物只贡献必需性项。");
    }
    return t;
  });

  /* 交叉项 */
  let gamma = 1, epiHit = null;
  const tg = traces.filter((t) => t.target).map((t) => t.target);
  outer: for (let i = 0; i < tg.length; i++) {
    for (let j = i + 1; j < tg.length; j++) {
      const h = getEpi(tg[i], tg[j]);
      if (h) { epiHit = h; gamma = epiGamma(h.epi); break outer; }
    }
  }

  const hasEss = traces.some((t) => t.dep);
  const hasPw = traces.some((t) => t.pwSrc !== "none");
  const layers = [];
  if (hasEss) layers.push("必需性骨干 VC-TOX");
  if (hasPw) layers.push("通路/表型层 VC-PHE·VC-PRT");
  if (epiHit) layers.push("实测上位性交叉项 VC-CMB");
  const tier = layers.length;
  let tierLabel = tier === 0 ? "T0 · 无任何层可装配（引擎拒绝定量）"
                             : "T" + tier + " · " + layers.join(" + ");
  if (tier > 0 && !hasEss) tierLabel += "　[无必需性骨干 → 不出剂量-响应与选择性]";

  const out = { spec, traces, tier, tierLabel, gamma, epiHit, led,
                hasEss, hasPw, viability: null, curve: null,
                selectivity: null, pathway: null, combination: null,
                confidence: null };

  /* 伪存活率 */
  if (!hasEss) {
    claim({ ev: "UNSUPPORTED", ttl: "伪存活率 / 剂量-响应", val: "不可用",
      basis: "没有任何化合物解析到 Chronos 依赖性数据",
      caveat: "引擎在这里返回「不可用」而不是 1.0。「算不出杀伤」与「无毒」是两回事。" });
  } else {
    const v = {};
    ESS_LINES.forEach((ln) => {
      v[ln] = viab(logViab(traces.map((t) => (t.press ? t.press[ln] : 0)),
                           spec.exposure, gamma));
    });
    out.viability = v;
    claim({ ev: "MODELED", ttl: "规格剂量下的伪存活率",
      val: ESS_LINES.map((ln) => LINE_LABEL[ln] + " " + v[ln].toFixed(3)).join("　/　"),
      basis: "Chronos gene effect（✅ 实测）经 Hill 占据率 + 指数杀伤（⚠️ 假设）组合",
      caveat: "λ 为约定锚定（ln2），未拟合任何存活率数据。绝对百分比无意义；" +
              "只有系间比值与剂量位移可解释。" });

    /* 剂量扫描 */
    const ms = [], N = 121, lo = Math.log(1e-3), hi = Math.log(1e3);
    for (let i = 0; i < N; i++) ms.push(Math.exp(lo + (hi - lo) * i / (N - 1)));
    const curves = {};
    ESS_LINES.forEach((ln) => {
      curves[ln] = ms.map((m) => viab(logViab(traces.map((t) =>
        t.dep ? killPressure(t.dep[ln], occupancy(t.conc * m, t.ic50, t.hill)) : 0
      ), spec.exposure, gamma)));
    });
    out.curve = { ms, curves };

    /* 选择性 */
    const byLevel = {};
    let iso = null, isoLv = null, mT = null, mN = null;
    [0.5, 0.3, 0.2, 0.1].forEach((lv) => {
      const a = crossing(ms, curves.hepg2, 1 - lv);
      const b = crossing(ms, curves.rpe1, 1 - lv);
      byLevel[Math.round(lv * 100) + "%"] = (a && b) ? b / a : null;
      if (a && b && iso === null) { iso = b / a; isoLv = lv; mT = a; mN = b; }
    });
    const m50t = crossing(ms, curves.hepg2, 0.5);
    const m90n = crossing(ms, curves.rpe1, 0.9);
    const cons = (m50t && m90n) ? m90n / m50t : null;

    const floor = {}, satP = {};
    ESS_LINES.forEach((ln) => {
      const p = traces.filter((t) => t.dep).map((t) => Math.max(0, -t.dep[ln]));
      satP[ln] = p.reduce((a, b) => a + b, 0);
      floor[ln] = viab(logViab(p, spec.exposure, gamma));
    });
    const ratio = satP.rpe1 > 0 ? satP.hepg2 / satP.rpe1 : null;

    let interp, why = null;
    if (iso !== null) {
      interp = iso > 1.2
        ? "等效剂量比 " + iso.toFixed(2) + "× —— 正常系要 " + iso.toFixed(2) +
          " 倍剂量才达到同样的 " + Math.round(isoLv * 100) + "% 抑制，存在选择性"
        : iso < 0.83
        ? "等效剂量比 " + iso.toFixed(2) + "× < 1 —— 正常系比肝癌系更敏感，选择性方向是反的"
        : "等效剂量比 " + iso.toFixed(2) + "× ≈ 1 —— 两系敏感性基本相同，无选择性（广谱杀伤）";
    } else if (floor.hepg2 > 0.5) {
      interp = "等效剂量比无法定义";
      why = "即使占据率饱和（剂量 →∞），靶点必需性最多把 HepG2 伪存活率压到 " +
        floor.hepg2.toFixed(3) + "，到不了 50%。含义：这些靶点的 on-target 必需性" +
        "不足以解释 50% 的细胞杀伤。若实验里确实看到 ≥50% 杀伤，那部分毒性来自" +
        "必需性以外的机制 —— 脱靶、活性代谢物、非依赖性应激 —— 必需性轴对它没有预测力。";
    } else {
      interp = "等效剂量比无法定义";
      why = "在 10⁻³–10³ × 规格剂量内两条曲线没有共同可比的抑制水平。";
    }
    out.selectivity = { iso, isoLv, mT, mN, byLevel, cons, floor, satP, ratio, interp, why };
    claim({ ev: "MODELED",
      ttl: "等效剂量比" + (isoLv ? "（正常系 ÷ 肝癌系，同为 " + Math.round(isoLv * 100) + "% 抑制）" : ""),
      val: iso === null ? "无法定义" : iso.toFixed(3) + " ×",
      basis: "HepG2 vs RPE1 Chronos 差值（✅ 实测）经组合律 L1–L4 传播",
      caveat: "对称指标：两系敏感性相同时恰为 1.0。它是比值，是曲线类量里最可解释的；" +
              "但它不是治疗指数，不能替代体内 PK/PD 与毒理。" });
    claim({ ev: "ANCHORED", ttl: "饱和杀伤压力比（HepG2 ÷ RPE1）",
      val: ratio === null ? "RPE1 无压力，比值发散" : ratio.toFixed(3) + " ×",
      basis: "两个 Chronos gene effect 的直接比值，不经过任何剂量假设",
      caveat: "本卡上唯一完全不依赖建模假设的选择性量 —— 它与浓度、IC50、Hill、暴露时长都无关。" });
  }

  /* 通路层 */
  const terms = [], srcs = [];
  traces.forEach((t) => {
    const p = getPhe(t.name);
    if (p) {
      terms.push([t.theta, p.pw_names, p.pw_vals]);
      srcs.push(t.name + " → VC-PHE（LOO r=" + p.loo_r.toFixed(3) + "）");
      return;
    }
    if (t.target) {
      const q = getPrt(t.target, spec.context);
      if (q) {
        terms.push([t.theta, q.pw_names, q.pw_vals]);
        srcs.push(t.name + " → VC-PRT " + t.target + "@" + spec.context +
          "（live r=" + q.live_r.toFixed(3) + " / 上限 " + q.ceiling.toFixed(3) + "）");
      }
    }
  });
  if (!terms.length) {
    out.pathway = { available: false, srcs: [], names: [], vals: [] };
    claim({ ev: "UNSUPPORTED", ttl: "通路响应谱", val: "不可用",
      basis: "无任何化合物落入 PhenoMap（8 化合物）或 PerturbLens（6 基因 × 4 系）覆盖",
      caveat: "引擎在这里拒绝外推，而不是给一个看起来合理的数。" });
  } else {
    const acc = Object.create(null);
    terms.forEach(([w, ns, vs]) => ns.forEach((n, i) => {
      acc[n] = (acc[n] || 0) + w * vs[i];
    }));
    const items = Object.keys(acc).map((k) => [k, acc[k]]).sort((a, b) => a[1] - b[1]);
    out.pathway = { available: true, srcs, names: items.map((x) => x[0]),
                    vals: items.map((x) => x[1]) };
    claim({ ev: "PREDICTED", ttl: "通路响应谱（占据率加权叠加）",
      val: items.length + " 条 Hallmark", basis: srcs.join("；"),
      caveat: "单位响应向量来自完全敲低 / 筛选浓度，按占据率线性缩放是 ⚠️ 假设；" +
              "多化合物按通路名相加没有考虑通路间串扰。" });
  }

  /* 组合项 */
  const nq = traces.filter((t) => t.dep).length;
  if (nq < 2 || !hasEss) {
    out.combination = { applicable: false,
      reason: nq < 2 ? "少于 2 个可定量化合物，无组合项" : "无必需性骨干" };
  } else {
    const press = {};
    ESS_LINES.forEach((ln) => { press[ln] = traces.map((t) => t.press ? t.press[ln] : 0); });
    const add = {}, cor = {};
    ESS_LINES.forEach((ln) => {
      add[ln] = viab(logViab(press[ln], spec.exposure, 1));
      cor[ln] = viab(logViab(press[ln], spec.exposure, gamma));
    });
    if (epiHit) {
      out.combination = { applicable: true, evidence: "measured", hit: epiHit,
                          gamma, add, cor };
      claim({ ev: "ANCHORED", ttl: "实测上位性交叉项（" + epiHit.pair + "）",
        val: "上位性 " + epiHit.epi.toFixed(3) + " → γ = " + gamma.toFixed(3),
        basis: "Norman 2019 双扰动 Perturb-seq，排名 " + epiHit.rank + "/" +
               D.cmb.pairs.length + "，" + epiHit.ncells + " 细胞",
        caveat: "原始测量在 K562 CRISPRa 遗传扰动上；迁移到你的细胞系与化学抑制是 ⚠️ 外推，" +
                "且单组合仅数百细胞。" });
    } else {
      const gLo = epiGamma(EPI_IQR[0]), gHi = epiGamma(EPI_IQR[1]);
      const band = {};
      ESS_LINES.forEach((ln) => {
        band[ln] = [viab(logViab(press[ln], spec.exposure, gLo)),
                    viab(logViab(press[ln], spec.exposure, gHi))].sort((a, b) => a - b);
      });
      out.combination = { applicable: true, evidence: "no_measurement_fallback_additive", gamma: 1,
                          add, cor: add, band };
      claim({ ev: "MODELED", ttl: "上位性交叉项", val: "无实测证据，回落到 Bliss 加和",
        basis: "该靶点对不在 ComboMap 的 " + D.cmb.pairs.length + " 对里",
        caveat: "Norman 2019 中 " + Math.round(FRAC_NOISE * 100) + "% 的组合新信号 >1× 噪声" +
                "（中位 " + NON_MED.toFixed(2) + "×），即加和预期系统性偏保守；" +
                "误差带按该经验分布的四分位给出，而非声称加和成立。" });
    }
  }

  /* 置信度 */
  const ctx = getTwn(spec.context) || {}, nrm = getTwn("rpe1") || {};
  const frac = ctx.frac || 0;
  const cov = traces.filter((t) => t.dep).length / Math.max(1, traces.length);
  const pwCov = traces.filter((t) => t.pwSrc !== "none").length / Math.max(1, traces.length);
  const ph = traces.some((t) => t.placeholder);
  let score = frac * (0.5 + 0.3 * cov + 0.2 * pwCov);
  if (ph) score *= 0.6;
  out.confidence = { ctx, nrm, frac, cov, pwCov, ph, score };
  claim({ ev: "MODELED", ttl: "综合置信度", val: score.toFixed(3),
    basis: "TwinCell " + spec.context + " 达 oracle 上限 " +
           (frac * 100).toFixed(1) + "%（✅ 实测）× 覆盖度",
    caveat: "用于排序的启发式分数，不是校准过的概率。" });

  /* 全局提示 */
  if (hasEss && traces.filter((t) => t.dep).every((t) =>
      ESS_LINES.every((ln) => t.dep[ln] >= 0))) {
    led.notes.push("所有已解析靶点在两系上的 Chronos 均为非负 —— 敲除它们不损伤（甚至促进）生长。" +
      "必需性轴在此规格下不产生任何杀伤信号，伪存活率恒为 1 并不代表无毒，" +
      "只代表这条轴对该组合无预测力。");
  }
  if (ESS_LINES.indexOf(spec.context) < 0) {
    led.notes.push("毒性主轴固定在 HepG2 ↔ RPE1 上 —— Chronos 依赖性数据只覆盖这两个系。" +
      "你选的上下文 " + spec.context + " 只参与通路层与置信度加权，不改变主轴数值。");
  }
  led.notes.push("力场类比的边界：与分子力场只有相对能量可解释一样，本引擎只有系间比值、" +
    "剂量位移、相对加和预期的偏离是可解释的。");
  led.notes.push("必需性与毒性通路激活仅弱耦合（Pearson 0.137）——「肝选择性必需」" +
    "不等于「特异性激活肝毒性通路」，两者须独立评估。");
  if (spec.context === "k562" && hasPw) {
    led.notes.push("K562 的跨系共识被原报告标为 EXPLORATORY，通路层外推需额外验证。");
  }
  return out;
}

/* ---------- 预设场景 ---------- */
const PRESETS = [
  { id: "A", label: "A · 选择性窗口", ctx: "hepg2", exp: 72,
    cmpds: [{ name: "Pemetrexed", conc: 0.5, ic50: 0.1, hill: 1 }] },
  { id: "A2", label: "A2 · 同条件广谱毒", ctx: "hepg2", exp: 72,
    cmpds: [{ name: "Bortezomib", conc: 0.5, ic50: 0.1, hill: 1 }] },
  { id: "B", label: "B · 同药加压", ctx: "hepg2", exp: 144,
    cmpds: [{ name: "Pemetrexed", conc: 10, ic50: 0.1, hill: 1 }] },
  { id: "C", label: "C · 组合·有实测", ctx: "hepg2", exp: 72,
    cmpds: [{ name: "cmpd-PLK4i", target: "PLK4", conc: 1, ic50: 0.2, hill: 1 },
            { name: "cmpd-STILi", target: "STIL", conc: 1, ic50: 0.3, hill: 1 }] },
  { id: "C2", label: "C2 · 组合·无实测", ctx: "hepg2", exp: 72,
    cmpds: [{ name: "Doxorubicin", conc: 0.5, ic50: 0.05, hill: 1 },
            { name: "Pemetrexed", conc: 0.5, ic50: 0.1, hill: 1 }] },
  { id: "D", label: "D · 覆盖不全", ctx: "hepg2", exp: 72,
    cmpds: [{ name: "Paclitaxel", conc: 0.1, ic50: 0.005, hill: 1 },
            { name: "Sorafenib", conc: 5, ic50: null, hill: 1 },
            { name: "内部代号-X271", conc: 1, ic50: 1, hill: 1 }] },
  { id: "E", label: "E · 换上下文", ctx: "jurkat", exp: 72,
    cmpds: [{ name: "cmpd-PLK4i", target: "PLK4", conc: 1, ic50: 0.2, hill: 1 },
            { name: "cmpd-STILi", target: "STIL", conc: 1, ic50: 0.3, hill: 1 }] }
];

/* ---------- 状态 ---------- */
let state = JSON.parse(JSON.stringify(PRESETS[0]));

/* ---------- 渲染工具 ---------- */
const esc = (s) => String(s).replace(/[&<>"]/g, (c) =>
  ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
const fmt = (v, n) => (v === null || v === undefined || !isFinite(v))
  ? "—" : Number(v).toFixed(n === undefined ? 3 : n);
const BADGE = {
  ANCHORED: ['anch', '✅ 实测'], PREDICTED: ['pred', '◐ 预测'],
  MODELED: ['mod', '⚠️ 假设'], UNSUPPORTED: ['uns', '⛔ 拒答']
};
const cssv = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

/* ---------- 规格面板 ---------- */
function buildStatic() {
  const ctx = $("ctx");
  Object.keys(CONTEXTS).forEach((k) => {
    const o = document.createElement("option");
    o.value = k; o.textContent = CONTEXTS[k]; ctx.appendChild(o);
  });
  const dl = $("names"), seen = Object.create(null);
  D.phe.drugs.concat(Object.keys(D.tox.drug2target))
    .concat(Object.keys(SUPP)).concat(D.prt.genes)
    .concat(["TYMS", "PSMB5", "PLK4", "STIL", "TOP2A", "TUBB", "HMGCR"])
    .forEach((n) => {
      if (seen[n]) return; seen[n] = 1;
      const o = document.createElement("option"); o.value = n; dl.appendChild(o);
    });
  $("presets").innerHTML = PRESETS.map((p) =>
    '<button class="preset" type="button" data-p="' + p.id + '">' +
    esc(p.label) + "</button>").join("");
  $("presets").addEventListener("click", (e) => {
    const b = e.target.closest("[data-p]"); if (!b) return;
    const p = PRESETS.find((x) => x.id === b.dataset.p);
    state = JSON.parse(JSON.stringify(p));
    syncSpec(); run();
  });
  $("addbtn").addEventListener("click", () => {
    state.cmpds.push({ name: "", conc: 1, ic50: null, hill: 1 });
    state.id = null; renderCmpds(); run();
  });
  $("ctx").addEventListener("change", () => { state.ctx = $("ctx").value; state.id = null; run(); });
  $("exp").addEventListener("input", () => {
    state.exp = Math.max(1, Number($("exp").value) || 72); state.id = null; run();
  });
  $("themebtn").addEventListener("click", () => {
    const cur = document.documentElement.getAttribute("data-theme");
    const dark = cur ? cur === "dark"
      : window.matchMedia("(prefers-color-scheme: dark)").matches;
    document.documentElement.setAttribute("data-theme", dark ? "light" : "dark");
    run();
  });
}

function syncSpec() {
  $("ctx").value = state.ctx; $("exp").value = state.exp;
  renderCmpds();
  document.querySelectorAll("[data-p]").forEach((b) =>
    b.setAttribute("aria-pressed", String(b.dataset.p === state.id)));
}

function renderCmpds() {
  const host = $("cmpds");
  host.innerHTML = state.cmpds.map((c, i) => {
    const hit = resolveEss(c.name, c.target);
    const rc = c.name ? (hit ? "hit" : "miss") : "";
    const rt = !c.name ? "填写化合物名或基因符号"
      : hit ? "→ " + hit.gene + "　HepG2 " + fmt(hit.dep.hepg2) +
              " / RPE1 " + fmt(hit.dep.rpe1)
      : "未命中 17,787 基因表与化合物映射";
    return '<div class="cmpd" data-i="' + i + '">' +
      '<div class="r1">' +
        '<input class="ctl" list="names" data-f="name" placeholder="化合物名或基因符号" value="' +
          esc(c.name || "") + '" aria-label="化合物名">' +
        '<button class="xbtn" type="button" data-del="' + i + '" aria-label="删除">×</button>' +
      "</div>" +
      '<div class="r2">' +
        '<label><span>浓度 µM</span><input class="ctl" data-f="conc" type="number" step="any" min="0" value="' + c.conc + '"></label>' +
        '<label><span>IC50 µM</span><input class="ctl" data-f="ic50" type="number" step="any" min="0" placeholder="占位 1.0" value="' + (c.ic50 === null || c.ic50 === undefined ? "" : c.ic50) + '"></label>' +
        '<label><span>Hill n</span><input class="ctl" data-f="hill" type="number" step="0.1" min="0.1" value="' + (c.hill || 1) + '"></label>' +
      "</div>" +
      '<div class="resolve ' + rc + '">' + esc(rt) + "</div>" +
    "</div>";
  }).join("");
}

$("cmpds").addEventListener("input", (e) => {
  const f = e.target.dataset.f; if (!f) return;
  const i = Number(e.target.closest("[data-i]").dataset.i);
  if (f === "name") { state.cmpds[i].name = e.target.value; state.cmpds[i].target = null; }
  else if (f === "ic50") state.cmpds[i].ic50 = e.target.value === "" ? null : Number(e.target.value);
  else state.cmpds[i][f] = Number(e.target.value);
  state.id = null;
  if (f === "name") {
    const box = e.target.closest("[data-i]").querySelector(".resolve");
    const hit = resolveEss(state.cmpds[i].name, null);
    box.className = "resolve " + (state.cmpds[i].name ? (hit ? "hit" : "miss") : "");
    box.textContent = !state.cmpds[i].name ? "填写化合物名或基因符号"
      : hit ? "→ " + hit.gene + "　HepG2 " + fmt(hit.dep.hepg2) + " / RPE1 " + fmt(hit.dep.rpe1)
      : "未命中 17,787 基因表与化合物映射";
  }
  run();
});
$("cmpds").addEventListener("click", (e) => {
  const b = e.target.closest("[data-del]"); if (!b) return;
  if (state.cmpds.length <= 1) return;
  state.cmpds.splice(Number(b.dataset.del), 1);
  state.id = null; renderCmpds(); run();
});

/* ---------- 装配链 ---------- */
const KERNELS = [
  { code: "VC-TOX", nm: "ToxSentinel", role: "必需性地形　17,787 靶点 × 2 系" },
  { code: "VC-PHE", nm: "PhenoMap", role: "化学表型向量　8 化合物" },
  { code: "VC-PRT", nm: "PerturbLens", role: "转录响应向量　6 基因 × 4 系" },
  { code: "VC-CMB", nm: "ComboMap", role: "上位性交叉项　126 组实测" },
  { code: "VC-TWN", nm: "TwinCell", role: "上下文置信权重　4 系" }
];
function renderChain(r) {
  const got = { "VC-TOX": null, "VC-PHE": null, "VC-PRT": null, "VC-CMB": null, "VC-TWN": null };
  const ess = r.traces.filter((t) => t.dep);
  if (ess.length) got["VC-TOX"] = ess.length + " 个靶点的 Chronos 值";
  const phe = r.traces.filter((t) => t.pwSrc === "VC-PHE");
  if (phe.length) got["VC-PHE"] = phe.map((t) => t.name).join("、") + " 的通路向量";
  const prt = r.traces.filter((t) => t.pwSrc === "VC-PRT");
  if (prt.length) got["VC-PRT"] = prt.map((t) => t.target).join("、") + "@" + r.spec.context;
  if (r.epiHit) got["VC-CMB"] = r.epiHit.pair + "　γ = " + r.gamma.toFixed(3);
  if (r.confidence && r.confidence.ctx.verdict) {
    got["VC-TWN"] = r.spec.context + " 达上限 " + (r.confidence.frac * 100).toFixed(1) + "%";
  }
  $("chain").innerHTML = KERNELS.map((k) => {
    const on = got[k.code];
    return '<div class="node ' + (on ? "on" : "off") + '">' +
      '<div class="code">' + k.code + "</div>" +
      '<div class="nm">' + k.nm + "</div>" +
      '<div class="role">' + esc(k.role) + "</div>" +
      '<div class="got">' + (on ? esc(on) : "未参与本次装配") + "</div>" +
    "</div>";
  }).join("");
  const laws = [
    ["L1", "θ = Cⁿ/(Cⁿ+IC50ⁿ)", "剂量→占据率"],
    ["L2", "响应 = θ × 单位向量", "占据率→表型"],
    ["L3", "log V = −λ Σ θᵢ·压力ᵢ", "Bliss 叠加"],
    ["L4", "log V ← γ(上位性) × log V", "交叉项"]
  ];
  $("laws").innerHTML = laws.map((l) =>
    '<span class="law"><b>' + l[0] + "</b>　" + esc(l[1]) + "　" + esc(l[2]) +
    "　⚠️ 未拟合</span>").join("") +
    '<span class="law">λ = ln2　约定锚定　t<sub>ref</sub> = 72 h</span>';
}

/* ---------- 图表 ---------- */
function doseChart(r) {
  if (!r.curve) return "";
  const W = 660, H = 260, ML = 46, MR = 62, MT = 14, MB = 40;
  const iw = W - ML - MR, ih = H - MT - MB;
  const lx = (m) => ML + (Math.log10(m) + 3) / 6 * iw;
  const ly = (v) => MT + (1 - v) * ih;
  const path = (ln) => r.curve.ms.map((m, i) =>
    (i ? "L" : "M") + lx(m).toFixed(1) + " " + ly(r.curve.curves[ln][i]).toFixed(1)).join(" ");
  const ticks = [-3, -2, -1, 0, 1, 2, 3];
  const yt = [0, 0.25, 0.5, 0.75, 1];
  const tCol = cssv("--target"), nCol = cssv("--normal");
  const ink3 = cssv("--ink-3"), line = cssv("--line"), ink2 = cssv("--ink-2");
  const endT = r.curve.curves.hepg2[r.curve.curves.hepg2.length - 1];
  const endN = r.curve.curves.rpe1[r.curve.curves.rpe1.length - 1];
  const specT = r.viability.hepg2, specN = r.viability.rpe1;

  let g = "";
  yt.forEach((v) => {
    g += '<line x1="' + ML + '" x2="' + (ML + iw) + '" y1="' + ly(v) + '" y2="' + ly(v) +
      '" stroke="' + line + '" stroke-width="1"/>' +
      '<text x="' + (ML - 8) + '" y="' + (ly(v) + 4) + '" text-anchor="end" font-size="10.5" ' +
      'fill="' + ink3 + '" font-family="' + cssv("--mono") + '">' + (v * 100) + "%</text>";
  });
  ticks.forEach((t) => {
    const x = lx(Math.pow(10, t));
    g += '<line x1="' + x + '" x2="' + x + '" y1="' + (MT + ih) + '" y2="' + (MT + ih + 4) +
      '" stroke="' + ink3 + '" stroke-width="1"/>' +
      '<text x="' + x + '" y="' + (MT + ih + 17) + '" text-anchor="middle" font-size="10.5" ' +
      'fill="' + ink3 + '" font-family="' + cssv("--mono") + '">' +
      (t === 0 ? "1×" : "10" + String(t).replace("-", "⁻")
        .replace(/[0-9]/g, (d) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[+d])) + "</text>";
  });
  const x1 = lx(1);
  g += '<line x1="' + x1 + '" x2="' + x1 + '" y1="' + MT + '" y2="' + (MT + ih) +
    '" stroke="' + ink2 + '" stroke-width="1" stroke-dasharray="3 3"/>' +
    '<text x="' + (x1 + 5) + '" y="' + (MT + 11) + '" font-size="10" fill="' + ink2 +
    '" font-family="' + cssv("--mono") + '">你的规格剂量</text>';
  g += '<path d="' + path("rpe1") + '" fill="none" stroke="' + nCol + '" stroke-width="2"/>';
  g += '<path d="' + path("hepg2") + '" fill="none" stroke="' + tCol + '" stroke-width="2"/>';
  g += '<circle cx="' + lx(1) + '" cy="' + ly(specN) + '" r="4.5" fill="' + nCol +
    '" stroke="' + cssv("--surface") + '" stroke-width="2"/>';
  g += '<circle cx="' + lx(1) + '" cy="' + ly(specT) + '" r="4.5" fill="' + tCol +
    '" stroke="' + cssv("--surface") + '" stroke-width="2"/>';
  g += '<text x="' + (ML + iw + 7) + '" y="' + (ly(endT) + 4) + '" font-size="10.5" fill="' +
    tCol + '" font-family="' + cssv("--mono") + '">HepG2</text>';
  g += '<text x="' + (ML + iw + 7) + '" y="' + (ly(endN) + 4) + '" font-size="10.5" fill="' +
    nCol + '" font-family="' + cssv("--mono") + '">RPE1</text>';
  g += '<text x="' + (ML + iw / 2) + '" y="' + (H - 3) + '" text-anchor="middle" font-size="10.5" ' +
    'fill="' + ink3 + '" font-family="' + cssv("--mono") + '">剂量倍数（× 你规格里的浓度组合，对数轴）</text>';

  return '<figure><div class="legend">' +
    '<span><i class="swatch" style="background:' + tCol + '"></i>HepG2 · 肝癌（目标）</span>' +
    '<span><i class="swatch" style="background:' + nCol + '"></i>RPE1 · 正常上皮（参照）</span></div>' +
    '<div class="chartbox scroll" id="dcbox"><svg viewBox="0 0 ' + W + " " + H +
    '" width="100%" style="min-width:560px;display:block" role="img" ' +
    'aria-label="剂量-响应曲线">' + g + "</svg>" +
    '<div class="tip" id="dctip"></div></div>' +
    "<figcaption>纵轴是<b>伪</b>存活率：λ 取 ln2，锚定「单个中位必需靶点在 100% 占据、72 h 恰好给 50%」。" +
    "这个锚定没有拟合任何存活率数据，所以绝对高度不可读，可读的是两条线的<b>横向距离</b>。</figcaption></figure>";
}

function pathwayChart(r) {
  if (!r.pathway || !r.pathway.available) return "";
  const names = r.pathway.names, vals = r.pathway.vals;
  const n = names.length, rowH = 19, ML = 176, MR = 54, W = 660;
  const H = n * rowH + 30, iw = W - ML - MR;
  const mx = Math.max.apply(null, vals.map(Math.abs)) || 1;
  const zx = ML + iw / 2;
  const sc = (v) => (v / mx) * (iw / 2);
  const slate = cssv("--slate"), ink3 = cssv("--ink-3"), line = cssv("--line"), ink2 = cssv("--ink-2");
  let g = '<line x1="' + zx + '" x2="' + zx + '" y1="14" y2="' + (H - 16) +
    '" stroke="' + ink2 + '" stroke-width="1"/>';
  names.forEach((nm, i) => {
    const y = 18 + i * rowH, w = sc(vals[i]);
    const x = w < 0 ? zx + w : zx;
    g += '<rect x="' + x.toFixed(1) + '" y="' + y + '" width="' + Math.max(1, Math.abs(w)).toFixed(1) +
      '" height="11" rx="2" fill="' + slate + '" opacity="' +
      (0.42 + 0.58 * Math.abs(vals[i]) / mx).toFixed(2) + '"/>';
    g += '<text x="' + (ML - 10) + '" y="' + (y + 9) + '" text-anchor="end" font-size="10.5" fill="' +
      ink3 + '">' + esc(nm) + "</text>";
  });
  g += '<text x="' + (zx - 8) + '" y="10" text-anchor="end" font-size="10" fill="' + ink3 +
    '" font-family="' + cssv("--mono") + '">← 受抑制</text>';
  g += '<text x="' + (zx + 8) + '" y="10" font-size="10" fill="' + ink3 +
    '" font-family="' + cssv("--mono") + '">被动员 →</text>';
  g += '<line x1="' + ML + '" x2="' + (ML + iw) + '" y1="' + (H - 16) + '" y2="' + (H - 16) +
    '" stroke="' + line + '"/>';
  g += '<text x="' + (ML + iw / 2) + '" y="' + (H - 3) + '" text-anchor="middle" font-size="10.5" fill="' +
    ink3 + '" font-family="' + cssv("--mono") + '">通路预测 Δ（占据率加权）　最大 |Δ| = ' +
    mx.toFixed(3) + "</text>";
  return '<figure><svg viewBox="0 0 ' + W + " " + H + '" width="100%" style="min-width:560px;display:block" ' +
    'role="img" aria-label="通路响应谱">' + g + "</svg>" +
    "<figcaption>零线两侧的位置就是方向，所以这里只用一种色调 —— 页面上的橙与蓝始终只表示细胞系。<br>" +
    "来源：" + r.pathway.srcs.map(esc).join("；") + "</figcaption></figure>";
}

function comboChart(r) {
  const c = r.combination;
  if (!c || !c.applicable) return "";
  const W = 660, H = 118, ML = 116, MR = 96;
  const iw = W - ML - MR;
  const all = [];
  ESS_LINES.forEach((ln) => { all.push(c.add[ln], c.cor[ln]);
    if (c.band) all.push(c.band[ln][0], c.band[ln][1]); });
  const lo = Math.max(0, Math.min.apply(null, all) - 0.08);
  const hi = Math.min(1, Math.max.apply(null, all) + 0.08);
  const sx = (v) => ML + (v - lo) / Math.max(1e-9, hi - lo) * iw;
  const ink3 = cssv("--ink-3"), line = cssv("--line"), surf = cssv("--surface");
  let g = "";
  ESS_LINES.forEach((ln, i) => {
    const y = 30 + i * 38;
    const col = ln === "hepg2" ? cssv("--target") : cssv("--normal");
    g += '<text x="' + (ML - 12) + '" y="' + (y + 4) + '" text-anchor="end" font-size="11" fill="' +
      ink3 + '">' + esc(LINE_LABEL[ln]) + "</text>";
    if (c.band) {
      g += '<rect x="' + sx(c.band[ln][0]) + '" y="' + (y - 7) + '" width="' +
        Math.max(2, sx(c.band[ln][1]) - sx(c.band[ln][0])) + '" height="14" rx="2" fill="' +
        col + '" opacity="0.16"/>';
    }
    g += '<line x1="' + sx(c.add[ln]) + '" x2="' + sx(c.cor[ln]) + '" y1="' + y + '" y2="' + y +
      '" stroke="' + col + '" stroke-width="2"/>';
    g += '<circle cx="' + sx(c.add[ln]) + '" cy="' + y + '" r="5" fill="' + surf +
      '" stroke="' + col + '" stroke-width="2"/>';
    g += '<circle cx="' + sx(c.cor[ln]) + '" cy="' + y + '" r="5" fill="' + col +
      '" stroke="' + surf + '" stroke-width="2"/>';
    g += '<text x="' + (ML + iw + 10) + '" y="' + (y + 4) + '" font-size="10.5" fill="' + ink3 +
      '" font-family="' + cssv("--mono") + '">' + c.add[ln].toFixed(3) + " → " +
      c.cor[ln].toFixed(3) + "</text>";
  });
  g += '<line x1="' + ML + '" x2="' + (ML + iw) + '" y1="' + (H - 24) + '" y2="' + (H - 24) +
    '" stroke="' + line + '"/>';
  [lo, (lo + hi) / 2, hi].forEach((v) => {
    g += '<text x="' + sx(v) + '" y="' + (H - 10) + '" text-anchor="middle" font-size="10" fill="' +
      ink3 + '" font-family="' + cssv("--mono") + '">' + v.toFixed(2) + "</text>";
  });
  g += '<text x="' + ML + '" y="12" font-size="10" fill="' + ink3 +
    '" font-family="' + cssv("--mono") + '">空心 = Bliss 加和预期　实心 = ' +
    (c.evidence === "measured" ? "上位性修正" : "同上（无实测证据）") +
    (c.band ? "　浅块 = 经验四分位区间" : "") + "</text>";
  return '<figure><svg viewBox="0 0 ' + W + " " + H + '" width="100%" style="min-width:560px;display:block" ' +
    'role="img" aria-label="组合项：加和预期 vs 上位性修正">' + g + "</svg></figure>";
}

/* ---------- 结果渲染 ---------- */
function renderResults(r) {
  const s = r.selectivity, out = [];
  out.push('<div class="tier"><b>装配层级　' + esc(r.tierLabel) + "</b></div>");

  /* 头条三格 */
  const cards = [];
  if (s && s.ratio !== null) {
    cards.push('<div><div class="k"><span class="badge anch">✅ 实测</span>饱和杀伤压力比</div>' +
      '<div class="v">' + fmt(s.ratio, 2) + " ×</div>" +
      '<div class="n">HepG2 ÷ RPE1，两个 Chronos 值的直接比值。<b>与浓度、IC50、Hill、时长全部无关</b> —— ' +
      "本页唯一不依赖任何建模假设的选择性量。</div></div>");
  } else {
    cards.push('<div><div class="k"><span class="badge uns">⛔ 拒答</span>饱和杀伤压力比</div>' +
      '<div class="v sm">不可用</div><div class="n">没有靶点解析到依赖性数据，或 RPE1 侧压力为零。</div></div>');
  }
  if (s) {
    cards.push('<div><div class="k"><span class="badge mod">⚠️ 假设</span>等效剂量比' +
      (s.isoLv ? "（同为 " + Math.round(s.isoLv * 100) + "% 抑制）" : "") + "</div>" +
      '<div class="v' + (s.iso === null ? " sm" : "") + '">' +
      (s.iso === null ? "无法定义" : fmt(s.iso, 2) + " ×") + "</div>" +
      '<div class="n">' + esc(s.interp) + (s.why ? "<br>" + esc(s.why) : "") + "</div></div>");
    cards.push('<div><div class="k"><span class="badge mod">⚠️ 假设</span>饱和下限（剂量 →∞）</div>' +
      '<div class="v sm">HepG2 ' + fmt(s.floor.hepg2) + "　/　RPE1 " + fmt(s.floor.rpe1) + "</div>" +
      '<div class="n">on-target 必需性最多能把伪存活率压到这里。压不下去的那部分杀伤，' +
      "必需性轴对它没有预测力。</div></div>");
  } else {
    cards.push('<div><div class="k"><span class="badge uns">⛔ 拒答</span>等效剂量比</div>' +
      '<div class="v sm">不可用</div><div class="n">无必需性骨干，引擎不出剂量-响应。</div></div>');
    cards.push('<div><div class="k"><span class="badge uns">⛔ 拒答</span>饱和下限</div>' +
      '<div class="v sm">不可用</div><div class="n">「算不出杀伤」不等于「无毒」。</div></div>');
  }
  out.push('<div class="headline">' + cards.join("") + "</div>");

  /* 化合物装配轨迹 */
  out.push('<section class="panel"><header><h2>化合物装配轨迹</h2>' +
    '<span class="hint">每一步用了哪个 kernel</span></header><div class="scroll">' +
    '<table class="dt"><thead><tr><th>化合物</th><th>解析</th><th class="num">占据率 θ</th>' +
    '<th class="num">HepG2 dep</th><th class="num">RPE1 dep</th>' +
    '<th class="num">杀伤压力 H / R</th><th>分类</th></tr></thead><tbody>' +
    r.traces.map((t) =>
      "<tr><td>" + esc(t.name || "（空）") + "</td><td style=\"font-size:11px;color:var(--ink-3)\">" +
      esc(t.via) + "</td>" +
      '<td class="num">' + fmt(t.theta) + "</td>" +
      '<td class="num">' + (t.dep ? fmt(t.dep.hepg2) : "—") + "</td>" +
      '<td class="num">' + (t.dep ? fmt(t.dep.rpe1) : "—") + "</td>" +
      '<td class="num">' + (t.press ? fmt(t.press.hepg2, 2) + " / " + fmt(t.press.rpe1, 2) : "—") + "</td>" +
      '<td style="font-size:11px">' + esc(t.klass || "—") + "</td></tr>"
    ).join("") + "</tbody></table></div>" +
    (r.traces.some((t) => t.warn.length)
      ? '<div class="pad" style="padding-top:12px;display:flex;flex-direction:column;gap:7px">' +
        r.traces.flatMap((t) => t.warn.map((w) =>
          '<div style="font-size:11.5px;color:var(--warn);line-height:1.5">⚠️ <b>' +
          esc(t.name || "（空）") + "</b>　" + esc(w) + "</div>")).join("") + "</div>"
      : "") + "</section>");

  /* 剂量响应 */
  if (r.curve) {
    out.push('<section class="panel"><header><h2>剂量-响应与选择性</h2>' +
      '<span class="hint">L1 → L3 → L4 的输出</span></header><div class="pad">' +
      doseChart(r) +
      '<div class="scroll" style="margin-top:14px"><table class="dt"><thead><tr>' +
      "<th>抑制水平</th><th class=\"num\">HepG2 所需剂量</th><th class=\"num\">RPE1 所需剂量</th>" +
      "<th class=\"num\">等效剂量比</th></tr></thead><tbody>" +
      [0.5, 0.3, 0.2, 0.1].map((lv) => {
        const key = Math.round(lv * 100) + "%";
        const a = crossing(r.curve.ms, r.curve.curves.hepg2, 1 - lv);
        const b = crossing(r.curve.ms, r.curve.curves.rpe1, 1 - lv);
        const rr = s.byLevel[key];
        return "<tr><td>" + key + "</td>" +
          '<td class="num">' + (a ? a.toPrecision(3) + " ×" : "不可达") + "</td>" +
          '<td class="num">' + (b ? b.toPrecision(3) + " ×" : "不可达") + "</td>" +
          '<td class="num"><b>' + (rr ? rr.toFixed(2) + " ×" : "—") + "</b></td></tr>";
      }).join("") + "</tbody></table></div>" +
      '<p style="font-size:11.5px;color:var(--ink-3);line-height:1.55;margin-top:10px">' +
      "比值随抑制水平变化不是数值毛病，是占据率趋饱和的真实后果：越接近饱和，" +
      "较不敏感的那一系需要越不成比例的剂量。要跨规格比较，请比同一行。</p>" +
      "</div></section>");
  }

  /* 组合项 */
  const c = r.combination;
  if (c && c.applicable) {
    const head = c.evidence === "measured"
      ? '<span class="badge anch">✅ 实测</span>　命中 ComboMap 实测对 <b>' + esc(c.hit.pair) +
        "</b>　排名 " + c.hit.rank + "/126　" + c.hit.ncells + " 细胞　" +
        (c.hit.non !== null ? "新信号 " + c.hit.non.toFixed(2) + "× 噪声" : "新信号缺测")
      : '<span class="badge mod">⚠️ 假设</span>　该靶点对不在 ComboMap 的 126 对里 → 回落 Bliss 加和。' +
        "Norman 2019 中 " + Math.round(FRAC_NOISE * 100) + "% 的实测组合新信号 >1× 噪声（中位 " +
        NON_MED.toFixed(2) + "×），所以加和预期系统性偏保守，浅色块是按经验四分位给的区间。";
    out.push('<section class="panel"><header><h2>组合项</h2>' +
      '<span class="hint">L4 交叉项</span></header><div class="pad">' +
      '<p style="font-size:12.5px;color:var(--ink-2);line-height:1.6;margin-bottom:14px">' +
      head + "</p>" + comboChart(r) + "</div></section>");
  }

  /* 通路 */
  if (r.pathway) {
    out.push('<section class="panel"><header><h2>通路响应谱</h2>' +
      '<span class="hint">L2 占据率加权</span></header><div class="pad">' +
      (r.pathway.available ? pathwayChart(r)
        : '<p style="font-size:12.5px;color:var(--ink-3);line-height:1.6">' +
          "⛔ 不可用 —— 没有化合物落入 PhenoMap（8 化合物）或 PerturbLens（6 基因 × 4 系）的覆盖。" +
          "引擎在这里拒绝外推，而不是给一个看起来合理的数。</p>") +
      "</div></section>");
  }

  /* 置信度 */
  const cf = r.confidence;
  out.push('<section class="panel"><header><h2>置信度与上下文</h2>' +
    '<span class="hint">VC-TWN 加权</span></header><div class="scroll">' +
    '<table class="dt"><tbody>' +
    ["<tr><th>细胞上下文</th><td>" + esc(CONTEXTS[r.spec.context]) + "</td></tr>",
     "<tr><th>数字孪生就绪度</th><td>" + esc(cf.ctx.verdict || "—") + "　达 oracle 上限 " +
       (cf.frac * 100).toFixed(1) + "%</td></tr>",
     "<tr><th>必需性覆盖</th><td>" + Math.round(cf.cov * 100) + "%　·　通路覆盖 " +
       Math.round(cf.pwCov * 100) + "%　·　IC50 占位 " + (cf.ph ? "是" : "否") + "</td></tr>",
     "<tr><th>综合分数</th><td><b>" + cf.score.toFixed(3) +
       "</b>　⚠️ 排序用启发式，不是校准过的概率</td></tr>",
     "<tr><th>细胞系备注</th><td style=\"font-size:11.5px;line-height:1.5\">" +
       esc(cf.ctx.note || "") + "</td></tr>"].join("") +
    "</tbody></table></div></section>");

  $("res").innerHTML = out.join("");
  if (r.curve) attachTip(r);
}

function attachTip(r) {
  const box = $("dcbox"), tip = $("dctip");
  if (!box || !tip) return;
  const svg = box.querySelector("svg");
  const W = 660, ML = 46, MR = 62, iw = W - ML - MR;
  svg.addEventListener("mousemove", (e) => {
    const rect = svg.getBoundingClientRect();
    const px = (e.clientX - rect.left) / rect.width * W;
    const t = (px - ML) / iw;
    if (t < 0 || t > 1) { tip.classList.remove("on"); return; }
    const m = Math.pow(10, -3 + 6 * t);
    let bi = 0, bd = Infinity;
    r.curve.ms.forEach((mm, i) => {
      const d = Math.abs(Math.log(mm) - Math.log(m));
      if (d < bd) { bd = d; bi = i; }
    });
    tip.innerHTML = "剂量 " + r.curve.ms[bi].toPrecision(3) + " ×<br>HepG2 " +
      (r.curve.curves.hepg2[bi] * 100).toFixed(1) + "%<br>RPE1 " +
      (r.curve.curves.rpe1[bi] * 100).toFixed(1) + "%";
    tip.style.left = Math.min(rect.width - 118, Math.max(0, e.clientX - rect.left + 12)) + "px";
    tip.style.top = Math.max(0, e.clientY - rect.top - 52) + "px";
    tip.classList.add("on");
  });
  svg.addEventListener("mouseleave", () => tip.classList.remove("on"));
}

/* ---------- 账本 ---------- */
function renderLedger(r) {
  $("claims").innerHTML = r.led.claims.map((c) => {
    const b = BADGE[c.ev];
    return '<div class="claim"><div class="lv"><span class="badge ' + b[0] + '">' + b[1] + "</span></div>" +
      "<div><div class=\"ttl\">" + esc(c.ttl) + "</div>" +
      '<div class="val">' + esc(c.val) + "</div>" +
      '<div class="meta"><b>依据</b>　' + esc(c.basis) +
      (c.caveat ? "<br><b>边界</b>　" + esc(c.caveat) : "") + "</div></div></div>";
  }).join("");
  $("notes").innerHTML = r.led.notes.map((n) =>
    '<p class="note">' + esc(n) + "</p>").join("");

  const card = {
    engine: "VC-FF v0.1.0",
    spec: { context: r.spec.ctx || r.spec.context, exposure_h: r.spec.exposure,
            compounds: r.spec.compounds },
    coverage_tier: r.tier, coverage_tier_label: r.tierLabel,
    compound_trace: r.traces.map((t) => ({
      name: t.name, target: t.target, resolved_via: t.via,
      occupancy: +t.theta.toFixed(4), ic50_is_placeholder: t.placeholder,
      dep: t.dep, pressure: t.press, pathway_source: t.pwSrc,
      kernels_used: t.kernels, warnings: t.warn })),
    pseudo_viability_at_spec_dose: r.viability,
    selectivity: r.selectivity && {
      iso_effect_ratio_x: r.selectivity.iso, iso_effect_level: r.selectivity.isoLv,
      iso_effect_ratio_by_level: r.selectivity.byLevel,
      conservative_window_x: r.selectivity.cons,
      saturation_floor: r.selectivity.floor,
      saturation_pressure_ratio: r.selectivity.ratio,
      interpretation: r.selectivity.interp, why_undefined: r.selectivity.why },
    combination: r.combination,
    pathway_response: r.pathway && { available: r.pathway.available,
      sources: r.pathway.srcs, names: r.pathway.names, values: r.pathway.vals },
    confidence: { context: r.spec.context, frac_of_oracle: r.confidence.frac,
      twin_readiness: r.confidence.ctx.verdict,
      essentiality_coverage: r.confidence.cov, pathway_coverage: r.confidence.pwCov,
      ic50_placeholder_used: r.confidence.ph, composite_score: r.confidence.score },
    honesty_ledger: { claims: r.led.claims, notes: r.led.notes,
      refused_outputs: ["NOAEL / LOAEL 数值",
        "绝对 IC50 或 EC50（效价是你的输入参数，不是本引擎的预测）",
        "绝对细胞存活率百分比（λ 为约定锚定，非拟合）",
        "临床安全窗、给药剂量、毒理放行结论"] }
  };
  $("jsonout").textContent = JSON.stringify(card, null, 2);
}

/* ---------- 主循环 ---------- */
function run() {
  const spec = { context: state.ctx, exposure: state.exp,
    compounds: state.cmpds.map((c) => ({ name: c.name, target: c.target || null,
      conc: c.conc, ic50: c.ic50, hill: c.hill })) };
  const r = evaluate(spec);
  r.spec.ctx = state.ctx;
  renderChain(r); renderResults(r); renderLedger(r);
  document.querySelectorAll("[data-p]").forEach((b) =>
    b.setAttribute("aria-pressed", String(b.dataset.p === state.id)));
}

buildStatic();
syncSpec();
run();
})();
