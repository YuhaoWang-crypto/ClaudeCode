# VC-FF · 虚拟细胞力场

把 ToxSentinel / PerturbLens / PhenoMap / ComboMap / TwinCell 这五个虚拟细胞模型，
从「查表式数据库」重构成一个**可装配的力场引擎**。

在线 demo：https://claude.ai/code/artifact/39537239-af69-4c9e-bbdb-da9442611cc3

## 为什么不是数据库

数据库对同一个基因永远返回同一行。你查 TYMS，它返回 HepG2 −0.8181、RPE1 −0.1708；
换浓度、加联用药、改暴露时长，返回值一个字都不会变。

力场的做法是把这五个模型当成**同一个泛函的五个项**：

| 谁提供什么 | 内容 |
|---|---|
| 模型提供**泛函形式** | 哪些通路会动、哪个细胞系脆弱、两个靶点怎么耦合 |
| 客户提供**标量参数** | 浓度、效价 IC50、Hill 系数、暴露时长、细胞系、组合设计 |
| 引擎按**组合律**装配 | L1–L4，见下 |

换一组参数，变的不只是数值，还有能装配上的**层数**，以及引擎该在**哪里拒绝回答**。

## 五个 kernel

每个 kernel 只声明三件事：提供力场的哪一项、输入接口是什么、覆盖不到时怎么降级。
引擎不硬编码任何一个 kernel 的内部，只按接口装配。少一个 kernel，输出就少一层，而不是崩掉。

| 代号 | 模型 | 提供的项 | 覆盖 |
|---|---|---|---|
| VC-TOX | ToxSentinel | 必需性地形（Chronos gene effect） | 17,787 靶点 × HepG2/RPE1 两系，最广的骨干 |
| VC-PHE | PhenoMap | 化学表型向量（Hallmark 通路签名） | 8 个预计算化合物 |
| VC-PRT | PerturbLens | 转录响应向量 | 6 基因 × 4 系 = 24 组合 |
| VC-CMB | ComboMap | 上位性交叉项 | Norman 2019 的 126 对双扰动 |
| VC-TWN | TwinCell | 上下文置信权重 | 4 个细胞系 |

化合物→靶点映射来自 Tahoe 的 260 条，另加 18 条人工补充（在结果里单独标注为「补充映射」，
多靶点化合物额外警告）。

## 四条组合律

全部是引擎自己引入的**建模假设，未在本数据上拟合或验证**：

```
L1  剂量→占据率    θ = C^h / (C^h + IC50^h)          Hill
L2  占据率→表型    响应 = θ × 完全敲低的单位响应向量    线性缩放
L3  多化合物叠加    log V = −λ Σ θᵢ · 压力ᵢ            Bliss 独立
L4  交叉项         log V ← γ(上位性) × log V          ComboMap 修正
```

λ 取 ln2，是**约定锚定**：使「单个中位必需靶点（Chronos = −1）在 100% 占据、72 h」
恰好给出 0.5 的伪存活率。它没有拟合任何存活率数据。

### 因此：只有相对量成立

和分子力场只有相对能量可解释完全同理。可解释的是系间比值、剂量位移、相对加和预期的偏离；
**不可解释的是绝对存活率百分比、绝对 IC50、NOAEL、安全窗数值**。

结果卡上最硬的一个数是**饱和杀伤压力比**（HepG2 ÷ RPE1）：它是两个 Chronos 值的直接比值，
与浓度、IC50、Hill、暴露时长全部无关，是唯一完全不依赖建模假设的选择性量。

## 装配分层与拒答

引擎按实际接上的层数报告 tier，并在覆盖不到的地方拒绝出数，而不是给一个看起来合理的值：

- 没有任何靶点解析到依赖性数据 → **不输出伪存活率**（返回「不可用」而不是 1.0，
  因为「算不出杀伤」与「无毒」是两回事）
- 没有化合物落入通路层覆盖 → 通路谱标 `available: false`，不外推
- 靶点对不在 ComboMap 的 126 对里 → 回落 Bliss 加和，并用 Norman 2019 的经验分布
  给区间（91% 的实测组合新信号 >1× 噪声，中位 2.72×，说明加和预期系统性偏保守）
- 所有靶点 Chronos 非负 → 明确说明「伪存活率恒为 1 不代表无毒，只代表这条轴无预测力」
- 曲线在饱和占据下仍到不了 50% 抑制 → 报告**饱和下限**，并指出多出来的杀伤来自
  必需性以外的机制（脱靶、活性代谢物、非依赖性应激）

每个输出量都带证据等级：`ANCHORED` 实测 / `PREDICTED` 模型预测 / `MODELED` 建模假设 /
`UNSUPPORTED` 拒答。

## 用法

```bash
python -m vcff.cli describe            # 引擎接口自述（kernel、组合律、覆盖度）
python -m vcff.cli scenarios           # 列出内置客户场景
python -m vcff.cli run A_selective_window
python -m vcff.cli run A_selective_window --json
echo '{"context":"hepg2","compounds":[{"name":"Bortezomib","conc_uM":0.5,"ic50_uM":0.1}]}' \
  | python -m vcff.cli stdin
python -m vcff.tests                   # 31 项不变量测试
```

作为库：

```python
from vcff import ForceField, AssaySpec, Compound

ff = ForceField()
card = ff.evaluate(AssaySpec(
    context="hepg2",
    exposure_h=72.0,
    compounds=[Compound(name="Pemetrexed", conc_uM=0.5, ic50_uM=0.1)],
))
print(card.to_json())
```

## 七个客户场景

同一批 kernel，不同 spec，不同输出、不同层数、不同可回答性。

| 场景 | 说明 | 结果 |
|---|---|---|
| A | Pemetrexed → TYMS，0.5 µM / 72 h | 压力比 4.79×，但饱和下限 0.567 —— on-target 必需性解释不了 50% 杀伤 |
| A2 | 同浓度同时长换 Bortezomib → PSMB5 | 压力比 0.99×，等效剂量比 0.97× —— 广谱毒性，无选择性 |
| B | 与 A 同药，浓度 ×20、时长 ×2 | 压力比不变（不依赖剂量），但曲线、下限、可达抑制水平全变 |
| C | PLK4 + STIL 双靶 | 命中 ComboMap 实测对，启用交叉项 γ = 1.207 |
| C2 | Doxorubicin + Pemetrexed | 未测过 → 回落加和 + 经验四分位区间 |
| D | 三药，缺 IC50 / 多靶点 / 完全未命中 | 逐条降级，置信度从 0.60 掉到 0.345 |
| E | 与 C 同组合换到 Jurkat | 主轴数值不变（Chronos 不覆盖该系），只有置信度变 |

## 不变量测试

测的不是「数值对不对」（没有金标准可对），而是**引擎的行为约束**：

- Hill 占据率的有界性、单调性、只依赖 C/IC50 之比
- 两系依赖性相同时等效剂量比恰为 1.0（选择性的正确零点）
- 饱和压力比与浓度 / IC50 / Hill / 时长全部无关
- IC50 放大 10 倍，达半效剂量倍数也放大 10 倍
- 无覆盖时返回「不可用」而不是 1.0
- 上位性查询与靶点顺序无关；未测过的对返回 None
- 与原报告一致：91% 的组合新信号超过噪声底
- Bliss 加和预期不高于任一单药
- 换上下文不改变毒性主轴数值

Python 与浏览器两套实现在全部 7 个场景上数值一致（差异 <0.03%，来自前端数据打包精度）。

## 文件

```
vcff/
  kernels.py      五个 kernel + 接口声明 + 补充靶点映射
  physics.py      四条组合律
  engine.py       装配器 + 结果卡
  spec.py         AssaySpec 输入适配
  honesty.py      证据等级与拒绝输出清单
  scenarios.py    七个客户场景
  tests.py        31 项不变量测试
  cli.py          命令行
  data/*.json     从五个模型 HTML 抽取的数据层
demo/
  vcff_demo.html  自包含交互式 demo
```

## 边界

毒性主轴固定在 HepG2 ↔ RPE1，因为 Chronos 依赖性数据只覆盖这两个系；换上下文只影响
通路层与置信度加权。组合律 L1–L4 全部是未拟合的建模假设。上位性的原始测量在 K562
CRISPRa 遗传扰动上，迁移到其他细胞系与化学抑制是外推，且单组合仅数百细胞。
必需性与毒性通路激活仅弱耦合（Pearson 0.137）——「肝选择性必需」不等于
「特异性激活肝毒性通路」。

研究级可行性演示，输出供实验设计与优先级排序参考，不构成临床或毒理放行结论。

## 数据来源

Replogle 2022 *Cell*（Perturb-seq）· Adduri 2025 *bioRxiv*（STATE）·
Zhang 2025 *bioRxiv*（Tahoe-100M）· Subramanian 2017 *Cell*（L1000）·
Norman 2019 *Science*（双扰动）· Corsello 2020 *Nat Cancer*（PRISM）·
Liberzon 2015 *Cell Syst*（Hallmark）· Tsherniak 2017 / Meyers 2017（DepMap Chronos）
