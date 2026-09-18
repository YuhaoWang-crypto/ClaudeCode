# `m23_virtual_tcell` 运行说明

本文只讲"怎么把 `grn_pipeline/m23_virtual_tcell.py` 跑起来、跑出来的东西怎么读"。
模型本身能替代哪些湿实验、不能替代哪些，见 [`VIRTUAL_TCELL_REPORT.md`](VIRTUAL_TCELL_REPORT.md)。

下文所有数字都是在本仓库当前代码上实跑得到的，不是照抄注释。

---

## 1. 需要哪些输入数据

**不需要任何输入数据文件。** 脚本不读本地文件、不连数据库、不走网络。
所有参数都硬编码在源码里，跑之前不用准备 IEDB、VDJdb、表达矩阵或 HLA 分型。

真正意义上的"输入"是三个可调的模型量，通过函数参数传：

| 量 | 含义 | 默认值 / 默认扫描范围 |
|---|---|---|
| `tau_s` | pMHC–TCR 复合物寿命（秒），即配体"质量" | L2 扫 0.2–60 s，L3 用 0.3–30 s 六个点 |
| `L` | 配体"数量"（任意单位） | L2 用 0.3 / 1 / 3 / 10 |
| `n_steps` | 校对链步数 N | L1 扫 N=0,1,2,4,6；L2/L3 固定 N=4 |

`VirtualTCell.__init__` 里另外九个参数（`alpha`、`beta`、`K_e`、`hill`、`gamma`、
`s_gain`、`s_decay`、`s_inhib`、`k_p`）是 **示意参数**：拓扑取自
Altan-Bonnet & Germain 2005，数值只是为了让定性区间出现，没有拟合过任何数据集。
改这些数会改 L2/L3 的阈值和排序间距，但不会改 L1 的结论（L1 只依赖 N 和 k_p 的比值）。

---

## 2. 依赖怎么装

```bash
pip install numpy scipy matplotlib
```

实测通过的版本：Python 3.11.15、numpy 2.4.6、scipy 1.17.1、matplotlib 3.11.2。

说明：

- 这个模块**不用** `networkx`（README 里那行 `pip install` 是给整条流水线装的）。
- `matplotlib` 只在 `main()` 内部 import。如果只想调 `VirtualTCell` / `kpr_output`
  这些函数、不出图，装 numpy + scipy 就够。
- 脚本自己设了 `matplotlib.use("Agg")`，无显示器的服务器、容器里可以直接跑。

---

## 3. 命令怎么跑

```bash
cd /path/to/ClaudeCode            # 必须在仓库根目录，-m 是按当前目录找 grn_pipeline 包的
python3 -m grn_pipeline.m23_virtual_tcell
```

在别的目录跑会报 `ModuleNotFoundError: No module named 'grn_pipeline'`；
要么 `cd` 到根目录，要么把根目录加进 `PYTHONPATH`。

- **耗时**：约 2.5 秒（单核，绝大部分花在 L2 的 ODE 扫描上）。
- **没有随机性**：没有任何随机种子，连跑两次标准输出逐行一致。
- 跑整条流水线时它是最后一个模块：`python3 -m grn_pipeline.run_all`。

当库用（`report` 就是 `main` 的别名，返回结果字典）：

```python
from grn_pipeline.m23_virtual_tcell import VirtualTCell, kpr_output, il2_readout

cell = VirtualTCell(n_steps=4)
cell.tau_threshold(L=1.0)        # 1.108  —— ERK 开启的配体寿命阈值（秒）
cell.steady_state(tau_s=30.0, L=1.0, e0=0.0)   # (0.981, ...) ERK 与 SHP-1 稳态
cell.bistable(tau_s=0.5, L=1.0)  # 同一配体下是否存在两个稳态
```

---

## 4. 输出是什么

**一、标准输出的三段文字报告**（L1 / L2 / L3），外加结尾的 HONESTY LINE。

**二、一张图**：`figures/m23_virtual_tcell.png`（三联图，150 dpi，约 180 KB）。
路径是从源码位置推出来的，跟你在哪个目录敲命令无关；同名文件会被直接覆盖。

**三、`main()` / `report()` 的返回字典**，六个键：

| 键 | 内容 |
|---|---|
| `exponents` | `[(N, 指数@k_off=1000, 指数@k_off=10, 指数@k_off=0.02), ...]` |
| `thresholds` | `{配体剂量 L: ERK 开启的 tau 阈值（秒）}`，范围内无解则为 `None` |
| `threshold_spread` | 阈值最大/最小之比 |
| `dose_span` | 剂量最大/最小之比 |
| `bistable_window` | `(tau_lo, tau_hi)`，L=1 下的双稳窗口；没找到则 `None` |
| `apl` | `[(配体名, tau, ERK 稳态, 相对 IL-2), ...]` |

---

## 5. 结果怎么解读

### L1（严格层）辨别指数

输出形如：

```
    N=4:   4.996 (k_off=1000)    4.636 (k_off=10)    1.078 (k_off=0.02)
```

读法：`d log C_N / d log tau`，即配体寿命翻一倍时，校对链末端产物放大多少个数量级的斜率。

- **快解离区（k_off ≫ k_p）**：指数 → N+1。N=4 时实测 4.996，也就是寿命差 2 倍
  被放大成输出差 2⁵ = 32 倍。这是校对链能区分配体的数学来源。
- **慢解离区（k_off ≪ k_p）**：指数塌回 1（N=4 时实测 1.078）。配体一旦全都能走完
  整条链，校对就不再区分任何东西了 —— 这是模型自己划的能力边界。

这一层是解析解的有限差分核验，**结论不依赖任何示意参数**，可以当硬结论用。
脚本里有 `assert` 守着：快解离区偏离 N+1 超过 0.05、或慢解离区偏离 1 超过 0.15 就会崩。
它崩了说明改坏了 `kpr_output`，不是数值噪声。

### L2（拓扑可信、参数示意）质量阈值与双稳

实测输出：

```
    ligand dose L=  0.3:  tau_threshold =   1.50 s
    ligand dose L= 10.0:  tau_threshold =   0.58 s
  -> a 33x change in ligand QUANTITY moves the quality threshold 2.58x.
    bistable window at L=1: tau in [0.02, 1.03] s (35/70 grid points)
```

读法：

- **该看的是比值，不是绝对秒数**。配体数量变 33 倍，质量阈值只挪 2.58 倍 ——
  这就是"T 细胞判的是配体质量，不是占据数"的定量表述。绝对的 1.11 s 只是示意参数的产物，
  不对应任何真实肽。
- **双稳窗口**指同一配体下，静息细胞和已激活细胞会停在不同的 ERK 态（滞后）。
  它的存在意味着"静息 vs 再刺激"不是同一条剂量-响应曲线，实验设计上不能混着比。
- 窗口边界 `[0.02, 1.03]` 是 70 点对数网格扫出来的，**精度受网格限制**，
  别把它当解析边界引用。

### L3（仅相对量）配体排序

实测输出：

```
    agonist          tau= 30.0s  ERK=0.981  IL-2(rel)=0.883
    antagonist       tau=  0.8s  ERK=0.063  IL-2(rel)=0.002
```

读法：**只有顺序是预测，纵轴不是 pg/mL**。

- 可以用的结论：agonist > weak agonist > partial > weak partial > antagonist > null，
  且 antagonist 与 weak partial 之间存在一个数量级以上的断崖（0.744 → 0.002），
  对应 ERK 开关的阈值位置。
- 不能用的结论：任何绝对分泌量；"agonist 是 antagonist 的 441 倍"这类比值也不行，
  因为它由 `il2_readout` 的 Hill 系数（n=3）决定，而这个数没有标定过。
- 表里的 `tau` 是**扫描取值**，不是哪个具名肽的实测寿命。
- 脚本用 `assert` 强制排序对 tau 单调；这条断言失败说明 L2 参数被改到了非单调区间。

### 这份输出不能拿去说的事

脚本结尾自己列了：绝对 IL-2 浓度（pg/mL）、某个具名肽的 tau、某个供体的响应、
MLR 同种反应性、ELISPOT 斑点数、前体频率、HLA 分型 —— 一个都不在里面。

---

## 6. 常见问题

| 现象 | 原因 / 处理 |
|---|---|
| `ModuleNotFoundError: No module named 'grn_pipeline'` | 不在仓库根目录跑。`cd` 过去，或设 `PYTHONPATH`。 |
| `ModuleNotFoundError: No module named 'scipy'` | 依赖没装，见 §2。 |
| 图没更新 | 图写的是仓库里的 `figures/`，不是当前目录；检查那个路径。 |
| 两次结果不一样 | 不应该发生。模型无随机性，出现差异说明代码或依赖版本被改动过。 |
| `AssertionError` | 见 §5 里对应层的说明，这是守门断言，不是偶发数值问题。 |

---

## 7. 相关文件

- `grn_pipeline/m23_virtual_tcell.py` —— 模型本体，296 行，三层都在里面。
- `VIRTUAL_TCELL_REPORT.md` —— 结论与训练数据现状（IEDB、VDJdb、CELLxGENE 全库统计）。
- `figures/m23_virtual_tcell.png` —— 本模块的输出图。
- `README.md` —— 整条流水线的模块表。

注：源码 docstring 末尾写的是 "See REPORT.md M23"，但 `REPORT.md` 里并没有 M23 小节，
对应的写作在 `VIRTUAL_TCELL_REPORT.md`。这里按实际情况指向后者。
