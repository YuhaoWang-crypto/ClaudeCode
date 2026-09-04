# grn-pipeline

A small, fully-runnable pipeline that applies four "irreducibility / symmetry"
mathematical tools to gene-regulatory and metabolic networks, on concrete
literature-grounded systems where every number is *computed*, not asserted.

| Module | Tool | System | Key result |
|---|---|---|---|
| `m1_symmetry` | graph automorphism → quotient | RTK/RAS/RAF/MEK/ERK | \|Aut\|=S₃; 9→7 nodes (3 RAS paralogues → 1 core) |
| `m2_crnt` | CRNT deficiency δ | A⇌B⇌C vs Schlögl | δ=0 monostable / δ=1 bistable switch |
| `m3_efm` | elementary flux modes | 4-metabolite network | 3 irreducible flux generators span the cone |
| `m4_dnb_lyapunov` | DNB / critical slowing / Lyapunov | 2-gene fold bifurcation | LLE→0, SD/autocorr/DNB rise at tipping point |
| `m5_kras_real` | symmetry breaking on a real target | KRAS G12C + covalent drugs (ChEMBL/Boltz/Inductive Bio) | covalent G12C drug breaks paralog symmetry S₃(6)→S₂(2) |
| `m6_integrate` | binding → network stability | sotorasib vs adagrasib | real ChEMBL+Boltz binding → engagement → DNB biomarker |
| `m7_screen` | Boltz-2.1 library screen | 10 G12C ligands (ChEMBL) | ranked by binding; 4 analogues out-rank sotorasib |
| `m8_clinical` | biomarker → trial endpoints | CodeBreaK 100 (NCT03600883) | layers mapped to ORR/DOR, PFS/OS, Cmax/AUC, QTc |
| `m9_occupancy` | PK occupancy → μ calibration | sotorasib (IC50=30 nM) | 98% occupancy at approved dose → network near tipping |
| `m10_validate` | Boltz ranking vs ChEMBL truth | 5 G12C ligands w/ measured IC50 | opt_score tracks potency (ρ=+0.6); binding_confidence doesn't (ρ=−0.2) |
| `m11_fibration` | input-tree fibration (Morone) | expanded MAPK paralogue graph | 27→11 fibers; generalises M1 automorphism to fiber representatives |
| `m12_dualphos` | real ERK double-phospho core | Markevich-style mass-action | CRNT deficiency δ=2; bistable ERK switch; EFM = 2 futile cycles |
| `m13_fim_sloppy` | FIM / sloppy / stiff axes | ERK dual-phospho ODE | sloppy spectrum (38 orders); flux-ratio observables load best |
| `m14_atlas` | 18-pathway systematic atlas | JAK-STAT…mevalonate | fibration compression + biomarker class per pathway; JAK-STAT top (3.0×) |
| `m15_markevich_mm` | exact Markevich 2004 MM ERK cycle | published parameters (JCB 2004) | reproduces bistable window [39.25, 57.38] nM + 3-state table to the decimal |
| `m16_erk_dnb` | DNB / critical slowing on the real switch | M15 saddle-nodes 39.25/57.38 nM | λ_max→0, τ≈4720 s at boundaries; SD/autocorr/DNB rise (early warning) |
| `m17_realdata` | validate M16 on real single-cell ERK imaging | Pertz-lab EKAR traces (FGF pulses + EGF dose) | lag-1 autocorr rises before ERK pulses (p≈0.004); variance flat; EGF all supra-threshold — partial validation |
| `m18_titration_benchmark` | positive control: MEKi titration across the real bifurcation | simulated from M15 Markevich switch | variance & lag-1 autocorr PEAK near threshold (≈5×), tracking τ — pipeline is sensitive, not blind |
| `m19_switch_library` | migrate the critical-slowing engine to many pathways | 10 canonical bistable switches (MAPK, Rb-E2F, apoptosis, Cdc2, CaMKII, Wnt, Cdc42, lac, Schlögl, master-TF) | 10/10 show variance+autocorr rising to their saddle-node — biomarker is universal to the bifurcation |
| `m20_literature_bistable` | multi-variable literature switches + hysteresis | Rb-E2F (Yao 2008), apoptosis (Eissing 2004) topology | bistability + hysteresis loops reproduced; eigenvalue→0 at folds (exact params not open-access-fetchable) |
| `m21_oscillators` | extend framework to oscillatory (Hopf) pathways | Goodwin (circadian), p53-Mdm2, Brusselator (glycolytic) | approaching Hopf: variance rises AND a spectral peak sharpens at the intrinsic frequency — distinct from saddle-node |
| `m20b_biomodels_exact` | fetch + simulate EXACT curated models (fills M20 gap) | Markevich2004 (BIOMD27), Legewie2006 apoptosis (BIOMD102) | download method = biomodels GitHub mirror + libRoadRunner; official Km5=78 confirms hand-coded M15 (states to the decimal); Legewie caspase switch bistable in XIAP synthesis |
| `m22_snic_mixed` | mixed bifurcation: saddle-node ON a limit cycle (SNIC) | θ / Ermentrout-Kopell normal form (cell-cycle / excitable) | finite-amplitude spikes whose period diverges (T~π/√I, log-log slope −0.50; frequency→0) — signature distinct from both Hopf and pure saddle-node; ISI mean+CV both grow |

## Run

```bash
pip install numpy scipy networkx matplotlib
python3 -m grn_pipeline.run_all       # full pipeline + figures
python3 -m grn_pipeline.m1_symmetry   # or any single module
```

Figures are written to `figures/`. A full write-up with numbers, rigour
labels, and the interpretation (including the Lyapunov-exponent biomarker
question) is in [`REPORT.md`](REPORT.md).

---

## `assaysim` — 湿实验数字孪生（已验证，非 mock）

把诊断/抗病毒湿实验抽象成可运行、可证伪的模型。全部验证：

```bash
pip install numpy scipy biopython rdkit
python3 -m assaysim.validate      # 47/47 通过
```

| 模块 | 模型 | 验证方式 | 结果 |
|---|---|---|---|
| `nn_thermo` | Allawi/SantaLucia 最近邻双链热力学 | 与 Biopython 独立实现交叉比对 | 200 条寡核苷酸 max\|ΔTm\| = 0.003 °C |
| `neutralization` | 多击中占据模型 Kd → NT50 | 数值解 vs 闭式解 `Kd·(2^(1/n)−1)` | 相对偏差 7e-14；n=100 时 NT50 = Kd/143.77 |
| `viral_dynamics` | 四状态靶细胞受限 ODE + 药物机制接口 | 守恒律、雅可比特征值、R₀ 阈值、KM 最终规模 | 全部通过 |
| `bridge` | 非细胞 → 细胞效力桥接 | ChEMBL_37 真实配对 + Murcko 骨架留出 | 见下 |
| `agid` | 琼脂免疫扩散反应-扩散 PDE | 求解器 vs erfc 解析解；Stokes-Einstein vs 实测 D | 扩散最大偏差 9.8e-6；R_h 复现实测值在 1.4% 内 |
| `agglutination` | Smoluchowski 聚集 + 抗体桥联 | 数值解 vs 常核闭式解 | 簇分布相对偏差 6.5e-9；质量守恒 1.3e-10 |

**真实数据结论**（ChEMBL_37 实时抓取，缓存于 `data/`）：

- HIV-1 RT 有 **2067** 个化合物同时具备酶法 IC50 与细胞法 EC50；流感 NA 有 **90** 个
- 线性桥接留出 RMSE **1.03 log** vs 直接当细胞值用的 **1.20 log**
- 细胞法终点自身的跨论文重现性约 **0.76 log** —— 预测误差已在噪声下限的 1.4 倍以内
- 在 HIV-1 RT 上标定的桥接迁移到流感 NA：RMSE **1.49**、系统性偏差 **+0.84 log**；
  就地重拟合则为 **1.01**。这把"结构可迁移、参数不可迁移"变成了可测量的量

AGID 与凝集的两个要点（都是模型的推论，不是外加假设）：

- **双向免疫扩散不会出现前带/后带。** 两孔是恒浓度储库，凝胶内浓度比从 ∞ 连续
  扫到 0，等价带总能找到位置——线只移动不消失（抗原越浓越靠近抗体孔）。真正的
  钟形前带/后带属于**试管沉淀**这种封闭体系：抗原过量时抗体被消耗进**不沉淀的
  可溶复合物**。`tube_precipitin()` 的沉淀峰精确落在化学计量等价点 A/B = 1/ν。
- **凝集的最优抗体浓度就是 Kd。** 桥联需要一端已结合抗体、另一端仍是空表位，
  效率 ∝ 4θ(1−θ)，θ=0.5 时取 1。所以定量凝集必须做系列稀释——单一浓度阴性
  无法区分抗体太少还是太多。

交互式工作台：`tools/twin_template.html`（用 `data/` 里的真实数据注入后发布）。
服务线的 in-silico 可行性分级见 [`BIOVENIC_INSILICO.md`](BIOVENIC_INSILICO.md)。

### ASFV p72 案例与跨实现交叉验证

```bash
python3 -m assaysim.case_asfv_p72                 # 全链路逐段跑，先判适用性
python3 -m assaysim.crossvalidate_av <bundle_dir> # 与 VC-NEU / VC-AGD 对拍
```

`case_asfv_p72` 用真实输入（UniProt P22776 · 646 aa · 73,179 Da；衣壳 T=277、
p72 8,280 拷贝，Liu 2019 Cell Res PMID 31649031）走完六段链路，并**先判定适用性
再计算**：2 段直接适用（AGID 浓度窗口、凝集稀释窗口），2 段结构适用但参数缺，
1 段明确不出预测（ChEMBL_37 里 ASFV 的 pchembl 活性为 **0 条**）。

`crossvalidate_av` 把 assaysim 与一套**独立编写**的浏览器端实现对拍
（不同语言、不同 lgamma、不同求根、不同单位约定）：

- 占据模型 NT50：11 组参数最大相对偏差 **1.05e-11**
- Stokes-Einstein 扩散层：R_h 与 D 均吻合到 **机器精度 (2e-16)**
- 靶细胞受限 ODE（Baccam 2006 参数）：峰值滴度最大相对偏差 **8.2e-4**

**吸附项与参数集是一个整体。** Baccam 2006 的方程 `dV/dt = p·I − c·V` 不含吸附损耗
`−βTV`，其 β 是在这个形式上拟合的。该参数集下 β·T0 = 12,800/day 而 c = 5.2/day，
相差 2462 倍——把它塞进含吸附项的方程，R0 会从 **21.8 掉到 0.009**，感染直接熄灭。
`CellVirusParams.absorption` 因此是显式开关，默认 True（我的形式），引用文献参数时
必须设为 False。本仓库用 `absorption=False` 复现了 Baccam 的 R0=21.78（文献 21.5）、
峰值第 2.33 天、以及"延迟给药"结论（第 1.5/2.0/2.5 天给药，峰值降低 98.9%/52.9%/0%）。

对拍过程中采纳了对方由偏比容 v̄=0.73 cm³/g **推导**无水球半径的做法，
替换掉原先四舍五入的经验常数 0.066，R_h 对实测值的复现从 10% 容差收紧到 1.4%。
