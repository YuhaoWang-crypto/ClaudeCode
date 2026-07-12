# 材料学计算 Demo 集

针对聊天截图里的几项材料学需求，用**本地轻量环境真实能跑**的方式做的概念验证（demo），
以及**本地跑不了、需上 HPC** 的部分对应的方法学方案。

## 诚实的能力边界

| | 内容 | 本环境能否运行 |
|---|---|---|
| **能做** | 解析/紧束缚模型、小分子 DFT（PySCF）、方法学梳理、文献检索 | ✅ 真实端到端跑通 |
| **不能做** | 周期性材料 DFT（VASP/QE）、声子/电子-声子、GW-BSE、NEB 能垒 | ❌ 无软件+无算力，需 HPC |

环境仅有 Python 3.11 + pip（numpy/matplotlib/pyscf 现装现用），无任何第一性原理引擎；
挂载的 MCP 工具是生物医药方向，无材料计算工具。

## Demo 清单

### Demo #1 — 量子点激子精细结构劈裂（FSS）解析模型 ✅
- 文件：`qd_fss_model.py` → 图 `qd_fss_demo.png`
- 复现 Gong et al., *PRL* **106**, 227401 (2011) 的两能级激子哈密顿量（Eq.2），
  给出「单轴应力 → 极化角 → FSS」关系，对应截图里问的**图2/图3**。
- 极限情形与论文解析结果精确吻合：β=0 时 FSS_min=2|κ|，α=0 时 FSS_min=2|δ|。
- **纯 numpy，本环境真实跑出图。**

### Demo #2 — 仿生酶催化金属中心 DFT 初筛 ✅
- 文件：`metal_center_screen.py` → 图 `metal_center_screen.png`
- 对 Cr/Mn/Fe/Co/Ni/Cu 六个金属中心，用最小化 [M–OH] 活性位点模型
  （PBE/def2-SVP，UKS，扫自旋态取最低能）算电子结构初筛描述符：
  自旋基态、HOMO-LUMO 能隙、金属 Mulliken 电荷与自旋布居。
- 对应王芮需求①「结构稳定性、电子结构、电荷分布初筛」。
- 结果：**Cr 同时具最强 Lewis 酸性 + 最强金属-氧氧化还原活性 + 较大能隙**，
  定性支撑 Cr 掺杂做仿生催化中心的思路。
- ⚠️ 是**分子近似**，非真实周期性 Cr-Co₃O₄，只作工作流演示。真实材料需 VASP/QE + NEB 能垒。

### Demo #3 — 审稿意见 → 补充计算方法学清单 ✅
- 文件：`demo3_method_checklist.md`
- 把闪烁体/STE 论文审稿人的每条意见（缺声子、缺电子-声子耦合、只做基态、STE 证据不足）
  逐条映射到**具体计算 → 方法 → 软件 → 量级**，并给出性价比排序的最小补充方案。
- 核心建议：先做**声子谱(Phonopy)** + **Huang-Rhys 因子**，把"强电子-声子耦合"定量化。

### Demo #4 — 熔盐宏观性质 + DeePMD 训练集打包工作流 ✅
- 文件：`moltensalt_mlp_pipeline.py` → 图 `moltensalt_demo.png`，数据集 `deepmd_dataset/`
- 对应客户（郭硕）：熔盐 AIMD → 训 DeepMD 势 → 算密度/粘度/扩散。
- 用经典 LJ 液体作透明代理，跑 **6 个温度点** NVT MD，出 g(r)、扩散系数-温度(Arrhenius)曲线，
  并把带能量/力标签的轨迹**打包成 DeePMD-kit `deepmd/npy` 格式**（客户指定的交付格式，已校验：
  coord/box/energy/force.npy + type.raw）。
- ⚠️ 力/能来自 LJ 而非 DFT。真实工作需 CP2K/VASP AIMD + DP 训练（HPC）；把 LJ 换成 AIMD 输出即可复用打包代码。

### Demo #5 — 蛋白 Cd²⁺ 结合位点几何预测 ✅
- 文件：`cd_binding_site_predictor.py` → 图 `cd_binding_sites.png`
- 对应客户（bizlikery）：找蛋白里结合 Cd²⁺ 的位点 → 设计突变体做 MST。
- 基于 HSAB 软硬酸碱（Cd²⁺ 偏好 Cys硫 > His氮 > 羧基氧），扫描配位原子做**团(clique)聚类**打分排序，
  直接给出突变靶点。**验证**：在碳酸酐酶(1CA2)上第1位点精确命中已知金属位点 His94/His96/His119。
- ⚠️ 是几何/启发式预测，非对接或 QM。真实工作用 AutoDock(Cd²⁺参数化)/MetalPDB + QM/MM 精修 top 位点。

### Demo #6 — L/D 手性氨基酸在 Au 上 SERS 差异的分子级 DFT 论证 ✅
- 文件：`sers_chirality_dft.py` → 图 `sers_chirality_demo.png`
- 对应客户（旺仔秋秋唐）：L/D-精氨酸在金上吸附的 SERS 差异，DFT 算构型/吸附能/振动频率。
- 用 L-丙氨酸作最小手性代理，DFT(PBE/6-31G) 优化 + 镜像构造 D-对映体，算两者振动频率，
  **证明孤立分子 L≡D（谱完全重合）**——所以 SERS 差异必来自它们在金表面上的**不同吸附几何**。
  附单原子 Au 吸附能计算演示吸附能工作流。
- ⚠️ 丙氨酸/单 Au 原子是代理。真实 L/D SERS 差异需周期性 Au(111) slab + 色散 + 溶剂 + 表面 Raman 张量（HPC）——
  本 demo 恰好从分子层面论证了"为什么必须做 slab 计算"。

## 运行方式

```bash
pip install numpy matplotlib pyscf ase dpdata biopython rdkit pyberny
python3 qd_fss_model.py             # demo #1，秒级，纯解析
python3 metal_center_screen.py      # demo #2，几分钟，金属中心 DFT 初筛
python3 moltensalt_mlp_pipeline.py  # demo #4，几分钟，MD + DeePMD 打包
python3 cd_binding_site_predictor.py [PDBID]   # demo #5，秒级，默认 1CA2
python3 sers_chirality_dft.py       # demo #6，约 10 分钟，DFT 优化+频率
```

## 三次需求的能力对照

| 需求方向 | 本环境可跑的 demo | 必须上 HPC 的真实工作 |
|---|---|---|
| 量子点 FSS（图2/3） | #1 解析模型 ✅ | 原子级赝势超胞 |
| 仿生酶催化金属筛选 | #2 分子 DFT 初筛 ✅ | 周期性 Cr-Co₃O₄ + NEB 能垒 |
| 闪烁体/STE 审稿补算 | #3 方法学清单 ✅ | 声子/电子-声子/GW-BSE |
| 熔盐 ML 势 | #4 MD+DeePMD 打包 ✅ | CP2K/VASP AIMD + DP 训练 |
| 蛋白-Cd²⁺ 结合位点 | #5 几何位点预测 ✅ | AutoDock/QM-MM 精修 |
| 手性 SERS on Au | #6 分子频率论证 ✅ | 周期性 Au slab + Raman |
