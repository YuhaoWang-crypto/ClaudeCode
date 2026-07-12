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

## 运行方式

```bash
pip install numpy matplotlib pyscf
python3 qd_fss_model.py          # demo #1，秒级
python3 metal_center_screen.py   # demo #2，约几分钟（6 金属 × 多自旋态 DFT）
```
