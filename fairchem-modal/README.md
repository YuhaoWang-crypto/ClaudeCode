# FAIRChem (UMA) on Modal — 材料 / 药物结晶 / MOF 吸附

用 Meta FAIR 的 **UMA** 通用原子间势（HuggingFace `facebook/UMA`）在 **Modal GPU** 上做：

| 任务 | UMA task head | 说明 |
|---|---|---|
| MOF 吸附（H₂O / CO₂） | `odac` | Open Direct Air Capture，专为 MOF + 气体吸附训练 |
| 药物多晶型能量排序 | `omc` | 有机分子晶体（organic molecular crystals） |
| 无机晶体材料弛豫 | `omat` | bulk 材料 |
| 分子单点/校验 | `omol` | 隔离分子 |

> Claude Code 的容器是 **CPU-only**，所以所有重计算都 offload 到 Modal 的 GPU。

---

## 一次性配置

### 1. 装 Modal CLI（本地容器里已经能装）
```bash
pip install modal
```

### 2. Modal 认证 —— 用**环境变量**（你选的方式，最安全）
在 **Claude Code on the web 的环境配置**里加两个变量，然后重启会话：
```
MODAL_TOKEN_ID=ak-...
MODAL_TOKEN_SECRET=as-...
```
（在 Modal 后台 Settings → API Tokens 新建一对**专用** token；不要用账号密码，不要贴进聊天。）

验证：
```bash
modal token verify   # 或任意 modal 命令能连上 workspace 即可
```

### 3. 放 HuggingFace token（UMA 是 gated 模型）
你的 HF token 需要**已接受 `facebook/UMA` 的 license**。存成 Modal secret：
```bash
modal secret create huggingface HF_TOKEN=hf_xxxxxxxx
```

---

## 跑起来

```bash
# 0) 冒烟测试：最便宜地验证「Modal GPU + gated UMA 权重 + fairchem 推理」整条链路
modal run modal_app.py::smoke_test

# 1) MOF 吸附能（默认 H2O，可 --adsorbate CO2）
modal run modal_app.py::mof_adsorption --cif your_mof.cif --adsorbate H2O

# 2) 药物多晶型排序（目录下放多个 .cif 候选晶型）
modal run modal_app.py::polymorph_rank --cif-dir ./polymorphs/

# 3) 无机晶体弛豫 + 能量
modal run modal_app.py::relax_material --cif your_crystal.cif
```

第一次跑先做 `smoke_test`：它只算一个水分子的能量，几秒钟、几分钱，用来确认
GPU、gated 权重下载、fairchem API 都通。通了再喂真实体系。

---

## 现实预期 / 坑（重要）

- **UMA 给的是能量和力，不是宏观"吸水率%"。** `mof_adsorption` 算的是单个吸附
  位点的吸附能（负值=有利吸附）。
- **要宏观吸附等温线 / uptake 曲线**：需要在 MLIP 之上再叠 **GCMC**（巨正则
  蒙特卡洛，典型工具 RASPA / GPU 版 gRASPA）。路线：UMA 提供能量 → GCMC 采样 →
  出等温线。这是下一步扩展，脚手架里先没含（GCMC 建议单独跑，或用 MLIP 拟合力场）。
- **`mof_adsorption` 目前是朴素实现**：把吸附质放在孔中心、固定骨架弛豫。做筛选要
  采样多个吸附位点/朝向再取最稳的。
- **药物多晶型**：多晶型能差常 <1 kJ/mol，逼近 MLIP 精度极限。当**初筛**用，
  排在前面的候选一定要用 **DFT-D（含色散校正）** 复核，别只信 MLIP 排序。
- **模型 tag / task 名**：代码里用的是 fairchem-core v2 API（`uma-s-1` +
  `get_predict_unit`）。首次运行如报模型名/task 名错误，按你 HF 上实际拿到的
  权重名改 `MODEL_TAG`。

---

## GPU / 费用

`modal_app.py` 默认用 `A10G`（便宜）。大 MOF 超胞或慢收敛时把函数装饰器里的
`gpu="A10G"` 改成 `"A100"`。每个任务开跑前会在终端打印用的卡型和体系大小。
