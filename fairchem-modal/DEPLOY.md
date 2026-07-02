# 常驻端点（Deployed Endpoints）使用说明

计算函数已 **`modal deploy`** 到你 Modal 账号（`wyh-58141`）的两个持久 app：

| App | 部署的函数 | 用途 |
|---|---|---|
| `fairchem-uma` | `_smoke` / `_mof_adsorption` / `_drug_loading` / `_polymorph_rank` / `_relax_material` | UMA 能量类 |
| `raspa-gcmc` | `_isotherm` / `_water_isotherm` | RASPA GCMC 等温线 |

**特性：** 空闲时缩到零、**不计费**；被调用时才拉起 GPU；部署是持久的——重启会话、换机器都在。管理台：`modal.com/apps/wyh-58141/main/deployed`。

---

## 怎么调用（从任何环境）

前提：该环境配好 `MODAL_TOKEN_ID` / `MODAL_TOKEN_SECRET`（你的 Modal 账号），并 `pip install modal`。**不需要**克隆计算代码——函数已在云端。

### 方式一：用封装好的 `client.py`（推荐）
```python
from client import mof_adsorption, drug_loading, polymorph_rank, water_isotherm, gas_isotherm

# MOF 吸附能
mof_adsorption("inputs/MOF-5.cif", adsorbate="H2O")

# MOF 载药亲和力
drug_loading("inputs/ZIF-8.cif", ["drugs/loading/ibuprofen.sdf",
                                  "drugs/loading/5-fluorouracil.sdf"])

# 吸水率曲线
water_isotherm("inputs/HKUST-1.cif", rh=(0.2, 0.5, 0.8))

# 气体吸附等温线
gas_isotherm("inputs/MOF-5.cif", molecule="CO2", pressures=(1e3, 1e4, 1e5))

# 药物多晶型排序
polymorph_rank(["a.cif", "b.cif"], atoms_per_molecule=20)
```
命令行快速自检：`python client.py smoke`

### 方式二：直接用 Modal SDK（连 client.py 都不需要）
```python
import modal
f = modal.Function.from_name("fairchem-uma", "_mof_adsorption")
result = f.remote(open("my_mof.cif").read(), "CO2")   # 传 CIF 文本 + 吸附质
```

### 方式三：异步 / 批量
```python
f = modal.Function.from_name("fairchem-uma", "_drug_loading")
call = f.spawn(mof_cif_text, drugs_dict)   # 立即返回，不阻塞
result = call.get()                        # 需要时再取结果
```

---

## 更新部署

改了 `modal_app.py` / `gcmc_raspa.py` 后重新部署：
```bash
modal deploy modal_app.py
modal deploy gcmc_raspa.py
```

## 停用
```bash
modal app stop fairchem-uma
modal app stop raspa-gcmc
```

---

## 可选：HTTP 端点

若想用 URL / curl / 网页前端调用（而非 Python SDK），可给函数加 `@modal.fastapi_endpoint`
包装成 HTTP 接口。因单次计算耗时数分钟，HTTP 场景建议用 `.spawn()` 返回任务 ID + 轮询，
而非同步等待。需要的话我可以再加这一层。
