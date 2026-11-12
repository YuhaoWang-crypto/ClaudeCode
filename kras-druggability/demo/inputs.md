# Boltz‑2 输入规格（可复现）

模型 `boltz-2.1`，工具 `structure_and_binding`，`binding = {type: ligand_protein_binding, binder_chain_id: "L"}`，`num_samples = 3`。

## 蛋白序列（KRAS 催化域 1–169，UniProt P01116）

仅第 12 位不同（P‑loop `…VVVG X GVGKS…`）：

- **WT**（X=G）
  ```
  MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK
  ```
- **G12C**（X=C）：第 12 位 `G→C` → `…VVVGACGVGKS…`
- **G12D**（X=D）：第 12 位 `G→D` → `…VVVGADGVGKS…`

## 配体 SMILES（来自 ChEMBL v34）

| 药物 | ChEMBL ID | SMILES |
|---|---|---|
| Sotorasib (AMG510) | CHEMBL4535757 | `C=CC(=O)N1CCN(c2nc(=O)n(-c3c(C)ccnc3C(C)C)c3nc(-c4c(O)cccc4F)c(F)cc23)[C@@H](C)C1` |
| Adagrasib (MRTX849) | CHEMBL4594350 | `C=C(F)C(=O)N1CCN(c2nc(OC[C@@H]3CCCN3C)nc3c2CCN(c2cccc4cccc(Cl)c24)C3)C[C@@H]1CC#N` |
| MRTX1133 | CHEMBL4858364 | `C#Cc1c(F)ccc2cc(O)cc(-c3ncc4c(N5CC6CCC(C5)N6)nc(OC[C@@]56CCCN5C[C@H](F)C6)nc4c3F)c12` |

## 体系 → prediction_id

| 体系 | prediction_id | idempotency_key |
|---|---|---|
| G12C + Sotorasib | `sab_pred_iuG44c2Kzi14J0ROeAA6` | `kras-demo-g12c-sotorasib-v1` |
| G12C + Adagrasib | `sab_pred_IVB9HLHv7cPQhFOakBUD` | `kras-demo-g12c-adagrasib-v1` |
| G12D + MRTX1133 | `sab_pred_aMxstjpy8MHiWcBSeqPE` | `kras-demo-g12d-mrtx1133-v1` |
| WT + MRTX1133 | `sab_pred_4tZO9yUiGvOCRAR0wCxS` | `kras-demo-wt-mrtx1133-v1` |

## ChEMBL 实验锚点

| 药物 | 靶点 / 细胞 | 活性 | 说明 |
|---|---|---|---|
| Sotorasib | GTPase KRas (CHEMBL2189121) | IC50 = 68 nM (pChEMBL 7.17) | 生化 |
| Sotorasib | MIA PaCa‑2 (KRAS **G12C**) | IC50 = 5 nM | 细胞，突变敏感 |
| Sotorasib | A549 (KRAS **G12S**，非 G12C) | IC50 = 36500 nM | 细胞，突变选择性对照 |
| MRTX1133 | KRAS **G12D** (1–169) TR‑FRET | IC50 = **0.4 nM** (pChEMBL 9.4) | 与本 demo 构建体一致 |
| MRTX1133 | AGS (KRAS **G12D**) | IC50 = 6 nM | 细胞 |
| MRTX1133 | MKN‑1 (KRAS **WT**) | IC50 **> 3000 nM** | 细胞，选择性对照 |

> 注：Adagrasib 的 ChEMBL 活性未在本次 demo 中单独拉取；它是 2022 年获批的第二个共价 KRAS **G12C** 抑制剂，此处仅作 G12C 共价结合的结构交叉验证。
