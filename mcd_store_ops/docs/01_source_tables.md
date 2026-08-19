# 01 · 一家门店每天在用的表（数据字典）

> 本文件由 `python -m mcdops.cli docs` 从 `src/mcdops/schema.py` 生成，不要手改。

共 30 张表，分 9 个业务域。
门店的现实是：这些表分散在 6~8 个互不相通的系统里，
把它们对齐到同一个 `store_id` + `biz_date`，是做任何报表的前提工作，
也是最容易出错的地方（尤其是跨零点的夜宵订单该算哪一天）。

## 交易

### `fct_pos_transaction` — POS 订单头

- **层级** fct · **来源系统** POS 收银系统 · **刷新** 近实时（15 分钟增量）
- **粒度** 1 行 = 1 张订单 · **主键** store_id, order_id · **日增行数** 1500 ~ 3000
- **为什么需要它**：所有生意类指标的地基：销售、交易笔数、客单价、渠道结构、时段曲线

| 字段 | 类型 | 说明 |
|---|---|---|
| `order_id` | string | 订单号 |
| `store_id` | string | 门店号 |
| `biz_date` | date | 营业日，跨零点的夜宵单仍算前一天 ← 最容易错的口径 |
| `order_ts` | timestamp | 下单时间 |
| `interval15` | string | 15 分钟时段键，如 '18:15'，看板时段矩阵的行 |
| `daypart` | string | 早餐/午市/下午茶/晚市/夜宵 |
| `channel` | string | FrontCounter/Kiosk/DriveThru/MobileApp/Delivery |
| `gross_amount` | decimal | 原价金额 |
| `discount_amount` | decimal | 折扣（券/会员/员工餐） |
| `refund_amount` | decimal | 退款 |
| `net_amount` | decimal | 净销售 = gross - discount - refund |
| `pay_type` | string | 现金/微信/支付宝/银行卡/平台代收 |
| `cashier_id` | string | 收银员，附加销售率下钻到人 |
| `status` | string | completed/voided/refunded |

### `fct_pos_line_item` — POS 订单行

- **层级** fct · **来源系统** POS 收银系统 · **刷新** 近实时
- **粒度** 1 行 = 订单里的 1 个商品 · **主键** store_id, order_id, line_no · **日增行数** 5000 ~ 10000
- **为什么需要它**：品类结构、单均件数、附加销售、理论食品成本都只能从这里算

| 字段 | 类型 | 说明 |
|---|---|---|
| `order_id` | string | 关联订单头 |
| `line_no` | int | 行号 |
| `sku_id` | string | 商品编码 |
| `category` | string | 汉堡/鸡类/薯条小食/饮料/麦咖啡/甜品派 |
| `qty` | int | 数量 |
| `unit_price` | decimal | 单价 |
| `line_net` | decimal | 行净额（拆分套餐折扣后） |
| `is_combo_component` | bool | 是否套餐组成部分 |
| `is_upsell` | bool | 是否附加销售识别出的加购/升级 |
| `promo_id` | string | 促销活动号 |

### `fct_payment` — 支付明细

- **层级** fct · **来源系统** POS / 聚合支付 · **刷新** 近实时
- **粒度** 1 行 = 1 笔支付（一单可拆多笔） · **主键** store_id, payment_id · **日增行数** 1600 ~ 3200
- **为什么需要它**：现金短溢、平台回款对账、支付渠道手续费

| 字段 | 类型 | 说明 |
|---|---|---|
| `payment_id` | string | 支付流水号 |
| `order_id` | string | 订单号 |
| `pay_type` | string | 支付方式 |
| `amount` | decimal | 金额 |
| `settle_status` | string | 已结算/在途/差异 |

### `dim_product` — 商品主数据

- **层级** dim · **来源系统** 总部主数据 · **刷新** 总部下发，周级
- **粒度** 1 行 = 1 个 SKU 的 1 个生效版本 · **主键** sku_id, valid_from · **日增行数** —（约 300 个在售 SKU）
- **为什么需要它**：标准成本 → 理论食品成本；有它才能把'实际 FC%'拆成'理论差'和'损耗差'

| 字段 | 类型 | 说明 |
|---|---|---|
| `sku_id` | string | 商品编码 |
| `name_cn` | string | 品名 |
| `category` | string | 品类 |
| `price` | decimal | 建议零售价 |
| `standard_cost` | decimal | 标准物料成本（BOM 展开） |
| `valid_from/valid_to` | date | 生效区间，做 SCD2 |

### `dim_promotion` — 促销主数据

- **层级** dim · **来源系统** 总部营销 · **刷新** 周级
- **粒度** 1 行 = 1 个活动 · **主键** promo_id · **日增行数** —
- **为什么需要它**：解释销售异动：涨了是活动还是本事，跌了是去年有活动今年没有

| 字段 | 类型 | 说明 |
|---|---|---|
| `promo_id` | string | 活动号 |
| `name_cn` | string | 活动名 |
| `mechanic` | string | 第二份半价/买赠/券 |
| `start_date/end_date` | date | 档期 |
| `funding` | string | 总部承担/门店承担 |

## 数字与外送

### `fct_delivery_order` — 外送订单

- **层级** fct · **来源系统** 美团/饿了么/麦乐送 中台 · **刷新** 近实时
- **粒度** 1 行 = 1 张外送订单 · **主键** platform, platform_order_id · **日增行数** 300 ~ 800
- **为什么需要它**：外送已是四分之一的生意，但它的时效、差错、评分和堂食完全不同源，必须单独看

| 字段 | 类型 | 说明 |
|---|---|---|
| `platform` | string | Meituan/Eleme/McDeliveryApp |
| `platform_order_id` | string | 平台单号 |
| `pos_order_id` | string | 回落到 POS 的单号，用于对账 |
| `accept_ts` | timestamp | 门店接单 |
| `ready_ts` | timestamp | 出餐完成 |
| `rider_pickup_ts` | timestamp | 骑手取走 |
| `delivered_ts` | timestamp | 送达 |
| `commission_amount` | decimal | 平台佣金 |
| `rating` | int | 顾客评分 1-5 |
| `issue_flag` | string | 漏餐/洒漏/超时/无 |

### `fct_kiosk_session` — 自助点餐机会话

- **层级** fct · **来源系统** 点餐机后台 · **刷新** 日批
- **粒度** 1 行 = 1 次触屏会话 · **主键** kiosk_id, session_id · **日增行数** 500 ~ 1200
- **为什么需要它**：点餐机放弃率是隐形流失：屏幕坏了/网络慢/菜单卡顿，POS 上完全看不出来

| 字段 | 类型 | 说明 |
|---|---|---|
| `session_id` | string | 会话号 |
| `kiosk_id` | string | 设备号 |
| `start_ts/end_ts` | timestamp | 起止 |
| `converted` | bool | 是否成单 |
| `abandon_step` | string | 放弃在哪一步：选餐/支付/等待 |

## 出品时效

### `fct_kds_order` — 厨房显示系统订单

- **层级** fct · **来源系统** KDS · **刷新** 近实时
- **粒度** 1 行 = 1 张订单在 1 个工作站的一段 · **主键** order_id, station, fire_ts · **日增行数** 2000 ~ 4000
- **为什么需要它**：把'慢'拆开：是点单慢、厨房慢、还是交付慢。没有 KDS 就只能靠感觉

| 字段 | 类型 | 说明 |
|---|---|---|
| `order_id` | string | 订单号 |
| `station` | string | 煎炉/炸炉/理货/饮料/打包 |
| `fire_ts` | timestamp | 进单 |
| `bump_ts` | timestamp | 敲单完成 |
| `prep_seconds` | int | bump - fire |
| `item_count` | int | 件数，慢单要先排除大单干扰 |
| `recall_flag` | bool | 是否被召回重做 |

### `fct_dt_timer` — 得来速计时

- **层级** fct · **来源系统** DT 计时器（地感线圈 + 计时主机） · **刷新** 近实时
- **粒度** 1 行 = 1 辆车的一次通过 · **主键** store_id, car_sequence, biz_date · **日增行数** 300 ~ 700
- **为什么需要它**：得来速是分段计时：菜单板→点单→缴费→取餐→离开，堵在哪一段，动作完全不同

| 字段 | 类型 | 说明 |
|---|---|---|
| `car_sequence` | int | 当日车序 |
| `lane` | string | 单车道/双车道 A/B |
| `arrive_ts` | timestamp | 进入车道 |
| `order_start_ts / order_end_ts` | timestamp | 点单段 |
| `window_arrive_ts` | timestamp | 到取餐窗 |
| `depart_ts` | timestamp | 驶离 |
| `total_seconds` | int | depart - arrive，看板主指标 |
| `window_seconds` | int | depart - window_arrive，最能反映后厨备餐 |
| `pull_forward_flag` | bool | 是否被请到等餐位（会人为压低窗口时长，看数要剔除） |

### `fct_order_issue` — 订单差错与投诉

- **层级** fct · **来源系统** POS 退换 + 客服工单 + 平台售后 · **刷新** 日批
- **粒度** 1 行 = 1 个差错事件 · **主键** issue_id · **日增行数** 5 ~ 40
- **为什么需要它**：订单准确率、千单投诉数的唯一来源；按 reason 归类才知道是配料错还是漏给

| 字段 | 类型 | 说明 |
|---|---|---|
| `issue_id` | string | 事件号 |
| `order_id` | string | 关联订单 |
| `source` | string | 店内退换/顾客热线/平台售后/线上评价 |
| `issue_type` | string | 漏餐/错餐/异物/温度/态度/超时 |
| `amount` | decimal | 补偿金额 |
| `channel` | string | 发生渠道 |

## 人员

### `dim_employee` — 员工主数据

- **层级** dim · **来源系统** 人事系统 · **刷新** 日批
- **粒度** 1 行 = 1 名员工 · **主键** emp_id · **日增行数** —（约 60~90 人在册）
- **为什么需要它**：工时成本要按人算时薪；流失率、认证覆盖都挂在这张表上

| 字段 | 类型 | 说明 |
|---|---|---|
| `emp_id` | string | 工号 |
| `role` | string | 店长/值班经理/组长/训练员/员工 |
| `hire_date / term_date` | date | 入离职，算 30 日流失 |
| `wage_per_hour` | decimal | 时薪 |
| `employment_type` | string | 全职/兼职/小时工 |

### `fct_schedule` — 排班计划

- **层级** fct · **来源系统** 排班系统 · **刷新** T-7 生成，日内可调
- **粒度** 1 行 = 1 人 1 个班次 · **主键** sched_id · **日增行数** 40 ~ 70
- **为什么需要它**：所有'该有多少人'的分母。排班 vs 实际的差，是人工成本失控的第一现场

| 字段 | 类型 | 说明 |
|---|---|---|
| `sched_id` | string | 排班号 |
| `emp_id` | string | 工号 |
| `biz_date` | date | 营业日 |
| `station` | string | 岗位 |
| `start_ts / end_ts` | timestamp | 班次起止 |
| `planned_hours` | decimal | 计划工时 |

### `fct_timecard` — 考勤打卡

- **层级** fct · **来源系统** 考勤机 / 排班系统 · **刷新** 近实时
- **粒度** 1 行 = 1 人 1 次上下班 · **主键** emp_id, clock_in_ts · **日增行数** 40 ~ 70
- **为什么需要它**：实际工时 → 人工成本率、SPLH；迟到与出勤 → 高峰为什么撑不住

| 字段 | 类型 | 说明 |
|---|---|---|
| `emp_id` | string | 工号 |
| `clock_in_ts / clock_out_ts` | timestamp | 打卡 |
| `actual_hours` | decimal | 实际工时（含餐休扣减） |
| `late_minutes` | int | 迟到分钟 |
| `station` | string | 实际所站岗位 |

### `fct_training_record` — 培训与岗位认证

- **层级** fct · **来源系统** 学习平台 · **刷新** 周批
- **粒度** 1 行 = 1 人 1 项认证 · **主键** emp_id, cert_code · **日增行数** —
- **为什么需要它**：高峰慢、废弃高、食安失分，往往先在'谁没认证却在站岗'这里露头

| 字段 | 类型 | 说明 |
|---|---|---|
| `emp_id` | string | 工号 |
| `cert_code` | string | 认证项：煎炉/炸炉/DT/食安基础 |
| `passed_at` | date | 通过日期 |
| `expires_at` | date | 有效期 |

## 库存

### `dim_item` — 原料主数据

- **层级** dim · **来源系统** 供应链系统 · **刷新** 周批
- **粒度** 1 行 = 1 个原料 · **主键** item_id · **日增行数** —（约 200 个原料）
- **为什么需要它**：盘点、订货、废弃全部以原料为单位；保质期字段决定先进先出怎么排

| 字段 | 类型 | 说明 |
|---|---|---|
| `item_id` | string | 原料码 |
| `name_cn` | string | 名称 |
| `uom` | string | 计量单位 |
| `standard_cost` | decimal | 标准单价 |
| `shelf_life_days` | int | 保质期 |
| `storage_type` | string | 常温/冷藏/冷冻 |

### `fct_goods_receipt` — 收货记录

- **层级** fct · **来源系统** 订货收货系统 · **刷新** 每次送货
- **粒度** 1 行 = 1 次收货的 1 个原料 · **主键** receipt_id, item_id · **日增行数** 0 ~ 120（隔日送货）
- **为什么需要它**：食品成本的进项；收货温度不合格是食安一票否决项

| 字段 | 类型 | 说明 |
|---|---|---|
| `receipt_id` | string | 收货单号 |
| `item_id` | string | 原料 |
| `qty / cost` | decimal | 数量与金额 |
| `receive_ts` | timestamp | 收货时间 |
| `temp_c` | decimal | 到货温度 |
| `temp_ok` | bool | 是否在标准区间 |

### `fct_inventory_count` — 盘点

- **层级** fct · **来源系统** 门店管理工作站 · **刷新** 日盘（关键品）+ 周盘（全品）
- **粒度** 1 行 = 1 次盘点的 1 个原料 · **主键** count_id, item_id · **日增行数** 20 ~ 200
- **为什么需要它**：实际食品成本 =（期初 + 进货 − 期末）／销售。没有盘点，FC% 就是猜的

| 字段 | 类型 | 说明 |
|---|---|---|
| `count_id` | string | 盘点单 |
| `count_type` | string | 日盘/周盘/月盘 |
| `item_id` | string | 原料 |
| `qty_counted` | decimal | 实盘量 |
| `qty_theoretical` | decimal | 理论量（由 BOM × 销量推） |
| `variance_cost` | decimal | 差异金额 ← 真正要追的数 |

### `fct_waste_log` — 废弃登记

- **层级** fct · **来源系统** 门店管理工作站 / 手工录入 · **刷新** 班次录入
- **粒度** 1 行 = 1 次废弃 · **主键** waste_id · **日增行数** 20 ~ 80
- **为什么需要它**：废弃分'原料废弃'（过期、掉地）和'成品废弃'（做多了、保温超时），两者的原因和对策完全不同

| 字段 | 类型 | 说明 |
|---|---|---|
| `waste_id` | string | 记录号 |
| `waste_type` | string | raw 原料 / complete 成品 / holding 保温超时 |
| `item_or_sku` | string | 原料或商品 |
| `qty / cost` | decimal | 数量与成本 |
| `reason` | string | 过期/做多/掉落/顾客退/设备故障 |
| `logged_ts` | timestamp | 登记时间 |
| `logged_by` | string | 登记人 ← 漏登记比废弃本身更可怕 |

### `fct_holding_batch` — 保温批次

- **层级** fct · **来源系统** 生产管理 / 保温柜计时 · **刷新** 近实时
- **粒度** 1 行 = 1 批出品 · **主键** batch_id · **日增行数** 150 ~ 400
- **为什么需要它**：在'新鲜'和'快'之间的那个平衡点，就写在这张表里

| 字段 | 类型 | 说明 |
|---|---|---|
| `batch_id` | string | 批次 |
| `sku_id` | string | 商品 |
| `cook_ts / expire_ts` | timestamp | 出品与到期 |
| `qty_made / qty_sold / qty_discard` | int | 做了/卖了/丢了 |
| `triggered_by` | string | 预测建议 / 人工判断 |

### `dim_sales_forecast` — 销售预测

- **层级** dim · **来源系统** 需求预测系统 · **刷新** T-1 生成，日内滚动
- **粒度** 1 行 = 1 个 15 分钟时段 × 渠道 · **主键** biz_date, interval15, channel · **日增行数** 约 360
- **为什么需要它**：排班、订货、备货全部由它驱动。预测不准 → 要么排多了亏人工，要么排少了亏速度

| 字段 | 类型 | 说明 |
|---|---|---|
| `biz_date` | date | 营业日 |
| `interval15` | string | 时段 |
| `channel` | string | 渠道 |
| `forecast_sales / forecast_gc` | decimal | 预测销售与笔数 |
| `forecast_version` | string | 版本，做事后准确度复盘 |

## 食安与品质

### `fct_food_safety_check` — 食安检查表

- **层级** fct · **来源系统** 食安检查 App · **刷新** 每班（开店/交接/闭店）
- **粒度** 1 行 = 1 次检查的 1 个检查项 · **主键** check_id, item_code · **日增行数** 60 ~ 150
- **为什么需要它**：唯一一个'红了就必须立刻停下手上一切事'的域。它不是 KPI，是准入条件

| 字段 | 类型 | 说明 |
|---|---|---|
| `check_id` | string | 检查单 |
| `shift` | string | 早班/中班/晚班 |
| `item_code` | string | 检查项（洗手、消毒液浓度、交叉污染…） |
| `result` | string | pass/fail/na |
| `is_critical` | bool | 是否关键项，关键项一票否决 |
| `checker_emp_id` | string | 检查人 |
| `corrective_action` | string | 整改动作 |

### `fct_temperature_log` — 温度记录

- **层级** fct · **来源系统** IoT 温感 + 手工探针 · **刷新** IoT 每 5 分钟 / 手工每 2 小时
- **粒度** 1 行 = 1 次测温 · **主键** log_id · **日增行数** 200 ~ 3000
- **为什么需要它**：冷链断了不会有人喊，只有这张表会说话

| 字段 | 类型 | 说明 |
|---|---|---|
| `log_id` | string | 记录号 |
| `equipment_id` | string | 设备：冷藏库/冷冻库/保温柜/油炸锅 |
| `ts` | timestamp | 时间 |
| `temp_c` | decimal | 温度 |
| `in_range` | bool | 是否达标 |
| `source` | string | iot/manual |

### `fct_travel_path_check` — 巡店检查（清洁与体验）

- **层级** fct · **来源系统** 巡店 App · **刷新** 每班 2~3 次
- **粒度** 1 行 = 1 次巡店的 1 个区域 · **主键** walk_id, zone · **日增行数** 8 ~ 20
- **为什么需要它**：顾客对'干净'的判断在洗手间和餐桌，而不在厨房；分区打分才知道派谁去

| 字段 | 类型 | 说明 |
|---|---|---|
| `walk_id` | string | 巡店号 |
| `zone` | string | 大堂/洗手间/厨房/外围/得来速车道 |
| `score` | int | 0-100 |
| `issue_note` | string | 问题描述 |
| `walk_ts` | timestamp | 时间 |

### `dim_equipment` — 设备台账

- **层级** dim · **来源系统** 资产/维保系统 · **刷新** 变更时
- **粒度** 1 行 = 1 台设备 · **主键** equipment_id · **日增行数** —（约 40 台）
- **为什么需要它**：把停机时长挂到具体设备上，才能判断是该修还是该换

| 字段 | 类型 | 说明 |
|---|---|---|
| `equipment_id` | string | 设备号 |
| `name_cn` | string | 名称 |
| `is_critical` | bool | 停了是否直接停售 |
| `temp_range_c` | array | 标准温度区间 |
| `install_date` | date | 启用日期 |
| `pm_cycle_days` | int | 保养周期 |

### `fct_equipment_ticket` — 设备工单

- **层级** fct · **来源系统** 维保工单系统 · **刷新** 近实时
- **粒度** 1 行 = 1 张工单 · **主键** ticket_id · **日增行数** 0 ~ 5
- **为什么需要它**：'今天薯条为什么慢'的答案，一半时候在这张表里（一台炸炉在修）

| 字段 | 类型 | 说明 |
|---|---|---|
| `ticket_id` | string | 工单号 |
| `equipment_id` | string | 设备 |
| `open_ts / close_ts` | timestamp | 报障与修复 |
| `downtime_minutes` | int | 停机时长 |
| `severity` | string | 停售/降级/不影响 |
| `vendor` | string | 服务商 |

## 顾客

### `fct_guest_survey` — 顾客调研与评价

- **层级** fct · **来源系统** 小票调研 + App + 平台评价 · **刷新** 日批
- **粒度** 1 行 = 1 份反馈 · **主键** survey_id · **日增行数** 10 ~ 60
- **为什么需要它**：唯一一个不来自门店自己系统的分数。它和内部指标背离时，通常是内部口径出了问题

| 字段 | 类型 | 说明 |
|---|---|---|
| `survey_id` | string | 编号 |
| `source` | string | receipt/app/meituan/eleme |
| `channel` | string | 对应消费渠道 |
| `score` | int | 1-5 |
| `nps` | int | 0-10（若采集） |
| `comment_tag` | string | 结构化标签：速度/口味/态度/干净/漏餐 |

## 现金与日结

### `fct_cash_declaration` — 现金交接

- **层级** fct · **来源系统** POS 现金管理 · **刷新** 每班
- **粒度** 1 行 = 1 次钱箱交接 · **主键** declaration_id · **日增行数** 3 ~ 8
- **为什么需要它**：现金短溢是内控红线，金额小但性质重

| 字段 | 类型 | 说明 |
|---|---|---|
| `declaration_id` | string | 交接号 |
| `shift / drawer_id` | string | 班次与钱箱 |
| `declared_amount` | decimal | 实盘 |
| `system_amount` | decimal | 系统应有 |
| `variance` | decimal | 短溢 |
| `cashier_id` | string | 责任人 |

### `fct_daily_sales_summary` — 日结汇总 (DSR)

- **层级** agg · **来源系统** 由上述 fct 表汇总 · **刷新** 日结后生成
- **粒度** 1 行 = 1 店 1 天 · **主键** store_id, biz_date · **日增行数** 1
- **为什么需要它**：对账口径的定稿：财务、总部、门店三方看同一行数

| 字段 | 类型 | 说明 |
|---|---|---|
| `store_id / biz_date` | - | 主键 |
| `net_sales / guest_count / avg_check` | decimal | 生意三件套 |
| `food_cost / labor_cost` | decimal | 两大可控成本 |
| `cash_variance` | decimal | 现金短溢 |
| `closed_by` | string | 日结人 |

## 聚合

### `agg_metric_interval` — 指标 × 15分钟时段

- **层级** agg · **来源系统** 本 package 的指标引擎 · **刷新** 15 分钟
- **粒度** 1 行 = 门店 × 营业日 × 15分钟 × 指标 · **主键** store_id, biz_date, interval15, metric_id · **日增行数** 约 2500
- **为什么需要它**：时段矩阵（看板上那张热力图）的物理落地。日总数正常但某两个时段崩了——只有这张表看得见

| 字段 | 类型 | 说明 |
|---|---|---|
| `metric_id` | string | 指标键，对应 config/metrics.json |
| `interval15` | string | 时段 |
| `value` | decimal | 指标值 |
| `target` | decimal | 该粒度目标 |
| `status` | string | green/amber/red |

### `agg_metric_dim` — 指标 × 下钻维度

- **层级** agg · **来源系统** 本 package 的指标引擎 · **刷新** 日批 + 日内滚动
- **粒度** 1 行 = 门店 × 营业日 × 指标 × 维度 × 维度值 · **主键** store_id, biz_date, metric_id, dim_name, dim_value · **日增行数** 约 800
- **为什么需要它**：归因矩阵（指标 × 渠道 / 品类 / 岗位 / 人）的物理落地

| 字段 | 类型 | 说明 |
|---|---|---|
| `metric_id` | string | 指标键 |
| `dim_name` | string | channel/daypart/category/station/cashier/platform |
| `dim_value` | string | 维度取值 |
| `value / target / status` | - | 同上 |
| `contribution` | decimal | 该维度对门店总偏差的贡献额 ← 排序用这个，不是用绝对值 |
