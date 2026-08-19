PY := PYTHONPATH=src python3 -m mcdops.cli

# LibreOffice 重算脚本（随 Claude 的 xlsx skill 分发）；没有它就跳过重算与校验
RECALC ?= /root/.claude/skills/synced/xlsx/scripts/recalc.py

.PHONY: all brief report dashboard excel reports docs validate verify dump clean

all: validate report dashboard excel brief   ## 校验 → 日报 → 看板 → Excel → 简报

brief:      ## 班前会简报（终端）
	@$(PY) brief

report:     ## output/daily_report.md
	@$(PY) report

dashboard:  ## output/dashboard.html
	@$(PY) dashboard

excel:      ## output/store_ops_package.xlsx（11 个 sheet，公式驱动）
	@$(PY) excel
	@if [ -f "$(RECALC)" ]; then \
		HOME=$${HOME:-/root} python3 "$(RECALC)" output/store_ops_package.xlsx 360; \
		$(MAKE) --no-print-directory verify; \
	else \
		echo "  ! 未找到 $(RECALC)：公式已写入但没有缓存值，用 Excel/LibreOffice 打开一次即可"; \
	fi

verify:     ## 核对 Excel 公式结果 == 指标引擎结果
	@PYTHONPATH=src python3 tests/verify_excel.py

reports:    ## 打印每日/每周报表清单
	@$(PY) reports

docs:       ## 从代码重新生成 docs/01 与 docs/02
	@$(PY) docs

validate:   ## 指标字典 ↔ 数据字典 ↔ 报表 三者自洽性检查
	@$(PY) validate

dump:       ## 导出模拟数据为 CSV
	@$(PY) dump

clean:
	@rm -rf output/data output/*.md output/*.html output/*.xlsx
	@find src -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
