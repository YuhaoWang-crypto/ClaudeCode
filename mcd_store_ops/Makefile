PY := PYTHONPATH=src python3 -m mcdops.cli

.PHONY: all brief report dashboard reports docs validate dump clean

all: validate report dashboard brief   ## 校验 → 日报 → 看板 → 简报

brief:      ## 班前会简报（终端）
	@$(PY) brief

report:     ## output/daily_report.md
	@$(PY) report

dashboard:  ## output/dashboard.html
	@$(PY) dashboard

reports:    ## 打印每日/每周报表清单
	@$(PY) reports

docs:       ## 从代码重新生成 docs/01 与 docs/02
	@$(PY) docs

validate:   ## 指标字典 ↔ 数据字典 ↔ 报表 三者自洽性检查
	@$(PY) validate

dump:       ## 导出模拟数据为 CSV
	@$(PY) dump

clean:
	@rm -rf output/data output/*.md output/*.html
	@find src -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
