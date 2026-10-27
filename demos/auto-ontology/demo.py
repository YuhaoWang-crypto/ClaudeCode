#!/usr/bin/env python3
"""NVIDIA Auto Ontology 离线演示 —— 不需要 NVIDIA API key。

前提：setup.sh 已经跑完，后端在 http://localhost:3001 运行，
示例库 pagila（DVD 租赁）和 chinook（音乐商店）已经摄入。

只用标准库 + psql 命令行，直接运行：
    python3 demo.py
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.error
import urllib.request

API = os.environ.get("AUTO_ONTOLOGY_BACKEND", "http://localhost:3001")
PAGILA = os.environ.get("PAGILA_URL", "postgresql://postgres:postgres@localhost:5432/pagila")
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
HEADERS = {"x-auto-ontology-user-id": "demo", "content-type": "application/json"}


def call(method: str, path: str, body: dict | None = None, raw: bool = False):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(API + path, data=data, method=method, headers=HEADERS)
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            text = resp.read().decode()
    except urllib.error.HTTPError as e:
        text = e.read().decode()
        try:
            return {"_error": e.code, **json.loads(text)}
        except json.JSONDecodeError:
            return {"_error": e.code, "detail": text}
    return text if raw else json.loads(text)


def run_sql(sql: str, limit: int = 5) -> str:
    wrapped = f"SELECT * FROM ({sql.rstrip(';')}) q LIMIT {limit}"
    res = subprocess.run(["psql", PAGILA, "-c", wrapped], capture_output=True, text=True)
    return (res.stdout or res.stderr).rstrip()


def section(n: int, title: str) -> None:
    print(f"\n{'=' * 72}\n{n}. {title}\n{'=' * 72}")


def main() -> int:
    section(0, "健康检查")
    health = call("GET", "/api/health")
    print(json.dumps(health, ensure_ascii=False))
    if health.get("status") != "ok":
        print("后端没起来，先运行 setup.sh")
        return 1

    section(1, "数据目录 (Catalog)：自动从源库抽取的结构")
    dbs = call("GET", "/api/datasources/dbs")["data"]
    for db in dbs:
        schemas = call("GET", f"/api/schemas/{db['id']}")["schemas"]
        parts = [f"{s['schema_name']}({s['tables_count']} 张表/视图)" for s in schemas]
        print(f"  • {db['name']:8s} → " + ", ".join(parts))

    section(2, "本体 (Ontology)：业务术语、同义词、认证状态")
    terms = call("GET", "/api/terms")["terms"]
    term_by_id = {t["id"]: t for t in terms}
    for t in terms:
        syn = ", ".join(t["synonyms"]) or "-"
        print(f"  • {t['name']:9s} 同义词[{syn}]  认证={t['certification']}")
        print(f"      {t['description']}")

    section(3, "语义关系：术语之间如何关联（本演示为手写；正式流程由 LLM 语义编译生成）")
    graph = call("GET", "/api/exploration/semantic-graph")["data"]
    names = {n["id"]: n["name"] for n in graph["nodes"]}
    for link in graph["links"]:
        kinds = ", ".join(link.get("relationship_types", []))
        print(f"  • {names.get(link['source'])} ──▶ {names.get(link['target'])}   ({kinds})")

    section(4, "同义词检索：用户说的词 → 本体里的正式概念")
    for word in ["movie", "song", "renter", "transaction"]:
        hits = call("POST", "/api/search/global-search", {"search_term": word}).get("data", [])
        found = ", ".join(f"{h['type']}:{h['name']}" for h in hits[:3]) or "(无)"
        print(f"  • “{word}” → {found}")

    section(5, "受治理的指标 (SQL Attributes)：挂在术语上的标准计算口径，直接在源库执行")
    attrs = call("GET", "/api/sql-attributes")["data"]
    for a in attrs:
        term = term_by_id.get(a["term_id"], {}).get("name", "?")
        print(f"\n  ▸ [{term}] {a['name']} — {a['description']}  (certified={a['certified']})")
        print(f"    SQL: {a['sql']}")
        if term in ("Track",):
            print("    (属于 chinook/SQLite 库，此处跳过执行)")
            continue
        print("    " + run_sql(a["sql"], 3).replace("\n", "\n    "))

    section(6, "护栏：校验 Agent 生成的 SQL 是否基于目录里真实存在的表")
    film_id = next(t["id"] for t in terms if t["name"] == "Film")
    cases = [
        ("合法：按影片统计收入", {"expression": (
            "SELECT f.title, sum(p.amount) FROM film f "
            "JOIN inventory i ON i.film_id = f.film_id "
            "JOIN rental r ON r.inventory_id = i.inventory_id "
            "JOIN payment p ON p.rental_id = r.rental_id GROUP BY f.title")}),
        ("幻觉：引用不存在的表 films", {"expression": "SELECT title FROM films"}),
        ("重复：和已有指标 film count 口径相同",
         {"expression": "SELECT count(film_id) FROM film", "term_id": film_id}),
    ]
    for label, body in cases:
        r = call("POST", "/api/sql-attributes/validate", body)
        verdict = "✅ 通过" if r.get("data", {}).get("valid") else f"❌ 拒绝：{r.get('detail')}"
        print(f"  • {label:28s} {verdict}")

    section(7, "自定义分析 (Custom Analyses)：保存下来的常用分析")
    for a in call("GET", "/api/custom-analyses")["data"]:
        print(f"\n  ▸ {a['name']} — {a['description']}")
        print("    " + run_sql(a["sql"], 5).replace("\n", "\n    "))

    section(8, "导出整个本体模型为 YAML（可版本管理 / 导入到别的实例）")
    yaml_text = call("POST", "/api/model/export", {}, raw=True)
    out = os.path.join(OUT_DIR, "exported_model.yaml")
    with open(out, "w") as f:
        f.write(yaml_text)
    print(f"  已写入 {out}  ({len(yaml_text.splitlines())} 行)")
    top = [ln for ln in yaml_text.splitlines() if ln and not ln.startswith(" ")]
    print("  顶层分区：", ", ".join(top))

    section(9, "自然语言问答（需要 LLM）")
    print("  问题：Which film category has the most films?")
    stream = call("POST", "/api/chat/completions",
                  {"question": "Which film category has the most films?"}, raw=True)
    events = [json.loads(ln[6:]) for ln in str(stream).splitlines()
              if ln.startswith("data: {")]
    for ev in events:
        if ev.get("type") == "step" and ev.get("phase") == "start":
            print(f"    … {ev.get('label')}")
        elif ev.get("type") == "result":
            print("  答案：", json.dumps(ev, ensure_ascii=False)[:800])
        elif ev.get("type") == "error":
            print(f"  ⚠ {ev.get('message')}")
            print("  → 这一步要调用大模型。在 .env 里填入 DEFAULT_MODELS_API_KEY"
                  "（build.nvidia.com 免费申请）后重启后端即可。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
