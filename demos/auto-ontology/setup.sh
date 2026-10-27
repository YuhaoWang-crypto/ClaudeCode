#!/usr/bin/env bash
# 在没有 Docker 的 Linux（Ubuntu 24.04 + PostgreSQL 16）上安装 NVIDIA Auto Ontology
# 并灌入演示数据。需要 root（或 sudo）、uv、Node 22 + pnpm。
#
#   bash setup.sh            # 安装 + 启动后端(:3001) 和前端(:3000)
#   python3 demo.py          # 跑命令行演示
#
# 可选：export DEFAULT_MODELS_API_KEY=nvapi-...  打开自然语言问答（见 README）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="${SRC_DIR:-$HOME/auto-ontology}"
PGPASS="${PGPASS:-postgres}"
ADMIN_EMAIL="${ADMIN_EMAIL:-admin@example.com}"
ADMIN_PASSWORD="${ADMIN_PASSWORD:-demo-password-123}"
SUDO=""; [ "$(id -u)" -ne 0 ] && SUDO="sudo"
as_pg() { if [ "$(id -u)" -eq 0 ]; then runuser -u postgres -- "$@"; else sudo -u postgres "$@"; fi; }

echo ">>> 1/7 获取源码 -> $SRC_DIR"
[ -d "$SRC_DIR/.git" ] || git clone --depth 1 https://github.com/NVIDIA/auto-ontology "$SRC_DIR"
cd "$SRC_DIR"

echo ">>> 2/7 PostgreSQL 16 + pgvector 0.8.1（apt 自带的 0.6 缺少 halfvec 类型）"
$SUDO apt-get install -y postgresql-16 postgresql-server-dev-16 build-essential >/dev/null
$SUDO service postgresql start
if ! as_pg psql -tAc "select 1 from pg_available_extension_versions where name='vector' and version='0.8.1'" 2>/dev/null | grep -q 1; then
  tmp=$(mktemp -d)
  git clone --depth 1 --branch v0.8.1 https://github.com/pgvector/pgvector.git "$tmp/pgvector"
  make -C "$tmp/pgvector" -j4 >/dev/null && $SUDO make -C "$tmp/pgvector" install >/dev/null
fi
$SUDO service postgresql restart  # 重启以加载新版 pgvector
as_pg psql -qc "ALTER USER postgres PASSWORD '$PGPASS';"
as_pg psql -tAc "select 1 from pg_database where datname='auto_ontology'" | grep -q 1 \
  || as_pg psql -qc "CREATE DATABASE auto_ontology;"

echo ">>> 3/7 写 .env"
CHINOOK="$SRC_DIR/dev_tools/fixtures/sql/chinook.sqlite"
cat > .env <<EOF
POSTGRES_HOST=localhost
POSTGRES_PORT=5432
POSTGRES_USER=postgres
POSTGRES_PASSWORD=$PGPASS
POSTGRES_DATABASE=auto_ontology
CONNECTION_STRINGS=postgresql://postgres:$PGPASS@localhost:5432/pagila,sqlite:///$CHINOOK
APP_URL=http://localhost:3000
INGESTION_SERVICE_URL=http://localhost:3002
AUTO_ONTOLOGY_ADMIN_EMAIL=$ADMIN_EMAIL
AUTO_ONTOLOGY_ADMIN_PASSWORD=$ADMIN_PASSWORD
AUTH_SECRET=$(openssl rand -hex 24)
EOF
if [ -n "${DEFAULT_MODELS_API_KEY:-}" ]; then
  echo "DEFAULT_MODELS_API_KEY=$DEFAULT_MODELS_API_KEY" >> .env
  echo "DEFAULT_MODELS_ENDPOINT=https://integrate.api.nvidia.com/v1" >> .env
fi
set -a; . ./.env; set +a

echo ">>> 4/7 Python 依赖 (uv sync)"
uv sync

echo ">>> 5/7 示例库 pagila + chinook，建表结构，生成目录 + 手写语义层"
# pagila.sql 由 PG17 导出，PG16 不认识 transaction_timeout
sed -i 's/^SET transaction_timeout = 0;/-- &/' dev_tools/fixtures/sql/pagila.sql
uv run --no-sync python -m dev_tools.fixtures.seed_fixtures
uv run --no-sync alembic upgrade head
uv run --no-sync python -m dev_tools.fixtures.seed_graph_fixture --reset

echo ">>> 6/7 启动后端 :3001"
nohup uv run --no-sync python -m auto_ontology.server > backend.log 2>&1 &
for _ in $(seq 1 60); do curl -sf localhost:3001/api/health >/dev/null && break; sleep 1; done

echo ">>> 7/7 前端依赖 + 启动前端 :3000"
( cd frontend && pnpm install --frozen-lockfile && sh prisma/migrate.sh && pnpm exec prisma generate )
( cd frontend && PYTHON_API_URL=http://127.0.0.1:3001 nohup pnpm exec next dev -p 3000 > ../frontend.log 2>&1 & )
for _ in $(seq 1 90); do curl -s -o /dev/null localhost:3000 && break; sleep 2; done

echo
echo "完成。UI: http://localhost:3000  （登录 $ADMIN_EMAIL / $ADMIN_PASSWORD）"
echo "      API: http://localhost:3001/docs"
echo "接着运行: python3 $SCRIPT_DIR/demo.py"
