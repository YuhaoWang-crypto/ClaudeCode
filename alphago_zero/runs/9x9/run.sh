#!/bin/bash
# 9x9 run. Re-run this script after an interruption: it resumes automatically.
cd "$(dirname "$0")/../../.."
R=""; [ -f alphago_zero/runs/9x9/latest.pt ] && R="--resume"
exec python3 -m alphago_zero.train --board-size 9 --komi 7.5 --blocks 5 --filters 48 \
  --value-hidden 64 --sims 100 --eval-sims 100 --temp-moves 12 --dir-alpha 0.15 \
  --window-games 3000 --time-budget-min 240 --out-dir alphago_zero/runs/9x9 $R
