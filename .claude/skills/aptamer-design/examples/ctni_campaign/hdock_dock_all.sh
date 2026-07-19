#!/bin/bash
# HDOCK needs SHORT file paths (fixed Fortran filename buffer). Run from /tmp/hd.
set -u
BASE=/tmp/claude-0/-home-user-ClaudeCode/38836661-4715-511c-bc80-0c9cc9a82492/scratchpad
HD=$BASE/tselex_env/hdock_raw/HDOCKlite-v1.1
export LD_LIBRARY_PATH=$BASE/tselex_env/libs
SRC=$BASE/tni/dock
RES=$SRC/hdock_results2.txt
: > "$RES"
for tgt in ctni tnni2 tnni1 lyso; do
  d=/tmp/hd_$tgt
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  cp "$SRC/tni1_rna.pdb" lig.pdb
  cp "$SRC/${tgt}_prot.pdb" rec.pdb
  echo "=== $tgt start $(date +%H:%M:%S) ===" >> "$RES"
  "$HD/hdock" rec.pdb lig.pdb -out out.dat > h.log 2>&1
  "$HD/createpl" out.dat top.pdb -nmax 10 -complex -models > c.log 2>&1
  best=$(grep -h "REMARK Score" model_*.pdb 2>/dev/null | awk '{print $3}' | sort -g | head -1)
  echo "$tgt SCORE $best" >> "$RES"
  echo "=== $tgt done $(date +%H:%M:%S) ===" >> "$RES"
done
echo "ALL_DONE" >> "$RES"
