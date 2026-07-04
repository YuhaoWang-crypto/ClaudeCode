cd /home/user/ClaudeCode/analysis/mdqm
for i in $(seq 1 90); do
  got=0
  modal volume get --force is621-md-vol WT.json ./WT.json 2>/dev/null && got=$((got+1))
  modal volume get --force is621-md-vol A248K.json ./A248K.json 2>/dev/null && got=$((got+1))
  if [ "$got" -eq 2 ]; then echo "MD_DONE both variants retrieved at $(date +%T)"; break; fi
  sleep 45
done
echo "poller exit; files:"; ls -la WT.json A248K.json 2>/dev/null
