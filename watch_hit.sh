#!/usr/bin/env bash
# Watch a run dir for a verified hit; on hit, record it and stop all workers.
# Usage: watch_hit.sh <run_dir>
set -u
run_dir="$1"
while true; do
  if grep -H "Verified: PASS" "$run_dir"/worker*.log > /tmp/p24_hit_lines.txt 2>/dev/null; then
    {
      date
      cat /tmp/p24_hit_lines.txt
      for log in "$run_dir"/worker*.log; do
        if grep -q "Verified: PASS" "$log"; then
          echo "--- $log tail ---"
          tail -n 8 "$log"
        fi
      done
    } > "$run_dir/HIT.txt"
    awk '{print $1}' "$run_dir/pids.txt" | xargs kill 2>/dev/null
    exit 0
  fi
  # stop if all workers exited without a hit
  alive=0
  while read -r pid _; do
    kill -0 "$pid" 2>/dev/null && alive=$((alive+1))
  done < "$run_dir/pids.txt"
  if [ "$alive" -eq 0 ]; then
    echo "all workers exited, no hit" > "$run_dir/EXHAUSTED.txt"
    date >> "$run_dir/EXHAUSTED.txt"
    exit 1
  fi
  sleep 30
done
