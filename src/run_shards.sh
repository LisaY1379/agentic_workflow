#!/usr/bin/env bash
# Launch N detached single-thread pomerance workers with distinct PRNG seeds.
# Usage: run_shards.sh <p> <mode> <workers> <trials_per_worker> <seed_base> <seed_step> <run_dir>
set -euo pipefail
p="$1"; mode="$2"; workers="$3"; trials="$4"; seed_base="$5"; seed_step="$6"; run_dir="$7"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bin="$here/pomerance"
mkdir -p "$run_dir"
: > "$run_dir/pids.txt"
for ((i=0; i<workers; i++)); do
  seed=$((seed_base + i * seed_step))
  log=$(printf "%s/worker%02d.log" "$run_dir" "$i")
  nohup "$bin" "$p" "$seed" "$trials" "$mode" > "$log" 2>&1 &
  pid=$!
  disown "$pid" 2>/dev/null || true
  echo "$pid $log $seed" >> "$run_dir/pids.txt"
done
echo "run_dir=$run_dir"
echo "launched=$workers workers for p=$p mode=$mode trials/worker=$trials"
cat "$run_dir/pids.txt"
