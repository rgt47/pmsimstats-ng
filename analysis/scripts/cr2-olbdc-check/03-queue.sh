#!/usr/bin/env bash
# Production queue, in priority order. Each cell writes its own rds.
set -uo pipefail
here="$(dirname "$0")"
log=analysis/data/cr2-olbdc-check/queue.log
cd "$here/../../.."
run() {
  echo "[$(date '+%H:%M:%S')] start $*" >> "$log"
  bash "$here/01-run.sh" "$@" >> "$log" 2>&1
  echo "[$(date '+%H:%M:%S')] end $*" >> "$log"
}
run p09 OLBDC 70 3000 8 0 ta0,g6form
run p02 OLBDC 70 3000 8 0 g6
run p09 OLBDC 35 1500 8 0 ta0
run p02 OLBDC 35 1500 8 0 g6
run p09 CO 70 1000 8 0 ta0
run p09 Hybrid 70 1000 8 0 ta0
