#!/usr/bin/env bash
# Top-up queue: brings the OL+BDC cells to 5000 datasets each with new
# seeds (replicate offset), after 03-queue.sh. N = 70 cells fit only
# the primary analysis (ta0 for p09, g6 for p02).
set -uo pipefail
here="$(dirname "$0")"
log=analysis/data/cr2-olbdc-check/queue.log
cd "$here/../../.."
run() {
  echo "[$(date '+%H:%M:%S')] start $*" >> "$log"
  bash "$here/01-run.sh" "$@" >> "$log" 2>&1
  echo "[$(date '+%H:%M:%S')] end $*" >> "$log"
}
run p09 OLBDC 35 3500 8 0 ta0 all 1500
run p02 OLBDC 35 3500 8 0 g6 all 1500
run p02 OLBDC 70 2000 8 0 g6 all 3000
run p09 OLBDC 70 2000 8 0 ta0 all 3000
