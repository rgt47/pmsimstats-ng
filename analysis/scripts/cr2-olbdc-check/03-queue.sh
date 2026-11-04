#!/usr/bin/env bash
# Production queue, in priority order. Each cell writes its own rds.
# Usage (from the repository root):
#   [PMSIM_PARAMS=NAME] bash analysis/scripts/cr2-olbdc-check/03-queue.sh
# With PMSIM_PARAMS set, every cell runs under that set into
# analysis/data/cr2-olbdc-check-NAME (or CHK_OUT); without it the
# queue runs in legacy mode into analysis/data/cr2-olbdc-check.
set -uo pipefail
here="$(dirname "$0")"
cd "$here/../../.."
if [ -n "${PMSIM_PARAMS:-}" ]; then
  export PMSIM_PARAMS
  export CHK_OUT="${CHK_OUT:-analysis/data/cr2-olbdc-check-$PMSIM_PARAMS}"
  # The log is written before the first cell, and guard_out_dir()
  # reads an unmarked folder holding files as 'simcore2026', so the
  # folder is marked (or refused) first.
  Rscript -e "source('analysis/scripts/carryover-sensitivity/simulation-core.R'); guard_out_dir('$CHK_OUT')" || exit 1
fi
out="${CHK_OUT:-analysis/data/cr2-olbdc-check}"
mkdir -p "$out"
log="$out/queue.log"
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
