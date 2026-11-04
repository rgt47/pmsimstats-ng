#!/usr/bin/env bash
# Top-up queue: brings the OL+BDC cells to 5000 datasets each with new
# seeds (replicate offset), after 03-queue.sh. N = 70 cells fit only
# the primary analysis (ta0 for p09, g6 for p02).
# Usage (from the repository root), with the same PMSIM_PARAMS (and
# CHK_OUT) as the 03-queue.sh run it tops up:
#   [PMSIM_PARAMS=NAME] bash analysis/scripts/cr2-olbdc-check/06-topup.sh
set -uo pipefail
here="$(dirname "$0")"
cd "$here/../../.."
if [ -n "${PMSIM_PARAMS:-}" ]; then
  export PMSIM_PARAMS
  export CHK_OUT="${CHK_OUT:-analysis/data/cr2-olbdc-check-$PMSIM_PARAMS}"
  # See 03-queue.sh: mark (or refuse) the folder before the log.
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
run p09 OLBDC 35 3500 8 0 ta0 all 1500
run p02 OLBDC 35 3500 8 0 g6 all 1500
run p02 OLBDC 70 2000 8 0 g6 all 3000
run p09 OLBDC 70 2000 8 0 ta0 all 3000
