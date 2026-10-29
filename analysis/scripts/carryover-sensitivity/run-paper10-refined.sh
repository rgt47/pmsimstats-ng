#!/bin/bash
# analysis/scripts/carryover-sensitivity/run-paper10-refined.sh
#
# Reruns every paper 10 driver under the refined DGP
# (paper10-common.R: c.cfct = 0.2, carryover decaying from the last
# on-drug mean, stop_on_nonpd = TRUE) and writes the replicate-level
# outputs to output/paper10-refined/. The manuscript
# analysis/report/10-interaction-test-calibration/report.Rmd switches
# each table from the pre-correction run to these outputs once the
# corresponding file holds the full replicate count.
#
# Estimated cost (2026-10-29, from CPU time per fit on an idle 8-core
# Apple Silicon host with 7 workers): about 10 CPU-hours in total, or
# 1.5 to 2 hours of wall time; Kenward-Roger with the unstructured
# covariance is the largest single item. On a shared host the wall
# time scales with the load. The drivers size their worker pool with
# parallelly::availableCores(), so they also run on the simulation host
# (simrun collects every file the run creates).
#
# Usage, from the repository root:
#   bash analysis/scripts/carryover-sensitivity/run-paper10-refined.sh
# A single driver can be run alone with the same thread settings.

set -euo pipefail
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1

D=analysis/scripts/carryover-sensitivity
mkdir -p "$D/output/paper10-refined"
LOG="$D/output/paper10-refined/run.log"

run() {
  echo "=== $* ($(date '+%Y-%m-%d %H:%M'))" | tee -a "$LOG"
  Rscript "$@" 2>&1 | tee -a "$LOG"
}

run $D/diagnose-kappa-analytic.R
run $D/check-null-architecture-equivalence.R
run $D/diagnose-null-type1.R --reps 5000
run $D/diagnose-null-type1-stages46.R --stage 4 --reps 5000
run $D/diagnose-null-type1-cr2.R --reps 3000
run $D/diagnose-alt-cell-cr2.R --reps 3000
run $D/diagnose-null-type1-gee.R --reps 5000
run $D/diagnose-null-type1-stages46.R --stage 6 --reps 1500
run $D/diagnose-null-type1-kr.R --reps 2000
echo "=== done ($(date '+%Y-%m-%d %H:%M'))" | tee -a "$LOG"
