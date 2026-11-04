#!/usr/bin/env bash
# Run one cell of 01-run-cr2-null.R with single-threaded BLAS.
# Usage: [PMSIM_PARAMS=NAME] [CHK_OUT=DIR] \
#   bash 01-run.sh SETTING DESIGN N REPS [CORES] [SMOKE] [VARIANTS] \
#     [ISO] [REP0]
# PMSIM_PARAMS (published | simcore2026 | fixture) is passed on as
# --params; the output folder is then analysis/data/cr2-olbdc-check-NAME
# unless CHK_OUT is set. Without it the cell runs in legacy mode (see
# 01-run-cr2-null.R).
set -euo pipefail
cd "$(dirname "$0")/../../.."
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export CHK_SETTING="$1" CHK_DESIGN="$2" CHK_N="$3" CHK_REPS="$4"
export CHK_CORES="${5:-8}" CHK_SMOKE="${6:-0}" CHK_VARIANTS="${7:-}"
export CHK_ISO="${8:-all}" CHK_REP0="${9:-0}"
params=()
if [ -n "${PMSIM_PARAMS:-}" ]; then
  params=(--params "$PMSIM_PARAMS")
fi
Rscript analysis/scripts/cr2-olbdc-check/01-run-cr2-null.R ${params[@]+"${params[@]}"}
