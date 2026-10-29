#!/usr/bin/env bash
# Run one cell of 01-run-cr2-null.R with single-threaded BLAS.
# Usage: bash 01-run.sh SETTING DESIGN N REPS [CORES] [SMOKE] [VARIANTS]
set -euo pipefail
cd "$(dirname "$0")/../../.."
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export CHK_SETTING="$1" CHK_DESIGN="$2" CHK_N="$3" CHK_REPS="$4"
export CHK_CORES="${5:-8}" CHK_SMOKE="${6:-0}" CHK_VARIANTS="${7:-}"
export CHK_ISO="${8:-all}" CHK_REP0="${9:-0}"
Rscript analysis/scripts/cr2-olbdc-check/01-run-cr2-null.R
