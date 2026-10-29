#!/bin/bash
## Runs the paper 11 driver with single-threaded BLAS, so the parallel
## workers do not oversubscribe the cores. Arguments are passed on to
## 01-run-combined-factorial.R (see its header).
set -euo pipefail
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "$repo_root"
Rscript analysis/scripts/dgp-combined/01-run-combined-factorial.R "$@"
