#!/usr/bin/env bash
# Run the paper 09 paired dropout driver with single-threaded BLAS.
# Usage (from the repository root):
#   bash analysis/scripts/informative-dropout-by-design/20-run-paired.sh
# Optional environment: P09_REPS, P09_NULL_REPS, P09_CORES, P09_OUT,
# P09_TAG (see 20-paired-dropout.R).
set -euo pipefail
export VECLIB_MAXIMUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript analysis/scripts/informative-dropout-by-design/20-paired-dropout.R
