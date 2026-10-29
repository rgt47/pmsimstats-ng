#!/bin/bash
# Run the paper 05 sweeps on the simulation host.
# The host runs the committed tree and has no nof1power 0.1.0 installed,
# so load the working-tree engine through pkgload (see 00-common.R).
export NOF1POWER_DEV=1
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript analysis/scripts/nof1-design-sensitivity/run-all.R || exit $?
Rscript analysis/scripts/nof1-design-sensitivity/16-summarize-all.R
