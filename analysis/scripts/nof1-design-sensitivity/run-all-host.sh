#!/bin/bash
# Run the paper 05 sweeps on the simulation host.
# The host runs the committed tree and has no nof1power 0.1.0 installed.
# furrr workers are separate R processes that load the installed package,
# so pkgload::load_all in the main session is not enough: install the
# working-tree engine into a private library that every worker inherits
# through R_LIBS.
set -e
LIB="$(pwd)/.nof1power-lib"
mkdir -p "$LIB"
R CMD INSTALL --no-test-load --library="$LIB" implementations/nof1power
export R_LIBS="$LIB${R_LIBS:+:$R_LIBS}"
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript -e "stopifnot(utils::packageVersion('nof1power') >= '0.1.0')"
Rscript analysis/scripts/nof1-design-sensitivity/run-all.R
Rscript analysis/scripts/nof1-design-sensitivity/16-summarize-all.R
