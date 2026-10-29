#!/bin/bash
# Run the paper 05 sweeps on the simulation host.
# The container's renv library holds an old nof1power (0.0.0.9000), and
# furrr workers are separate R processes that load the installed copy,
# so pkgload::load_all in the main session is not enough. Install the
# working-tree engine over it from an Rscript session, where renv is
# active and the engine's dependencies resolve. The job runs on a fresh
# unpacked copy of the project, so the install does not persist.
set -e
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
Rscript -e "install.packages('implementations/nof1power', repos = NULL, type = 'source')"
Rscript -e "stopifnot(utils::packageVersion('nof1power') >= '0.1.0'); print(.libPaths())"
Rscript analysis/scripts/nof1-design-sensitivity/run-all.R
Rscript analysis/scripts/nof1-design-sensitivity/16-summarize-all.R
