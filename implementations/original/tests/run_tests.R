#!/usr/bin/env Rscript
# run_tests.R
#
# Runs the tinytest suite for implementations/original/ by sourcing the
# R/ files directly into the global environment, then invoking
# tinytest::run_test_dir on the tests directory.
#
# Run from the project root:
#   Rscript implementations/original/tests/run_tests.R

suppressPackageStartupMessages({
  library(data.table)
  library(corpcor)
  library(MASS)
  library(nlme)
  library(Matrix)
  library(tinytest)
})

here <- function(...) {
  file.path('implementations', 'original', ...)
}

for (f in list.files(here('R'), full.names = TRUE)) {
  source(f)
}

results <- tinytest::run_test_dir(
  here('tests'),
  ignore.file = 'run_tests.R',
  verbose = 1
)

print(summary(results))

n_fail <- sum(!as.logical(results))
quit(status = if (n_fail > 0) 1 else 0)
