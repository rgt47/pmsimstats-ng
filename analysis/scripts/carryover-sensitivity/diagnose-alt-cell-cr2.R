## analysis/scripts/carryover-sensitivity/diagnose-alt-cell-cr2.R
##
## Paper 10 (analysis/report/10-interaction-test-calibration/),
## Section "The remedy at an alternative cell". Runs the CR2 driver
## (diagnose-null-type1-cr2.R) at the alternative c_bm = 0.45 of the
## covariance-moderation architecture: Hybrid design, N = 70,
## exponential decay, t1half = 1.0 week, under the refined DGP of
## paper10-common.R. Each replicate is analyzed by M3 for E1, E2, E3
## (model-based and CR2) and by G9 and F1, so power differences
## between procedures are paired.
##
## Corrected 2026-10-29 (referee report 2026-10-28): the 2026-05-20
## version hard-coded c.cfct = 0.1 (docs/61) and predates carryover
## that decays from the last on-drug mean.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-alt-cell-cr2.R [--reps N]

repo_root <- here::here()
args <- c('--cbm', '0.45', '--seed', '20260521',
          '--out', 'diag-alt-cell-cr2.rds',
          commandArgs(trailingOnly = TRUE))
status <- system2('Rscript',
  c(file.path(repo_root,
    'analysis/scripts/carryover-sensitivity/diagnose-null-type1-cr2.R'),
    args))
if (!identical(status, 0L)) stop('diagnose-null-type1-cr2.R failed')
