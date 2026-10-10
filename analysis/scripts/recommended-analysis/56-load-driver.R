## Load the data-generating process and the analyses of the stress-test
## driver (11-strawman-stress-test.R) without running it, so that the
## docs/56 examples use exactly the driver's code. Only top-level
## assignments of the names below, and the loop that reads the
## published construct from 04-power-simulation.R, are evaluated.
##
## carryover, allocation: the driver's --carryover and --allocation.
## 'recursive' and 'fixed' reproduce the construct of the stored
## stress-test results (run 2026-10-07); the driver now defaults to
## 'anchored' and 'random'.
load_driver <- function(carryover = 'recursive', allocation = 'fixed') {
  suppressMessages({
    library(data.table); library(MASS); library(nlme)
  })
  drv <- 'analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R'
  stopifnot(file.exists(drv))
  env <- new.env()
  env$CARRYOVER <- carryover
  env$ALLOCATION <- allocation
  want <- c('here04', 'keep', 'dgp', 'D', 'PIS', 'TV_MAX', 'base', 'cell',
            'decay', 'build', 'draw', 'dropout_visit', 'TH_A', 'GRID', 'dbc',
            'na3', 'lctrl', 'fit_lme', 'cr2', 'mb', 'F_S', 'F_F1', 'F_F2',
            'F_E', 'with_x', 'make_trial', 'cell_setup', 'alloc')
  for (e in parse(drv)) {
    is_assign <- is.call(e) && identical(e[[1]], as.name('<-')) &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% want
    is_dgp_loop <- is.call(e) && identical(e[[1]], as.name('for')) &&
      grepl('here04', paste(deparse(e), collapse = ' '))
    if (is_assign || is_dgp_loop) eval(e, env)
  }
  stopifnot(all(want %in% ls(env)), all(env$keep %in% ls(env$dgp)))
  env
}
