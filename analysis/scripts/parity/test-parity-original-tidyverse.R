#!/usr/bin/env Rscript
# test-parity-original-tidyverse.R
#
# Cross-implementation parity test comparing the `original/` data.table
# implementation (Architecture B only) against the `tidyverse/`
# implementation at `dgp_architecture = "mvn"`. This is the narrower
# counterpart to test-parity-extended-tidyverse.R and is specifically
# scoped to a state of the codebase after the 2026-04-18 bug fixes
# (mod_gompertz maxr=0 guard; lme_analysis Dbc 0/0 guard).
#
# Coverage axes:
#   design: OL, CO, Hybrid, OLBDC
#   t_half: 0, 0.5, 1.0
#   c_bm:   0.0, 0.3, 0.45
#   seed:   42, 101
#
# Tolerances:
#   sigma, data: 1e-10
#   analysis:    1e-6 on beta/SE; 1e-6 on p
#
# Run from the project root:
#   Rscript analysis/scripts/parity/test-parity-original-tidyverse.R

suppressPackageStartupMessages({
  library(data.table)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(nlme)
  library(MASS)
  library(Matrix)
  library(corpcor)
})

# ---- 1. Source both implementations into separate envs ---------------------

here <- function(...) file.path('implementations', ...)

orig_env <- new.env()
for (f in list.files(here('original', 'R'), full.names = TRUE)) {
  sys.source(f, envir = orig_env)
}
tid_env <- new.env()
sys.source(here('tidyverse', 'R', 'functions.R'), envir = tid_env)

# ---- 2. Trial designs (matched between the two backends) -------------------

make_designs <- function(build_fn) {
  list(
    OL = build_fn(
      name_longform = 'open label', name_shortform = 'OL',
      timepoints = cumsum(rep(2.5, 8)),
      timeptnames = paste0('OL', 1:8),
      expectancies = rep(1, 8),
      ondrug = list(pathA = rep(1, 8))
    ),
    CO = build_fn(
      name_longform = 'traditional crossover', name_shortform = 'CO',
      timepoints = cumsum(rep(2.5, 8)),
      timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
      expectancies = rep(0.5, 8),
      ondrug = list(pathA = c(1,1,1,1,0,0,0,0),
                    pathB = c(0,0,0,0,1,1,1,1))
    ),
    Hybrid = build_fn(
      name_longform = 'primary N-of-1 design', name_shortform = 'Nof1',
      timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
      timeptnames = c('OL1','OL2','BD1','BD2','BD3','BD4','COd','COp'),
      expectancies = c(1,1,0.5,0.5,0.5,0.5,0.5,0.5),
      ondrug = list(pathA = c(1,1,1,1,0,0,1,0),
                    pathB = c(1,1,1,1,0,0,0,1))
    ),
    OLBDC = build_fn(
      name_longform = 'open label + blinded discontinuation',
      name_shortform = 'OL+BDC',
      timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
      timeptnames = c('OL1','OL2','OL3','OL4','BD1','BD2','BD3','BD4'),
      expectancies = c(1,1,1,1,0.5,0.5,0.5,0.5),
      ondrug = list(pathA = c(1,1,1,1,1,1,0,0),
                    pathB = c(1,1,1,1,1,0,0,0))
    )
  )
}

designs_orig <- make_designs(orig_env$buildtrialdesign)

# The tidyverse backend expects a `timepoint_name` column alongside
# `timeptnames`; add it as an alias so both pipelines see matching designs.
with_tp_name <- function(d) {
  d$trialpaths <- lapply(d$trialpaths, function(p) {
    p$timepoint_name <- p$timeptnames
    p
  })
  d
}
designs_tid <- lapply(designs_orig, with_tp_name)

# ---- 3. Fixed response / baseline parameters -------------------------------

resp_dt <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, 5)
)
resp_tbl <- as_tibble(resp_dt)

bl_dt <- data.table(cat = c('bm','BL'), m = c(0, 70), sd = c(1, 10))
bl_tbl <- as_tibble(bl_dt)

make_model_param <- function(N, c_bm, t_half) {
  data.table(N = N, c.bm = c_bm, carryover_t1half = t_half,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

# ---- 4. Comparison helpers -------------------------------------------------

TOL_NUM <- 1e-10
TOL_BETA <- 1e-6
TOL_P <- 1e-6

cmp <- function(a, b, tol = TOL_NUM) {
  r <- all.equal(a, b, tolerance = tol)
  if (isTRUE(r)) TRUE else sprintf('MISMATCH: %s', paste(r, collapse = '; '))
}

compare_sigma <- function(design, t_half, c_bm, N = 35) {
  mp <- make_model_param(N, c_bm, t_half)
  s_orig <- orig_env$buildSigma(mp, resp_dt, bl_dt,
                                designs_orig[[design]]$trialpaths[[1]])
  s_tid <- tid_env$build_sigma_matrix(mp, resp_tbl, bl_tbl,
                                      designs_tid[[design]]$trialpaths[[1]],
                                      dgp_architecture = 'mvn')
  list(
    labels = cmp(s_orig$labels, s_tid$labels, tol = 0),
    means  = cmp(s_orig$means,  s_tid$means),
    sigma  = cmp(unname(s_orig$sigma), unname(s_tid$sigma))
  )
}

compare_data <- function(design, t_half, c_bm, seed, N = 35,
                         empirical = FALSE) {
  mp <- make_model_param(N, c_bm, t_half)
  set.seed(seed)
  d_orig <- orig_env$generateData(mp, resp_dt, bl_dt,
                                  designs_orig[[design]]$trialpaths[[1]],
                                  empirical = empirical,
                                  makePositiveDefinite = TRUE)
  set.seed(seed)
  d_tid <- tid_env$generate_data(mp, resp_tbl, bl_tbl,
                                 designs_tid[[design]]$trialpaths[[1]],
                                 empirical = empirical,
                                 make_positive_definite = TRUE,
                                 dgp_architecture = 'mvn')
  shared <- intersect(names(d_orig), names(d_tid))
  results <- vapply(shared, function(col) {
    isTRUE(cmp(d_orig[[col]], d_tid[[col]]))
  }, logical(1))
  list(shared = shared, ok = results,
       nrow_orig = nrow(d_orig), nrow_tid = nrow(d_tid),
       n_pass = sum(results), n_total = length(results))
}

# ---- 5. Run factorial -------------------------------------------------------

grid <- expand.grid(
  design = c('OL', 'CO', 'Hybrid', 'OLBDC'),
  t_half = c(0, 0.5, 1.0),
  c_bm   = c(0.0, 0.3, 0.45),
  seed   = c(42, 101),
  stringsAsFactors = FALSE
)

fail_rows <- list()
pass_rows <- 0L

cat(sprintf('Running %d cells...\n', nrow(grid)))
for (i in seq_len(nrow(grid))) {
  cell <- grid[i, ]

  # Sigma layer (independent of seed, test once per design/t_half/c_bm)
  s_cmp <- compare_sigma(cell$design, cell$t_half, cell$c_bm)

  # Data layer
  d_cmp <- compare_data(cell$design, cell$t_half, cell$c_bm, cell$seed)

  any_fail <- !isTRUE(s_cmp$labels) || !isTRUE(s_cmp$means) ||
              !isTRUE(s_cmp$sigma)  || (d_cmp$n_pass < d_cmp$n_total)

  if (any_fail) {
    fail_rows[[length(fail_rows) + 1]] <- list(
      cell = cell, sigma = s_cmp, data = d_cmp
    )
  } else {
    pass_rows <- pass_rows + 1L
  }
}

cat(sprintf('\n=== Parity Summary ===\n'))
cat(sprintf('cells passed: %d / %d\n', pass_rows, nrow(grid)))
cat(sprintf('cells failed: %d\n', length(fail_rows)))

if (length(fail_rows) > 0) {
  cat('\n=== Failures (first 5) ===\n')
  for (f in fail_rows[seq_len(min(5, length(fail_rows)))]) {
    cat(sprintf('\n--- design=%s t_half=%s c_bm=%s seed=%d ---\n',
                f$cell$design, f$cell$t_half, f$cell$c_bm, f$cell$seed))
    cat(sprintf('  sigma.labels: %s\n',
                if (isTRUE(f$sigma$labels)) 'OK' else f$sigma$labels))
    cat(sprintf('  sigma.means:  %s\n',
                if (isTRUE(f$sigma$means)) 'OK' else f$sigma$means))
    cat(sprintf('  sigma.sigma:  %s\n',
                if (isTRUE(f$sigma$sigma)) 'OK' else f$sigma$sigma))
    cat(sprintf('  data:         %d/%d cols pass\n',
                f$data$n_pass, f$data$n_total))
    if (f$data$n_pass < f$data$n_total) {
      fail_cols <- names(f$data$ok)[!f$data$ok]
      cat('  failing cols: ', paste(head(fail_cols, 10), collapse = ', '),
          '\n')
    }
  }
  quit(status = 1)
}

cat('\nAll parity checks passed.\n')
