## analysis/scripts/dgp-combined/01-run-combined-factorial.R
##
## Factorial simulation driver for Architecture C (combined DGP).
## Extends the Section 3 grid from paper 01 with a 3x3 c.bm_a x c.bm_b
## parameter panel, holding all other design choices identical to the
## existing Architecture A and B comparisons.
##
## Usage:
##   Rscript analysis/scripts/dgp-combined/01-run-combined-factorial.R
##           [--smoke] [--dev] [--reps N]
##
## --smoke : 2 reps, one cell (pipeline sanity only)
## --dev   : 50 reps, trimmed grid (carryover_t1half {0, 1}, CO only,
##           boundary cells only: a=0.45/b=0 and a=0/b=0.45 plus full)
## default : 1000 reps, full grid
##
## Validation gate (embedded): the (c_bm_a=0.45, c_bm_b=0) and
## (c_bm_a=0, c_bm_b=0.45) cells are compared against the existing
## Architecture A and B baselines from analysis/data/. A warning
## is emitted if mean power deviates by more than 3 pp.
##
## Outputs:
##   analysis/data/dgp-combined/01-combined-replicates.rds
##   analysis/data/dgp-combined/01-combined-summary.rds

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(furrr)
})

repo_root <- here::here()
source(file.path(repo_root,
  'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

## Load the prazosin/PTSD reference calibration from the package data.
## These match the Section 3.1 production run exactly so that the
## Architecture C boundary cells (c_bm_a=0.45, c_bm_b=0) and
## (c_bm_a=0, c_bm_b=0.45) reproduce the Architecture A and B
## results from that run.
local({
  e <- new.env()
  load(file.path(repo_root, 'data/extracted_rp.rda'), envir = e)
  load(file.path(repo_root, 'data/extracted_bp.rda'), envir = e)
  assign('pkg_rp', e$extracted_rp, envir = .GlobalEnv)
  assign('pkg_bp', e$extracted_bp, envir = .GlobalEnv)
})
cat('Reference calibration loaded.\n')
cat('  resp_param:\n'); print(pkg_rp)
cat('  baseline_param:\n'); print(pkg_bp)

args        <- commandArgs(trailingOnly = TRUE)
smoke_mode  <- '--smoke' %in% args
dev_mode    <- '--dev'   %in% args
reps_idx    <- which(args == '--reps')
n_reps_override <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else NA_integer_
n_reps <- if (!is.na(n_reps_override)) n_reps_override else
          if (smoke_mode) 2L else if (dev_mode) 50L else 1000L

seed <- 20260610L
set.seed(seed)

## -----------------------------------------------------------------
## Parameter grid
## -----------------------------------------------------------------

## The 3x3 c.bm_a x c.bm_b panel. Boundary row/column cells
## (one parameter = 0) reproduce pure Architecture A and B,
## providing the internal validation gate.
cbm_levels <- c(0, 0.22, 0.45)

if (smoke_mode) {
  grid_combined <- tidyr::expand_grid(
    carryover_form = 'exponential',
    weibull_shape  = 1.0,
    dgp_arch       = 'combined',
    t1half         = 0.5,
    design         = 'CO',
    N              = 35L,
    c_bm_a         = 0.45,
    c_bm_b         = 0.45
  )
} else if (dev_mode) {
  grid_combined <- tidyr::expand_grid(
    carryover_form = 'exponential',
    weibull_shape  = 1.0,
    dgp_arch       = 'combined',
    t1half         = c(0.0, 1.0),
    design         = 'CO',
    N              = 35L,
    c_bm_a         = cbm_levels,
    c_bm_b         = cbm_levels
  )
} else {
  grid_combined <- tidyr::expand_grid(
    carryover_form = 'exponential',
    weibull_shape  = 1.0,
    dgp_arch       = 'combined',
    t1half         = c(0.0, 0.5, 1.0),
    design         = c('CO', 'Hybrid', 'OLBDC'),
    N              = c(35L, 70L),
    c_bm_a         = cbm_levels,
    c_bm_b         = cbm_levels
  )
}

stopifnot(nrow(grid_combined) > 0)
cat(sprintf('Grid cells: %d, n_reps per cell: %d\n',
            nrow(grid_combined), n_reps))

plan(
  if (smoke_mode) sequential else multicore,
  workers = if (smoke_mode) 1L
            else max(1L, parallel::detectCores() - 1L)
)

t_start <- Sys.time()

results <- future_map_dfr(
  seq_len(nrow(grid_combined)),
  function(i) {
    cell        <- grid_combined[i, ]
    cell_result <- simulate_cell(cell, n_reps, robust = FALSE,
                                 resp_param     = pkg_rp,
                                 baseline_param = pkg_bp)
    bind_cols(cell[rep(1L, nrow(cell_result)), ], cell_result)
  },
  .options = furrr_options(seed = seed),
  .progress = !smoke_mode
)

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('Simulation complete. Elapsed: %.1f s.\n', elapsed))

## -----------------------------------------------------------------
## Save replicates
## -----------------------------------------------------------------

out_dir <- file.path(repo_root, 'analysis/data/dgp-combined')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(results,
  file.path(out_dir, '01-combined-replicates.rds'))
cat('Replicates saved to', file.path(out_dir,
  '01-combined-replicates.rds'), '\n')

## -----------------------------------------------------------------
## Summary
## -----------------------------------------------------------------

summary_dt <- results |>
  group_by(dgp_arch, t1half, design, N, c_bm_a, c_bm_b,
           carryover_form, spec) |>
  summarise(
    n_reps     = n(),
    n_conv     = sum(converged, na.rm = TRUE),
    power      = mean(p_value[converged] < 0.05, na.rm = TRUE),
    mcse_power = sqrt(power * (1 - power) / max(n_conv, 1L)),
    mean_beta  = mean(estimate[converged], na.rm = TRUE),
    sd_beta    = sd(estimate[converged], na.rm = TRUE),
    conv_rate  = n_conv / n(),
    .groups = 'drop'
  )

saveRDS(summary_dt,
  file.path(out_dir, '01-combined-summary.rds'))
cat('Summary saved to', file.path(out_dir,
  '01-combined-summary.rds'), '\n')

## -----------------------------------------------------------------
## Print summary table
## -----------------------------------------------------------------

cat('\nPower summary (A2 spec, N=35, primary carryover cells):\n')
summary_dt |>
  dplyr::filter(spec == 'A2', N == 35L,
                design == 'CO',
                t1half %in% c(0.0, 1.0)) |>
  dplyr::select(t1half, c_bm_a, c_bm_b, power, mcse_power,
                conv_rate) |>
  dplyr::arrange(t1half, c_bm_b, c_bm_a) |>
  print(n = Inf)

invisible(NULL)
