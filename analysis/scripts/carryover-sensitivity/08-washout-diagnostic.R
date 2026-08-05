## analysis/scripts/carryover-sensitivity/08-washout-diagnostic.R
##
## Reference-cell diagnostic for the E4 (Washout-adjusted) analysis
## specification, added to simulation-core.R as the replacement middle
## arm for manuscript 02. This is not a production run: it fits the
## reference cells only, so the swap decision can be checked cheaply
## before committing to a full factorial rerun.
##
## Two questions it must answer:
##   1. Does E4 recover power over E1, and does it do so through the
##      variance channel (same point estimate, smaller standard error)
##      as the reference implementation's output suggests?
##   2. Does E4 hold nominal Type I error? The reference
##      implementation's analogous arm reached 0.053, so a power gain
##      is not creditable until size is verified.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/08-washout-diagnostic.R [--reps N]
##
## Writes output/08-washout-diagnostic.rds

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

args <- commandArgs(trailingOnly = TRUE)
reps_idx <- which(args == '--reps')
n_reps <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else 500

## Seed matches the production factorial so the E1/E2/E3 columns here
## are directly comparable with 01-factorial.rds.
seed <- 20260415L
set.seed(seed)

grid <- tidyr::expand_grid(
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  dgp_arch       = 'mvn',
  t1half         = c(0.5, 1.0),
  design         = c('CO', 'OLBDC', 'Hybrid'),
  N              = 70,
  c_bm           = c(0, 0.45)
)

n_cores <- max(1, future::availableCores() - 2)
plan(multisession, workers = n_cores)
cat(sprintf('Cells: %d | reps/cell: %d | workers: %d\n',
            nrow(grid), n_reps, n_cores))

t_start <- Sys.time()
results <- future_map_dfr(
  seq_len(nrow(grid)),
  function(i) {
    cell <- grid[i, ]
    cell_result <- simulate_cell(cell, n_reps, robust = FALSE)
    bind_cols(cell[rep(1, nrow(cell_result)), ], cell_result)
  },
  .options = furrr_options(seed = seed),
  .progress = interactive()
)
t_elapsed <- as.numeric(Sys.time() - t_start, units = 'secs')
cat(sprintf('Completed in %.1f seconds\n', t_elapsed))

## Cell-level summary. true_value follows the docs/19 calibration,
## theta_true = -c_bm * sigma_BR with sigma_BR = 8.
summ <- results |>
  filter(!is.na(estimate)) |>
  group_by(t1half, design, c_bm, spec) |>
  summarise(
    n_conv        = n(),
    power         = mean(p_value < 0.05),
    mc_se_power   = sqrt(mean(p_value < 0.05) *
                         (1 - mean(p_value < 0.05)) / n()),
    mean_estimate = mean(estimate),
    empirical_se  = sd(estimate),
    mean_model_se = mean(std_error),
    kappa         = mean(std_error) / sd(estimate),
    .groups = 'drop'
  ) |>
  mutate(true_value = -c_bm * 8)

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(results = results, summary = summ, grid = grid,
       meta = list(script = '08-washout-diagnostic.R',
                   date = Sys.time(), seed = seed, n_reps = n_reps,
                   elapsed_secs = t_elapsed,
                   r_version = R.version.string)),
  file.path(out_dir, '08-washout-diagnostic.rds')
)

options(width = 200)
cat('\n===== TYPE I ERROR (c_bm = 0) =====\n')
print(as.data.frame(
  summ |> filter(c_bm == 0) |>
    select(t1half, design, spec, power, mc_se_power, kappa) |>
    arrange(t1half, design, spec)))

cat('\n===== POWER (c_bm = 0.45) =====\n')
print(as.data.frame(
  summ |> filter(c_bm == 0.45) |>
    select(t1half, design, spec, power, mean_estimate,
           empirical_se, mean_model_se) |>
    arrange(t1half, design, spec)))

message(sprintf('Wrote %s',
  file.path(out_dir, '08-washout-diagnostic.rds')))
