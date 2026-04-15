## analysis/02-carryover-sensitivity/01-run-factorial.R
##
## Factorial simulation driver for manuscript 02. Crosses DGP decay
## form x analysis-model carryover specification x architecture, and
## writes replicate-level results to output/01-factorial.rds.
##
## Usage:
##   Rscript analysis/02-carryover-sensitivity/01-run-factorial.R [--dev]
##
## --dev reduces replicates to 50 per cell for development timing;
## production omits the flag and runs 500 replicates.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(purrr)
  library(furrr)
  library(nlme)
})

repo_root <- here::here()
source(file.path(repo_root,
  'implementations/tidyverse/R/functions.R'))

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- '--dev' %in% args
n_reps <- if (dev_mode) 50 else 500

seed <- 20260415L
set.seed(seed)

grid <- tidyr::expand_grid(
  carryover_form  = c('linear', 'exponential', 'weibull'),
  weibull_shape   = c(0.7, 1.0, 1.5),
  analysis_spec   = c('A1', 'A2', 'A3'),
  dgp_arch        = c('mean_moderation', 'mvn'),
  t1half          = c(0.0, 0.5, 1.0),
  design          = c('CO', 'Hybrid', 'OLBDC'),
  N               = c(35, 70),
  c_bm            = c(0.0, 0.30, 0.45)
) |>
  filter(!(carryover_form != 'weibull' & weibull_shape != 1.0))

stopifnot(nrow(grid) > 0)

fit_one_cell <- function(row, n_reps) {
  ## TODO: replace with actual simulation body that:
  ##   1. builds trial design with row$design and row$t1half
  ##   2. calls generate_data() with model_param carrying
  ##      carryover_form, weibull_shape, dgp_architecture, c.bm
  ##   3. fits three analysis models (A1 binary, A2 Dbc,
  ##      A3 binary + lagged) and extracts the bm-by-treatment
  ##      interaction p-value and point estimate
  ##   4. returns a tibble of replicate-level results
  tibble::tibble(
    rep = seq_len(n_reps),
    p_value = NA_real_,
    estimate = NA_real_,
    converged = NA
  )
}

plan(multicore, workers = max(1, parallel::detectCores() - 1))

results <- grid |>
  mutate(row_id = row_number()) |>
  (\(g) future_map_dfr(
    seq_len(nrow(g)),
    \(i) bind_cols(g[i, ], fit_one_cell(g[i, ], n_reps)),
    .options = furrr_options(seed = seed),
    .progress = interactive()
  ))()

out_dir <- file.path(repo_root,
  'analysis/02-carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '01-run-factorial.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

saveRDS(
  list(results = results, grid = grid, meta = meta),
  file.path(out_dir, '01-factorial.rds')
)

message('Wrote ', file.path(out_dir, '01-factorial.rds'),
  ' (', nrow(results), ' replicate rows across ',
  nrow(grid), ' cells)')
