## analysis/scripts/carryover-sensitivity/01-run-factorial.R
##
## Factorial simulation driver for manuscript 02. Crosses DGP decay
## form x analysis-model carryover specification x architecture, and
## writes replicate-level results to output/01-factorial.rds.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/01-run-factorial.R [--dev] [--smoke]
##
## --dev   : 50 reps, reduced grid (96 cells). Keeps all DGP forms and
##           both architectures; trims t1half to {0, 1}, design to
##           {CO, Hybrid}, c_bm to {0, 0.45}, weibull_shape to 0.7.
##           Writes output/01-dev.rds (does not overwrite production).
## --smoke : 2 replicates, one cell per design (pipeline sanity only)
## default : 500 replicates, full 540-cell grid (production)

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
dev_mode <- '--dev' %in% args
smoke_mode <- '--smoke' %in% args
n_reps <- if (smoke_mode) 2 else if (dev_mode) 50 else 500

seed <- 20260415L
set.seed(seed)

if (smoke_mode) {
  grid <- tidyr::expand_grid(
    carryover_form  = 'exponential',
    weibull_shape   = 1.0,
    dgp_arch        = 'mvn',
    t1half          = 0.5,
    design          = c('OL', 'CO', 'OLBDC', 'Hybrid'),
    N               = 35,
    c_bm            = 0.45
  )
} else if (dev_mode) {
  ## Reduced grid for development. All DGP forms and both architectures
  ## are retained (the primary axes); secondary axes are trimmed to two
  ## representative values each. t1half = 0 provides the no-carryover
  ## power ceiling; t1half = 1.0 is the reference cell. c_bm = 0 gives
  ## the type-I error check; c_bm = 0.45 is the reference effect size.
  ## Weibull uses shape = 0.7 (slow-elimination tail, most distinct from
  ## exponential).
  grid <- tidyr::expand_grid(
    carryover_form  = c('linear', 'exponential', 'weibull'),
    weibull_shape   = c(0.7, 1.0),
    dgp_arch        = c('mean_moderation', 'mvn'),
    t1half          = c(0.0, 1.0),
    design          = c('CO', 'Hybrid'),
    N               = c(35, 70),
    c_bm            = c(0.0, 0.45)
  ) |>
    dplyr::filter(
      (carryover_form == 'weibull'  & weibull_shape == 0.7) |
      (carryover_form != 'weibull'  & weibull_shape == 1.0)
    )
} else {
  grid <- tidyr::expand_grid(
    carryover_form  = c('linear', 'exponential', 'weibull'),
    weibull_shape   = c(0.7, 1.0, 1.5),
    dgp_arch        = c('mean_moderation', 'mvn'),
    t1half          = c(0.0, 0.5, 1.0),
    design          = c('CO', 'Hybrid', 'OLBDC'),
    N               = c(35, 70),
    c_bm            = c(0.0, 0.30, 0.45)
  ) |>
    dplyr::filter(
      !(carryover_form != 'weibull' & weibull_shape != 1.0)
    )
}

stopifnot(nrow(grid) > 0)
cat(sprintf('Grid cells: %d, n_reps per cell: %d\n',
            nrow(grid), n_reps))

plan(
  if (smoke_mode) sequential else multicore,
  workers = if (smoke_mode) 1
            else max(1, parallel::detectCores() - 1)
)

t_start <- Sys.time()

results <- future_map_dfr(
  seq_len(nrow(grid)),
  function(i) {
    cell <- grid[i, ]
    cell_result <- simulate_cell(cell, n_reps)
    bind_cols(
      cell[rep(1, nrow(cell_result)), ],
      cell_result
    )
  },
  .options = furrr_options(seed = seed),
  .progress = interactive()
)

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f %s\n',
            as.numeric(t_elapsed, units = 'secs'), 'seconds'))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '01-run-factorial.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  smoke_mode = smoke_mode,
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

out_file <- if (smoke_mode) '01-smoke.rds' else if (dev_mode) '01-dev.rds' else '01-factorial.rds'
saveRDS(
  list(results = results, grid = grid, meta = meta),
  file.path(out_dir, out_file)
)

message(sprintf(
  'Wrote %s (%d rows across %d cells)',
  file.path(out_dir, out_file), nrow(results), nrow(grid)
))
