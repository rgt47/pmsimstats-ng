## analysis/02-carryover-sensitivity/04-run-sensitivity-blocks.R
##
## Tier 2 sensitivity blocks S1-S5, anchored at the reference
## configuration specified in simulation-grid-plan.md.
##
## Usage:
##   Rscript analysis/02-carryover-sensitivity/04-run-sensitivity-blocks.R [--dev]
##
## --dev : 50 reps/cell (development); default is 500 (production).

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
  'analysis/02-carryover-sensitivity/simulation-core.R'))

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- '--dev' %in% args
n_reps <- if (dev_mode) 50 else 500

seed <- 20260415L
set.seed(seed)

## Reference configuration (must match simulation-grid-plan.md)
ref <- tibble::tibble(
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  dgp_arch       = 'mvn',
  t1half         = 1.0,
  design         = 'Hybrid',
  N              = 70,
  c_bm           = 0.45
)

## --------------------------------------------------------------
## Block S1: AR(1) autocorrelation sensitivity
## --------------------------------------------------------------
S1 <- ref |>
  dplyr::slice(rep(1, 4)) |>
  dplyr::mutate(rho = c(0.5, 0.7, 0.8, 0.9),
                block = 'S1')

## --------------------------------------------------------------
## Block S2: analyst-truth half-life mismatch
## 2 true x 4 assumed = 8 cells; analysis_form/shape tied to
## the analyst's assumption (exponential).
## --------------------------------------------------------------
S2 <- tidyr::expand_grid(
  t1half = c(0.5, 1.0),
  analysis_t1half = c(0.25, 0.5, 1.0, 2.0)
) |>
  dplyr::mutate(
    carryover_form  = 'exponential',
    weibull_shape   = 1.0,
    analysis_form   = 'exponential',
    analysis_shape  = 1.0,
    dgp_arch        = ref$dgp_arch,
    design          = ref$design,
    N               = ref$N,
    c_bm            = ref$c_bm,
    block           = 'S2'
  )

## --------------------------------------------------------------
## Block S3: dropout rate x mechanism
## --------------------------------------------------------------
S3 <- tidyr::expand_grid(
  dropout_rate = c(0.0, 0.1, 0.2, 0.3),
  dropout_mech = c('MCAR', 'MAR')
) |>
  dplyr::filter(!(dropout_rate == 0 & dropout_mech == 'MAR')) |>
  dplyr::mutate(
    carryover_form = ref$carryover_form,
    weibull_shape  = ref$weibull_shape,
    dgp_arch       = ref$dgp_arch,
    t1half         = ref$t1half,
    design         = ref$design,
    N              = ref$N,
    c_bm           = ref$c_bm,
    block          = 'S3'
  )

## --------------------------------------------------------------
## Block S4: effect-size curve
## --------------------------------------------------------------
S4 <- tibble::tibble(
  c_bm = c(0.0, 0.10, 0.20, 0.30, 0.45, 0.60)
) |>
  dplyr::mutate(
    carryover_form = ref$carryover_form,
    weibull_shape  = ref$weibull_shape,
    dgp_arch       = ref$dgp_arch,
    t1half         = ref$t1half,
    design         = ref$design,
    N              = ref$N,
    block          = 'S4'
  )

## --------------------------------------------------------------
## Block S5: rho x carryover (exploratory)
## --------------------------------------------------------------
S5 <- tidyr::expand_grid(
  rho    = c(0.5, 0.8),
  t1half = c(0.0, 0.5, 1.0)
) |>
  dplyr::mutate(
    carryover_form = ref$carryover_form,
    weibull_shape  = ref$weibull_shape,
    dgp_arch       = ref$dgp_arch,
    design         = ref$design,
    N              = ref$N,
    c_bm           = ref$c_bm,
    block          = 'S5'
  )

## Bind with a common column set; missing columns -> NA
all_cols <- unique(c(names(S1), names(S2), names(S3),
                     names(S4), names(S5)))
fill_cols <- function(d) {
  for (col in setdiff(all_cols, names(d))) d[[col]] <- NA
  d[, all_cols]
}
grid <- dplyr::bind_rows(
  fill_cols(S1), fill_cols(S2), fill_cols(S3),
  fill_cols(S4), fill_cols(S5)
)

cat(sprintf('Sensitivity blocks: %d cells total\n', nrow(grid)))
cat(sprintf('  S1: %d  S2: %d  S3: %d  S4: %d  S5: %d\n',
            sum(grid$block == 'S1'), sum(grid$block == 'S2'),
            sum(grid$block == 'S3'), sum(grid$block == 'S4'),
            sum(grid$block == 'S5')))
cat(sprintf('n_reps per cell: %d\n', n_reps))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- future_map_dfr(
  seq_len(nrow(grid)),
  function(i) {
    cell <- grid[i, ]
    cell_result <- simulate_cell(cell, n_reps)
    dplyr::bind_cols(
      cell[rep(1, nrow(cell_result)), ],
      cell_result
    )
  },
  .options = furrr_options(seed = seed),
  .progress = interactive()
)

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root,
  'analysis/02-carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '04-run-sensitivity-blocks.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  reference = as.list(ref),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

out_file <- if (dev_mode) '04-sensitivity-dev.rds' else '04-sensitivity.rds'
saveRDS(list(results = results, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
