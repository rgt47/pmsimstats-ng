## analysis/scripts/carryover-sensitivity/17-run-sensitivity-s1234-g9.R
##
## Extends Blocks S1-S4 (originally run by 04-run-sensitivity-blocks.R
## with only G1-G3, via the Tier 1 fit_all_specs() pipeline) to the
## full nine-specification G1-G9 set, via simulate_cell_s7()/
## fit_spec_s7() instead. Produces a second, independent summary
## (`analysis/data/02-sensitivity-summary-g9.rds`) alongside the
## original three-spec `02-sensitivity-summary.rds`; the existing
## S1-S4 line plots (04-sensitivity-figures-extra-slim.R) continue to
## read the original file and are unaffected. Only the S1-S4
## Hendrickson-style heatmaps (13-hendrickson-heatmap-tier2.R) are
## updated to read this new file, since a 9-line line plot would be
## visually unreadable but a 9-column heatmap is not.
##
## Same four cell grids as 04-run-sensitivity-blocks.R (S1
## autocorrelation, S2 half-life mismatch, S3 dropout, S4 effect
## size); S5 (exploratory, archived-drafts-only) is not included.
## G8 (E2cr2) is fit directly here, at every cell, sharing common
## random numbers with G1-G7/G9 (unlike Block S10, which spliced in
## an independently-seeded G8 from Block S6).
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/17-run-sensitivity-s1234-g9.R [--dev] [--reps N]
##
## --dev    : 50 reps/cell (development); default is 500 (production).
## --reps N : override the rep count for the current mode.

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
reps_idx <- which(args == '--reps')
n_reps_override <- if (length(reps_idx) && reps_idx < length(args))
  as.integer(args[reps_idx + 1]) else NA_integer_
n_reps <- if (!is.na(n_reps_override)) n_reps_override else
          if (dev_mode) 50 else 500

seed <- 20260415L
set.seed(seed)

specs_g9 <- c('E1', 'E3', 'E2', 'E7', 'E9',
             'E1cr2', 'E3cr2', 'E2cr2', 'E7cr2')

## Reference configuration (must match simulation-grid-plan.md /
## 04-run-sensitivity-blocks.R)
ref <- tibble::tibble(
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  dgp_arch       = 'mvn',
  t1half         = 1.0,
  design         = 'Hybrid',
  N              = 70,
  c_bm           = 0.45
)

S1 <- ref |>
  dplyr::slice(rep(1, 4)) |>
  dplyr::mutate(rho = c(0.5, 0.7, 0.8, 0.9), block = 'S1')

S2 <- tidyr::expand_grid(
  t1half = c(0.5, 1.0),
  analysis_t1half = c(0.25, 0.5, 1.0, 2.0)
) |>
  dplyr::mutate(
    carryover_form  = 'exponential', weibull_shape = 1.0,
    analysis_form   = 'exponential', analysis_shape = 1.0,
    dgp_arch = ref$dgp_arch, design = ref$design,
    N = ref$N, c_bm = ref$c_bm, block = 'S2'
  )

S3 <- tidyr::expand_grid(
  dropout_rate = c(0.0, 0.1, 0.2, 0.3),
  dropout_mech = c('MCAR', 'MAR')
) |>
  dplyr::filter(!(dropout_rate == 0 & dropout_mech == 'MAR')) |>
  dplyr::mutate(
    carryover_form = ref$carryover_form, weibull_shape = ref$weibull_shape,
    dgp_arch = ref$dgp_arch, t1half = ref$t1half, design = ref$design,
    N = ref$N, c_bm = ref$c_bm, block = 'S3'
  )

S4 <- tibble::tibble(c_bm = c(0.0, 0.10, 0.20, 0.30, 0.45, 0.60)) |>
  dplyr::mutate(
    carryover_form = ref$carryover_form, weibull_shape = ref$weibull_shape,
    dgp_arch = ref$dgp_arch, t1half = ref$t1half, design = ref$design,
    N = ref$N, block = 'S4'
  )

all_cols <- unique(c(names(S1), names(S2), names(S3), names(S4)))
fill_cols <- function(d) {
  for (col in setdiff(all_cols, names(d))) d[[col]] <- NA
  d[, all_cols]
}
grid <- dplyr::bind_rows(fill_cols(S1), fill_cols(S2),
                         fill_cols(S3), fill_cols(S4))

cat(sprintf('S1-S4 (G1-G9): %d cells total\n', nrow(grid)))
cat(sprintf('  S1: %d  S2: %d  S3: %d  S4: %d\n',
            sum(grid$block == 'S1'), sum(grid$block == 'S2'),
            sum(grid$block == 'S3'), sum(grid$block == 'S4')))
cat(sprintf('n_reps per cell: %d, specs: %s\n',
            n_reps, paste(specs_g9, collapse = ', ')))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  cell <- grid[i, ]
  cell_result <- simulate_cell_s7(cell, n_reps, specs = specs_g9)
  dplyr::bind_cols(cell[rep(1, nrow(cell_result)), ], cell_result)
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root, 'analysis/data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '17-run-sensitivity-s1234-g9.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_g9,
  reference = as.list(ref),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_g9 <- results |>
  group_by(across(all_of(names(grid))), spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    mean_best_t_half = mean(best_t_half, na.rm = TRUE),
    .groups = 'drop'
  )

print(as.data.frame(summary_g9) |> dplyr::select(block, spec, power) |>
       dplyr::arrange(block, spec) |> head(40))

out_file <- if (dev_mode) '02-sensitivity-summary-g9-dev.rds' else '02-sensitivity-summary-g9.rds'
saveRDS(list(results = results, summary = summary_g9, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
