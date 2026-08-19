## analysis/scripts/carryover-sensitivity/24-run-decay-shape-sensitivity-archA.R
##
## Architecture A (mean moderation) counterpart to
## 23-run-decay-shape-sensitivity.R, for Supplement Section S12
## Panel C (the Architecture A counterpart to Figure 3). Same
## 30-cell design and more extreme Weibull range (k = 0.25, 0.5,
## 2.0, 4.0; k = 1.0 dropped as numerically identical to
## Exponential) as the Architecture B version, only dgp_arch
## switched from 'mvn' to 'mean_moderation'. This keeps S12 Panel C
## internally consistent with the main manuscript's Figure 3, rather
## than reading the unmodified shared production grid (still at the
## earlier k = 0.7, 1.0, 1.5) as it did previously.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/24-run-decay-shape-sensitivity-archA.R [--dev] [--reps N]
##
## --dev    : 3 reps/cell (smoke test); default is 250 (matching
##            S12 Panels A/B's precision).
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
          if (dev_mode) 3 else 250

seed <- 20260415L
set.seed(seed)

specs_tier1 <- c('E1', 'E3', 'E2')

forms <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70)
) |>
  dplyr::mutate(c_bm = 0.45, t1half = 1.0, dgp_arch = 'mean_moderation')

grid <- dplyr::bind_rows(
  forms |> dplyr::mutate(carryover_form = 'exponential', weibull_shape = 1.0),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 0.25),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 0.5),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 2.0),
  forms |> dplyr::mutate(carryover_form = 'weibull', weibull_shape = 4.0)
)

cat(sprintf('Decay-shape sensitivity, Architecture A (k=0.25, 0.5, 2.0, 4.0; G1-G3): %d cells, n_reps = %d\n',
            nrow(grid), n_reps))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  cell <- grid[i, ]
  cell_result <- simulate_cell_s7(cell, n_reps, specs = specs_tier1)
  dplyr::bind_cols(cell[rep(1, nrow(cell_result)), ], cell_result)
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root, 'analysis/data')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '24-run-decay-shape-sensitivity-archA.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_tier1,
  weibull_shapes = c(0.25, 0.5, 2.0, 4.0),
  dgp_arch = 'mean_moderation',
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_tab <- results |>
  group_by(design, N, carryover_form, weibull_shape, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    .groups = 'drop'
  )

print(as.data.frame(summary_tab) |>
       dplyr::select(design, N, carryover_form, weibull_shape, spec, power) |>
       dplyr::arrange(design, N, carryover_form, weibull_shape, spec))

out_file <- if (dev_mode) '02-decay-shape-sensitivity-archA-dev.rds' else '02-decay-shape-sensitivity-archA.rds'
saveRDS(list(results = results, summary = summary_tab, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
