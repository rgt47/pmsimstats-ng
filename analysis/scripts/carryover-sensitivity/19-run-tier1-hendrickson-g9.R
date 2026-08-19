## analysis/scripts/carryover-sensitivity/19-run-tier1-hendrickson-g9.R
##
## Extends the Tier 1 Hendrickson-style grid figures (Figures 1-2,
## fig-hendrickson-a/b, script 12-hendrickson-heatmap.R Panels A/B)
## from the original three specifications (G1-G3) to the full
## nine-specification G1-G9 set, via simulate_cell_s7()/fit_spec_s7()
## (the same pipeline used for the S1-S4/S6/S7-S11 extensions).
##
## Panel A varies c_bm in {0, 0.30, 0.45} at t1half = 1.0. Panel B
## varies t1half in {0, 0.5, 1.0} at c_bm = 0.45. Both are restricted
## to the exponential DGP and Architecture B (mvn), matching what
## 12-hendrickson-heatmap.R currently plots for Panels A/B, and cross
## design in {CO, Hybrid, OLBDC} with N in {35, 70}. The shared cell
## (t1half = 1.0, c_bm = 0.45) is fit once, not twice.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/19-run-tier1-hendrickson-g9.R [--dev] [--reps N]
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

panel_a <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70),
  c_bm   = c(0, 0.30, 0.45)
) |>
  dplyr::mutate(t1half = 1.0, panel = 'A')

panel_b <- tidyr::expand_grid(
  design = c('CO', 'Hybrid', 'OLBDC'),
  N      = c(35, 70),
  t1half = c(0, 0.5)
) |>
  dplyr::mutate(c_bm = 0.45, panel = 'B')

grid <- dplyr::bind_rows(panel_a, panel_b) |>
  dplyr::mutate(carryover_form = 'exponential', weibull_shape = 1.0,
               dgp_arch = 'mvn')

cat(sprintf('Tier 1 Hendrickson A/B (G1-G9): %d unique cells, n_reps = %d\n',
            nrow(grid), n_reps))

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
  script = '19-run-tier1-hendrickson-g9.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_g9,
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_g9 <- results |>
  group_by(design, N, c_bm, t1half, panel, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    .groups = 'drop'
  )

print(as.data.frame(summary_g9) |> dplyr::select(panel, design, N, c_bm, t1half, spec, power) |>
       dplyr::arrange(panel, design, N, spec) |> head(40))

out_file <- if (dev_mode) '02-grid-summary-hendrickson-g9-dev.rds' else '02-grid-summary-hendrickson-g9.rds'
saveRDS(list(results = results, summary = summary_g9, grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
