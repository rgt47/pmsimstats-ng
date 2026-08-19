## analysis/scripts/carryover-sensitivity/14-run-sensitivity-s9.R
##
## Tier 2 sensitivity block S9: does the CO design's longer, graded
## washout (off-drug occasions span 2.5-20 weeks, versus 1-2 weeks
## under Hybrid/OL+BDC) change how the alternative-paradigm
## specifications (E7, E9; see 09-run-sensitivity-s7.R) perform
## relative to E1, or is the Hybrid reference-cell ranking
## design-general?
##
## An earlier version of this block tested E5 (Lag x bm) / E6
## (Washout x bm) here instead, asking whether CO's graded t_sd
## range relieved E6's collinearity problem. It did (Section 4.5),
## but the joint-test power gap did not close and was slightly wider
## under CO than under Hybrid, consistent with the structural
## over-parameterization diagnosis that led to replacing E5/E6 with
## E7/E9 (simulation-core.R). This version re-asks the analogous
## design-transferability question for the replacement
## specifications.
##
## Two-cell grid: design = CO (all else as the S7/S8 reference cell)
## crossed with t1/2 in {1.0, 0.5}. Full five-specification G1-G5 set,
## matching S7/S8, so the CO ranking is directly comparable to the
## Hybrid one at every specification, not just E1/E7/E9.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/14-run-sensitivity-s9.R [--dev] [--reps N]
##
## --dev    : 50 reps/cell (development); default is 500 (production).
## --reps N : override the rep count for the current mode.

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
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

specs_s9 <- c('E1', 'E3', 'E2', 'E7', 'E9')

grid <- tibble::tibble(t1half = c(1.0, 0.5)) |>
  dplyr::mutate(
    carryover_form = 'exponential',
    weibull_shape  = 1.0,
    dgp_arch       = 'mvn',
    design         = 'CO',
    N              = 70,
    c_bm           = 0.45,
    block          = 'S9'
  )

cat(sprintf('S9: %d cells (CO design, t1/2), n_reps = %d, specs = %s\n',
            nrow(grid), n_reps, paste(specs_s9, collapse = ', ')))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  cell <- grid[i, ]
  simulate_cell_s7(cell, n_reps, specs = specs_s9) |>
    dplyr::mutate(t1half = cell$t1half, .before = 1)
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '14-run-sensitivity-s9.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_s9,
  grid = as.list(grid),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_s9 <- results |>
  group_by(t1half, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    mean_estimate_bmDb = mean(estimate_bmDb, na.rm = TRUE),
    mean_p_bmDb_reject = mean(p_value_bmDb < 0.05, na.rm = TRUE),
    .groups = 'drop'
  )

print(as.data.frame(summary_s9))

out_file <- if (dev_mode) '14-sensitivity-s9-dev.rds' else '14-sensitivity-s9.rds'
saveRDS(list(results = results, summary = summary_s9,
             grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
