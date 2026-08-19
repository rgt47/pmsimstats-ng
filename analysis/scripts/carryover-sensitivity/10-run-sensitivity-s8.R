## analysis/scripts/carryover-sensitivity/10-run-sensitivity-s8.R
##
## Tier 2 sensitivity block S8: DGP architecture (A, mean
## moderation, vs B, MVN differential correlation) crossed against
## five analysis specifications and three carryover half-lives
## (t1/2 in {0, 0.5, 1.0}), at the same reference cell used by
## S1-S4/S7 otherwise (Hybrid, N = 70, exponential DGP,
## c_bm = 0.45).
##
## This is NOT a re-run of paper 01's architecture comparison (power
## loss under carryover, by architecture). It asks a narrower
## question specific to this manuscript's scope: does the ranking
## among analysis specifications established under Architecture B
## (Sections 3.3-3.4) hold under Architecture A, or is the choice of
## carryover representation itself architecture-dependent? The rest
## of the manuscript's grid remains restricted to Architecture B per
## its stated scope; this block is the one architecture-sensitivity
## check.
##
## Five specifications (E4 excluded; it was a one-off diagnostic
## probe in Section 4.1, not a formal arm):
##   E1  Unadjusted        Sx ~ bm + t + Db  + bm:Db
##   E3  Lag-adjusted      Sx ~ bm + t + Db  + bm:Db + L
##   E2  Exposure-weighted Sx ~ bm + t + Dbc + bm:Dbc
##   E5  Lag x bm          Sx ~ bm + t + Db  + bm:Db + L   + bm:L
##   E6  Washout x bm      Sx ~ bm + t + Db  + bm:Db + tsd + bm:tsd
##
## E5/E6 report the joint 2-df Wald test (see simulation-core.R,
## "Biomarker-interacted carryover terms (S7)"); E1/E3/E2 report
## their single-coefficient test as in the Tier 1 grid.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/10-run-sensitivity-s8.R [--dev] [--reps N]
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

specs_s8 <- c('E1', 'E3', 'E2', 'E5', 'E6')

## Reference configuration (must match the S1-S4/S7 reference cell,
## except t1half, which is now a 3-level grid axis rather than fixed)
ref <- tibble::tibble(
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  design         = 'Hybrid',
  N              = 70,
  c_bm           = 0.45
)

grid <- tidyr::expand_grid(
    dgp_arch = c('mvn', 'mean_moderation'),
    t1half   = c(0, 0.5, 1.0)
  ) |>
  dplyr::cross_join(ref) |>
  dplyr::mutate(block = 'S8')

cat(sprintf('S8: %d cells (architecture x half-life), n_reps = %d, specs = %s\n',
            nrow(grid), n_reps, paste(specs_s8, collapse = ', ')))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- purrr::map_dfr(seq_len(nrow(grid)), function(i) {
  cell <- grid[i, ]
  simulate_cell_s7(cell, n_reps, specs = specs_s8) |>
    dplyr::mutate(dgp_arch = cell$dgp_arch, t1half = cell$t1half,
                 .before = 1)
})

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '10-run-sensitivity-s8.R',
  date = Sys.time(),
  seed = seed,
  n_reps = n_reps,
  dev_mode = dev_mode,
  specs = specs_s8,
  reference = as.list(ref),
  elapsed_secs = as.numeric(t_elapsed, units = 'secs'),
  r_version = R.version.string,
  git_sha = tryCatch(
    system('git rev-parse --short HEAD', intern = TRUE),
    error = function(e) NA_character_
  )
)

summary_s8 <- results |>
  group_by(dgp_arch, t1half, spec) |>
  summarise(
    n = n(),
    n_converged = sum(converged, na.rm = TRUE),
    non_convergence_rate = 1 - n_converged / n,
    power = mean(p_value < 0.05, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_converged),
    mean_estimate_bmDb = mean(estimate_bmDb, na.rm = TRUE),
    .groups = 'drop'
  )

print(as.data.frame(summary_s8))

out_file <- if (dev_mode) '10-sensitivity-s8-dev.rds' else '10-sensitivity-s8.rds'
saveRDS(list(results = results, summary = summary_s8,
             grid = grid, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows across %d cells)',
        file.path(out_dir, out_file), nrow(results), nrow(grid)))
