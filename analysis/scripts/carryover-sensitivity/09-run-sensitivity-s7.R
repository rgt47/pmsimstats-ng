## analysis/scripts/carryover-sensitivity/09-run-sensitivity-s7.R
##
## Tier 2 sensitivity block S7: biomarker-interacted carryover
## terms, at the same reference cell used by S1-S4
## (Hybrid, N = 70, exponential DGP, t1/2 = 1.0, c_bm = 0.45).
##
## Compares four analysis specifications, all targeting bm:Db:
##   E1  Unadjusted        Sx ~ bm + t + Db + bm:Db
##   E4  Washout-adjusted  Sx ~ bm + t + Db + bm:Db + tsd
##   E5  Lag x bm          Sx ~ bm + t + Db + bm:Db + L   + bm:L
##   E6  Washout x bm      Sx ~ bm + t + Db + bm:Db + tsd + bm:tsd
##
## E5/E6 add a second interaction coefficient (bm:L or bm:tsd), so
## "power" for those two specifications is the rejection rate of a
## joint 2-df Wald test of bm:Db and the second coefficient together
## (nlme::anova.lme(Terms = ...)), not a single-coefficient test as
## for E1/E4. See simulation-core.R, "Biomarker-interacted carryover
## terms (S7)", for the rationale.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/09-run-sensitivity-s7.R [--dev] [--reps N]
##
## --dev    : 50 reps (development); default is 500 (production).
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

## Reference configuration (must match simulation-grid-plan.md /
## the S1-S4 reference cell in 04-run-sensitivity-blocks.R)
ref <- tibble::tibble(
  carryover_form = 'exponential',
  weibull_shape  = 1.0,
  dgp_arch       = 'mvn',
  t1half         = 1.0,
  design         = 'Hybrid',
  N              = 70,
  c_bm           = 0.45,
  block          = 'S7'
)

cat(sprintf('S7: 1 cell, n_reps = %d\n', n_reps))

plan(multicore, workers = max(1, parallel::detectCores() - 1))

t_start <- Sys.time()

results <- simulate_cell_s7(ref, n_reps)

t_elapsed <- Sys.time() - t_start
cat(sprintf('Completed in %.1f seconds\n',
            as.numeric(t_elapsed, units = 'secs')))

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- list(
  script = '09-run-sensitivity-s7.R',
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

summary_s7 <- results |>
  group_by(spec) |>
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

print(summary_s7)

out_file <- if (dev_mode) '09-sensitivity-s7-dev.rds' else '09-sensitivity-s7.rds'
saveRDS(list(results = results, summary = summary_s7,
             ref = ref, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows)',
        file.path(out_dir, out_file), nrow(results)))
