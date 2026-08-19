## analysis/scripts/carryover-sensitivity/09-run-sensitivity-s7.R
##
## Tier 2 sensitivity block S7: alternative-paradigm specifications
## addressing this manuscript's two open methodological questions
## (unknown decay half-life; is a mixed model even necessary), at
## the same reference cell used by S1-S4 (Hybrid, N = 70,
## exponential DGP, t1/2 = 1.0, c_bm = 0.45).
##
## Compares three analysis specifications:
##   E1  Unadjusted           Sx ~ bm + t + Db + bm:Db
##   E7  AIC-selected t1/2    E2's formula, half-life chosen by AIC
##                            over {0.25, 0.5, 1.0, 2.0} rather than
##                            assumed (targets bm:Dbc)
##   E9  Paired-difference    Per-patient mean on-drug minus
##                            off-drug Sx, regressed on bm across
##                            patients (diff ~ bm, plain OLS); a
##                            biomarker-interaction extension of the
##                            paired t-test Senn, Julious & Araujo
##                            (2014) found outperforms
##                            carryover-adjusted mixed models for
##                            ordinary N-of-1 treatment-effect
##                            estimation
##
## An earlier version of this block tested E5 (Lag x bm) and E6
## (Washout x bm), carryover terms crossed with the biomarker; both
## underperformed E1 at every cell examined across S7-S9 and were
## replaced by E7/E9, which probe the analyst's unknown decay
## half-life and inferential paradigm choice directly, rather than
## adding a second, structurally redundant interaction coefficient
## (Section 4.5). See
## simulation-core.R, "Biomarker-interacted carryover terms and
## their replacements (S7)", for the full rationale.
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
    mean_best_t_half = mean(best_t_half, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_s7)

out_file <- if (dev_mode) '09-sensitivity-s7-dev.rds' else '09-sensitivity-s7.rds'
saveRDS(list(results = results, summary = summary_s7,
             ref = ref, meta = meta),
        file.path(out_dir, out_file))

message(sprintf('Wrote %s (%d rows)',
        file.path(out_dir, out_file), nrow(results)))
