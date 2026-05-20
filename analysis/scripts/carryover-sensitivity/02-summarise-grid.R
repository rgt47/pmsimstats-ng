## analysis/scripts/carryover-sensitivity/02-summarise-grid.R
##
## Aggregate replicate-level results from 01-factorial.rds into
## cell-level Morris-style performance measures with Monte Carlo
## standard errors. Writes both an analysis-side copy and a
## manuscript-side mirror.
##
## Usage:
##   Rscript analysis/scripts/carryover-sensitivity/02-summarise-grid.R [--dev]
##
## --dev : read from 01-dev.rds / 04-sensitivity-dev.rds and write
##         the summaries to the standard manuscript data paths
##         (analysis/data/02-grid-summary.rds etc.), overwriting
##         production summaries. Use for review renders only; restore
##         production summaries with git checkout analysis/data/ when
##         done.
##
## Reporting follows Morris, White, and Crowther (2019), 'Using
## simulation studies to evaluate statistical methods', Statistics
## in Medicine 38(11):2074-2102, doi 10.1002/sim.8086. Per-cell
## columns and their MC SE formulas correspond to Table 6 of that
## paper.
##
## Coverage of nominal 95% Wald confidence intervals is computed
## from the per-replicate `std_error` column emitted by
## simulation-core.R (added 2026-05-06). For older raw outputs
## that lack `std_error`, the coverage column is set to NA.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
dev_mode <- '--dev' %in% args

repo_root <- here::here()

tier1_file <- if (dev_mode) '01-dev.rds' else '01-factorial.rds'
if (dev_mode) message('Dev mode: reading ', tier1_file,
                      ' (production summaries will be overwritten)')

raw <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output', tier1_file))

alpha <- 0.05

## Cell-level summarisation. Morris (2019) Table 6 column names in
## comments below.
summarise_cells <- function(results, grid_names) {
  ## Handle legacy raw outputs that lack `std_error`: insert a
  ## column of NAs so that the coverage summarise step still
  ## produces a numeric (NA) value rather than erroring.
  if (!('std_error' %in% colnames(results))) {
    results <- results |> mutate(std_error = NA_real_)
  }
  results |>
    group_by(across(all_of(c(grid_names, 'spec')))) |>
    summarise(
      ## Counts
      n_reps = n(),
      n_converged = sum(converged, na.rm = TRUE),

      ## Power and Type I error
      power = mean(p_value < alpha, na.rm = TRUE),
      mc_se_power = sqrt(power * (1 - power) / n_reps),

      ## Calibrated true value on the analysis-coefficient scale.
      ## Per docs/19-mean-moderation-implementation-notes.tex the
      ## population effect size is c_bm * sigma_BR units of BR per
      ## SD of biomarker. The analysis model fits Sx = ... +
      ## bm:X + ..., where Sx is the symptom score (BL - factors)
      ## so an increase in BR (better drug response in higher-
      ## biomarker participants) corresponds to a NEGATIVE
      ## analysis-model coefficient on bm:X. The calibrated
      ## estimand on the analysis-coefficient scale is therefore
      ##   true_value = - c_bm * sigma_BR
      ## with sigma_BR = 8 in default_resp_param() (the Hendrickson
      ## (2020) calibration target). Architecture A and Architecture
      ## B share this estimand under the calibration, so bias is
      ## commensurable across DGMs in this column.
      sigma_BR = 8,
      true_value = -first(c_bm) * sigma_BR,
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean_estimate - true_value,
      empirical_se = sd(estimate, na.rm = TRUE),
      mc_se_bias = empirical_se / sqrt(n_reps),

      ## MSE = bias^2 + empirical variance.
      ## MC SE for MSE follows Morris Table 6:
      ##   MC SE(MSE) ~ sqrt(Var((est - truth)^2) / n_reps)
      ## computed here as the empirical SD of squared deviations
      ## divided by sqrt(n_reps).
      mse = bias^2 + empirical_se^2,
      mc_se_mse = sd((estimate - true_value)^2, na.rm = TRUE) /
        sqrt(n_reps),

      ## Coverage of nominal 95% Wald CI for the interaction
      ## coefficient: |est - true_value| < z_{0.975} * std_error.
      ## Requires std_error in the raw data (added 2026-05-06).
      coverage = mean(
        abs(estimate - true_value) <
          qnorm(0.975) * std_error,
        na.rm = TRUE),
      mc_se_coverage = sqrt(
        coverage * (1 - coverage) / n_reps),

      ## Non-convergence
      converged_frac = mean(converged, na.rm = TRUE),
      non_convergence_rate = 1 - converged_frac,
      mc_se_non_convergence = sqrt(
        non_convergence_rate * (1 - non_convergence_rate) / n_reps),

      .groups = 'drop'
    )
}

grid_names <- names(raw$grid)
summary_grid <- summarise_cells(raw$results, grid_names)

## Sanity checks: every cell has the configured n_reps; bias is
## finite where any replicate converged; MC SE columns are non-
## negative.
stopifnot(all(summary_grid$n_reps == raw$meta$n_reps))
stopifnot(all(is.finite(summary_grid$bias) |
              summary_grid$n_converged == 0))
stopifnot(all(summary_grid$mc_se_power >= 0,
              summary_grid$mc_se_bias >= 0,
              summary_grid$mc_se_mse >= 0,
              summary_grid$mc_se_non_convergence >= 0))

out_analysis <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/02-grid-summary.rds')
out_manuscript <- file.path(repo_root,
  'analysis/data/02-grid-summary.rds')

saveRDS(
  list(summary = summary_grid, meta = raw$meta),
  out_analysis
)
saveRDS(
  list(summary = summary_grid, meta = raw$meta),
  out_manuscript
)

message('Wrote ', out_analysis)
message('Wrote ', out_manuscript)

## -----------------------------------------------------------------
## Tier 2 sensitivity blocks (optional)
##
## If output/04-sensitivity.rds exists, summarise it analogously
## and mirror to analysis/data/02-sensitivity-summary.rds.
## -----------------------------------------------------------------

tier2_file <- if (dev_mode) '04-sensitivity-dev.rds' else '04-sensitivity.rds'
tier2_path <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output', tier2_file)

if (file.exists(tier2_path)) {
  raw2 <- readRDS(tier2_path)
  summary_grid2 <- summarise_cells(raw2$results, names(raw2$grid))

  out2_analysis <- file.path(repo_root,
    'analysis/scripts/carryover-sensitivity/output/04-sensitivity-summary.rds')
  out2_manuscript <- file.path(repo_root,
    'analysis/data/02-sensitivity-summary.rds')

  saveRDS(list(summary = summary_grid2, meta = raw2$meta),
          out2_analysis)
  saveRDS(list(summary = summary_grid2, meta = raw2$meta),
          out2_manuscript)

  message('Wrote ', out2_analysis)
  message('Wrote ', out2_manuscript)
} else {
  message('Skipped Tier 2: ', tier2_path, ' not found.')
}
