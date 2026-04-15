## analysis/02-carryover-sensitivity/02-summarise-grid.R
##
## Aggregate replicate-level results from 01-factorial.rds into
## cell-level power, type-I-error, bias, and Monte Carlo SE. Writes
## both an analysis-side copy and a manuscript-side mirror.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

repo_root <- here::here()

raw <- readRDS(file.path(repo_root,
  'analysis/02-carryover-sensitivity/output/01-factorial.rds'))

alpha <- 0.05

summary_grid <- raw$results |>
  group_by(across(all_of(c(names(raw$grid), 'spec')))) |>
  summarise(
    n_reps = n(),
    n_converged = sum(converged, na.rm = TRUE),
    power = mean(p_value < alpha, na.rm = TRUE),
    mc_se_power = sqrt(power * (1 - power) / n_reps),
    mean_estimate = mean(estimate, na.rm = TRUE),
    sd_estimate = sd(estimate, na.rm = TRUE),
    converged_frac = mean(converged, na.rm = TRUE),
    .groups = 'drop'
  )

out_analysis <- file.path(repo_root,
  'analysis/02-carryover-sensitivity/output/02-grid-summary.rds')
out_manuscript <- file.path(repo_root,
  'manuscripts/data/02-grid-summary.rds')

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
## and mirror to manuscripts/data/02-sensitivity-summary.rds. When
## the Tier 2 production run has not yet been executed, this
## section is a no-op.
## -----------------------------------------------------------------

tier2_path <- file.path(repo_root,
  'analysis/02-carryover-sensitivity/output/04-sensitivity.rds')

if (file.exists(tier2_path)) {
  raw2 <- readRDS(tier2_path)

  summary_grid2 <- raw2$results |>
    group_by(across(all_of(c(names(raw2$grid), 'spec')))) |>
    summarise(
      n_reps = n(),
      n_converged = sum(converged, na.rm = TRUE),
      power = mean(p_value < alpha, na.rm = TRUE),
      mc_se_power = sqrt(power * (1 - power) / n_reps),
      mean_estimate = mean(estimate, na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      converged_frac = mean(converged, na.rm = TRUE),
      .groups = 'drop'
    )

  out2_analysis <- file.path(repo_root,
    'analysis/02-carryover-sensitivity/output/04-sensitivity-summary.rds')
  out2_manuscript <- file.path(repo_root,
    'manuscripts/data/02-sensitivity-summary.rds')

  saveRDS(list(summary = summary_grid2, meta = raw2$meta),
          out2_analysis)
  saveRDS(list(summary = summary_grid2, meta = raw2$meta),
          out2_manuscript)

  message('Wrote ', out2_analysis)
  message('Wrote ', out2_manuscript)
} else {
  message('Skipped Tier 2: ', tier2_path, ' not found.')
}
