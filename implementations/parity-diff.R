#!/usr/bin/env Rscript
# parity-diff.R
#
# Compares two parity reports (default: parity_baseline.rds vs
# parity_report.rds) and prints the set of cells that differ. Exits 0
# if the reports are identical, 1 otherwise.
#
# Usage:
#   Rscript implementations/parity-diff.R
#   Rscript implementations/parity-diff.R <baseline.rds> <current.rds>

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
baseline_path <- if (length(args) >= 1) {
  args[1]
} else {
  'implementations/parity_baseline.rds'
}
current_path <- if (length(args) >= 2) {
  args[2]
} else {
  'implementations/parity_report.rds'
}

for (p in c(baseline_path, current_path)) {
  if (!file.exists(p)) {
    stop(sprintf('Report file not found: %s', p), call. = FALSE)
  }
}

a <- readRDS(baseline_path)
b <- readRDS(current_path)

cat(sprintf('Baseline:  %s  (run %s, mode=%s)\n',
            baseline_path, format(a$timestamp), a$mode))
cat(sprintf('Current :  %s  (run %s, mode=%s)\n',
            current_path, format(b$timestamp), b$mode))

# ---- Grid (sigma + data) diff --------------------------------------------

if (!identical(dim(a$grid), dim(b$grid)) ||
    !identical(sort(colnames(a$grid)), sort(colnames(b$grid)))) {
  cat('\nGrid shape changed between baseline and current.\n')
  cat(sprintf('  baseline: %d rows x %d cols\n',
              nrow(a$grid), ncol(a$grid)))
  cat(sprintf('  current : %d rows x %d cols\n',
              nrow(b$grid), ncol(b$grid)))
  quit(status = 1)
}

join_keys <- c('design', 'arch', 't_half', 'c_bm', 'seed')
joined <- dplyr::inner_join(
  a$grid |> dplyr::rename_with(~ paste0('a_', .), -dplyr::all_of(join_keys)),
  b$grid |> dplyr::rename_with(~ paste0('b_', .), -dplyr::all_of(join_keys)),
  by = join_keys)

grid_diff <- joined |>
  dplyr::rowwise() |>
  dplyr::mutate(
    changed = !isTRUE(a_sigma_labels == b_sigma_labels) ||
              !isTRUE(a_sigma_means  == b_sigma_means)  ||
              !isTRUE(a_sigma_sigma  == b_sigma_sigma)  ||
              !isTRUE(a_data_pass    == b_data_pass)    ||
              !isTRUE(a_data_total   == b_data_total)) |>
  dplyr::ungroup() |>
  dplyr::filter(changed)

# ---- Analysis diff --------------------------------------------------------

if (is.null(a$analysis) || is.null(b$analysis) ||
    !identical(sort(colnames(a$analysis)), sort(colnames(b$analysis)))) {
  cat('\nAnalysis report shape changed or missing.\n')
  an_diff <- tibble::tibble()
} else {
  an_joined <- dplyr::inner_join(
    a$analysis |> dplyr::rename_with(~ paste0('a_', .),
                                      -dplyr::all_of(join_keys)),
    b$analysis |> dplyr::rename_with(~ paste0('b_', .),
                                      -dplyr::all_of(join_keys)),
    by = join_keys)
  an_diff <- an_joined |>
    dplyr::rowwise() |>
    dplyr::mutate(
      changed = !isTRUE(a_beta_ok == b_beta_ok) ||
                !isTRUE(a_se_ok   == b_se_ok)   ||
                !isTRUE(a_p_ok    == b_p_ok)) |>
    dplyr::ungroup() |>
    dplyr::filter(changed)
}

# ---- Report ---------------------------------------------------------------

cat(sprintf('\nSigma/Data cells with changed outcome : %d\n',
            nrow(grid_diff)))
cat(sprintf('Analysis cells with changed outcome   : %d\n',
            nrow(an_diff)))

if (nrow(grid_diff) > 0) {
  cat('\nSigma/Data changes:\n')
  print(grid_diff |> dplyr::select(dplyr::all_of(join_keys),
    a_sigma_labels, b_sigma_labels,
    a_sigma_means, b_sigma_means, a_sigma_sigma, b_sigma_sigma,
    a_data_pass, b_data_pass, a_data_total, b_data_total))
}
if (nrow(an_diff) > 0) {
  cat('\nAnalysis changes:\n')
  print(an_diff |> dplyr::select(dplyr::all_of(join_keys),
    a_beta_ok, b_beta_ok, a_se_ok, b_se_ok, a_p_ok, b_p_ok))
}

if (nrow(grid_diff) == 0 && nrow(an_diff) == 0) {
  cat('\nNO CHANGES: current report matches baseline.\n')
  quit(status = 0)
} else {
  cat('\nReports differ.\n')
  quit(status = 1)
}
