## analysis/scripts/dgp-combined/02-summarise-combined.R
##
## Summaries for paper 11 from the output of
## 01-run-combined-factorial.R (2026-10-29 rebuild). All cells are at a
## total of N = 70.
##
## Inputs (analysis/data/dgp-combined/):
##   01-combined-replicates.rds       main pass (E2, E2cr2), nine cells
##   01-combined-replicates-g9.rds    optional G9 pass
##   01-combined-replicates-hl02.rds  optional t1/2 = 0.2 pass
##   03-analytic.rds                  from 03-analytic-projection.R
##
## Outputs (analysis/data/dgp-combined/):
##   02-power-grid-co.csv, 02-power-grid-hybrid.csv,
##   02-power-grid-olbdc.csv   E2 power and MCSE, all half-lives
##   02-power-wide.csv         E2 power, t1/2 = 1, 3 x 3 panels
##   02-benchmarks.rds         additivity benchmarks with paired MCSE,
##                             null rejection, equal-total contrasts,
##                             carryover loss, design contrasts
##
## Because the nine weight cells of a replicate share their random
## numbers, contrasts between cells of the same design and half-life
## are paired. Their MCSE comes from a bootstrap over replicates
## (B = 2000). Contrasts across designs or half-lives use independent
## streams and combine the cell MCSEs in quadrature.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

repo_root <- here::here()
args <- commandArgs(trailingOnly = TRUE)
dir_idx <- which(args == '--dir')
out_dir <- if (length(dir_idx) && dir_idx < length(args))
  args[dir_idx + 1] else file.path(repo_root, 'analysis/data/dgp-combined')
read_opt <- function(f) {
  p <- file.path(out_dir, f)
  if (file.exists(p)) readRDS(p) else NULL
}

reps <- bind_rows(readRDS(file.path(out_dir, '01-combined-replicates.rds')),
                  read_opt('01-combined-replicates-hl02.rds'),
                  read_opt('01-combined-replicates-g9.rds'))
stopifnot(nrow(reps) > 0, !anyNA(reps$design))

cell_summary <- reps |>
  group_by(design, t1half, c_bm_a, c_bm_b, method) |>
  summarise(
    n_reps     = n(),
    n_conv     = sum(converged),
    power      = mean(p_value[converged] < 0.05),
    mcse_power = sqrt(power * (1 - power) / max(n_conv, 1L)),
    mean_est   = mean(estimate[converged]),
    mcse_est   = sd(estimate[converged]) / sqrt(n_conv),
    sd_est     = sd(estimate[converged]),
    mean_se    = mean(std_error[converged]),
    kappa      = mean_se / sd_est,
    .groups = 'drop'
  )

## -----------------------------------------------------------------
## CSV artifacts (E2, N = 70)
## -----------------------------------------------------------------

e2 <- dplyr::filter(cell_summary, method == 'E2')
for (dn in c('CO', 'Hybrid', 'OLBDC')) {
  e2 |>
    dplyr::filter(design == dn) |>
    dplyr::transmute(N = 70L, t1half, c_bm_a, c_bm_b, power, mcse_power,
                     mean_est, sd_est) |>
    dplyr::arrange(t1half, c_bm_a, c_bm_b) |>
    write.csv(file.path(out_dir, sprintf('02-power-grid-%s.csv',
                                         tolower(dn))),
              row.names = FALSE)
}
e2 |>
  dplyr::filter(t1half == 1) |>
  dplyr::select(design, c_bm_a, c_bm_b, power) |>
  tidyr::pivot_wider(names_from = c_bm_b, names_prefix = 'cbmb_',
                     values_from = power) |>
  dplyr::mutate(N = 70L, .before = 1) |>
  dplyr::arrange(design, c_bm_a) |>
  write.csv(file.path(out_dir, '02-power-wide.csv'), row.names = FALSE)

## -----------------------------------------------------------------
## Paired contrasts within a design and half-life
## -----------------------------------------------------------------

key <- function(a, b) sprintf('%.3f_%.3f', a, b)
clamp <- function(p, n) pmin(pmax(p, 0.5 / n), 1 - 0.5 / n)

## Rejection matrix: replicates x cells, for one design, half-life
## and method. Non-converged fits are NA and dropped per cell.
reject_matrix <- function(d) {
  d |>
    dplyr::mutate(cell = key(c_bm_a, c_bm_b),
                  rej = ifelse(converged, as.numeric(p_value < 0.05),
                               NA_real_)) |>
    dplyr::select(rep, cell, rej) |>
    tidyr::pivot_wider(names_from = cell, values_from = rej) |>
    dplyr::arrange(rep)
}

contrasts_from_power <- function(p, n) {
  g <- function(a, b) p[[key(a, b)]]
  out <- list()
  z <- qnorm(0.975)
  c_emp <- qnorm(1 - clamp(g(0, 0), n) / 2)
  for (a in c(0.225, 0.45)) for (b in c(0.225, 0.45)) {
    lam <- function(pp, crit) qnorm(clamp(pp, n)) + crit
    out[[sprintf('prob_excess_%s', key(a, b))]] <-
      g(a, b) - (g(a, 0) + g(0, b) - g(0, 0))
    out[[sprintf('probit_excess_%s', key(a, b))]] <-
      g(a, b) - pnorm(lam(g(a, 0), z) + lam(g(0, b), z) - z)
    out[[sprintf('probitcal_excess_%s', key(a, b))]] <-
      g(a, b) - pnorm(lam(g(a, 0), c_emp) + lam(g(0, b), c_emp) - c_emp)
  }
  out$double_A <- g(0.45, 0) - (2 * g(0.225, 0) - g(0, 0))
  out$double_B <- g(0, 0.45) - (2 * g(0, 0.225) - g(0, 0))
  out$probit_double_A <- g(0.45, 0) -
    pnorm(2 * (qnorm(clamp(g(0.225, 0), n)) + z) - z)
  out$probit_double_B <- g(0, 0.45) -
    pnorm(2 * (qnorm(clamp(g(0, 0.225), n)) + z) - z)
  out$eqsum_A_minus_B <- g(0.45, 0) - g(0, 0.45)
  out$eqsum_mix_minus_B <- g(0.225, 0.225) - g(0, 0.45)
  out$eqsum_A_minus_mix <- g(0.45, 0) - g(0.225, 0.225)
  unlist(out)
}

paired_contrasts <- function(d, B = 2000L) {
  m <- reject_matrix(d)
  cols <- setdiff(names(m), 'rep')
  pw <- function(mm) setNames(as.list(colMeans(mm[, cols], na.rm = TRUE)),
                              cols)
  n <- nrow(m)
  est <- contrasts_from_power(pw(m), n)
  boot <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    contrasts_from_power(pw(m[idx, ]), n)
  })
  tibble(contrast = names(est), value = unname(est),
         mcse = apply(boot, 1, sd))
}

set.seed(20261029L)
main_e2 <- reps |> dplyr::filter(method %in% c('E2', 'E2cr2'))
contrasts <- main_e2 |>
  group_by(design, t1half, method) |>
  group_modify(~ paired_contrasts(.x)) |>
  ungroup()

## -----------------------------------------------------------------
## Paired contrasts of the fitted estimate (E2)
## -----------------------------------------------------------------

est_pairs <- reps |>
  dplyr::filter(method == 'E2', converged) |>
  dplyr::mutate(cell = key(c_bm_a, c_bm_b)) |>
  dplyr::select(design, t1half, rep, cell, estimate) |>
  tidyr::pivot_wider(names_from = cell, values_from = estimate) |>
  group_by(design, t1half) |>
  summarise(
    additivity_dev = mean(`0.450_0.450` - `0.450_0.000` - `0.000_0.450` +
                            `0.000_0.000`, na.rm = TRUE),
    additivity_dev_mcse = sd(`0.450_0.450` - `0.450_0.000` -
                               `0.000_0.450` + `0.000_0.000`,
                             na.rm = TRUE) / sqrt(n()),
    A_minus_B = mean(`0.450_0.000` - `0.000_0.450`, na.rm = TRUE),
    A_minus_B_mcse = sd(`0.450_0.000` - `0.000_0.450`, na.rm = TRUE) /
      sqrt(n()),
    .groups = 'drop'
  )

## -----------------------------------------------------------------
## Null rejection, carryover loss, design contrasts
## -----------------------------------------------------------------

null_rej <- cell_summary |>
  dplyr::filter(c_bm_a == 0, c_bm_b == 0) |>
  dplyr::select(design, t1half, method, n_conv, power, mcse_power)

carry_loss <- e2 |>
  dplyr::filter((c_bm_a == 0.45 & c_bm_b == 0) |
                  (c_bm_a == 0 & c_bm_b == 0.45) |
                  (c_bm_a == 0.225 & c_bm_b == 0.225)) |>
  dplyr::select(design, t1half, c_bm_a, c_bm_b, power, mcse_power) |>
  group_by(design, c_bm_a, c_bm_b) |>
  summarise(
    p0 = power[t1half == 0], p1 = power[t1half == 1],
    loss = p0 - p1,
    loss_mcse = sqrt(mcse_power[t1half == 0]^2 +
                       mcse_power[t1half == 1]^2),
    rel_loss = loss / p0,
    .groups = 'drop'
  )

design_contr <- e2 |>
  dplyr::select(design, t1half, c_bm_a, c_bm_b, power, mcse_power) |>
  tidyr::pivot_wider(names_from = design,
                     values_from = c(power, mcse_power)) |>
  dplyr::mutate(
    CO_minus_OLBDC = power_CO - power_OLBDC,
    CO_minus_OLBDC_mcse = sqrt(mcse_power_CO^2 + mcse_power_OLBDC^2),
    Hybrid_minus_OLBDC = power_Hybrid - power_OLBDC,
    Hybrid_minus_OLBDC_mcse = sqrt(mcse_power_Hybrid^2 +
                                     mcse_power_OLBDC^2)
  )

## -----------------------------------------------------------------
## Simulated mean estimate against the analytic GLS projection
## -----------------------------------------------------------------

analytic <- read_opt('03-analytic.rds')
est_vs_analytic <- if (is.null(analytic)) NULL else
  e2 |>
    dplyr::select(design, t1half, c_bm_a, c_bm_b, mean_est, mcse_est,
                  sd_est, power) |>
    dplyr::inner_join(analytic$cells,
                      by = c('design', 't1half', 'c_bm_a', 'c_bm_b')) |>
    dplyr::mutate(gls_power = pnorm(Lambda - qnorm(0.975)) +
                    pnorm(-Lambda - qnorm(0.975)))

## -----------------------------------------------------------------
## G9 against E2cr2 on the same trials (the G9 pass analyzes the first
## replicates of the main pass, under the same seeds)
## -----------------------------------------------------------------

g9_compare <- NULL
if (any(reps$method == 'G9')) {
  g9 <- dplyr::filter(reps, method == 'G9')
  n_g9 <- max(g9$rep)
  g9_compare <- reps |>
    dplyr::filter(method %in% c('E2', 'E2cr2', 'G9'), rep <= n_g9) |>
    dplyr::semi_join(dplyr::distinct(g9, design, t1half, c_bm_a, c_bm_b),
                     by = c('design', 't1half', 'c_bm_a', 'c_bm_b')) |>
    group_by(design, t1half, c_bm_a, c_bm_b, method) |>
    summarise(n_reps = n(), n_conv = sum(converged),
              power = mean(p_value[converged] < 0.05),
              mcse_power = sqrt(power * (1 - power) / max(n_conv, 1L)),
              mean_h = mean(h[converged]),
              share_h_true = mean(abs(h[converged] - t1half) < 1e-9),
              .groups = 'drop')
}

bench <- list(cell_summary = cell_summary, contrasts = contrasts,
              g9_compare = g9_compare,
              est_pairs = est_pairs, null_rej = null_rej,
              carry_loss = carry_loss, design_contr = design_contr,
              est_vs_analytic = est_vs_analytic,
              created = format(Sys.time(), '%Y-%m-%d %H:%M %Z'))
saveRDS(bench, file.path(out_dir, '02-benchmarks.rds'))

## -----------------------------------------------------------------
## Console summary
## -----------------------------------------------------------------

cat('\nNull rejection (c_bm_a = c_bm_b = 0):\n')
print(as.data.frame(null_rej), digits = 3)
cat('\nE2 power panels:\n')
e2 |>
  dplyr::select(design, t1half, c_bm_a, c_bm_b, power) |>
  tidyr::pivot_wider(names_from = c_bm_b, values_from = power) |>
  as.data.frame() |> print(digits = 3)
cat('\nPaired contrasts (E2):\n')
contrasts |> dplyr::filter(method == 'E2') |>
  tidyr::pivot_wider(names_from = contrast,
                     values_from = c(value, mcse)) |>
  as.data.frame() |> print(digits = 3)
cat('\nEstimate additivity and equal-total difference (E2):\n')
print(as.data.frame(est_pairs), digits = 3)
cat('\nCarryover loss, t1/2 0 to 1 (E2):\n')
print(as.data.frame(carry_loss), digits = 3)
if (!is.null(est_vs_analytic)) {
  cat('\nMean estimate against the GLS projection (E2):\n')
  print(as.data.frame(est_vs_analytic |>
    dplyr::select(design, t1half, c_bm_a, c_bm_b, mean_est, theta,
                  sd_est, se_gls, power, gls_power)), digits = 3)
}
invisible(NULL)
