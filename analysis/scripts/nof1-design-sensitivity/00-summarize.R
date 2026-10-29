#' 00-summarize.R
#'
#' Performance measures per cell, with Monte Carlo standard errors
#' (Morris, White and Crowther 2019). Sourced by plot-all.R, by the
#' check scripts and by report.Rmd, so that the tables, the figures
#' and the prose all use one definition.
#'
#' Definitions, per cell (group of replicate rows):
#'   power       share of replicates with p < alpha (Type I error when
#'               the true effect is zero)
#'   target      ETbeta, the coefficient the analysis estimates in a
#'               pooled fit of every replicate's participants (the
#'               large-sample value of the estimator)
#'   bias        mean(beta) - target
#'   bias_estimand  mean(beta) minus delta_on, the design-independent
#'               estimand (the mean drug response on drug); it is the
#'               zero-mean check at the null and equals bias there
#'   emp_se      SD of beta across replicates
#'   model_se    mean model-based SE
#'   kappa       model_se / emp_se (1 is calibrated)
#'   coverage    share of replicates whose 95% interval around beta
#'               contains target (t interval when a df is stored)

suppressPackageStartupMessages({
  library(dplyr)
})

summarize_cells <- function(df, by, alpha = 0.05) {
  df |>
    mutate(
      crit = ifelse(is.na(.data$betaDF), stats::qnorm(0.975),
                    stats::qt(0.975, .data$betaDF))
    ) |>
    group_by(across(all_of(by))) |>
    summarise(
      n_rep = n(),
      n_fit = sum(!is.na(.data$p)),
      nonconv = mean(is.na(.data$p)),
      power = mean(.data$p < alpha, na.rm = TRUE),
      mean_beta = mean(.data$beta, na.rm = TRUE),
      target = mean(.data$ETbeta, na.rm = TRUE),
      emp_se = stats::sd(.data$beta, na.rm = TRUE),
      model_se = mean(.data$betaSE, na.rm = TRUE),
      coverage = mean(
        abs(.data$beta - .data$ETbeta) <= .data$crit * .data$betaSE,
        na.rm = TRUE),
      delta_on_cell = mean(.data$delta_on, na.rm = TRUE),
      singular = mean(.data$issingular, na.rm = TRUE),
      .groups = 'drop'
    ) |>
    mutate(
      mcse_power = sqrt(.data$power * (1 - .data$power) / .data$n_fit),
      bias = .data$mean_beta - .data$target,
      mcse_bias = .data$emp_se / sqrt(.data$n_fit),
      mcse_emp_se = .data$emp_se / sqrt(2 * (.data$n_fit - 1)),
      bias_estimand = .data$mean_beta - .data$delta_on_cell,
      kappa = .data$model_se / .data$emp_se,
      mcse_coverage = sqrt(.data$coverage * (1 - .data$coverage) /
                             .data$n_fit)
    )
}

#' Normal-approximation 95% interval for a rate, from its MCSE.
rate_ci <- function(rate, mcse) {
  c(lower = rate - stats::qnorm(0.975) * mcse,
    upper = rate + stats::qnorm(0.975) * mcse)
}

#' TRUE when a rate is within two MCSE of a nominal value, with the
#' MCSE computed at the nominal value (the band used for 'nominal').
within_nominal <- function(rate, n, nominal = 0.05) {
  abs(rate - nominal) <= 2 * sqrt(nominal * (1 - nominal) / n)
}
