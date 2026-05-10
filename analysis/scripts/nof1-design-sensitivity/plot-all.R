#' plot-all.R
#'
#' Render per-sweep figures from the RDS artifacts. Intended to
#' run after run-all.R completes. Produces one PDF per sweep
#' under analysis/figures/sensitivity/ and returns a named list
#' of ggplot objects for downstream assembly.

suppressPackageStartupMessages({
  library(nof1power)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})

sens_data_dir <- file.path(
  "analysis", "data", "derived_data", "sensitivity"
)
sens_fig_dir <- file.path("analysis", "figures", "sensitivity")
if (!dir.exists(sens_fig_dir)) {
  dir.create(sens_fig_dir, recursive = TRUE)
}

#' Generic Morris metric helper (mirrors report.Rmd).
morris <- function(df, true_beta = NA_real_, alpha = 0.05) {
  n <- sum(!is.na(df$p))
  power <- mean(df$p < alpha, na.rm = TRUE)
  mcse_power <- sqrt(power * (1 - power) / max(n, 1))
  emp_se <- stats::sd(df$beta, na.rm = TRUE)
  tibble(
    n_sim = n, power = power, mcse_power = mcse_power,
    mean_beta = mean(df$beta, na.rm = TRUE),
    empirical_SE = emp_se,
    frac_singular = mean(df$issingular, na.rm = TRUE)
  )
}

#=====================================================================
# Dispatcher
#=====================================================================

plot_sweep <- function(id, name) {
  path <- file.path(
    sens_data_dir,
    sprintf("sensitivity-%02d-%s.rds", id, name)
  )
  if (!file.exists(path)) {
    message("Skipping (missing artifact): ", path)
    return(invisible(NULL))
  }
  obj <- readRDS(path)
  res <- if (is.list(obj) && "results" %in% names(obj))
    obj$results else obj

  plotter <- switch(
    sprintf("%02d", id),
    "01" = plot_s1, "02" = plot_s2, "03" = plot_s3,
    "04" = plot_s4, "05" = plot_s5, "06" = plot_s6,
    "07" = plot_s7, "08" = plot_s8, "09" = plot_s9,
    "10" = plot_s10, "11" = plot_s11, "12" = plot_s12
  )
  if (is.null(plotter)) return(invisible(NULL))

  p <- plotter(res)
  out <- file.path(
    sens_fig_dir,
    sprintf("sens-%02d-%s.pdf", id, name)
  )
  ggsave(out, p, width = 7, height = 5, device = "pdf")
  message("Rendered: ", out)
  invisible(p)
}

#=====================================================================
# Per-sweep plotters
#=====================================================================

plot_s1 <- function(res) {
  agg <- res |>
    group_by(design_name, delta) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = delta, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    labs(title = "S1. Power vs effect size",
         x = "True delta (nightmares/week)", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s2 <- function(res) {
  agg <- res |>
    group_by(design_name, halflife_days) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = halflife_days, y = power,
                  color = design_name, group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S2. Power vs carryover half-life",
         x = "Half-life (days)", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s3 <- function(res) {
  agg <- res |>
    group_by(design_name, decay_form) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = decay_form, y = power,
                  fill = design_name)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7) +
    labs(title = "S3. Power by DGP decay form (analysis = exponential)",
         x = "DGP form", y = "Power", fill = "Design") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_s4 <- function(res) {
  agg <- res |>
    group_by(design_name, tau) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = tau, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S4. Power vs between-patient SD",
         x = "Baseline (BL) SD (tau)", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s5 <- function(res) {
  agg <- res |>
    group_by(design_name, sigma) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = sigma, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S5. Power vs within-patient SD",
         x = "Bio-response (BR) SD (sigma)", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s6 <- function(res) {
  agg <- res |>
    group_by(design_name, k) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = k, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S6. Power vs timepoints per patient (k)",
         x = "Cycles per patient k", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s7 <- function(res) {
  agg <- res |>
    group_by(design_name, n_total) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = n_total, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S7. Power vs total participants",
         x = "N (total)", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s8 <- function(res) {
  agg <- res |>
    group_by(design_name, rho) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = rho, y = power, color = design_name,
                  group = design_name)) +
    geom_line() + geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S8. Power vs AR(1) rho (analysis uses corCAR1)",
         x = "AR(1) rho", y = "Power",
         color = "Design") +
    theme_minimal()
}

plot_s9 <- function(res) {
  agg <- res |>
    group_by(period_days, halflife_days) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = factor(period_days),
                  y = factor(halflife_days), fill = power)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", power)), size = 3) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1,
                         limits = c(0, 1)) +
    labs(title = "S9. Power heatmap (N-of-1) vs period and half-life",
         x = "Period length (days)",
         y = "Carryover half-life (days)",
         fill = "Power") +
    theme_minimal()
}

plot_s10 <- function(res) {
  agg <- res |>
    group_by(design_name, dgp_architecture, c_bm,
             halflife_days) |>
    group_modify(~ morris(.x))
  ggplot(agg, aes(x = c_bm, y = power,
                  color = design_name, linetype = factor(halflife_days),
                  group = interaction(design_name, halflife_days))) +
    geom_line() + geom_point(size = 2) +
    facet_wrap(~ dgp_architecture) +
    geom_hline(yintercept = 0.8, linetype = "dashed") +
    labs(title = "S10. Biomarker-interaction power",
         x = "c.bm", y = "Power",
         color = "Design", linetype = "Half-life (days)") +
    theme_minimal()
}

plot_s11 <- function(res) {
  agg <- res |>
    group_by(design_name, delta_label, carryover, ar1) |>
    group_modify(~ morris(.x))
  agg$condition <- paste0("carryover=", agg$carryover,
                          ", ar1=", agg$ar1)
  ggplot(agg,
         aes(x = condition, y = power, fill = design_name)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7) +
    facet_wrap(~ delta_label,
               labeller = labeller(
                 delta_label = c(null = "Null (Type I)",
                                 alt = "Alternative (Power)"))) +
    geom_hline(yintercept = 0.05, linetype = "dotted") +
    labs(title = "S11. Misspecification factorial",
         x = "DGP condition", y = "Rejection rate",
         fill = "Design") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_s12 <- function(res) {
  agg <- res |>
    group_by(design_name, config, tau, sigma, halflife_days) |>
    group_modify(~ morris(.x, true_beta = 0))
  agg$ci_lo <- agg$power - 1.96 * agg$mcse_power
  agg$ci_hi <- agg$power + 1.96 * agg$mcse_power
  agg$config_label <- sprintf(
    "cfg%d (tau=%g, sigma=%g, hl=%gd)",
    agg$config, agg$tau, agg$sigma, agg$halflife_days
  )
  ggplot(agg,
         aes(x = config_label, y = power, color = design_name,
             group = design_name)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.2,
                  position = position_dodge(width = 0.4)) +
    geom_point(size = 2.5,
               position = position_dodge(width = 0.4)) +
    labs(title = "S12. Null-effect Type I error with 95% CI",
         x = "Configuration", y = "Empirical Type I",
         color = "Design") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

#=====================================================================
# Main
#=====================================================================

if (sys.nframe() == 0L) {
  sweeps <- list(
    c(1, "effect-size"),
    c(2, "carryover-halflife"),
    c(3, "decay-form"),
    c(4, "tau2"),
    c(5, "sigma2"),
    c(6, "cycles"),
    c(7, "patients"),
    c(8, "ar1"),
    c(9, "period"),
    c(10, "biomarker"),
    c(11, "misspec"),
    c(12, "morris-null")
  )

  for (sw in sweeps) {
    plot_sweep(as.integer(sw[1]), sw[2])
  }
  message("plot-all.R complete.")
}
