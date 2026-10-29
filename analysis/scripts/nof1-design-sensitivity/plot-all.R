#' plot-all.R
#'
#' Render one figure per sweep from the RDS artifacts. Run after
#' run-all.R. Output: analysis/figures/sensitivity/sens-NN-name.pdf.
#' Every point carries a 95% Monte Carlo interval (1.96 MCSE). Sweeps
#' whose artifact lacks the schema-2 metadata (the 2026-05 runs) are
#' skipped.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
})

source(here::here(
  'analysis', 'scripts', 'nof1-design-sensitivity', '00-summarize.R'
))

sens_data_dir <- Sys.getenv(
  'SENS_OUT_DIR',
  here::here('analysis', 'data', 'derived_data', 'sensitivity')
)
sens_fig_dir <- Sys.getenv(
  'SENS_FIG_DIR',
  here::here('analysis', 'figures', 'sensitivity')
)
if (!dir.exists(sens_fig_dir)) {
  dir.create(sens_fig_dir, recursive = TRUE)
}

with_ci <- function(sm, y = 'power', se = 'mcse_power') {
  sm$lo <- sm[[y]] - 1.96 * sm[[se]]
  sm$hi <- sm[[y]] + 1.96 * sm[[se]]
  sm
}

#' Power (or Type I error) against a numeric axis, one line per design.
line_plot <- function(sm, x, color = 'design_name', title, xlab,
                      ylab = 'Rejection rate', facet = NULL,
                      hline = c(0.05, 0.8)) {
  sm <- with_ci(sm)
  p <- ggplot(sm, aes(x = .data[[x]], y = power,
                      color = .data[[color]],
                      group = .data[[color]])) +
    geom_hline(yintercept = hline[1], linetype = 'dotted',
               linewidth = 0.3) +
    geom_hline(yintercept = hline[2], linetype = 'dashed',
               linewidth = 0.3) +
    geom_line() +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3) +
    labs(title = title, x = xlab, y = ylab, color = 'Design') +
    theme_minimal()
  if (!is.null(facet)) p <- p + facet_wrap(facet, labeller = label_both)
  p
}

load_art <- function(id, name) {
  path <- file.path(sens_data_dir,
                    sprintf('sensitivity-%02d-%s.rds', id, name))
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!identical(obj$run_meta$schema, 2L)) {
    message('Skipping (pre-schema-2 artifact): ', path)
    return(NULL)
  }
  obj
}

plot_s1 <- function(r) {
  sm <- summarize_cells(r, c('design_name', 'delta_on'))
  line_plot(sm, 'delta_on', title = 'S1. Rejection rate against the effect',
            xlab = 'Mean drug response on drug (outcome units)')
}

plot_s2 <- function(r) {
  r$halflife_days <- round(r$carryover_t1half * 7, 4)
  sm <- summarize_cells(r, c('design_name', 'halflife_days',
                             'analysis_t1half'))
  line_plot(sm, 'halflife_days',
            title = 'S2. Power against carryover half-life',
            xlab = 'Data-generating half-life (days)',
            facet = 'analysis_t1half')
}

plot_s3 <- function(r) {
  r$decay_form <- ifelse(r$carryover_form == 'weibull',
                         sprintf('weibull_k=%g', r$weibull_shape),
                         r$carryover_form)
  sm <- with_ci(summarize_cells(r, c('design_name', 'decay_form',
                                     'analysis_t1half')))
  ggplot(sm, aes(x = decay_form, y = power, color = design_name)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3,
                    position = position_dodge(width = 0.5)) +
    facet_wrap('analysis_t1half', labeller = label_both) +
    labs(title = 'S3. Power by data-generating decay form',
         x = 'Decay form', y = 'Power', color = 'Design') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_s4 <- function(r) {
  sm <- summarize_cells(r, c('axis', 'design_name', 'tau', 'bl_sd'))
  a <- line_plot(filter(sm, axis == 'tau'), 'tau',
                 title = 'S4a. Power against treatment heterogeneity',
                 xlab = 'Heterogeneity SD (tau, outcome units)')
  a
}

plot_s4b <- function(r) {
  sm <- summarize_cells(r, c('axis', 'design_name', 'tau', 'bl_sd'))
  line_plot(filter(sm, axis == 'bl_sd'), 'bl_sd',
            title = 'S4b. Power against baseline SD',
            xlab = 'Baseline SD (outcome units)')
}

plot_s5 <- function(r) {
  r$sd_scale <- round(r$br_sd / 1.0, 4)
  sm <- summarize_cells(r, c('design_name', 'sd_scale'))
  line_plot(sm, 'sd_scale',
            title = 'S5. Power against within-participant SD',
            xlab = 'Multiplier on all component SDs')
}

plot_s6 <- function(r) {
  r$periods <- r$n_switch + 1
  sm <- summarize_cells(r, c('design_name', 'periods'))
  line_plot(sm, 'periods',
            title = 'S6. Power against number of periods',
            xlab = 'Periods in 16 weekly visits (RCT: 1)')
}

plot_s7 <- function(r) {
  sm <- summarize_cells(r, c('design_name', 'N'))
  line_plot(sm, 'N', title = 'S7. Power against sample size',
            xlab = 'Total participants N')
}

plot_s8 <- function(r) {
  sm <- summarize_cells(r, c('design_name', 'c.br'))
  line_plot(sm, 'c.br', title = 'S8. Power against AR(1) correlation',
            xlab = 'AR(1) correlation per week (rho)')
}

plot_s9 <- function(r) {
  r$period_days <- as.numeric(sub('^Alt-(\\d+)d$', '\\1', r$design_name))
  r$halflife_days <- round(r$carryover_t1half * 7, 4)
  sm <- summarize_cells(r, c('period_days', 'halflife_days'))
  ggplot(sm, aes(x = factor(period_days), y = factor(halflife_days),
                 fill = power)) +
    geom_tile(color = 'white') +
    geom_text(aes(label = sprintf('%.2f', power)), size = 3) +
    scale_fill_distiller(palette = 'RdYlBu', direction = -1,
                         limits = c(0, 1)) +
    labs(title = 'S9. Power of the alternating design',
         x = 'Period length (days)', y = 'Carryover half-life (days)',
         fill = 'Power') +
    theme_minimal()
}

plot_s10 <- function(r) {
  r$halflife_days <- round(r$carryover_t1half * 7, 4)
  sm <- summarize_cells(r, c('design_name', 'dgp_architecture', 'c.bm',
                             'halflife_days'))
  line_plot(sm, 'c.bm',
            title = 'S10. Main-effect power against biomarker correlation',
            xlab = 'c_bm', facet = c('dgp_architecture', 'halflife_days'))
}

plot_s11 <- function(r) {
  sm <- with_ci(summarize_cells(
    r, c('design_name', 'delta_label', 'carryover', 'ar1',
         'analysis_t1half')))
  sm$condition <- paste0('carryover ', sm$carryover, ', AR(1) ', sm$ar1)
  ggplot(sm, aes(x = condition, y = power, color = design_name,
                 shape = factor(analysis_t1half))) +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3,
                    position = position_dodge(width = 0.6)) +
    geom_hline(yintercept = 0.05, linetype = 'dotted') +
    facet_wrap('delta_label', scales = 'free_y') +
    labs(title = 'S11. Misspecification factorial', x = NULL,
         y = 'Rejection rate', color = 'Design',
         shape = 'Analysis half-life (wk)') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

plot_s12 <- function(r) {
  sm <- with_ci(summarize_cells(
    r, c('design_name', 'config', 'analysis_t1half')))
  ggplot(sm, aes(x = factor(config), y = power, color = design_name,
                 shape = factor(analysis_t1half))) +
    geom_hline(yintercept = 0.05, linetype = 'dashed') +
    geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.3,
                    position = position_dodge(width = 0.6)) +
    labs(title = 'S12. Type I error with 95% Monte Carlo interval',
         x = 'Configuration', y = 'Type I error', color = 'Design',
         shape = 'Analysis half-life (wk)') +
    theme_minimal()
}

sweeps <- list(
  list(1, 'effect-size', plot_s1),
  list(2, 'carryover-halflife', plot_s2),
  list(3, 'decay-form', plot_s3),
  list(4, 'tau2', plot_s4),
  list(4, 'tau2', plot_s4b, 'b'),
  list(5, 'sigma2', plot_s5),
  list(6, 'cycles', plot_s6),
  list(7, 'patients', plot_s7),
  list(8, 'ar1', plot_s8),
  list(9, 'period', plot_s9),
  list(10, 'biomarker', plot_s10),
  list(11, 'misspec', plot_s11),
  list(12, 'morris-null', plot_s12)
)

for (sw in sweeps) {
  art <- load_art(sw[[1]], sw[[2]])
  if (is.null(art)) next
  suffix <- if (length(sw) > 3) sw[[4]] else ''
  out <- file.path(sens_fig_dir,
                   sprintf('sens-%02d%s-%s.pdf', sw[[1]], suffix, sw[[2]]))
  ggsave(out, sw[[3]](art$results), width = 7, height = 5,
         device = 'pdf')
  message('Rendered: ', out)
}
message('plot-all.R complete.')
