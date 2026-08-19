## analysis/scripts/carryover-sensitivity/20-decay-curve-figure.R
##
## Illustrates the five DGP decay-shape levels used throughout this
## manuscript's Tier 1 grid (exponential; Weibull
## k = 0.25, 0.5, 2.0, 4.0) as closed-form curves of the carryover
## multiplier phi(t_sd) against time since discontinuation, at the
## reference half-life t1half = 1.0 weeks. No simulation: phi_exp/
## phi_weib are the exact formulas from Section 2.3 (### Exponential
## decay, ### Weibull decay). Weibull k = 1.0 is not included: it
## duplicates Exponential exactly (phi_weib(k=1) = phi_exp), so it is
## not a distinct decay-shape condition
## (23-run-decay-shape-sensitivity.R). Writes
## analysis/figures/02xs-decay-curves.pdf.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

repo_root <- here::here()
fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

t1half <- 1.0
t_sd <- seq(0, 4, length.out = 200)

phi_exp <- function(t_sd, t1half) {
  lambda <- log(2) / t1half
  exp(-lambda * t_sd)
}

phi_weib <- function(t_sd, t1half, k) {
  lambda_w <- (log(2))^(1 / k) / t1half
  exp(-(lambda_w * t_sd)^k)
}

curves <- bind_rows(
  tibble(t_sd = t_sd, phi = phi_exp(t_sd, t1half),
        form = 'Exponential'),
  tibble(t_sd = t_sd, phi = phi_weib(t_sd, t1half, 0.25),
        form = 'Weibull k=0.25'),
  tibble(t_sd = t_sd, phi = phi_weib(t_sd, t1half, 0.5),
        form = 'Weibull k=0.5'),
  tibble(t_sd = t_sd, phi = phi_weib(t_sd, t1half, 2.0),
        form = 'Weibull k=2.0'),
  tibble(t_sd = t_sd, phi = phi_weib(t_sd, t1half, 4.0),
        form = 'Weibull k=4.0')
) |>
  mutate(form = factor(form,
                       levels = c('Exponential', 'Weibull k=0.25',
                                 'Weibull k=0.5', 'Weibull k=2.0',
                                 'Weibull k=4.0')))

form_colours <- c('Exponential' = '#1f78b4', 'Weibull k=0.25' = '#b2182b',
                  'Weibull k=0.5' = '#e31a1c', 'Weibull k=2.0' = '#33a02c',
                  'Weibull k=4.0' = '#006d2c')
form_linetypes <- c('Exponential' = 'solid', 'Weibull k=0.25' = 'dashed',
                    'Weibull k=0.5' = 'longdash', 'Weibull k=2.0' = 'dotted',
                    'Weibull k=4.0' = 'dotdash')

p <- ggplot(curves, aes(x = t_sd, y = phi, colour = form, linetype = form)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0.5, linetype = 'dotted', colour = 'grey60',
            linewidth = 0.4) +
  geom_vline(xintercept = t1half, linetype = 'dotted', colour = 'grey60',
            linewidth = 0.4) +
  scale_colour_manual(values = form_colours) +
  scale_linetype_manual(values = form_linetypes) +
  labs(x = expression('Time since discontinuation ('*t[sd]*', weeks)'),
      y = expression(phi(t[sd])),
      colour = 'Decay form', linetype = 'Decay form',
      title = expression('Carryover decay multiplier at'~t['1/2']==1.0~'week')) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'right')

ggsave(file.path(fig_dir, '02xs-decay-curves.pdf'), p, width = 5.4, height = 3.2)

message('Wrote ', file.path(fig_dir, '02xs-decay-curves.pdf'))
