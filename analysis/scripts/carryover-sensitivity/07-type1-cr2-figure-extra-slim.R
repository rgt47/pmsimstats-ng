## analysis/scripts/carryover-sensitivity/07-type1-cr2-figure-extra-slim.R
##
## Filtered S6 type-I figure for report-extra-slim.Rmd: drops the A1
## analysis specification (S6 is already Architecture B only). Reads
## the same diag-s6-cr2.rds as 07-type1-cr2-figure.R and writes
## analysis/figures/02xs-type1-cr2.pdf, leaving the full-size figure
## untouched.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

repo_root <- here::here()
fig_dir   <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
out_path  <- file.path(fig_dir, '02xs-type1-cr2.pdf')

s6_path <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/diag-s6-cr2.rds')

placeholder <- function() {
  p <- ggplot() +
    annotate('text', x = 0.5, y = 0.5,
             label = 'Pending S6 CR2 type-I run.', size = 4) +
    xlim(0, 1) + ylim(0, 1) + theme_void()
  ggsave(out_path, p, width = 6, height = 3)
}

if (!file.exists(s6_path)) { placeholder(); quit(save = 'no') }
s6 <- readRDS(s6_path)
if (is.null(s6$type1)) { placeholder(); quit(save = 'no') }

n <- s6$n_reps
## binomial 95% Monte Carlo band for nominal 0.05
band_lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / n)
band_hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / n)

cell_levels <- c('Reference (N=70)', 'Small N (N=35)',
                 'High autocorr (rho=0.9)',
                 'Half-life mis-spec (assume 0.25)', '30% MCAR dropout')
cell_short  <- c('Reference', 'Small N', 'High rho',
                 'HL mis-spec', '30% dropout')

d <- s6$type1 |>
  dplyr::filter(spec %in% c('A2', 'A3')) |>
  mutate(
    cell = factor(cell, levels = cell_levels, labels = cell_short),
    spec = factor(spec, levels = c('A2', 'A3'),
                  labels = c('A2 Dbc (matched)', 'A3 lagged'))) |>
  pivot_longer(c(type1_mod, type1_cr2),
               names_to = 'inference', values_to = 'type1') |>
  mutate(inference = factor(inference,
                            levels = c('type1_mod', 'type1_cr2'),
                            labels = c('Model-based',
                                       'CR2 cluster-robust')))

dodge <- position_dodge(width = 0.55)

p <- ggplot(d, aes(cell, type1, colour = inference, shape = inference,
                   group = inference)) +
  annotate('rect', xmin = -Inf, xmax = Inf,
           ymin = band_lo, ymax = band_hi, fill = 'grey88') +
  geom_hline(yintercept = 0.05, linetype = 'dashed', colour = 'grey40') +
  ## Model-based is an open circle (shape 1); against the light grey
  ## MC band a thin light-grey outline nearly disappears. Use a dark
  ## outline colour, a thicker stroke, and a slightly larger size so
  ## the open marker reads clearly at the same visual weight as the
  ## filled CR2 marker.
  geom_point(position = dodge, size = 2.6, stroke = 1.4) +
  facet_wrap(~ spec) +
  scale_colour_manual(name = 'Inference',
    values = c('Model-based' = '#4d4d4d',
               'CR2 cluster-robust' = '#0072B2')) +
  scale_shape_manual(name = 'Inference',
    values = c('Model-based' = 1, 'CR2 cluster-robust' = 16)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = 'Empirical type-I error',
       title = paste0('S6: type-I error at the null (c_bm = 0), ',
                      'model-based vs CR2'),
       subtitle = paste0('Dashed line = nominal 0.05; grey band = ',
                         'binomial 95% MC interval at ', n, ' reps')) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'grey92', colour = NA),
        axis.text.x = element_text(size = 7, angle = 35, hjust = 1),
        legend.position = 'top')

ggsave(out_path, p, width = 6, height = 3.6)
message('Wrote ', out_path)
