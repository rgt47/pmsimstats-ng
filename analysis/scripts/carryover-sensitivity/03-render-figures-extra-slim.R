## analysis/scripts/carryover-sensitivity/03-render-figures-extra-slim.R
##
## Filtered Tier 1 figures for report.Rmd: covariance-moderation
## architecture only, linear DGP decay dropped, all three analysis
## specifications retained. Display names and ordering come from
## spec-labels.R (stored E1/E3/E2 shown as Unadjusted, Lag-adjusted,
## Exposure-weighted). Reads the same 02-grid-summary.rds as
## 03-render-figures.R and writes to analysis/figures/ with a 02xs-
## prefix, so the full-scope figures used by the archived drafts are
## untouched.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

repo_root <- here::here()

source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))

grid <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary.rds'))$summary |>
  dplyr::filter(dgp_arch == 'mvn', carryover_form != 'linear',
                spec %in% spec_order) |>
  dplyr::mutate(spec = spec_factor(spec))

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_paper <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey92', colour = NA),
    legend.position = 'top'
  )

spec_col <- scale_colour_manual(name = 'Analysis spec',
                                values = spec_colours)
spec_lty <- scale_linetype_manual(name = 'Analysis spec',
                                  values = spec_linetypes)
spec_shp <- scale_shape_manual(name = 'Analysis spec',
                               values = spec_shapes)

## -----------------------------------------------------------------
## Figure 02xs-a: power vs carryover half-life by analysis spec
## Facet: design. Fixed: c_bm = 0.45, N = 70, exponential DGP.
## -----------------------------------------------------------------

d_a <- grid |>
  dplyr::filter(carryover_form == 'exponential', c_bm == 0.45, N == 70) |>
  dplyr::mutate(design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'),
                                labels = c('CO', 'OL+BDC', 'Hybrid')))

p_a <- ggplot(d_a, aes(t1half, power, colour = spec,
                       linetype = spec, shape = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ design, nrow = 1) +
  scale_y_continuous(limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  spec_col + spec_lty + spec_shp +
  labs(
    x = expression('Carryover half-life'~t['1/2']~'(weeks)'),
    y = 'Power',
    title = 'Power vs carryover under matched exponential DGP',
    subtitle = expression(
      'Covariance-moderation architecture,'~N==70*','~c[bm]==0.45*','~
      '500 reps/cell')
  ) +
  theme_paper

ggsave(file.path(fig_dir, '02xs-power-by-spec.pdf'),
  p_a, width = 7, height = 3.4)

## -----------------------------------------------------------------
## Figure 02xs-b: matched-vs-mismatched heatmap
## Cells: DGP decay form x analysis spec; fill = power
## Fixed: design = Hybrid, N = 70, c_bm = 0.45, t1half = 1.0
## Reads the dedicated decay-shape sensitivity data
## (23-run-decay-shape-sensitivity.R, k = 0.25, 0.5, 2.0, 4.0), not
## the shared 02-grid-summary.rds grid (still at the earlier
## k = 0.7, 1.0, 1.5, and retained there for Panel C / Figure 3,
## which is fit separately; see report.Rmd Section 2.6 for why this
## axis is not part of the shared production grid this manuscript
## otherwise reads).
## -----------------------------------------------------------------

d_b <- readRDS(file.path(repo_root,
  'analysis/data/02-decay-shape-sensitivity.rds'))$summary |>
  dplyr::filter(design == 'Hybrid', N == 70, spec %in% spec_order) |>
  dplyr::mutate(
    spec = spec_factor(spec),
    dgp_label = dplyr::case_when(
      carryover_form == 'exponential' ~ 'Exponential',
      carryover_form == 'weibull' & weibull_shape == 0.25 ~ 'Weibull (k=0.25)',
      carryover_form == 'weibull' & weibull_shape == 0.5 ~ 'Weibull (k=0.5)',
      carryover_form == 'weibull' & weibull_shape == 2.0 ~ 'Weibull (k=2.0)',
      carryover_form == 'weibull' & weibull_shape == 4.0 ~ 'Weibull (k=4.0)'
    ),
    dgp_label = factor(dgp_label, levels = c(
      'Exponential', 'Weibull (k=0.25)', 'Weibull (k=0.5)',
      'Weibull (k=2.0)', 'Weibull (k=4.0)'))
  )

## The twelve cells fall in a narrow high-power band, so a fill scale
## fixed to the full [0, 1] range renders them as one dark hue with no
## visible contrast. Stretch the sequential scale to the observed
## range instead, with a small pad.
pad <- diff(range(d_b$power)) * 0.12
fill_lim <- range(d_b$power) + c(-pad, pad)

p_b <- ggplot(d_b, aes(spec, dgp_label, fill = power)) +
  geom_tile(colour = 'white') +
  geom_text(aes(label = sprintf('%.2f', power),
                colour = power < mean(fill_lim)), size = 3,
            show.legend = FALSE) +
  scale_colour_manual(values = c('TRUE' = 'grey15', 'FALSE' = 'white')) +
  scale_fill_gradient(low = '#eff3ff', high = '#08306b',
    limits = fill_lim, breaks = scales::pretty_breaks(4),
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_colorbar(barwidth = grid::unit(4.2, 'cm'),
                            barheight = grid::unit(0.35, 'cm'))) +
  labs(
    x = 'Analysis specification',
    y = 'DGP decay form',
    fill = 'Power',
    title = 'Decay-form x analysis-spec sensitivity',
    subtitle = expression(
      'Covariance-moderation architecture, Hybrid design,'~N==70*','~
      c[bm]==0.45*','~t['1/2']==1.0)
  ) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 0, size = 8))

ggsave(file.path(fig_dir, '02xs-heatmap-matched-vs-mismatched.pdf'),
  p_b, width = 6.8, height = 4.6)

## -----------------------------------------------------------------
## Figure 02xs-c: type-I error check
## Boxplot of power (= type I error) at c_bm = 0 across spec/design
## -----------------------------------------------------------------

d_c <- grid |>
  dplyr::filter(c_bm == 0, N == 70) |>
  dplyr::mutate(design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'),
                                labels = c('CO', 'OL+BDC', 'Hybrid')))

p_c <- ggplot(d_c, aes(spec, power, fill = spec)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  facet_wrap(~ design) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = spec_colours) +
  labs(
    x = 'Analysis specification',
    y = 'Empirical type-I error rate',
    title = 'Type-I error control under null',
    subtitle = expression(
      'Covariance-moderation architecture,'~c[bm]==0*','~N==70*','~
      'pooled across DGP decay forms and carryover levels')
  ) +
  theme_paper +
  theme(axis.text.x = element_text(size = 7, angle = 20, hjust = 1))

ggsave(file.path(fig_dir, '02xs-type1-boxplots.pdf'),
  p_c, width = 6.8, height = 3.4)

message('Wrote three extra-slim figures to ', fig_dir, ':')
message('  02xs-power-by-spec.pdf')
message('  02xs-heatmap-matched-vs-mismatched.pdf')
message('  02xs-type1-boxplots.pdf')
