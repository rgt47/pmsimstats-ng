## analysis/scripts/carryover-sensitivity/03-render-figures-extra-slim.R
##
## Filtered Tier 1 figures for report-extra-slim.Rmd: Architecture B
## only, linear DGP decay dropped, A1 analysis specification dropped.
## Reads the same 02-grid-summary.rds as 03-render-figures.R and
## writes to analysis/figures/ with a 02xs- prefix so the full-size
## report.Rmd and report-slim.Rmd figures are untouched.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

repo_root <- here::here()

grid <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary.rds'))$summary |>
  dplyr::filter(dgp_arch == 'mvn', carryover_form != 'linear',
                spec %in% c('A2', 'A3'))

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_paper <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey92', colour = NA),
    legend.position = 'top'
  )

spec_labels <- c(A2 = 'A2 Dbc (matched)', A3 = 'A3 lagged')

## -----------------------------------------------------------------
## Figure 02xs-a: power vs carryover half-life by analysis spec
## Facet: design. Fixed: c_bm = 0.45, N = 70, exponential DGP.
## -----------------------------------------------------------------

d_a <- grid |>
  dplyr::filter(carryover_form == 'exponential', c_bm == 0.45, N == 70) |>
  dplyr::mutate(
    spec = factor(spec, levels = c('A2', 'A3'), labels = spec_labels),
    design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'))
  )

p_a <- ggplot(d_a, aes(t1half, power, colour = spec,
                       linetype = spec, shape = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_wrap(~ design, nrow = 1) +
  scale_y_continuous(limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_linetype_manual(name = 'Analysis spec',
    values = c('solid', 'longdash')) +
  scale_shape_manual(name = 'Analysis spec', values = c(17, 4)) +
  labs(
    x = expression('Carryover half-life'~t['1/2']~'(weeks)'),
    y = 'Power',
    colour = 'Analysis spec',
    title = 'Power vs carryover under matched exponential DGP',
    subtitle = expression(
      'Architecture B,'~N==70*','~c[bm]==0.45*','~'500 reps/cell')
  ) +
  theme_paper

ggsave(file.path(fig_dir, '02xs-power-by-spec.pdf'),
  p_a, width = 7, height = 3.2)

## -----------------------------------------------------------------
## Figure 02xs-b: matched-vs-mismatched heatmap
## Cells: DGP decay form x analysis spec; fill = power
## Fixed: design = Hybrid, N = 70, c_bm = 0.45, t1half = 1.0
## -----------------------------------------------------------------

d_b <- grid |>
  dplyr::filter(design == 'Hybrid', c_bm == 0.45, N == 70, t1half == 1.0) |>
  dplyr::mutate(
    dgp_label = dplyr::case_when(
      carryover_form == 'exponential' ~ 'Exponential',
      carryover_form == 'weibull' & weibull_shape == 0.7 ~ 'Weibull (k=0.7)',
      carryover_form == 'weibull' & weibull_shape == 1.0 ~ 'Weibull (k=1.0)',
      carryover_form == 'weibull' & weibull_shape == 1.5 ~ 'Weibull (k=1.5)'
    ),
    dgp_label = factor(dgp_label, levels = c(
      'Exponential', 'Weibull (k=0.7)', 'Weibull (k=1.0)', 'Weibull (k=1.5)')),
    spec = factor(spec, levels = c('A2', 'A3'), labels = spec_labels)
  )

## All eight cells fall in a narrow high-power band (~0.66-0.87), so
## a fill scale fixed to the full [0, 1] power range renders them as
## a single dark hue with almost no visible contrast. Stretch the
## sequential (single-hue) scale to the observed range instead, with
## a small pad, so the modest cell-to-cell differences that matter
## for this figure's claim are visible.
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
      'Architecture B, Hybrid design,'~N==70*','~
      c[bm]==0.45*','~t['1/2']==1.0)
  ) +
  theme_paper +
  theme(axis.text.x = element_text(angle = 0))

ggsave(file.path(fig_dir, '02xs-heatmap-matched-vs-mismatched.pdf'),
  p_b, width = 6, height = 3.6)

## -----------------------------------------------------------------
## Figure 02xs-c: type-I error check
## Boxplot of power (= type I error) at c_bm = 0 across spec/design
## -----------------------------------------------------------------

d_c <- grid |>
  dplyr::filter(c_bm == 0, N == 70) |>
  dplyr::mutate(
    spec = factor(spec, levels = c('A2', 'A3'), labels = spec_labels),
    design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'))
  )

p_c <- ggplot(d_c, aes(spec, power, fill = spec)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  facet_wrap(~ design) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_fill_manual(values = c('A2 Dbc (matched)' = '#33a02c',
                                'A3 lagged' = '#e31a1c')) +
  labs(
    x = 'Analysis specification',
    y = 'Empirical type-I error rate',
    title = 'Type-I error control under null',
    subtitle = expression(
      'Architecture B,'~c[bm]==0*','~N==70*','~
      'pooled across DGP decay forms and carryover levels')
  ) +
  theme_paper

ggsave(file.path(fig_dir, '02xs-type1-boxplots.pdf'),
  p_c, width = 6, height = 3.2)

message('Wrote three extra-slim figures to ', fig_dir, ':')
message('  02xs-power-by-spec.pdf')
message('  02xs-heatmap-matched-vs-mismatched.pdf')
message('  02xs-type1-boxplots.pdf')
