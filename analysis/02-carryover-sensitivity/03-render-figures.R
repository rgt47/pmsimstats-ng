## analysis/02-carryover-sensitivity/03-render-figures.R
##
## Render manuscript-02 figures from the summarised grid. Outputs
## PDF files into manuscripts/figures/.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

repo_root <- here::here()

grid <- readRDS(file.path(repo_root,
  'manuscripts/data/02-grid-summary.rds'))$summary

fig_dir <- file.path(repo_root, 'manuscripts/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_paper <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey92', colour = NA),
    legend.position = 'top'
  )

## -----------------------------------------------------------------
## Figure 02a: power vs carryover half-life by analysis spec
## Facets: design x architecture
## Fixed: c_bm = 0.45, N = 70, exponential DGP (weibull shape = 1)
## -----------------------------------------------------------------

d_a <- grid |>
  dplyr::filter(
    carryover_form == 'exponential',
    c_bm == 0.45, N == 70
  ) |>
  dplyr::mutate(
    dgp_arch = factor(dgp_arch,
      levels = c('mean_moderation', 'mvn'),
      labels = c('Arch A (mean mod)', 'Arch B (MVN)')),
    spec = factor(spec, levels = c('A1', 'A2', 'A3'),
      labels = c('A1 binary', 'A2 Dbc (matched)', 'A3 lagged')),
    design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'))
  )

p_a <- ggplot(d_a, aes(t1half, power, colour = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  facet_grid(dgp_arch ~ design) +
  scale_y_continuous(limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    x = expression('Carryover half-life'~t['1/2']~'(weeks)'),
    y = 'Power',
    colour = 'Analysis spec',
    title = 'Power vs carryover under matched exponential DGP',
    subtitle = expression(
      N==70*','~c[bm]==0.45*','~'50 reps/cell (dev)')
  ) +
  theme_paper

ggsave(file.path(fig_dir, '02-power-by-spec.pdf'),
  p_a, width = 7, height = 4.5)

## -----------------------------------------------------------------
## Figure 02b: matched-vs-mismatched heatmap
## Cells: DGP decay form x analysis spec; fill = power
## Fixed: dgp_arch = mvn, design = Hybrid, N = 70, c_bm = 0.45,
##        t1half = 1.0
## -----------------------------------------------------------------

d_b <- grid |>
  dplyr::filter(
    dgp_arch == 'mvn', design == 'Hybrid',
    c_bm == 0.45, N == 70, t1half == 1.0
  ) |>
  dplyr::mutate(
    dgp_label = dplyr::case_when(
      carryover_form == 'linear'      ~ 'Linear',
      carryover_form == 'exponential' ~ 'Exponential',
      carryover_form == 'weibull' & weibull_shape == 0.7 ~ 'Weibull (k=0.7)',
      carryover_form == 'weibull' & weibull_shape == 1.0 ~ 'Weibull (k=1.0)',
      carryover_form == 'weibull' & weibull_shape == 1.5 ~ 'Weibull (k=1.5)'
    ),
    dgp_label = factor(dgp_label, levels = c(
      'Linear', 'Exponential',
      'Weibull (k=0.7)', 'Weibull (k=1.0)', 'Weibull (k=1.5)')),
    spec = factor(spec, levels = c('A1', 'A2', 'A3'),
      labels = c('A1 binary', 'A2 Dbc', 'A3 lagged'))
  )

p_b <- ggplot(d_b, aes(spec, dgp_label, fill = power)) +
  geom_tile(colour = 'white') +
  geom_text(aes(label = sprintf('%.2f', power)), size = 3) +
  scale_fill_gradient2(low = 'white', mid = '#fdd',
    high = '#08519c', midpoint = 0.4, limits = c(0, 1)) +
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

ggsave(file.path(fig_dir, '02-heatmap-matched-vs-mismatched.pdf'),
  p_b, width = 6, height = 4)

## -----------------------------------------------------------------
## Figure 02c: type-I error check
## Bar of power (= type I error) at c_bm = 0 across all specs/designs
## -----------------------------------------------------------------

d_c <- grid |>
  dplyr::filter(c_bm == 0, N == 70) |>
  dplyr::mutate(
    spec = factor(spec, levels = c('A1', 'A2', 'A3')),
    design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid')),
    dgp_arch = factor(dgp_arch,
      levels = c('mean_moderation', 'mvn'),
      labels = c('Arch A', 'Arch B'))
  )

p_c <- ggplot(d_c, aes(spec, power, fill = dgp_arch)) +
  geom_boxplot(outlier.size = 0.6) +
  geom_hline(yintercept = 0.05, linetype = 'dashed') +
  facet_wrap(~ design) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(
    x = 'Analysis specification',
    y = 'Empirical type-I error rate',
    fill = NULL,
    title = 'Type-I error control under null',
    subtitle = expression(
      c[bm]==0*','~N==70*','~
      'pooled across DGP decay forms and carryover levels')
  ) +
  theme_paper

ggsave(file.path(fig_dir, '02-type1-boxplots.pdf'),
  p_c, width = 7, height = 3.5)

message('Wrote three figures to ', fig_dir, ':')
message('  02-power-by-spec.pdf')
message('  02-heatmap-matched-vs-mismatched.pdf')
message('  02-type1-boxplots.pdf')
