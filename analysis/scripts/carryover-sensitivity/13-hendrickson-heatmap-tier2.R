## analysis/scripts/carryover-sensitivity/13-hendrickson-heatmap-tier2.R
##
## Tier 2 sensitivity blocks (S1-S4, S7, S8, S9), in the same
## Hendrickson-style annotated-heatmap format as
## 12-hendrickson-heatmap.R: a blue (low power)-to-red (high power)
## fill scale with the power value printed in each cell. Unlike the
## Tier 1 heatmaps, none of these blocks vary trial design or N (all
## are anchored at the single Hybrid, N = 70 reference cell), so
## there is no design x N facet grid here; each panel instead facets
## on whatever secondary axis that block varies (S2: true half-life;
## S3: dropout mechanism; S8: DGP architecture), with analysis
## specification as the shared y-axis throughout. These are
## additional to, not replacements for, the existing line-plot
## figures for S1-S4 (04-sensitivity-figures-extra-slim.R) and the
## S7/S8 kable tables in report.Rmd; the line plots read half-life
## and dropout as continuous trends better than a heatmap can, and
## the heatmaps read the specification comparison at a glance better
## than the line plots can. Writes to analysis/figures/ with a
## 02xs-heatmap-sens- prefix.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

repo_root <- here::here()

source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

hendrickson_fill <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  limits = c(0, 1), breaks = c(0.25, 0.50, 0.75, 1.00),
  guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(2.6, 'cm')))

theme_heat <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'grey85', colour = NA),
    legend.position = 'right',
    axis.text.x = element_text(size = 8)
  )

heat_panel <- function(d, x_var, x_lab, title, facet_var = NULL) {
  p <- ggplot(d, aes(x = factor(!!rlang::sym(x_var)), y = spec, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', power)),
              colour = 'grey10', size = 2.8) +
    hendrickson_fill +
    labs(x = x_lab, y = 'Analysis specification', fill = 'Power',
        title = title) +
    theme_heat
  if (!is.null(facet_var))
    p <- p + facet_wrap(as.formula(paste('~', facet_var)))
  p
}

## -----------------------------------------------------------------
## S1: autocorrelation
## -----------------------------------------------------------------

s2 <- readRDS(file.path(repo_root,
  'analysis/data/02-sensitivity-summary.rds'))$summary |>
  dplyr::filter(spec %in% spec_order) |>
  dplyr::mutate(spec = spec_factor(spec))

d1 <- s2 |> dplyr::filter(block == 'S1')
p1 <- heat_panel(d1, 'rho', expression('AR(1) autocorrelation'~rho),
  'S1: Sensitivity to within-factor autocorrelation')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S1.pdf'), p1,
      width = 5.2, height = 2.6)

## -----------------------------------------------------------------
## S2: half-life mismatch (facet: true half-life)
## -----------------------------------------------------------------

d2 <- s2 |> dplyr::filter(block == 'S2') |>
  dplyr::mutate(true_lab = factor(sprintf('True t[1/2] = %.1f', t1half)))
p2 <- heat_panel(d2, 'analysis_t1half',
  expression('Analyst-assumed half-life ('*t['1/2']*', weeks)'),
  'S2: Cost of analyst-truth half-life mismatch', 'true_lab')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S2.pdf'), p2,
      width = 7.2, height = 2.8)

## -----------------------------------------------------------------
## S3: dropout (facet: mechanism)
## -----------------------------------------------------------------

d3 <- s2 |> dplyr::filter(block == 'S3')
p3 <- heat_panel(d3, 'dropout_rate', 'Dropout rate',
  'S3: Sensitivity to dropout', 'dropout_mech')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S3.pdf'), p3,
      width = 6.4, height = 2.8)

## -----------------------------------------------------------------
## S4: effect-size curve
## -----------------------------------------------------------------

d4 <- s2 |> dplyr::filter(block == 'S4')
p4 <- heat_panel(d4, 'c_bm', expression(c[bm]),
  'S4: Biomarker-moderation effect-size curve')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S4.pdf'), p4,
      width = 5.8, height = 2.6)

## -----------------------------------------------------------------
## S1-S4 (G1-G9 extension): the same four blocks, all nine
## specifications fit directly (common random numbers within each
## cell), from 17-run-sensitivity-s1234-g9.R. Preliminary n_sim = 100
## (medium smoke test); additional to, not a replacement for, the
## three-spec panels above.
## -----------------------------------------------------------------

s1234g9_path <- file.path(repo_root,
  'analysis/data/02-sensitivity-summary-g9.rds')
if (file.exists(s1234g9_path)) {
  s1234g9 <- readRDS(s1234g9_path)$summary |>
    dplyr::mutate(spec = factor(g_labels_short[spec],
                                levels = unname(g_labels_short[g_order])))

  d1g9 <- s1234g9 |> dplyr::filter(block == 'S1')
  p1g9 <- heat_panel(d1g9, 'rho', expression('AR(1) autocorrelation'~rho),
    'S1 (G1-G9): sensitivity to within-factor autocorrelation')
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S1-g9.pdf'), p1g9,
        width = 7.2, height = 4.2)

  d2g9 <- s1234g9 |> dplyr::filter(block == 'S2') |>
    dplyr::mutate(true_lab = factor(sprintf('True t[1/2] = %.1f', t1half)))
  p2g9 <- heat_panel(d2g9, 'analysis_t1half',
    expression('Analyst-assumed half-life ('*t['1/2']*', weeks)'),
    'S2 (G1-G9): cost of analyst-truth half-life mismatch', 'true_lab')
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S2-g9.pdf'), p2g9,
        width = 8.4, height = 4.4)

  d3g9 <- s1234g9 |> dplyr::filter(block == 'S3')
  p3g9 <- heat_panel(d3g9, 'dropout_rate', 'Dropout rate',
    'S3 (G1-G9): sensitivity to dropout', 'dropout_mech')
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S3-g9.pdf'), p3g9,
        width = 7.6, height = 4.4)

  d4g9 <- s1234g9 |> dplyr::filter(block == 'S4')
  p4g9 <- heat_panel(d4g9, 'c_bm', expression(c[bm]),
    'S4 (G1-G9): biomarker-moderation effect-size curve')
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S4-g9.pdf'), p4g9,
        width = 7.2, height = 4.2)
} else {
  message('Skipping S1-S4 G1-G9 panels: ', s1234g9_path, ' not found.')
}

## -----------------------------------------------------------------
## S7: alternative-paradigm specifications (AIC-selected half-life,
## paired-difference regression) vs Unadjusted. All single-coefficient
## tests; see report.Rmd Section 3.6.
## -----------------------------------------------------------------

s7 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/09-sensitivity-s7.rds'))
lab_g5 <- c(E1 = 'G1: Unadjusted', E3 = 'G2: Lag-\nadjusted',
           E2 = 'G3: Exposure-\nweighted', E7 = 'G4: AIC-selected\nt1/2',
           E9 = 'G5: Paired-\ndifference')
d7 <- s7$summary |>
  dplyr::mutate(spec = factor(lab_g5[spec], levels = unname(lab_g5)),
               dummy = 'Reference cell')
p7 <- ggplot(d7, aes(x = spec, y = dummy, fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
           size = 2.8) +
  hendrickson_fill +
  labs(x = NULL, y = NULL, fill = 'Power',
      title = 'S7: G1-G5 at the Hybrid reference cell') +
  theme_heat +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
       axis.text.x = element_text(size = 7))
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S7.pdf'), p7,
      width = 7.2, height = 2.0)

## -----------------------------------------------------------------
## S8: architecture sensitivity of the specification ranking
## -----------------------------------------------------------------

s8 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/10-sensitivity-s8.rds'))
## S8 is now a 3-way grid (architecture x t1/2 x spec). Panel A is
## Architecture A, Panel B is Architecture B: within each, rows are
## half-life and columns are analysis specification, so the panel
## structure directly answers "does the specification ranking change
## with carryover severity, separately in each architecture?" rather
## than collapsing half-life into a single reference value.
arch8 <- c(mvn = 'B (covariance)', mean_moderation = 'A (mean)')
d8 <- s8$summary |>
  dplyr::mutate(spec = factor(lab_g5[spec], levels = unname(lab_g5)),
               arch = factor(arch8[dgp_arch], levels = unname(arch8)),
               t1half_lab = factor(sprintf('t[1/2] = %.1f', t1half),
                                   levels = sprintf('t[1/2] = %.1f',
                                                    c(0, 0.5, 1.0))))

s8_panel <- function(d, title) {
  ggplot(d, aes(x = spec, y = t1half_lab, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
             size = 2.8) +
    hendrickson_fill +
    labs(x = 'Analysis specification', y = expression(t['1/2']~'(weeks)'),
        fill = 'Power', title = title) +
    theme_heat +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
}

p8a <- s8_panel(d8 |> dplyr::filter(arch == 'A (mean)'),
  'S8 Panel A: Architecture A (mean moderation)')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S8-a.pdf'), p8a,
      width = 6.4, height = 3.2)

p8b <- s8_panel(d8 |> dplyr::filter(arch == 'B (covariance)'),
  'S8 Panel B: Architecture B (covariance moderation)')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S8-b.pdf'), p8b,
      width = 6.4, height = 3.2)

## -----------------------------------------------------------------
## S9: does the ranking among E1/E7/E9 established at the Hybrid
## reference cell (S7) hold under the CO design?
## -----------------------------------------------------------------

s9 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/14-sensitivity-s9.rds'))
d9 <- s9$summary |>
  dplyr::mutate(spec = factor(lab_g5[spec], levels = unname(lab_g5)),
               t1half_lab = factor(sprintf('t[1/2] = %.1f', t1half),
                                   levels = sprintf('t[1/2] = %.1f',
                                                    c(0.5, 1.0))))
p9 <- ggplot(d9, aes(x = spec, y = t1half_lab, fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
           size = 2.8) +
  hendrickson_fill +
  labs(x = NULL, y = expression(t['1/2']~'(weeks)'), fill = 'Power',
      title = 'S9: G1-G5 under the CO design') +
  theme_heat +
  theme(axis.text.x = element_text(size = 7))
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S9.pdf'), p9,
      width = 7.2, height = 2.4)

## -----------------------------------------------------------------
## S10: G1-G9 (G1-G5 base specifications plus CR2 cluster-robust
## variants of G1, G2, G4; G8 = G3+CR2 already exists from Block S6;
## G10 = G5+CR2 does not apply, no cluster structure left in G5's
## one-row-per-patient design). G8's number comes from an
## independent simulation (Block S6, different seed), not shared
## random numbers with the rest of this panel.
## -----------------------------------------------------------------

s10 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/15-sensitivity-s10-g.rds'))
s6 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/diag-s6-cr2.rds'))
g8_power <- as.data.frame(s6$summary) |>
  dplyr::filter(cell == 'Reference (N=70)', spec == 'E2') |>
  dplyr::pull(power_cr2)

lab10 <- c(E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4', E9 = 'G5',
          E1cr2 = 'G6', E3cr2 = 'G7', E7cr2 = 'G9')
d10 <- s10$summary |>
  dplyr::transmute(spec = lab10[spec], power) |>
  dplyr::bind_rows(tibble::tibble(spec = 'G8', power = g8_power)) |>
  dplyr::mutate(spec = factor(spec,
                              levels = c('G1','G2','G3','G4','G5',
                                        'G6','G7','G8','G9')),
               dummy = 'Reference cell')

p10 <- ggplot(d10, aes(x = spec, y = dummy, fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
           size = 2.8) +
  hendrickson_fill +
  labs(x = NULL, y = NULL, fill = 'Power',
      title = 'S10: G1-G9, Hybrid reference cell') +
  theme_heat +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S10.pdf'), p10,
      width = 7.2, height = 1.9)

## -----------------------------------------------------------------
## S6: empirical Type I error at the null, model-based vs CR2, across
## the five S6 stress cells and G1/G2/G3. Replaces the
## 02xs-type1-cr2.pdf dot plot (07-type1-cr2-figure-extra-slim.R)
## with two annotated heatmap panels sharing one diverging fill scale
## centered at nominal 0.05 (blue = conservative, red =
## anti-conservative), matching Panel D's convention in
## 12-hendrickson-heatmap.R.
## -----------------------------------------------------------------

s6_diag <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/diag-s6-cr2.rds'))
cell_levels <- c('Reference (N=70)', 'Small N (N=35)',
                 'High autocorr (rho=0.9)',
                 'Half-life mis-spec (assume 0.25)', '30% MCAR dropout')
cell_short  <- c('Reference', 'Small N', 'High rho',
                 'HL mis-spec', '30% dropout')
lab_g3 <- c(E1 = 'G1: Unadjusted', E3 = 'G2: Lag-adjusted',
           E2 = 'G3: Exposure-weighted')

d6 <- s6_diag$type1 |>
  dplyr::filter(spec %in% names(lab_g3)) |>
  dplyr::mutate(cell = factor(cell, levels = cell_levels, labels = cell_short),
               spec = factor(lab_g3[spec], levels = unname(lab_g3)))

type1_fill_s6 <- scale_fill_gradient2(
  low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0.05,
  limits = range(c(d6$type1_mod, d6$type1_cr2, 0.05)),
  breaks = scales::pretty_breaks(4),
  guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(2.4, 'cm')))

s6_panel <- function(var, title) {
  ggplot(d6, aes(x = cell, y = spec, fill = .data[[var]])) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.3f', .data[[var]])),
             colour = 'grey10', size = 2.6) +
    type1_fill_s6 +
    labs(x = NULL, y = NULL, fill = 'Type I error', title = title) +
    theme_heat +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
}

p6a <- s6_panel('type1_mod', 'S6 Panel A: model-based Type I error')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S6-model.pdf'), p6a,
      width = 6.4, height = 2.8)

p6b <- s6_panel('type1_cr2', 'S6 Panel B: CR2 cluster-robust Type I error')
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S6-cr2.pdf'), p6b,
      width = 6.4, height = 2.8)

## -----------------------------------------------------------------
## S6 (G1-G9 extension): power and Type I error across the same five
## stress cells, all nine specifications fit directly (common random
## numbers within each cell/arm), from 18-run-sensitivity-s6-g9.R.
## Preliminary n_sim = 100 (medium smoke test); supersedes the
## S6-model/S6-cr2 panels above for the full spec set, does not
## replace them (those stay as the original 3-spec, higher-precision
## view).
## -----------------------------------------------------------------

s6g9_path <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/18-sensitivity-s6-g9.rds')
if (file.exists(s6g9_path)) {
  s6g9 <- readRDS(s6g9_path)
  d6g9 <- s6g9$summary |>
    dplyr::mutate(cell = factor(cell, levels = cell_levels, labels = cell_short),
                 spec = factor(g_labels_short[spec],
                               levels = unname(g_labels_short[g_order])))

  p6g9_power <- ggplot(d6g9 |> dplyr::filter(arm == 'power'),
                       aes(x = cell, y = spec, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
             size = 2.6) +
    hendrickson_fill +
    labs(x = NULL, y = NULL, fill = 'Power',
        title = 'S6 (G1-G9): power across five stress cells') +
    theme_heat +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S6-g9-power.pdf'), p6g9_power,
        width = 7.2, height = 3.4)

  d6g9_null <- d6g9 |> dplyr::filter(arm == 'null')
  type1_fill_s6g9 <- scale_fill_gradient2(
    low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0.05,
    limits = range(c(d6g9_null$power, 0.05)),
    breaks = scales::pretty_breaks(4),
    guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                           barheight = grid::unit(2.8, 'cm')))
  p6g9_null <- ggplot(d6g9_null, aes(x = cell, y = spec, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.3f', power)), colour = 'grey10',
             size = 2.6) +
    type1_fill_s6g9 +
    labs(x = NULL, y = NULL, fill = 'Type I error',
        title = 'S6 (G1-G9): Type I error across five stress cells') +
    theme_heat +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))
  ggsave(file.path(fig_dir, '02xs-heatmap-sens-S6-g9-null.pdf'), p6g9_null,
        width = 7.2, height = 3.4)
} else {
  message('Skipping S6 G1-G9 panels: ', s6g9_path, ' not found.')
}

## -----------------------------------------------------------------
## S11: Type I error at the null (c_bm = 0) for G1-G5 plus G9, at
## the Hybrid reference cell, checking whether G4's AIC-selection
## step (post-selection inference risk) inflates Type I error.
## -----------------------------------------------------------------

s11 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/16-sensitivity-s11.rds'))
lab11 <- c(E1 = 'G1', E3 = 'G2', E2 = 'G3', E7 = 'G4', E9 = 'G5',
          E7cr2 = 'G9')
d11 <- s11$summary |>
  dplyr::transmute(spec = factor(lab11[spec],
                                 levels = c('G1','G2','G3','G4','G5','G9')),
                   type1, dummy = 'Null cell (Hybrid)')

type1_fill_s11 <- scale_fill_gradient2(
  low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0.05,
  limits = range(c(d11$type1, 0.05)),
  breaks = scales::pretty_breaks(4),
  guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(2.4, 'cm')))

p11 <- ggplot(d11, aes(x = spec, y = dummy, fill = type1)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.3f', type1)), colour = 'grey10',
           size = 2.8) +
  type1_fill_s11 +
  labs(x = NULL, y = NULL, fill = 'Type I error',
      title = 'S11: Type I error, G1-G5 and G9 at the null') +
  theme_heat +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S11.pdf'), p11,
      width = 6.4, height = 1.9)

message('Wrote Tier 2 Hendrickson-style heatmaps to ', fig_dir, ':')
message('  02xs-heatmap-sens-S1.pdf .. S4.pdf, S1-g9.pdf .. S4-g9.pdf,')
message('  S6-model/cr2.pdf, S6-g9-power/null.pdf, S7.pdf, S8-a/b.pdf, S9.pdf, S10.pdf, S11.pdf')
