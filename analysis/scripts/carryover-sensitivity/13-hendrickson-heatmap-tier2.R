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
## S7: alternative-paradigm specifications (AIC-selected half-life,
## paired-difference regression) vs Unadjusted. All single-coefficient
## tests; see report.Rmd Section 3.6.
## -----------------------------------------------------------------

s7 <- readRDS(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/09-sensitivity-s7.rds'))
lab7 <- c(E1 = 'Unadjusted', E7 = 'AIC-selected\nt1/2',
         E9 = 'Paired-\ndifference')
d7 <- s7$summary |>
  dplyr::mutate(spec = factor(lab7[spec], levels = unname(lab7)),
               dummy = 'Reference cell')
p7 <- ggplot(d7, aes(x = spec, y = dummy, fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
           size = 2.8) +
  hendrickson_fill +
  labs(x = NULL, y = NULL, fill = 'Power',
      title = 'S7: Alternative-paradigm specifications') +
  theme_heat +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S7.pdf'), p7,
      width = 5.4, height = 1.9)

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
lab8 <- c(E1 = 'Unadjusted', E3 = 'Lag-adjusted',
         E2 = 'Exposure-\nweighted', E7 = 'AIC-selected\nt1/2',
         E9 = 'Paired-\ndifference')
arch8 <- c(mvn = 'B (covariance)', mean_moderation = 'A (mean)')
d8 <- s8$summary |>
  dplyr::mutate(spec = factor(lab8[spec], levels = unname(lab8)),
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
lab9 <- c(E1 = 'Unadjusted', E7 = 'AIC-selected\nt1/2',
         E9 = 'Paired-\ndifference')
d9 <- s9$summary |>
  dplyr::mutate(spec = factor(lab9[spec], levels = unname(lab9)),
               t1half_lab = factor(sprintf('t[1/2] = %.1f', t1half),
                                   levels = sprintf('t[1/2] = %.1f',
                                                    c(0.5, 1.0))))
p9 <- ggplot(d9, aes(x = spec, y = t1half_lab, fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
           size = 2.8) +
  hendrickson_fill +
  labs(x = NULL, y = expression(t['1/2']~'(weeks)'), fill = 'Power',
      title = 'S9: Alternative-paradigm specifications under CO') +
  theme_heat
ggsave(file.path(fig_dir, '02xs-heatmap-sens-S9.pdf'), p9,
      width = 6.0, height = 2.4)

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

message('Wrote Tier 2 Hendrickson-style heatmaps to ', fig_dir, ':')
message('  02xs-heatmap-sens-S1.pdf .. S4.pdf, S7.pdf, S8-a/b.pdf, S9.pdf, S10.pdf')
