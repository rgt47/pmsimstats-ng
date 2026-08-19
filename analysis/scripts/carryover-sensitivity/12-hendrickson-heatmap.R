## analysis/scripts/carryover-sensitivity/12-hendrickson-heatmap.R
##
## Tier 1 grid, presented in the style of Hendrickson et al. (2020,
## Figure 4): an annotated heatmap with trial design as column
## facets, N as row-block facets (35 top, 70 bottom), a blue
## (low)-yellow (mid)-red (high) fill scale for power, and the power
## value printed in each cell. Hendrickson's row axis within each N
## block was "censoring parameters"; there is no censoring axis in
## this manuscript's grid, so it is replaced with analysis
## specification, the axis this manuscript's Tier 1 grid actually
## varies. Reads the same 02-grid-summary.rds as
## 03-render-figures-extra-slim.R and writes to analysis/figures/
## with a 02xs- prefix, so the full-scope figures used by the
## archived drafts are untouched.
##
## Panel A mirrors Hendrickson Fig. 4A (biomarker effect on the
## x-axis, carryover fixed): t1/2 fixed at 1.0 weeks (nonzero
## carryover, the manuscript's reference half-life), exponential
## DGP, c_bm on the x-axis.
##
## Panel B mirrors Hendrickson Fig. 4B (carryover on the x-axis,
## biomarker effect fixed): c_bm fixed at 0.45, exponential DGP,
## t1/2 on the x-axis.
##
## Panel C has no Hendrickson analogue (their grid had one decay
## form). It varies the manuscript's own decay-shape axis (Section
## 2.5): exponential and the three Weibull shapes, at the reference
## cell (c_bm = 0.45, t1/2 = 1.0).

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

repo_root <- here::here()

source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))

decay_levels <- c('Exponential', 'Weibull k=0.7', 'Weibull k=1.0',
                  'Weibull k=1.5')

grid <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary.rds'))$summary |>
  dplyr::filter(dgp_arch == 'mvn',
                carryover_form %in% c('exponential', 'weibull'),
                spec %in% spec_order) |>
  dplyr::mutate(
    spec = spec_factor(spec),
    design = factor(design, levels = c('CO', 'Hybrid', 'OLBDC'),
                    labels = c('CO', 'Hybrid', 'OL+BDC')),
    N_label = factor(paste0('N = ', N), levels = c('N = 35', 'N = 70')),
    decay_label = factor(
      dplyr::case_when(
        carryover_form == 'exponential' ~ 'Exponential',
        weibull_shape == 0.7 ~ 'Weibull k=0.7',
        weibull_shape == 1.0 ~ 'Weibull k=1.0',
        weibull_shape == 1.5 ~ 'Weibull k=1.5'
      ), levels = decay_levels)
  )

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Hendrickson's fill scale reads as a reversed RdYlBu (blue = low
## power, red = high power), with black cell labels throughout.
hendrickson_fill <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  limits = c(0, 1), breaks = c(0.25, 0.50, 0.75, 1.00),
  guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(3.2, 'cm')))

theme_heat <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'grey85', colour = NA),
    strip.text.y = element_text(angle = 0),
    legend.position = 'right',
    axis.text.x = element_text(size = 8)
  )

hendrickson_heatmap <- function(d, x_var, x_lab, title) {
  ggplot(d, aes(x = factor(!!rlang::sym(x_var)), y = spec, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', power)),
              colour = 'grey10', size = 2.6) +
    facet_grid(N_label ~ design, switch = 'y') +
    hendrickson_fill +
    labs(x = x_lab, y = 'Analysis specification', fill = 'Power',
        title = title) +
    theme_heat
}

## -----------------------------------------------------------------
## Panel A: biomarker effect (t1/2 = 1.0, exponential DGP)
## -----------------------------------------------------------------

d_a <- grid |> dplyr::filter(carryover_form == 'exponential', t1half == 1.0)

p_a <- hendrickson_heatmap(d_a, 'c_bm', expression('Biomarker moderation ('*c[bm]*')'),
  expression('A: Effect of trial design, analysis specification,'~
             'and biomarker effect on power'~(t['1/2']==1.0)))

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-a.pdf'),
  p_a, width = 7.4, height = 4.6)

## -----------------------------------------------------------------
## Panel B: carryover half-life (c_bm = 0.45, exponential DGP)
## -----------------------------------------------------------------

d_b <- grid |> dplyr::filter(carryover_form == 'exponential', c_bm == 0.45)

p_b <- hendrickson_heatmap(d_b, 't1half', expression('Carryover half-life ('*t['1/2']*', weeks)'),
  expression('B: Effect of carryover half-life on power'~(c[bm]==0.45)))

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-b.pdf'),
  p_b, width = 7.4, height = 4.6)

## -----------------------------------------------------------------
## Panel C: DGP decay shape (c_bm = 0.45, t1/2 = 1.0)
## -----------------------------------------------------------------

d_c <- grid |> dplyr::filter(c_bm == 0.45, t1half == 1.0)

p_c <- hendrickson_heatmap(d_c, 'decay_label', 'DGP decay shape',
  expression('C: Effect of DGP decay shape on power'~
             (c[bm]==0.45*','~t['1/2']==1.0))) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-c.pdf'),
  p_c, width = 7.4, height = 4.6)

## -----------------------------------------------------------------
## Panel D: Type I error under the null (c_bm = 0, N = 70), pooled
## across decay forms and half-lives, by design and specification.
## Replaces the 02xs-type1-boxplots.pdf boxplot
## (03-render-figures-extra-slim.R) with an annotated heatmap of the
## mean rejection rate per design x spec cell. Type I error values
## cluster in a narrow band near nominal 0.05, so the power heatmaps'
## fixed 0-1 fill scale would show almost no contrast here; a
## diverging scale centered at 0.05 (blue = conservative, red =
## anti-conservative) is used instead, matching the calibration
## framing of Section 3.5's kappa diagnostic.
## -----------------------------------------------------------------

d_d <- grid |>
  dplyr::filter(c_bm == 0, N == 70) |>
  dplyr::group_by(design, spec) |>
  dplyr::summarise(type1 = mean(power), .groups = 'drop')

type1_fill <- scale_fill_gradient2(
  low = '#2166ac', mid = 'white', high = '#b2182b', midpoint = 0.05,
  limits = range(c(d_d$type1, 0.05)),
  breaks = scales::pretty_breaks(4),
  guide = guide_colorbar(barwidth  = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(3.2, 'cm')))

p_d <- ggplot(d_d, aes(x = design, y = spec, fill = type1)) +
  geom_tile(colour = 'white', linewidth = 0.4) +
  geom_text(aes(label = sprintf('%.3f', type1)),
           colour = 'grey10', size = 2.8) +
  type1_fill +
  labs(x = 'Trial design', y = 'Analysis specification', fill = 'Type I error',
      title = 'D: Empirical Type I error under the null',
      subtitle = expression(c[bm]==0*','~N==70*','~
                            'pooled across decay forms and half-lives'~
                            '(nominal'~alpha==0.05*')')) +
  theme_heat

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-d.pdf'),
  p_d, width = 6.0, height = 3.2)

message('Wrote four Hendrickson-style heatmaps to ', fig_dir, ':')
message('  02xs-heatmap-hendrickson-a.pdf (biomarker effect)')
message('  02xs-heatmap-hendrickson-b.pdf (carryover half-life)')
message('  02xs-heatmap-hendrickson-c.pdf (DGP decay shape)')
message('  02xs-heatmap-hendrickson-d.pdf (Type I error)')
