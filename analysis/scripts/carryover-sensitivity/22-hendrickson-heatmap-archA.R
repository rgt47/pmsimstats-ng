## analysis/scripts/carryover-sensitivity/22-hendrickson-heatmap-archA.R
##
## Architecture A (mean moderation) counterpart to Figures 1-3 of the
## main manuscript (12-hendrickson-heatmap.R Panels A/B/C), for the
## Supplement's exploratory look at whether the Tier 1 grid's shape
## changes under mean moderation. The main manuscript's scope is
## Architecture B only (Section 2.3); these figures are supplementary
## material, not a replacement for Figures 1-3.
##
## Panels A/B (G1-G9) are from a dedicated 30-cell simulation
## restricted to Architecture A
## (21-run-tier1-hendrickson-g9-archA.R,
## 02-grid-summary-hendrickson-g9-archA.rds, n_sim = 250, after a
## 100-rep smoke test passed sanity checks), mirroring
## 19-run-tier1-hendrickson-g9.R's
## Architecture B design exactly except for dgp_arch.
##
## Panel C (G1-G3, decay-shape axis) needs no new simulation: the
## production Tier 1 grid (02-grid-summary.rds, from
## 01-run-factorial.R) already crosses both architectures at
## n_sim = 500, so this panel is read directly from existing data at
## full production precision, unlike Panels A/B.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

repo_root <- here::here()

source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/spec-labels.R'))

decay_levels <- c('Exponential', 'Weibull k=0.7', 'Weibull k=1.0',
                  'Weibull k=1.5')

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

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
## Panels A/B: G1-G9, Architecture A
## -----------------------------------------------------------------

grid_g9 <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary-hendrickson-g9-archA.rds'))$summary |>
  dplyr::mutate(
    spec = factor(g_labels_short[spec], levels = unname(g_labels_short[g_order])),
    design = factor(design, levels = c('CO', 'Hybrid', 'OLBDC'),
                    labels = c('CO', 'Hybrid', 'OL+BDC')),
    N_label = factor(paste0('N = ', N), levels = c('N = 35', 'N = 70'))
  )

d_a <- grid_g9 |> dplyr::filter(panel == 'A')

p_a <- hendrickson_heatmap(d_a, 'c_bm', expression('Biomarker moderation ('*c[bm]*')'),
  expression('A (Architecture A): trial design, analysis specification,'~
             'and biomarker effect on power'~(t['1/2']==1.0))) +
  labs(y = NULL)

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-a-archA.pdf'),
  p_a, width = 7.6, height = 6.4)

d_b <- dplyr::bind_rows(
  grid_g9 |> dplyr::filter(panel == 'B'),
  grid_g9 |> dplyr::filter(panel == 'A', c_bm == 0.45)
)

p_b <- hendrickson_heatmap(d_b, 't1half', expression('Carryover half-life ('*t['1/2']*', weeks)'),
  expression('B (Architecture A): effect of carryover half-life on power'~(c[bm]==0.45))) +
  labs(y = NULL)

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-b-archA.pdf'),
  p_b, width = 7.6, height = 6.4)

## -----------------------------------------------------------------
## Panel C: DGP decay shape, Architecture A (G1-G3, full production
## precision, from the existing 540-cell grid)
## -----------------------------------------------------------------

grid_c <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary.rds'))$summary |>
  dplyr::filter(dgp_arch == 'mean_moderation',
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
  ) |>
  dplyr::filter(c_bm == 0.45, t1half == 1.0)

p_c <- hendrickson_heatmap(grid_c, 'decay_label', 'DGP decay shape',
  expression('C (Architecture A): effect of DGP decay shape on power'~
             (c[bm]==0.45*','~t['1/2']==1.0))) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 7))

ggsave(file.path(fig_dir, '02xs-heatmap-hendrickson-c-archA.pdf'),
  p_c, width = 7.4, height = 4.6)

message('Wrote three Architecture A Hendrickson-style heatmaps to ', fig_dir, ':')
message('  02xs-heatmap-hendrickson-a-archA.pdf (biomarker effect, G1-G9, n_sim=250)')
message('  02xs-heatmap-hendrickson-b-archA.pdf (carryover half-life, G1-G9, n_sim=250)')
message('  02xs-heatmap-hendrickson-c-archA.pdf (DGP decay shape, G1-G3, n_sim=500)')
