## analysis/scripts/carryover-sensitivity/45-ladder-heatmaps.R
##
## Hendrickson-style power heatmaps (Hendrickson et al. 2020, Figure 4)
## for the five rungs of paper 02's analysis ladder, A0 to A4, on the
## corrected 30-cell grid (simrun job 27). Panel A varies the biomarker
## effect at t1/2 = 1 week; panel B varies the half-life at
## c_bm = 0.45. A4 (stored E7null_cr2) was fitted afterward to the same
## trials by 41-run-g10-backfill.R; replicates whose G9 refit failed the
## pairing guard are dropped. The nine-specification versions written
## by 12-hendrickson-heatmap.R are left untouched. Reads stored data
## only. Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/45-ladder-heatmaps.R
## Writes analysis/figures/02-heatmap-ladder-{a,b}.pdf (or to
## $HEATMAP_FIG_DIR).

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

root <- here::here()
source(file.path(root, 'analysis/scripts/carryover-sensitivity',
                 'spec-labels.R'))
dat_dir <- file.path(root, 'analysis/data/carryover-sensitivity')
s <- readRDS(file.path(dat_dir, 'tier1-g-vs-f-e0', 'summary.rds'))
bf <- readRDS(file.path(dat_dir, 'g10-backfill', 'all.rds'))
guard <- readRDS(file.path(dat_dir, 'g10-backfill', 'guard.rds'))

meta <- s |> distinct(cell_id, design, N, c_bm, t1half, panel)
a4 <- bf |>
  filter(source == 't1', method == 'E7null_cr2', converged,
         !is.na(p_value)) |>
  semi_join(guard |> filter(match), by = c('source', 'cell_id', 'rep')) |>
  group_by(cell_id, method) |>
  summarise(reject = mean(p_value < 0.05), .groups = 'drop') |>
  left_join(meta, by = 'cell_id')
stopifnot(nrow(a4) == 30, !anyNA(a4$panel))

d <- s |>
  ## A0 (E0) is liberal, so its power is not comparable; omitted.
  filter(method %in% c('E1', 'E1cr2', 'E2cr2')) |>
  select(cell_id, design, N, c_bm, t1half, panel, method, reject) |>
  bind_rows(a4) |>
  transmute(design = factor(design, levels = c('CO', 'Hybrid', 'OLBDC'),
                            labels = c('CO', 'Hybrid', 'OL+BDC')),
            N_label = factor(paste0('N = ', N),
                             levels = c('N = 35', 'N = 70')),
            c_bm, t1half, panel,
            spec = factor(to_a(method), levels = rev(setdiff(a_order, 'A0'))),
            power = reject)
stopifnot(!anyNA(d$spec))

fill_scale <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, 'RdYlBu')),
  limits = c(0, 1), breaks = c(0.25, 0.50, 0.75, 1.00),
  guide = guide_colorbar(barwidth = grid::unit(0.35, 'cm'),
                         barheight = grid::unit(3.2, 'cm')))

heat <- function(dd, x_var, x_lab, title) {
  ggplot(dd, aes(x = factor(.data[[x_var]]), y = spec, fill = power)) +
    geom_tile(colour = 'white', linewidth = 0.4) +
    geom_text(aes(label = sprintf('%.2f', power)), colour = 'grey10',
              size = 2.6) +
    facet_grid(N_label ~ design, switch = 'y') +
    fill_scale +
    labs(x = x_lab, y = NULL, fill = 'Power', title = title) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = 'grey85', colour = NA),
          strip.text.y = element_text(angle = 0),
          axis.text.x = element_text(size = 8))
}

fig_dir <- Sys.getenv('HEATMAP_FIG_DIR', file.path(root, 'analysis/figures'))

p_a <- heat(d |> filter(panel == 'A'), 'c_bm',
  expression('Biomarker moderation ('*c[bm]*')'),
  expression('A: Power by design, analysis and biomarker effect'~
             (t['1/2']==1.0)))
ggsave(file.path(fig_dir, '02-heatmap-ladder-a.pdf'), p_a,
       width = 7.0, height = 4.4)

## The t1/2 = 1 week column of panel B is the c_bm = 0.45 column of
## panel A (shared cell).
d_b <- bind_rows(d |> filter(panel == 'B'),
                 d |> filter(panel == 'A', c_bm == 0.45))
p_b <- heat(d_b, 't1half',
  expression('Carryover half-life ('*t['1/2']*', weeks)'),
  expression('B: Power by design, analysis and half-life'~
             (c[bm]==0.45)))
ggsave(file.path(fig_dir, '02-heatmap-ladder-b.pdf'), p_b,
       width = 7.0, height = 4.4)
