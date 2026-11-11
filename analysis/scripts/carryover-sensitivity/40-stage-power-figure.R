## analysis/scripts/carryover-sensitivity/40-stage-power-figure.R
##
## Power against carryover half-life for the five rungs of paper 02's
## analysis ladder (A0 to A4) on the corrected 30-cell grid (simrun
## job 27), one panel per design. A4 (stored E7null_cr2, G10 in the
## companion papers) was fitted afterward to the same trials by
## 41-run-g10-backfill.R; replicates whose G9 refit failed the pairing
## guard are dropped. It replaces the job 26 figure 02xs-power-by-spec
## for paper 02. It reads the per-replicate output and runs no
## simulation. Run from the repository root:
##   Rscript analysis/scripts/carryover-sensitivity/40-stage-power-figure.R
## Writes analysis/figures/02-power-by-half-life.{pdf,png}.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

root <- here::here()
source(file.path(root, 'analysis/scripts/carryover-sensitivity',
                 'spec-labels.R'))
dat_dir <- file.path(root, 'analysis/data/carryover-sensitivity')
a <- readRDS(file.path(dat_dir, 'tier1-g-vs-f-e0', 'all.rds'))
s <- readRDS(file.path(dat_dir, 'tier1-g-vs-f-e0', 'summary.rds'))
bf <- readRDS(file.path(dat_dir, 'g10-backfill', 'all.rds'))
guard <- readRDS(file.path(dat_dir, 'g10-backfill', 'guard.rds'))

## A0 (E0) is liberal, so its power is not comparable; it is omitted.
stored <- c('E1', 'E1cr2', 'E2cr2')
meta <- a |> distinct(cell_id, design, N, c_bm, t1half)
a4 <- bf |>
  filter(source == 't1', method == 'E7null_cr2') |>
  semi_join(guard |> filter(match), by = c('source', 'cell_id', 'rep')) |>
  select(cell_id, rep, method, p_value, converged) |>
  left_join(meta, by = 'cell_id')
stopifnot(!anyNA(a4$t1half))

a <- a |>
  filter(method %in% stored) |>
  select(cell_id, rep, method, p_value, converged, design, N, c_bm,
         t1half) |>
  bind_rows(a4) |>
  filter(c_bm == 0.45, N == 70)

cell <- a |>
  filter(converged, !is.na(p_value)) |>
  group_by(design, t1half, method) |>
  summarise(n = n(), power = mean(p_value < 0.05), .groups = 'drop') |>
  mutate(mcse = sqrt(power * (1 - power) / n))

## Reject rates of the stored analyses must reproduce the stored summary.
chk <- cell |>
  filter(method %in% stored) |>
  inner_join(s |> filter(c_bm == 0.45, N == 70) |>
               select(design, t1half, method, reject),
             by = c('design', 't1half', 'method'))
stopifnot(nrow(chk) == sum(cell$method %in% stored),
          max(abs(chk$power - chk$reject)) < 1e-9,
          sum(cell$method == 'E7null_cr2') == 9)

d <- cell |>
  mutate(Analysis = factor(to_a(method), levels = setdiff(a_order, 'A0')),
         design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'),
                         labels = c('CO', 'OL+BDC', 'Hybrid')))

cols <- c(A0 = '#7F7F7F', A1 = '#56B4E9', A2 = '#0072B2',
          A3 = '#009E73', A4 = '#D55E00')
shp <- c(A0 = 15, A1 = 1, A2 = 16, A3 = 16, A4 = 17)
lty <- c(A0 = 'dotted', A1 = 'dashed', A2 = 'solid', A3 = 'solid',
         A4 = 'solid')
dodge <- position_dodge(width = 0.06)

p <- ggplot(d, aes(t1half, power, colour = Analysis, shape = Analysis,
                   linetype = Analysis, group = Analysis)) +
  geom_line(position = dodge, linewidth = 0.5) +
  geom_errorbar(aes(ymin = power - 2 * mcse, ymax = power + 2 * mcse),
                position = dodge, width = 0, linewidth = 0.3) +
  geom_point(position = dodge, size = 1.8, stroke = 0.7) +
  facet_wrap(~ design, nrow = 1) +
  scale_colour_manual(values = cols) +
  scale_shape_manual(values = shp) +
  scale_linetype_manual(values = lty) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = 'Carryover half-life (weeks)', y = 'Power') +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = 'grey92', colour = NA),
        legend.position = 'top', legend.title = element_blank())

fig_dir <- Sys.getenv('POWER_FIG_DIR', file.path(root, 'analysis/figures'))
ggsave(file.path(fig_dir, '02-power-by-half-life.pdf'), p,
       width = 6.5, height = 3.2)
ggsave(file.path(fig_dir, '02-power-by-half-life.png'), p,
       width = 6.5, height = 3.2, dpi = 300)
