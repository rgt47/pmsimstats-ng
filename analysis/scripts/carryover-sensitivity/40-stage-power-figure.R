## analysis/scripts/carryover-sensitivity/40-stage-power-figure.R
##
## Power against carryover half-life for the G family on the corrected
## 30-cell grid (simrun job 27), one panel per design. It replaces the
## job 26 figure 02xs-power-by-spec for paper 02. It reads the
## per-replicate output and runs no simulation. Run from the
## repository root:
##   Rscript analysis/scripts/carryover-sensitivity/40-stage-power-figure.R
## Writes analysis/figures/02-power-by-half-life.{pdf,png}.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

root <- here::here()
all_f <- file.path(root, 'analysis/data/carryover-sensitivity',
                   'tier1-g-vs-f-e0', 'all.rds')
sum_f <- file.path(root, 'analysis/data/carryover-sensitivity',
                   'tier1-g-vs-f-e0', 'summary.rds')
a <- readRDS(all_f)
s <- readRDS(sum_f)

lab <- c(E0 = 'E0', E1 = 'G1', E1cr2 = 'G6', E2 = 'G3', E2cr2 = 'G8',
         E7 = 'G4', E7cr2 = 'G9')
a <- a |> filter(method %in% names(lab), c_bm == 0.45, N == 70)

cell <- a |>
  filter(converged, !is.na(p_value)) |>
  group_by(design, t1half, method) |>
  summarise(n = n(), power = mean(p_value < 0.05), .groups = 'drop') |>
  mutate(mcse = sqrt(power * (1 - power) / n))

## Reject rates must reproduce the stored summary.
chk <- cell |>
  inner_join(s |> filter(c_bm == 0.45, N == 70) |>
               select(design, t1half, method, reject),
             by = c('design', 't1half', 'method'))
stopifnot(nrow(chk) == nrow(cell), max(abs(chk$power - chk$reject)) < 1e-9)

d <- cell |>
  mutate(Analysis = factor(lab[method], levels = unname(lab)),
         design = factor(design, levels = c('CO', 'OLBDC', 'Hybrid'),
                         labels = c('CO', 'OL+BDC', 'Hybrid')))

cols <- c(E0 = '#7F7F7F', G1 = '#56B4E9', G6 = '#0072B2',
          G3 = '#88CCAA', G8 = '#009E73', G4 = '#E69F00', G9 = '#D55E00')
shp <- c(E0 = 15, G1 = 1, G6 = 16, G3 = 1, G8 = 16, G4 = 1, G9 = 16)
lty <- c(E0 = 'dotted', G1 = 'dashed', G6 = 'solid', G3 = 'dashed',
         G8 = 'solid', G4 = 'dashed', G9 = 'solid')
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

fig_dir <- file.path(root, 'analysis/figures')
ggsave(file.path(fig_dir, '02-power-by-half-life.pdf'), p,
       width = 6.5, height = 3.2)
ggsave(file.path(fig_dir, '02-power-by-half-life.png'), p,
       width = 6.5, height = 3.2, dpi = 300)
