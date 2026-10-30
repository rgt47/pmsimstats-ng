# Size and power of the one-cycle crossover against the Hybrid design,
# drawn from two stored simulation grids, for docs/68 Section 9.
#
# Source 1: paper 01 G9 planning grid (19-g9-planning-and-fig4b.R,
#   block P, simrun job 40). Extended DGP, G9 only, N = 70; 2,000
#   null and 1,000 alternative replicates per cell.
# Source 2: paper 02 analysis-ladder grid (27-run-tier1-g-vs-f-e0.R,
#   simrun job 27). Extended DGP, E0 and G1 to G9 plus E9, N = 35 and
#   70; 500 replicates per cell.
#
# No simulation is run here. Outputs: three PNG heatmaps in
# docs/figures/ and the underlying tables in analysis/data/
# tutorial-design/.
#
# Usage: Rscript analysis/scripts/tutorial-design/01-co-vs-hybrid-heatmaps.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

root <- here::here()
out_dat <- file.path(root, 'analysis', 'data', 'tutorial-design')
out_fig <- file.path(root, 'docs', 'figures')
dir.create(out_dat, recursive = TRUE, showWarnings = FALSE)

# Reference palette (dataviz skill, references/palette.md): blue ramp
# for magnitude; blue <-> red with a gray midpoint for polarity.
seq_lo <- '#cde2fb'
seq_hi <- '#0d366b'
div_lo <- '#2a78d6'
div_mid <- '#f0efec'
div_hi <- '#e34948'
ink <- '#0b0b0b'
ink_muted <- '#52514e'

theme_heat <- function() {
  theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = ink_muted),
          axis.title = element_text(colour = ink_muted),
          strip.text = element_text(colour = ink, face = 'bold'),
          plot.title = element_text(colour = ink, face = 'bold',
                                    size = 10),
          legend.title = element_text(colour = ink_muted, size = 9),
          plot.background = element_rect(fill = '#fcfcfb',
                                         colour = NA))
}

# Text on a tile: white on the dark half of the blue ramp, else ink.
label_ink <- function(x, cut) ifelse(x > cut, '#ffffff', ink)

# Every simulated power lies between 0.35 and 1, so the power ramp
# starts at 0.3 to keep the cells distinguishable; values below it
# would be drawn at the lightest step.
power_scale <- function() {
  scale_fill_gradient(low = seq_lo, high = seq_hi, limits = c(0.3, 1),
                      oob = scales::squish, name = 'Power')
}

mcse <- function(p, n) sqrt(p * (1 - p) / n)

## Source 1: G9 at N = 70 ---------------------------------------------

p1 <- readRDS(file.path(root, 'analysis', 'data', 'quick-sim',
                        'carryover-closed-form', 'g9-planning',
                        'summary.rds'))$P |>
  filter(design %in% c('Hybrid', 'CO'), N == 70) |>
  transmute(design = factor(design, c('Hybrid', 'CO'),
                            c('Hybrid', 'Crossover')),
            c_bm, t_half, n = n_p, rate = power, mcse,
            e9_power)
stopifnot(nrow(p1) == 24, all(!is.na(p1$rate)))

size_g9 <- p1 |>
  filter(c_bm == 0) |>
  select(design, t_half, size = rate, mcse, n)
power_g9 <- p1 |>
  filter(c_bm > 0) |>
  select(design, c_bm, t_half, power = rate, mcse, e9_power, n)
write.csv(size_g9, file.path(out_dat, 'g9-size-n70.csv'),
          row.names = FALSE)
write.csv(power_g9, file.path(out_dat, 'g9-power-n70.csv'),
          row.names = FALSE)

diff_g9 <- power_g9 |>
  select(design, c_bm, t_half, power, mcse) |>
  pivot_wider(names_from = design, values_from = c(power, mcse)) |>
  mutate(diff = power_Hybrid - power_Crossover,
         diff_se = sqrt(mcse_Hybrid^2 + mcse_Crossover^2))
write.csv(diff_g9, file.path(out_dat, 'g9-power-difference-n70.csv'),
          row.names = FALSE)

hl_lab <- function(x) factor(x, sort(unique(x)),
                             paste(sort(unique(x)), 'wk'))

f1a <- power_g9 |>
  mutate(t_half = hl_lab(t_half), c_bm = factor(c_bm)) |>
  ggplot(aes(t_half, c_bm, fill = power)) +
  geom_tile(colour = '#fcfcfb', linewidth = 1) +
  geom_text(aes(label = sprintf('%.2f', power),
                colour = label_ink(power, 0.7)), size = 3.2) +
  scale_colour_identity() +
  power_scale() +
  facet_wrap(~ design) +
  labs(x = 'Carryover half-life', y = expression(c[bm]),
       title = 'A. G9 power, N = 70') +
  theme_heat()

f1b <- diff_g9 |>
  mutate(t_half = hl_lab(t_half), c_bm = factor(c_bm)) |>
  ggplot(aes(t_half, c_bm, fill = diff)) +
  geom_tile(colour = '#fcfcfb', linewidth = 1) +
  geom_text(aes(label = sprintf('%+.2f', diff)), colour = ink,
            size = 3.2) +
  scale_fill_gradient2(low = div_hi, mid = div_mid, high = div_lo,
                       midpoint = 0, limits = c(-0.35, 0.35),
                       name = 'Hybrid\nminus CO') +
  labs(x = 'Carryover half-life', y = expression(c[bm]),
       title = 'B. Hybrid minus crossover') +
  theme_heat()

ggsave(file.path(out_fig, '68-fig1-g9-power-hybrid-vs-co.png'),
       f1a + f1b + plot_layout(widths = c(2, 1)),
       width = 9, height = 3.2, dpi = 200, bg = '#fcfcfb')

## Source 2: the analysis ladder, N = 35 and 70 -----------------------

ladder <- c(E0 = 'E0', E1 = 'G1', E1cr2 = 'G6', E2 = 'G3',
            E2cr2 = 'G8', E7 = 'G4', E7cr2 = 'G9', E9 = 'G5 (E9)')

s2 <- as.data.frame(readRDS(file.path(
  root, 'analysis', 'data', 'carryover-sensitivity',
  'tier1-g-vs-f-e0', 'summary.rds'))) |>
  filter(design %in% c('Hybrid', 'CO'), method %in% names(ladder)) |>
  transmute(design = factor(design, c('Hybrid', 'CO'),
                            c('Hybrid', 'Crossover')),
            N, c_bm, t_half = t1half,
            analysis = factor(ladder[method], rev(ladder)),
            rate = reject, n, mcse = mcse(reject, n))
stopifnot(!anyNA(s2$rate), all(s2$n >= 490))

size_ladder <- s2 |> filter(c_bm == 0)
power_ladder <- s2 |> filter(c_bm == 0.45)
write.csv(size_ladder, file.path(out_dat, 'ladder-size.csv'),
          row.names = FALSE)
write.csv(power_ladder, file.path(out_dat, 'ladder-power-cbm045.csv'),
          row.names = FALSE)

f2 <- size_ladder |>
  mutate(col = paste0(design, '\nN = ', N),
         col = factor(col, unique(col[order(design, N)]))) |>
  ggplot(aes(col, analysis, fill = pmin(rate, 0.10))) +
  geom_tile(colour = '#fcfcfb', linewidth = 1) +
  geom_text(aes(label = sprintf('%.3f', rate)), colour = ink,
            size = 3.2) +
  scale_fill_gradient2(low = div_lo, mid = div_mid, high = div_hi,
                       midpoint = 0.05, limits = c(0, 0.10),
                       breaks = c(0, 0.025, 0.05, 0.075, 0.10),
                       labels = c('0', '0.025', '0.05', '0.075',
                                  '>= 0.10'),
                       name = 'Type I\nerror') +
  labs(x = NULL, y = NULL,
       title = paste('Type I error by analysis, half-life 1 week,',
                     '500 replicates per cell')) +
  theme_heat()

ggsave(file.path(out_fig, '68-fig2-ladder-size-hybrid-vs-co.png'), f2,
       width = 6, height = 3.6, dpi = 200, bg = '#fcfcfb')

f3 <- power_ladder |>
  mutate(t_half = hl_lab(t_half),
         design_n = factor(paste0(design, ', N = ', N),
                           c('Hybrid, N = 35', 'Crossover, N = 35',
                             'Hybrid, N = 70', 'Crossover, N = 70'))) |>
  ggplot(aes(t_half, analysis, fill = rate)) +
  geom_tile(colour = '#fcfcfb', linewidth = 1) +
  geom_text(aes(label = sprintf('%.2f', rate),
                colour = label_ink(rate, 0.7)), size = 3) +
  scale_colour_identity() +
  power_scale() +
  facet_wrap(~ design_n, nrow = 1) +
  labs(x = 'Carryover half-life', y = NULL,
       title = expression(bold('Power by analysis, ') *
                            bolditalic(c)[bold(bm)] * bold(' = 0.45'))) +
  theme_heat()

ggsave(file.path(out_fig, '68-fig3-ladder-power-hybrid-vs-co.png'), f3,
       width = 9, height = 3.6, dpi = 200, bg = '#fcfcfb')

## Markdown tables for docs/68 ----------------------------------------

fmt <- function(p, se) sprintf('%.3f (%.3f)', p, se)

cat('\n### G9 size, N = 70\n\n')
size_g9 |>
  mutate(cell = fmt(size, mcse)) |>
  select(design, t_half, cell) |>
  pivot_wider(names_from = t_half, values_from = cell) |>
  knitr::kable(format = 'pipe') |>
  print()

cat('\n### G9 power, N = 70\n\n')
power_g9 |>
  select(design, c_bm, t_half, power, mcse) |>
  pivot_wider(names_from = design, values_from = c(power, mcse)) |>
  left_join(select(diff_g9, c_bm, t_half, diff, diff_se),
            by = c('c_bm', 't_half')) |>
  transmute(c_bm, t_half,
            Hybrid = fmt(power_Hybrid, mcse_Hybrid),
            Crossover = fmt(power_Crossover, mcse_Crossover),
            Difference = fmt(diff, diff_se)) |>
  knitr::kable(format = 'pipe') |>
  print()

cat('\n### Ladder size, half-life 1 week\n\n')
size_ladder |>
  mutate(col = paste0(design, ' N=', N)) |>
  select(analysis, col, rate) |>
  mutate(rate = sprintf('%.3f', rate)) |>
  pivot_wider(names_from = col, values_from = rate) |>
  arrange(desc(analysis)) |>
  knitr::kable(format = 'pipe') |>
  print()

cat('\n### Ladder power, c_bm = 0.45, N = 70\n\n')
power_ladder |>
  filter(N == 70) |>
  mutate(col = paste0(design, ' ', t_half)) |>
  select(analysis, col, rate) |>
  mutate(rate = sprintf('%.3f', rate)) |>
  pivot_wider(names_from = col, values_from = rate) |>
  arrange(desc(analysis)) |>
  knitr::kable(format = 'pipe') |>
  print()
