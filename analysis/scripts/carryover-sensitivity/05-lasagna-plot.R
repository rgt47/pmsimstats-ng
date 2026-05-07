## analysis/scripts/carryover-sensitivity/05-lasagna-plot.R
##
## Nested-loop / lasagna plot for the 540-cell Tier 1 factorial,
## following the Morris, White, and Crowther (2019) recommendation
## for high-dimensional simulation visualisation
## (doi 10.1002/sim.8086, Section 7).
##
## Layout: one column per analysis specification (A1, A2, A3); rows
## are the 180 unique combinations of (decay form, weibull shape,
## architecture, design, N, c_bm), grouped by carryover half-life
## along the x axis. Colour encodes power.
##
## Output: analysis/figures/02-lasagna-power.pdf

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

repo_root <- here::here()

grid <- readRDS(file.path(repo_root,
  'analysis/data/02-grid-summary.rds'))$summary

fig_dir <- file.path(repo_root, 'analysis/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Build a stable cell ordering. Within a panel, we group cells by
## (decay form, architecture, design, N, c_bm); within each group,
## the t1half axis runs along x.
g <- grid |>
  mutate(
    decay_label = case_when(
      carryover_form == 'linear'      ~ 'linear',
      carryover_form == 'exponential' ~ 'exponential',
      carryover_form == 'weibull'     ~ sprintf('Weibull k=%.1f',
                                                 weibull_shape),
      TRUE ~ carryover_form
    ),
    arch_label = factor(dgp_arch,
      levels = c('mean_moderation', 'mvn'),
      labels = c('Architecture A', 'Architecture B')),
    cell_id = paste(decay_label, arch_label, design,
                    sprintf('N=%d', N),
                    sprintf('c_bm=%.2f', c_bm),
                    sep = ' | ')
  ) |>
  arrange(decay_label, arch_label, design, N, c_bm, t1half)

cell_levels <- unique(g$cell_id)
g <- g |>
  mutate(cell_id = factor(cell_id, levels = rev(cell_levels)))

p_lasagna <- ggplot(g, aes(x = factor(t1half), y = cell_id,
                           fill = power)) +
  geom_tile(colour = 'white', linewidth = 0.05) +
  scale_fill_gradientn(
    colours = c('#3b4cc0', '#7ea1f3', '#dddddd', '#f3a07e',
                '#b40426'),
    limits = c(0, 1),
    name   = 'Power'
  ) +
  facet_wrap(~ spec, ncol = 3) +
  labs(
    x = expression(paste('Carryover half-life ', t[1/2],
                         ' (weeks)')),
    y = NULL,
    title = 'Lasagna plot: power across the 540-cell factorial',
    subtitle = paste('Each row is a unique cell (decay form,',
                     'architecture, design, N, c_bm); columns',
                     'within a panel are carryover half-lives;',
                     'panels are analysis specifications.')
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.y = element_text(size = 4),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = 'grey92',
                                     colour = NA),
    legend.position = 'top'
  )

out <- file.path(fig_dir, '02-lasagna-power.pdf')
ggsave(out, p_lasagna, width = 11, height = 14)
message('Wrote ', out)
