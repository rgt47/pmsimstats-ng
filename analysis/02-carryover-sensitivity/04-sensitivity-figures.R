## analysis/02-carryover-sensitivity/04-sensitivity-figures.R
##
## Render manuscript-02 sensitivity-block figures (Tier 2). One PDF
## per block written to manuscripts/figures/02-sens-{S1..S5}.pdf.
##
## Behaviour when Tier 2 results are not yet available:
## manuscripts/data/02-sensitivity-summary.rds is missing -> each
## block writes a minimal placeholder PDF with the message
## "Pending Tier 2 production run." so that knitr::include_graphics
## calls in the manuscript do not fail.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

repo_root <- here::here()
fig_dir <- file.path(repo_root, 'manuscripts/figures')
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

summary_path <- file.path(repo_root,
  'manuscripts/data/02-sensitivity-summary.rds')

theme_paper <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = 'grey92', colour = NA),
    legend.position = 'top'
  )

placeholder <- function(block, path) {
  p <- ggplot() +
    annotate('text', x = 0.5, y = 0.55,
             label = sprintf('Block %s', block),
             fontface = 'bold', size = 5) +
    annotate('text', x = 0.5, y = 0.45,
             label = 'Pending Tier 2 production run.',
             size = 3.5) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void()
  ggsave(path, p, width = 5, height = 3.2)
}

spec_scale <- scale_colour_manual(
  values = c(A1 = '#1f78b4', A2 = '#33a02c', A3 = '#e31a1c'),
  labels = c(A1 = 'A1 binary',
             A2 = 'A2 Dbc (matched)',
             A3 = 'A3 lagged')
)

if (!file.exists(summary_path)) {
  message('Tier 2 summary not found: ', summary_path)
  message('Writing placeholders for S1-S5.')
  for (blk in c('S1', 'S2', 'S3', 'S4', 'S5')) {
    placeholder(blk,
      file.path(fig_dir, sprintf('02-sens-%s.pdf', blk)))
  }
  quit(save = 'no', status = 0)
}

s2 <- readRDS(summary_path)$summary |>
  dplyr::mutate(spec = factor(spec, levels = c('A1', 'A2', 'A3')))

## -----------------------------------------------------------------
## S1: power vs rho by spec
## -----------------------------------------------------------------
d1 <- s2 |> dplyr::filter(block == 'S1')
p1 <- ggplot(d1, aes(rho, power, colour = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  spec_scale +
  labs(x = expression('AR(1) autocorrelation'~rho),
       y = 'Power', colour = 'Analysis spec',
       title = 'S1: Sensitivity to within-factor autocorrelation',
       subtitle = 'Reference: Hybrid, Arch B, N=70, c_bm=0.45, exp DGP, t_{1/2}=1.0') +
  theme_paper

ggsave(file.path(fig_dir, '02-sens-S1.pdf'), p1,
       width = 5.5, height = 3.8)

## -----------------------------------------------------------------
## S2: analyst-truth half-life mismatch
## -----------------------------------------------------------------
d2 <- s2 |> dplyr::filter(block == 'S2') |>
  dplyr::mutate(true = factor(paste0('true t_{1/2} = ', t1half)))
p2 <- ggplot(d2, aes(analysis_t1half, power,
                     colour = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ true) +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2)) +
  scale_y_continuous(limits = c(0, 1)) +
  spec_scale +
  labs(x = expression('Analyst-assumed half-life (log scale, weeks)'),
       y = 'Power', colour = 'Analysis spec',
       title = 'S2: Cost of analyst-truth half-life mismatch',
       subtitle = 'A1 and A3 do not depend on assumed half-life; A2 does') +
  theme_paper

ggsave(file.path(fig_dir, '02-sens-S2.pdf'), p2,
       width = 7, height = 3.8)

## -----------------------------------------------------------------
## S3: dropout
## -----------------------------------------------------------------
d3 <- s2 |> dplyr::filter(block == 'S3')
p3 <- ggplot(d3, aes(dropout_rate, power,
                     colour = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ dropout_mech) +
  scale_y_continuous(limits = c(0, 1)) +
  spec_scale +
  labs(x = 'Dropout rate',
       y = 'Power', colour = 'Analysis spec',
       title = 'S3: Sensitivity to dropout',
       subtitle = 'MCAR and MAR (severity-biased); reference cell as above') +
  theme_paper

ggsave(file.path(fig_dir, '02-sens-S3.pdf'), p3,
       width = 7, height = 3.8)

## -----------------------------------------------------------------
## S4: effect-size curve
## -----------------------------------------------------------------
d4 <- s2 |> dplyr::filter(block == 'S4')
p4 <- ggplot(d4, aes(c_bm, power, colour = spec, group = spec)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed',
             colour = 'grey50') +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  scale_y_continuous(limits = c(0, 1)) +
  spec_scale +
  labs(x = expression(c[bm]),
       y = 'Power (c_bm > 0) or type-I error (c_bm = 0)',
       colour = 'Analysis spec',
       title = 'S4: Biomarker-moderation effect-size curve',
       subtitle = 'Dashed reference line at nominal alpha = 0.05') +
  theme_paper

ggsave(file.path(fig_dir, '02-sens-S4.pdf'), p4,
       width = 5.5, height = 3.8)

## -----------------------------------------------------------------
## S5: rho x carryover interaction (exploratory)
## -----------------------------------------------------------------
d5 <- s2 |> dplyr::filter(block == 'S5')
p5 <- ggplot(d5, aes(t1half, power, colour = spec, group = spec)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ rho, labeller = label_bquote(rho == .(rho))) +
  scale_y_continuous(limits = c(0, 1)) +
  spec_scale +
  labs(x = expression('Carryover half-life'~t['1/2']~'(weeks)'),
       y = 'Power', colour = 'Analysis spec',
       title = 'S5: Autocorrelation x carryover interaction',
       subtitle = 'Exploratory: is the spec ranking rho-sensitive?') +
  theme_paper

ggsave(file.path(fig_dir, '02-sens-S5.pdf'), p5,
       width = 7, height = 3.8)

message('Wrote Tier 2 figures to ', fig_dir)
