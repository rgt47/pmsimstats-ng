## Figures for the power-simulation results in
## docs/36-hendrickson-pd-and-carryover.md.
##
## Reads the summaries written by 05-summarize-power.R for both analysis
## models and writes PNGs to docs/figures/. Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hendrickson-problems/07-power-figures.R
suppressMessages({library(data.table); library(ggplot2)})
dir <- 'analysis/data/quick-sim/hendrickson-problems/power-sim'
fig_dir <- 'docs/figures'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## Categorical slots 1-5 of the reference palette (CVD checks pass;
## three slots are below 3:1 contrast, so every arm also has its own
## point shape, and the whitepaper carries the values in tables).
## Slot 6 (green) is used for F, the current covar form. C and F are the
## separability comparison; their colors pass the validator directly
## (normal-vision dE 15.6), and F is also dashed with its own shape so
## the pair never relies on color alone.
arm_lab <- c(A = 'A published', B = 'B AR(1)', C = 'C AR(1) + graded, separable',
             D = 'D mean moderation', E = 'E all three',
             F = 'F AR(1) + graded, covar')
pal <- setNames(c('#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4',
                  '#008300'), arm_lab)
shp <- setNames(c(16, 17, 15, 18, 25, 1), arm_lab)
lty <- setNames(c('solid', 'solid', 'solid', 'solid', 'solid', 'dashed'), arm_lab)
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'

theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
          axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
          plot.title = element_text(colour = ink, face = 'bold'),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, hjust = 0),
          plot.title.position = 'plot', plot.caption.position = 'plot',
          legend.position = 'top', legend.justification = 'left',
          legend.title = element_blank(),
          legend.text = element_text(colour = ink))
}

read_both <- function(stem) {
  out <- list()
  for (a in c('published', 'car1')) {
    f <- file.path(dir, sprintf('%s-reps500%s.csv', stem,
                                if (a == 'published') '' else '-car1'))
    if (file.exists(f)) out[[a]] <- fread(f)[, analysis := a]
  }
  d <- rbindlist(out)
  d[, model := factor(ifelse(analysis == 'published',
    'Published analysis', 'corCAR1 analysis'),
    levels = c('Published analysis', 'corCAR1 analysis'))]
  d[, arm := factor(arm_lab[arm], levels = arm_lab)]
  d[, design := factor(design, levels = c('Hybrid', 'OL+BDC', 'CO'))]
  d[, t_lab := factor(format(t_half, drop0trailing = TRUE),
                      levels = c('0', '0.1', '0.2', '0.5', '1'))]
  d
}
pw <- read_both('power-summary')
t1 <- read_both('type1')

models_note <- paste0('Rows: published analysis = random intercept only; ',
                      'corCAR1 analysis = the same plus continuous-time AR(1) residuals.')
lines_fig <- function(d, y, ylab, title, subtitle, caption, ylim, href) {
  caption <- paste(caption, models_note, sep = '\n')
  ggplot(d, aes(x = t_lab, y = .data[[y]], colour = arm, shape = arm,
                linetype = arm, group = arm)) +
    href +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.4, fill = 'white', stroke = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = shp) +
    scale_linetype_manual(values = lty) +
    scale_y_continuous(limits = ylim) +
    facet_grid(model ~ design) +
    labs(x = 'Carryover half-life (weeks)', y = ylab, title = title,
         subtitle = subtitle, caption = caption) +
    theme_fig() +
    guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2),
           linetype = guide_legend(nrow = 2))
}

## Figure 1: nominal power at c_bm = 0.25, feasible in every arm.
f1 <- lines_fig(pw[c.bm == 0.25], 'power', 'Power (nominal, alpha = 0.05)',
  'Graded coupling and mean moderation remove the collapse at small carryover',
  'N = 70, c_bm = 0.25 (feasible in every arm), 500 replicates per cell',
  'Steps on the x-axis are evenly spaced for legibility; the half-lives are not.',
  c(0, 1), NULL)
ggsave(file.path(fig_dir, '36-fig2-power.png'), f1, width = 9, height = 6.2,
       dpi = 200, bg = surface)

## Figure 2: Type I error, with a +/- 2 Monte Carlo SE band around 0.05.
band <- list(
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = 0.05 - 2 * sqrt(0.05 * 0.95 / 500),
           ymax = 0.05 + 2 * sqrt(0.05 * 0.95 / 500), fill = '#e6e5e1', alpha = 0.7),
  geom_hline(yintercept = 0.05, colour = ink2, linewidth = 0.4, linetype = 'dashed'))
f2 <- lines_fig(t1, 'type1', 'Type I error',
  'With AR(1) data, the published analysis rejects too often',
  'N = 70, c_bm = 0, 500 replicates per cell',
  'Dashed line: nominal 0.05; shaded band: +/- 2 Monte Carlo standard errors.',
  c(0, 0.25), band)
ggsave(file.path(fig_dir, '36-fig1-type1.png'), f2, width = 9, height = 6.2,
       dpi = 200, bg = surface)

## Figure 3: mean interaction estimate at c_bm = 0.25, against the
## expected slope c_bm * sigma_BR / sigma_bm (sign: symptom reduction).
target <- -0.25 * 8 / 15.36159
f3 <- lines_fig(pw[c.bm == 0.25], 'mean_beta', 'Mean interaction estimate',
  'Under the step rule the estimated interaction collapses toward zero',
  'N = 70, c_bm = 0.25, 500 replicates per cell',
  sprintf('Dashed line: expected slope %.3f = -c_bm * sigma_BR / sigma_bm.', target),
  c(-0.16, 0.01),
  geom_hline(yintercept = target, colour = ink2, linewidth = 0.4,
             linetype = 'dashed'))
ggsave(file.path(fig_dir, '36-fig3-estimate.png'), f3, width = 9, height = 6.2,
       dpi = 200, bg = surface)
cat('Wrote 3 figures; analyses present:',
    paste(levels(droplevels(pw$model)), collapse = ' | '), '\n')
