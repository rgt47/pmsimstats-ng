## 16-covariance-rule-figures.R
##
## Redraws Figures 1 and 3 of the closed-form study (docs/37) with the two
## covariance rules only (published step coupling and graded coupling), for
## paper 01, whose scope is the covariance-moderation specification. The
## data are the stored outputs of 01-e9-closed-form.R, so nothing is
## recomputed; the full four-rule figures (37-fig1, 37-fig3) are left in
## place because paper 04 and docs/37 use them.
##
## Usage: Rscript analysis/scripts/quick-sim/carryover-closed-form/16-covariance-rule-figures.R
suppressMessages({library(data.table); library(ggplot2)})
dat_dir <- 'analysis/data/quick-sim/carryover-closed-form'
fig_dir <- 'docs/figures'
CBM <- 0.25
ALPHA <- 0.05
mc_reps <- 5000L

ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
rule_lab <- c(step = 'Step coupling (published)', graded = 'Graded coupling')
pal <- setNames(c('#2a78d6', '#1baf7a'), rule_lab)
lty <- setNames(c('solid', 'solid'), rule_lab)
var_lev <- c(E9b = 'E9b: baseline counted off drug (as in the published lmer)',
             E9 = 'E9: post-baseline visits only')

theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
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

x0 <- 0.018
prep <- function(dt) {
  dt <- dt[rule %in% names(rule_lab)]
  dt[, rule := factor(rule_lab[rule], levels = rule_lab)]
  dt[, variant := factor(var_lev[variant], levels = var_lev)]
  dt[, x := ifelse(t_half == 0, x0, t_half)]
  dt
}

curves <- fread(file.path(dat_dir, 'e9-closed-form-curves.csv'))
val <- fread(file.path(dat_dir, 'e9-closed-form-vs-mc.csv'))
fd <- prep(curves[construct == 'CS'])
pts <- prep(val[c.bm == CBM & design == 'Hybrid' & construct == 'CS'])

g <- ggplot(fd[t_half > 0], aes(x, power, colour = rule, linetype = rule)) +
  geom_hline(yintercept = ALPHA, colour = ink2, linewidth = 0.4,
             linetype = 'dashed') +
  geom_line(linewidth = 0.9) +
  geom_point(data = fd[t_half == 0], size = 2.6, show.legend = FALSE) +
  geom_point(data = pts, aes(y = mc_power), shape = 21, fill = 'white',
             size = 2.2, stroke = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  scale_linetype_manual(values = lty) +
  scale_x_log10(breaks = c(x0, 0.03, 0.1, 1, 4),
                labels = c('0', '0.03', '0.1', '1', '4')) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ variant, nrow = 1) +
  labs(x = 'Carryover half-life (weeks, log scale; 0 at the left edge)',
       y = 'Power (alpha = 0.05)',
       title = 'Power of the paired-difference statistic against carryover half-life',
       subtitle = sprintf(
         'Hybrid design, N = 70, c_bm = %.2f, compound-symmetry construct', CBM),
       caption = paste('Lines: closed form. Filled points: closed form at t1/2 = 0.',
                       sprintf('Open points: Monte Carlo, %d replicates. Dashed line: alpha.',
                               mc_reps), sep = '\n')) +
  theme_fig() +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))
out1 <- file.path(fig_dir, '37-fig1-power-closed-form-cov.png')
ggsave(out1, g, width = 10, height = 5.1, dpi = 200, bg = surface)
cat('Wrote', out1, '\n')

hm <- fread(file.path(dat_dir, 'e9-power-heatmap.csv'))
hd <- hm[rule %in% names(rule_lab)]
hd[, rule := factor(rule_lab[rule], levels = rule_lab)]
hd[, variant := factor(variant, levels = c('E9b', 'E9'),
                       labels = c('E9b (as published)', 'E9'))]
hd[, t_lab := factor(format(t_half, drop0trailing = TRUE),
                     levels = c('0', '0.1', '0.2', '0.5', '1', '2'))]
hd[, c_lab := factor(format(c.bm, nsmall = 1),
                     levels = c('0.0', '0.1', '0.2', '0.3'))]
hd[, label := ifelse(feasible, sprintf('%.2f', power), 'not PD')]
hd[, ink_on := ifelse(!feasible | power < 0.55, ink, '#ffffff')]
g3 <- ggplot(hd, aes(t_lab, c_lab)) +
  geom_tile(data = hd[feasible == TRUE], aes(fill = power), colour = surface,
            linewidth = 1) +
  geom_tile(data = hd[feasible == FALSE], fill = '#d9d8d4', colour = surface,
            linewidth = 1) +
  geom_text(aes(label = label, colour = ink_on), size = 2.8) +
  scale_colour_identity() +
  scale_fill_gradient(low = '#eaf1fb', high = '#1a4f9c', limits = c(0, 1),
                      name = 'Power') +
  facet_grid(variant ~ rule) +
  labs(x = 'Carryover half-life (weeks)',
       y = expression(paste('Strength: ', c[bm])),
       title = 'Power by interaction strength and carryover half-life, per rule',
       subtitle = 'Hybrid design, N = 70, compound-symmetry construct, closed form',
       caption = paste('Gray cells: c_bm at or above the CS ceiling, so the correlation matrix',
                       'is not positive definite and no distribution exists.',
                       'Zero strength gives the test size.')) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = surface, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_text(colour = ink2),
        axis.title = element_text(colour = ink2),
        strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
        plot.title = element_text(colour = ink, face = 'bold'),
        plot.subtitle = element_text(colour = ink2),
        plot.caption = element_text(colour = ink2, hjust = 0),
        plot.title.position = 'plot', plot.caption.position = 'plot',
        legend.title = element_text(colour = ink2),
        legend.text = element_text(colour = ink2))
out3 <- file.path(fig_dir, '37-fig3-power-heatmap-cov.png')
ggsave(out3, g3, width = 9, height = 5.6, dpi = 200, bg = surface)
cat('Wrote', out3, '\n')
