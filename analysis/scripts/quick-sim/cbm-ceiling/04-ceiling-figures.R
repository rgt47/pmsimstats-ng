## Design-decision figures for docs/34-cbm-feasible-range.md.
##
## Every figure answers one question: for a given design choice, what
## is the largest biomarker-response correlation c_bm the simulation
## can test? Reads the outputs of 01-cbm-ceiling.R,
## 02-ceiling-determinants.R and 03-decay-shape.R and writes PNGs to
## docs/figures/.
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/cbm-ceiling/04-ceiling-figures.R

library(data.table)
library(ggplot2)

in_dir <- 'analysis/data/quick-sim/cbm-ceiling'
fig_dir <- 'docs/figures'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## Categorical slots 1 and 2 of the reference palette, validated for
## CVD separation (protan dE 24.7) and contrast on the light surface.
pal <- c('AR(1), covar' = '#2a78d6', 'CS, orig' = '#eb6834')
shp <- c('AR(1), covar' = 16, 'CS, orig' = 17)
ink <- '#0b0b0b'
ink2 <- '#52514e'
surface <- '#fcfcfb'
ref_c <- 0.45
y_lab <- expression('Largest testable ' * c[bm])

theme_ceiling <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.background = element_rect(fill = surface, colour = NA),
      panel.background = element_rect(fill = surface, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
      axis.text = element_text(colour = ink2),
      axis.title = element_text(colour = ink2),
      strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
      plot.title = element_text(colour = ink, face = 'bold'),
      plot.subtitle = element_text(colour = ink2),
      plot.caption = element_text(colour = ink2, hjust = 0),
      plot.title.position = 'plot',
      plot.caption.position = 'plot',
      legend.position = 'top',
      legend.justification = 'left',
      legend.title = element_blank(),
      legend.text = element_text(colour = ink)
    )
}

## Named in the caption rather than in the panel, where an in-plot
## label collides with points that sit near 0.45.
ref_line <- function() {
  geom_hline(yintercept = ref_c, colour = ink2, linewidth = 0.4,
             linetype = 'dashed')
}
ref_caption <- paste0(
  'Values above a curve cannot be simulated under that structure. ',
  'Dashed line: paper 01 reference value c_bm = 0.45.')

## Per-panel log breaks at the data values, so neighboring panels with
## different maxima do not print crowded labels.
log_breaks_from <- function(vals) {
  function(lims) vals[vals >= lims[1] * 0.999 & vals <= lims[2] * 1.001]
}

save <- function(p, name, w = 7.5, h = 4.3) {
  ggsave(file.path(fig_dir, name), p, width = w, height = h, dpi = 200,
         bg = surface)
}

line_panel <- function(d, x, xlab, title, subtitle, facet = NULL,
                       logx = FALSE, free_x = FALSE) {
  p <- ggplot(d, aes(x = .data[[x]], y = c_star, colour = structure,
                     shape = structure)) +
    ref_line() +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_point(size = 2.6, stroke = 0.8, colour = surface,
               aes(shape = structure), na.rm = TRUE,
               show.legend = FALSE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    scale_colour_manual(values = pal) +
    scale_shape_manual(values = shp) +
    scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.15)) +
    labs(x = xlab, y = y_lab, title = title, subtitle = subtitle,
         caption = ref_caption) +
    theme_ceiling()
  if (logx) {
    vals <- sort(unique(d[[x]]))
    brk <- if (free_x) log_breaks_from(vals) else vals
    p <- p + scale_x_continuous(trans = 'log2', breaks = brk,
      labels = function(b) format(b, drop0trailing = TRUE))
  }
  if (!is.null(facet)) {
    p <- p + facet_wrap(facet, scales = if (free_x) 'free_x' else 'fixed')
  }
  p
}

det <- fread(file.path(in_dir, 'determinants.csv'))
long <- melt(det, measure.vars = c('ar1', 'cs'), variable.name = 'structure',
             value.name = 'c_star')
long[, structure := factor(ifelse(structure == 'ar1', 'AR(1), covar',
                                  'CS, orig'), levels = names(pal))]

## Figure 1: the paper 01 designs, by carryover half-life.
des <- fread(file.path(in_dir, 'design.csv'))
des <- des[(within == 'ar1' & cross == 'decay' & coupling == 'graded') |
           (within == 'cs' & cross == 'const' & coupling == 'step')]
des[, structure := factor(ifelse(within == 'ar1', 'AR(1), covar',
                                 'CS, orig'), levels = names(pal))]
des[, design := factor(design, levels = c('CO', 'Hybrid', 'OL+BDC'))]
p1 <- line_panel(des, 't_half', 'Carryover half-life (weeks)',
  'How large a biomarker correlation can each paper 01 design test?',
  paste0('Architecture B (AR(1)) and the published specification (CS, ',
         'step coupling). The published limit of\n0.84 with carryover is ',
         'not usable: there its coupling carries no on/off contrast'),
  facet = 'design') +
  scale_x_continuous(breaks = c(0, 0.5, 1))
save(p1, '34-fig1-paper01-designs.png', w = 8, h = 4.3)

## Figure 5: number of visits.
d2 <- long[experiment == 'E1 occasions']
d2[, label := sub('spacing ([0-9.]+) wk', '\\1-week spacing', label)]
save(line_panel(d2, 'n', 'Number of visits (log scale)',
  'Adding visits lowers the largest testable correlation',
  'Equal spacing, first half of visits on drug, one switch',
  facet = 'label', logx = TRUE), '34-fig5-visits.png')

## Figure 6: number of on/off switches.
d3 <- long[experiment == 'E2 switches']
d3[, label := sub('spacing ([0-9.]+) wk', '\\1-week spacing', label)]
save(line_panel(d3, 'n_switch', 'Number of on/off switches (log scale)',
  'Under AR(1), each additional switch lowers the testable correlation',
  '16 visits, 8 on drug, equal spacing. Compound symmetry is unaffected',
  facet = 'label', logx = TRUE), '34-fig6-switches.png')

## Figure 7: visit spacing, uniform and at the switch.
d4a <- long[experiment == 'E3 spacing']
d4a[, x := as.numeric(sub('spacing ([0-9.]+) wk', '\\1', label))]
d4a[, panel := 'Uniform spacing (study length varies)']
d4b <- long[experiment == 'E4 switch gap']
d4b[, x := as.numeric(sub('switch gap ([0-9.]+) wk', '\\1', label))]
d4b[, panel := 'Gap across the switch (study length 17.5 wk)']
d4 <- rbind(d4a, d4b)
d4[, panel := factor(panel, levels = unique(panel))]
save(line_panel(d4, 'x', 'Weeks (log scale)',
  'Under AR(1), a short gap at the switch limits the testable correlation',
  paste0('8 visits, first 4 on drug. No AR(1) point at 0.5-week uniform ',
         'spacing: the response block\nis not positive definite. ',
         'Compound symmetry is unaffected by spacing'),
  facet = 'panel', logx = TRUE, free_x = TRUE),
  '34-fig7-spacing.png')

## Figure 8: on-drug fraction.
d5 <- long[experiment == 'E6 on fraction']
save(line_panel(d5, 'n_on', 'On-drug visits (of 8)',
  'More on-drug visits lower the testable correlation under AR(1)',
  paste0('8 visits at 2.5-week spacing, one switch. Under compound ',
         'symmetry the limit is\nlowest near half on drug')) +
  scale_x_continuous(breaks = 1:7),
  '34-fig8-on-fraction.png', w = 6.5)

## Figure 9: lookup grid of spacing by switches (AR(1)).
g <- fread(file.path(in_dir, 'determinants-grid.csv'))
g[, spacing_f := factor(format(spacing, drop0trailing = TRUE),
                        levels = format(sort(unique(spacing)),
                                        drop0trailing = TRUE))]
g[, switch_f := factor(n_switch, levels = sort(unique(n_switch)))]
g[, below := ar1 < ref_c]
p6 <- ggplot(g, aes(x = spacing_f, y = switch_f, fill = ar1)) +
  geom_tile(colour = surface, linewidth = 1.2) +
  geom_tile(data = g[below == FALSE], fill = NA, colour = ink,
            linewidth = 0.9, width = 0.94, height = 0.94) +
  geom_text(aes(label = sprintf('%.2f', ar1),
                colour = ar1 > 0.40), size = 3.4, show.legend = FALSE) +
  scale_fill_gradient(low = '#cde2fb', high = '#104281',
                      limits = c(0.1, 0.55),
                      name = expression('Largest testable ' * c[bm])) +
  scale_colour_manual(values = c(`TRUE` = '#ffffff', `FALSE` = ink)) +
  labs(x = 'Weeks between visits', y = 'On/off switches',
       title = 'Design lookup: largest testable correlation under AR(1)',
       subtitle = paste0('16 visits, 8 on drug, equal spacing, no ',
                         'carryover. Outlined cells can test the paper 01 ',
                         'value 0.45'),
       caption = paste0('Compound symmetry gives 0.25 in every cell. ',
                        'rho = 0.7, c1 = 0.2, cx = 0.1.')) +
  theme_ceiling() +
  theme(panel.grid.major = element_blank(),
        legend.position = 'right',
        legend.title = element_text(colour = ink2, size = 9))
save(p6, '34-fig9-lookup.png', w = 7.5, h = 3.8)

## Figure 10: real schedules against equal spacing (AR(1) only).
d7 <- det[experiment == 'E7 real vs equal']
d7[, path := sub(' (actual|equalized)$', '', label)]
d7[, version := sub('^.* (actual|equalized)$', '\\1', label)]
w7 <- dcast(d7, path ~ version, value.var = 'ar1')
w7[, path := factor(path, levels = path[order(actual)])]
p7 <- ggplot(w7, aes(y = path)) +
  geom_vline(xintercept = ref_c, colour = ink2, linewidth = 0.4,
             linetype = 'dashed') +
  geom_segment(aes(x = actual, xend = equalized, yend = path),
               colour = '#b7d3f6', linewidth = 1.6) +
  geom_point(aes(x = actual, shape = 'Actual schedule'), size = 3,
             colour = pal[[1]], fill = pal[[1]]) +
  geom_point(aes(x = equalized, shape = 'Equal spacing, same study length'),
             size = 3, colour = pal[[1]], fill = surface, stroke = 1.1) +
  geom_text(data = w7[abs(actual - equalized) > 1e-6],
            aes(x = actual, label = sprintf('%.3f', actual)),
            vjust = 2.1, size = 3, colour = ink2) +
  geom_text(aes(x = equalized, label = ifelse(
              abs(actual - equalized) > 1e-6, sprintf('%.3f', equalized),
              sprintf('%.3f (already evenly spaced)', equalized))),
            hjust = -0.25, size = 3, colour = ink2) +
  scale_shape_manual(values = c('Actual schedule' = 21,
                                'Equal spacing, same study length' = 21),
                     guide = guide_legend(override.aes = list(
                       fill = c(pal[[1]], surface)))) +
  scale_x_continuous(limits = c(0.40, 0.80),
                     breaks = seq(0.40, 0.80, 0.05)) +
  labs(x = expression('Largest testable ' * c[bm] * ' under AR(1)'),
       y = NULL,
       title = 'Evening out visit spacing would raise the testable correlation',
       subtitle = paste0('Binding paths of the paper 01 designs, no ',
                         'carryover. Dashed line: reference value 0.45')) +
  theme_ceiling()
save(p7, '34-fig10-real-vs-equal.png', w = 7, h = 3.4)

## Figure 2: sensitivity to rho at t1/2 = 0.
rho <- fread(file.path(in_dir, 'rho.csv'))
rho <- rho[t_half == 0 & config %in% c('ar1/decay/graded',
                                       'cs/const/graded')]
rho[, structure := factor(ifelse(config == 'ar1/decay/graded',
  'AR(1), covar', 'CS, orig'), levels = names(pal))]
rho[, design := factor(design, levels = c('CO', 'Hybrid', 'OL+BDC'))]
p8 <- line_panel(rho, 'rho', expression('Within-factor correlation ' * rho),
  'The assumed serial correlation changes which structure allows more',
  paste0('Paper 01 designs, no carryover, c1 = 0.2, cx = 0.1. ',
         'Missing points: response block not positive definite'),
  facet = 'design') +
  scale_x_continuous(breaks = c(0.3, 0.5, 0.7, 0.9))
save(p8, '34-fig2-rho.png', w = 8, h = 4)

## Figure 3: every 8-visit, 4-on-drug schedule in the analysis. Under
## CS the limit is identical across all of them.
inv <- det[n == 8 & n_on == 4]
inv[, feature := fcase(
  experiment == 'E3 spacing', 'Uniform spacing',
  experiment == 'E4 switch gap', 'Gap at the switch',
  experiment == 'E5 within-block gaps', 'Gaps within blocks',
  experiment == 'E7 real vs equal', 'Real schedules',
  experiment == 'E8 arrangement', 'Order of on-drug visits')]
inv <- inv[!is.na(feature)]
inv[, variant := fcase(
  feature == 'Uniform spacing', sub('spacing (.*) wk', '\\1 weeks apart', label),
  feature == 'Gap at the switch', sub('switch gap (.*) wk', '\\1-week gap', label),
  feature == 'Order of on-drug visits', sub(' Hybrid weeks', ', Hybrid weeks', label),
  default = label)]
inv[, feature := factor(feature, levels = c('Uniform spacing',
  'Gap at the switch', 'Gaps within blocks', 'Order of on-drug visits',
  'Real schedules'))]
inv[, variant := factor(variant, levels = rev(unique(variant)))]
invl <- melt(inv, id.vars = c('feature', 'variant'),
             measure.vars = c('ar1', 'cs'), variable.name = 'structure',
             value.name = 'c_star')
invl[, structure := factor(ifelse(structure == 'ar1', 'AR(1), covar',
                                  'CS, orig'), levels = names(pal))]
p9 <- ggplot(invl, aes(x = c_star, y = variant, colour = structure,
                       shape = structure)) +
  geom_vline(xintercept = ref_c, colour = ink2, linewidth = 0.4,
             linetype = 'dashed') +
  geom_point(size = 2.6, na.rm = TRUE) +
  scale_colour_manual(values = pal) +
  scale_shape_manual(values = shp) +
  scale_x_continuous(limits = c(0.3, 0.75), breaks = seq(0.3, 0.75, 0.05)) +
  facet_grid(feature ~ ., scales = 'free_y', space = 'free_y',
             switch = 'y') +
  labs(x = y_lab, y = NULL,
       title = 'Under compound symmetry, no schedule choice moves the limit',
       subtitle = paste0('Every schedule in the analysis with 8 visits, 4 on ',
                         'drug. CS gives 0.346 in all of them;\nAR(1) ranges ',
                         'from 0.33 to 0.74'),
       caption = paste0('No AR(1) point for uniform 0.5-week spacing: the ',
                        'response block is not positive definite.\n',
                        'Dashed line: paper 01 reference value c_bm = 0.45.')) +
  theme_ceiling() +
  theme(strip.placement = 'outside',
        strip.text.y.left = element_text(angle = 0, hjust = 1,
                                         colour = ink, face = 'bold'),
        panel.spacing.y = unit(0.6, 'lines'))
save(p9, '34-fig3-cs-invariance.png', w = 8.5, h = 6.8)

## Figure 4: CS lookup. Only visit count and on-drug share matter.
csg <- fread(file.path(in_dir, 'determinants-cs-grid.csv'))
csg[, frac_f := factor(sprintf('%d/8', as.integer(round(frac * 8))),
                       levels = sprintf('%d/8', 1:7))]
csg[, n_f := factor(n, levels = sort(unique(n)))]
p10 <- ggplot(csg, aes(x = frac_f, y = n_f, fill = cs)) +
  geom_tile(colour = surface, linewidth = 1.2) +
  geom_tile(data = csg[cs >= ref_c], fill = NA, colour = ink,
            linewidth = 0.9, width = 0.94, height = 0.94) +
  geom_text(aes(label = sprintf('%.2f', cs), colour = cs > 0.30),
            size = 3.4, show.legend = FALSE) +
  scale_fill_gradient(low = '#fde0d2', high = '#8a3210',
                      limits = c(0.1, 0.55),
                      name = expression('Largest testable ' * c[bm])) +
  scale_colour_manual(values = c(`TRUE` = '#ffffff', `FALSE` = ink)) +
  labs(x = 'Share of visits with the coupling on',
       y = 'Number of visits',
       title = 'Compound-symmetry lookup: only visit count and on-drug share matter',
       subtitle = paste0('Exact for any spacing, switch count or order. ',
                         'Outlined cells can test the paper 01 value 0.45'),
       caption = paste0('rho = 0.7, c1 = 0.2, cx = 0.1. Under the step ',
                        'coupling with carryover, every post-exposure ',
                        'visit counts as on.')) +
  theme_ceiling() +
  theme(panel.grid.major = element_blank(),
        legend.position = 'right',
        legend.title = element_text(colour = ink2, size = 9))
save(p10, '34-fig4-cs-lookup.png', w = 8, h = 3.8)

## Figure 11: decay shape (Weibull k, paper 02 family) under AR(1).
## Half-life is ordinal, so it takes a single-hue ramp (validated
## light-to-dark, adjacent dL >= 0.06) with shape as second encoding.
ds <- fread(file.path(in_dir, 'decay-shape.csv'))
k_of <- c(weibull_0.25 = 0.25, weibull_0.5 = 0.5, exponential = 1,
          weibull_2 = 2, weibull_4 = 4)
dsa <- ds[cfg == 'ar1/graded' & shape != 'none']
dsa[, k := k_of[shape]]
dsa[, t_lab := factor(sprintf('t1/2 = %g wk', t_half),
                      levels = sprintf('t1/2 = %g wk', c(0.5, 1, 2)))]
lev <- c('Hybrid', 'OL+BDC', 'Multi-cycle', 'CO')
dsa[, design := factor(design, levels = lev)]
base0 <- ds[cfg == 'ar1/graded' & shape == 'none']
base0[, design := factor(design, levels = lev)]
hl_pal <- setNames(c('#86b6ef', '#2a78d6', '#104281'), levels(dsa$t_lab))
hl_shp <- setNames(c(16, 17, 15), levels(dsa$t_lab))
p11 <- ggplot(dsa, aes(x = k, y = c_star, colour = t_lab, shape = t_lab)) +
  geom_hline(yintercept = ref_c, colour = ink2, linewidth = 0.4,
             linetype = 'dashed') +
  geom_hline(data = base0, aes(yintercept = c_star), colour = ink2,
             linewidth = 0.5, linetype = 'dotted') +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.6, stroke = 0.8, colour = surface,
             show.legend = FALSE) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = hl_pal) +
  scale_shape_manual(values = hl_shp) +
  scale_x_continuous(trans = 'log2', breaks = unname(k_of),
                     labels = c('0.25', '0.5', '1\n(exp.)', '2', '4')) +
  scale_y_continuous(limits = c(0.2, 0.65), breaks = seq(0.2, 0.65, 0.1)) +
  facet_wrap(~ design, nrow = 2) +
  labs(x = 'Weibull shape k (log scale): heavier tail to the left, faster washout to the right',
       y = y_lab,
       title = 'Faster washout lowers the testable correlation toward its no-carryover value',
       subtitle = paste0('AR(1) (covar), graded coupling, paper 02 Weibull ',
                         'family. Multi-cycle: 16 weekly visits, pattern ',
                         '1100 repeated'),
       caption = paste0('Dotted line: limit without carryover. ',
                        'Dashed line: paper 01 reference value c_bm = 0.45. ',
                        'Counterfactual: the package decays the\n',
                        'biomarker correlation exponentially; these ',
                        'curves show the limit if it followed each ',
                        'Weibull shape.')) +
  theme_ceiling()
save(p11, '34-fig11-decay-shape.png', w = 8.5, h = 6)

cat('Wrote 11 figures to', fig_dir, '\n')
