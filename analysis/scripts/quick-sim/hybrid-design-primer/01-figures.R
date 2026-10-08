## Figures and supporting numbers for docs/40 (the Hendrickson Hybrid
## design: response components, the modified Gompertz model, and the
## analysis model). Everything is computed from the published construct
## as reproduced in 04-power-simulation.R (58b32a9 parameters), plus the
## Monte Carlo summaries of carryover-closed-form/07-time-adjustment.R.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hybrid-design-primer/01-figures.R
suppressMessages({library(data.table); library(ggplot2); library(splines)})
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
ta_dir <- 'analysis/data/quick-sim/carryover-closed-form/time-adjustment'
fig_dir <- 'docs/figures'
out_dir <- 'analysis/data/quick-sim/hybrid-design-primer'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
gomp <- dgp$gomp
GP <- list(BR = c(m = 10.98604, d = 5, r = 0.42), TV = c(m = 6.50647, d = 5, r = 0.35),
           PB = c(m = 6.50647, d = 5, r = 0.35))

ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'; grid_col <- '#e6e5e1'
## Validated categorical slots (validate_palette.js: CVD and normal-vision
## checks pass; two slots are below 3:1 contrast, so every series is
## direct-labeled or has its own line type, and the values are tabulated).
comp_pal <- c(BR = '#2a78d6', PB = '#eda100', TV = '#1baf7a', Total = ink)
comp_lty <- c(BR = 'solid', PB = 'dashed', TV = 'dotdash', Total = 'solid')
ramp3 <- c('#9ec5f0', '#2a78d6', '#1a4f9c')
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
          plot.title = element_text(colour = ink, face = 'bold'),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, hjust = 0),
          plot.title.position = 'plot', plot.caption.position = 'plot',
          legend.position = 'top', legend.justification = 'left',
          legend.title = element_blank(), legend.text = element_text(colour = ink))
}
save_fig <- function(g, name, w, h) {
  f <- file.path(fig_dir, name)
  ggsave(f, g, width = w, height = h, dpi = 200, bg = surface)
  cat('Wrote', f, '\n')
}

## ---------------- Gompertz landmarks ----------------
## Modified Gompertz g(t) = m (exp(-d e^{-rt}) - e^{-d}) / (1 - e^{-d}).
landmarks <- rbindlist(lapply(names(GP), function(k) {
  p <- GP[[k]]
  t_inf <- log(p[['d']]) / p[['r']]
  t_half <- uniroot(function(t) gomp(t, p[['m']], p[['d']], p[['r']]) - p[['m']] / 2,
                    c(0, 50), tol = 1e-10)$root
  t_90 <- uniroot(function(t) gomp(t, p[['m']], p[['d']], p[['r']]) - 0.9 * p[['m']],
                  c(0, 80), tol = 1e-10)$root
  data.table(component = k, m = p[['m']], d = p[['d']], r = p[['r']],
             t_inflection = t_inf,
             frac_at_inflection = gomp(t_inf, p[['m']], p[['d']], p[['r']]) / p[['m']],
             t_half_max = t_half, t_90 = t_90,
             at_wk4 = gomp(4, p[['m']], p[['d']], p[['r']]) / p[['m']],
             at_wk8 = gomp(8, p[['m']], p[['d']], p[['r']]) / p[['m']],
             at_wk1 = gomp(1, p[['m']], p[['d']], p[['r']]) / p[['m']])
}))
cat('\nGompertz landmarks (weeks; fractions of m)\n')
print(landmarks[component != 'PB'], digits = 3)
fwrite(landmarks, file.path(out_dir, 'gompertz-landmarks.csv'))

## ---------------- Figure 1: the design ----------------
hyb <- dgp$designs$Hybrid
wk <- cumsum(hyb$t_wk)
path_lab <- c('Path 1: stop after wk 10, drug then placebo',
              'Path 2: stop after wk 10, placebo then drug',
              'Path 3: stop after wk 9, drug then placebo',
              'Path 4: stop after wk 9, placebo then drug')
sched <- rbindlist(lapply(seq_along(hyb$paths), function(p) {
  od <- hyb$paths[[p]]
  ## each visit closes an interval; the drug state applies to the interval
  data.table(path = factor(path_lab[p], levels = rev(path_lab)),
             start = c(0, head(wk, -1)), end = wk, week = wk,
             state = ifelse(od == 1, 'On prazosin', 'Placebo'),
             blind = ifelse(hyb$e == 1, 'Open label', 'Blinded'))
}))
phase <- data.table(xmin = c(0, 8, 12), xmax = c(8, 12, 20),
                    lab = c('Open label (e = 1)', 'Blinded discontinuation (e = 0.5)',
                            'Blinded crossover (e = 0.5)'))
g1 <- ggplot(sched) +
  geom_rect(data = phase, aes(xmin = xmin, xmax = xmax, ymin = 0.4, ymax = 4.6),
            fill = c('#f1efe9', '#e9eef6', '#eef4ee'), colour = NA, inherit.aes = FALSE) +
  geom_text(data = phase, aes(x = (xmin + xmax) / 2, y = 4.85, label = lab),
            colour = ink2, size = 3.3, inherit.aes = FALSE) +
  geom_rect(aes(xmin = start + 0.08, xmax = end - 0.08,
                ymin = as.numeric(path) - 0.32, ymax = as.numeric(path) + 0.32,
                fill = state), colour = NA) +
  geom_point(aes(x = week, y = as.numeric(path)), shape = 21, size = 2.6,
             fill = surface, colour = ink, stroke = 0.6) +
  geom_point(data = data.table(week = 0, y = 1:4), aes(x = week, y = y), shape = 21,
             size = 2.6, fill = ink, colour = ink) +
  scale_fill_manual(values = c('On prazosin' = '#2a78d6', 'Placebo' = '#d9d8d4')) +
  scale_x_continuous(breaks = c(0, wk), minor_breaks = NULL, expand = expansion(add = 0.4)) +
  scale_y_continuous(breaks = 1:4, labels = rev(path_lab), limits = c(0.4, 5.05)) +
  labs(x = 'Week (circles: assessment visits; filled circle: baseline)', y = NULL,
       title = 'The Hybrid design: four paths from two randomizations',
       subtitle = paste('Every path is on drug through week 9. Randomization 1 sets when the drug',
                        'stops (after week 10 or week 9);\nrandomization 2 sets the crossover',
                        'order. N = 70 is split 18, 18, 17, 17 in the published code.')) +
  theme_fig() + theme(panel.grid.major.y = element_blank())
save_fig(g1, '40-fig1-hybrid-design.png', 10, 4.4)

## ---------------- Figure 2: Gompertz anatomy ----------------
tt <- seq(0, 20, by = 0.05)
p <- GP$BR
anat <- data.table(t = tt, y = gomp(tt, p[['m']], p[['d']], p[['r']]))
lmB <- landmarks[component == 'BR']
gA <- ggplot(anat, aes(t, y)) +
  geom_hline(yintercept = p[['m']], colour = ink2, linetype = 'dashed', linewidth = 0.4) +
  geom_line(colour = comp_pal[['BR']], linewidth = 1) +
  geom_point(data = data.table(t = c(lmB$t_inflection, lmB$t_half_max, lmB$t_90),
                               y = gomp(c(lmB$t_inflection, lmB$t_half_max, lmB$t_90),
                                        p[['m']], p[['d']], p[['r']])),
             colour = ink, size = 2) +
  annotate('text', x = 19.8, y = p[['m']] + 0.45, hjust = 1, size = 3.2, colour = ink2,
           label = sprintf('ceiling m = %.2f', p[['m']])) +
  annotate('text', x = lmB$t_inflection + 0.4, y = lmB$frac_at_inflection * p[['m']] - 0.2,
           hjust = 0, size = 3.2, colour = ink,
           label = sprintf('inflection: ln(d)/r = %.1f wk, %.0f%% of m',
                           lmB$t_inflection, 100 * lmB$frac_at_inflection)) +
  annotate('text', x = lmB$t_half_max + 0.4, y = p[['m']] / 2 + 0.1, hjust = 0, size = 3.2,
           colour = ink, label = sprintf('half of m at %.1f wk', lmB$t_half_max)) +
  annotate('text', x = lmB$t_90 + 0.4, y = 0.9 * p[['m']] - 0.7, hjust = 0, size = 3.2,
           colour = ink, label = sprintf('90%% of m at %.1f wk', lmB$t_90)) +
  labs(x = 'Weeks on drug', y = 'Improvement (points)',
       title = 'A. The published BR curve',
       subtitle = 'm = 10.99, d = 5, r = 0.42 per week') +
  theme_fig()
vary <- function(which, vals, lab) {
  rbindlist(lapply(seq_along(vals), function(k) {
    q <- p; q[[which]] <- vals[k]
    data.table(t = tt, y = gomp(tt, q[['m']], q[['d']], q[['r']]),
               level = factor(sprintf('%s = %s', lab, as.character(vals[k])),
                              levels = sprintf('%s = %s', lab, as.character(vals))))
  }))
}
## The curves of panels C and D converge at the right edge, so identity is
## carried by an in-panel legend rather than end labels.
panel <- function(dt, title) {
  ggplot(dt, aes(t, y, colour = level)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = ramp3) +
    coord_cartesian(ylim = c(0, 12)) +
    labs(x = 'Weeks on drug', y = 'Improvement (points)', title = title) +
    theme_fig() +
    theme(legend.position = 'inside', legend.position.inside = c(0.98, 0.04),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = surface, colour = grid_col))
}
gB <- panel(vary('m', c(5, 8, 11), 'm'), 'B. Ceiling m: scales the curve')
gC <- panel(vary('d', c(2, 5, 10), 'd'), 'C. Displacement d: delays the rise')
gD <- panel(vary('r', c(0.2, 0.35, 0.5), 'r'), 'D. Rate r: compresses time')
g2 <- gridExtra::arrangeGrob(gA, gB, gC, gD, ncol = 2,
                             top = grid::textGrob('The modified Gompertz curve',
                                                  x = 0.01, hjust = 0,
                                                  gp = grid::gpar(fontface = 'bold', fontsize = 14)))
ggsave(file.path(fig_dir, '40-fig2-gompertz-anatomy.png'), g2, width = 11, height = 8,
       dpi = 200, bg = surface)
cat('Wrote', file.path(fig_dir, '40-fig2-gompertz-anatomy.png'), '\n')

## ---------------- component means by path ----------------
comp_means <- function(th) {
  rbindlist(lapply(seq_along(hyb$paths), function(k) {
    pi <- dgp$path_info(hyb$t_wk, hyb$e, hyb$paths[[k]])
    b <- dgp$build(pi, 'A', 0, th)
    n <- b$n
    data.table(path = k, week = c(0, pi$w), on = c(FALSE, pi$on), e = c(0, pi$e),
               tod = c(0, pi$tod), tsd = c(0, pi$tsd), tpb = c(0, pi$tpb),
               TV = c(0, b$mu[3:(2 + n)]), PB = c(0, b$mu[(3 + n):(2 + 2 * n)]),
               BR = c(0, b$mu[(3 + 2 * n):(2 + 3 * n)]))
  }))
}
cm <- rbindlist(lapply(c(0, 1), function(th) cbind(t_half = th, comp_means(th))))
cm[, Total := TV + PB + BR]
fwrite(cm, file.path(out_dir, 'hybrid-component-means.csv'))
cat('\nHybrid component means, no carryover\n')
print(dcast(melt(cm[t_half == 0], id.vars = c('t_half', 'path', 'week'),
                 measure.vars = c('TV', 'PB', 'BR', 'Total')),
            variable + path ~ week, value.var = 'value'), digits = 3)
cat('\nBR means with t1/2 = 1 week\n')
print(dcast(cm[t_half == 1], path ~ week, value.var = 'BR'), digits = 3)

## ---------------- Figure 3: component trajectories ----------------
plab3 <- c('Path 1', 'Path 2', 'Path 3', 'Path 4')
long3 <- melt(cm, id.vars = c('t_half', 'path', 'week', 'on'),
              measure.vars = c('BR', 'PB', 'TV', 'Total'), variable.name = 'component')
long3[, component := factor(component, levels = c('Total', 'BR', 'PB', 'TV'))]
long3[, path_f := factor(plab3[path], levels = plab3)]
long3[, row := factor(ifelse(t_half == 0, 'No carryover', 'Carryover half-life 1 week'),
                      levels = c('No carryover', 'Carryover half-life 1 week'))]
bars <- unique(cm[, .(t_half, path, week, on)])
bars[, `:=`(start = shift(week, fill = 0), end = week), by = .(t_half, path)]
bars <- bars[on == TRUE & week > 0]
bars[, path_f := factor(plab3[path], levels = plab3)]
bars[, row := factor(ifelse(t_half == 0, 'No carryover', 'Carryover half-life 1 week'),
                     levels = levels(long3$row))]
lab3 <- long3[path == 1 & week == 20 & t_half == 0]
g3 <- ggplot(long3, aes(week, value, colour = component, linetype = component)) +
  geom_rect(data = bars, aes(xmin = start, xmax = end, ymin = -1.6, ymax = -0.8),
            fill = ink2, inherit.aes = FALSE) +
  geom_line(aes(linewidth = component == 'Total')) +
  geom_point(size = 1.3) +
  scale_linewidth_manual(values = c(`TRUE` = 1.1, `FALSE` = 0.7), guide = 'none') +
  scale_colour_manual(values = comp_pal, labels = c(Total = 'Total improvement',
                                                    BR = 'BR (drug)', PB = 'PB (expectancy)',
                                                    TV = 'TV (natural course)')) +
  scale_linetype_manual(values = comp_lty, labels = c(Total = 'Total improvement',
                                                      BR = 'BR (drug)', PB = 'PB (expectancy)',
                                                      TV = 'TV (natural course)')) +
  facet_grid(row ~ path_f) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20)) +
  labs(x = 'Week', y = 'Mean improvement from baseline (points)',
       title = 'Expected component trajectories in the four Hybrid paths',
       subtitle = 'Published parameters; dark bars mark on-drug intervals. TV and PB are identical in every path; only BR differs.') +
  theme_fig()
save_fig(g3, '40-fig3-component-trajectories.png', 11, 6.6)

## ---------------- the interaction on the trajectory template ----------------
## Conditional on the biomarker, only BR moves. Under covariance
## moderation, Cor(B, BR_t) = c_bm g_t, so
##   E(BR_t | b) = mu_t + c_bm g_t sd_BR b      (b = standardized biomarker);
## under mean moderation the shift beta_bm sd_BR b h_t is added directly.
## The vertical gap between b = +1 and b = -1 is 2 sd_BR (strength) g_t,
## and the E9 interaction slope of a path is
##   beta_bm:D = -(mean gap on drug - mean gap off drug) / (2 sd_bm).
STR <- 0.25; sd_BR <- 8; sd_bm <- dgp$bm_sd
prof <- function(pi, rule, th) {
  br_mu <- dgp$build(pi, 'A', 0, th)$mu[(3 + 2 * length(pi$w)):(2 + 3 * length(pi$w))]
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, length(pi$w))
  phi[pi$on | pi$tsd <= 0] <- 0
  switch(rule,
    step = as.numeric(br_mu != 0),
    graded = ifelse(pi$on, 1, phi),
    mean = as.numeric(pi$on))
}
scen <- data.table(
  sid = c('step0', 'step01', 'graded1', 'mean1'),
  rule = c('step', 'step', 'graded', 'mean'), th = c(0, 0.1, 1, 1),
  lab = c('Step coupling, no carryover', 'Step coupling, half-life 0.1 week',
          'Graded coupling, half-life 1 week', 'Mean moderation, half-life 1 week'))
cond <- rbindlist(lapply(seq_len(nrow(scen)), function(s) {
  sc <- scen[s]
  rbindlist(lapply(seq_along(hyb$paths), function(k) {
    pi <- dgp$path_info(hyb$t_wk, hyb$e, hyb$paths[[k]])
    b <- dgp$build(pi, 'A', 0, sc$th)
    n <- b$n
    g <- prof(pi, sc$rule, sc$th)
    base <- data.table(week = pi$w, on = pi$on, g = g,
                       TV = b$mu[3:(2 + n)], PB = b$mu[(3 + n):(2 + 2 * n)],
                       BR = b$mu[(3 + 2 * n):(2 + 3 * n)])
    cbind(scenario = sc$sid, path = k, rbind(
      data.table(week = 0, on = FALSE, g = 0, TV = 0, PB = 0, BR = 0), base))
  }))
}))
cond[, `:=`(hi = TV + PB + BR + STR * sd_BR * g, lo = TV + PB + BR - STR * sd_BR * g,
            gap = 2 * STR * sd_BR * g)]
est <- cond[week > 0, .(gap_on = mean(gap[on]), gap_off = mean(gap[!on])),
            by = .(scenario, path)]
est[, beta := -(gap_on - gap_off) / (2 * sd_bm)]
cat('\nE9 interaction slope by path, from the biomarker gap (strength 0.25)\n')
print(dcast(est, scenario ~ path, value.var = 'beta'), digits = 3)
fwrite(est, file.path(out_dir, 'interaction-gap-by-path.csv'))

## Figure 4: the biomarker fan on the template (published rule, no carryover)
f4 <- cond[scenario == 'step0']
f4[, path_f := factor(plab3[path], levels = plab3)]
fan_long <- melt(f4, id.vars = c('path_f', 'week'), measure.vars = c('hi', 'lo'),
                 variable.name = 'grp')
fan_long[, grp := factor(ifelse(grp == 'hi', 'Biomarker 1 SD above the mean',
                                'Biomarker 1 SD below the mean'),
                         levels = c('Biomarker 1 SD above the mean',
                                    'Biomarker 1 SD below the mean'))]
nd4 <- melt(f4, id.vars = c('path_f', 'week'), measure.vars = c('TV', 'PB'),
            variable.name = 'component')
bars4 <- bars[t_half == 0]
g4f <- ggplot() +
  geom_rect(data = bars4, aes(xmin = start, xmax = end, ymin = -1.6, ymax = -0.8),
            fill = ink2) +
  geom_ribbon(data = f4, aes(week, ymin = lo, ymax = hi), fill = '#2a78d6', alpha = 0.18) +
  geom_line(data = nd4, aes(week, value, linetype = component), colour = '#9b9a96',
            linewidth = 0.5) +
  geom_line(data = fan_long, aes(week, value, colour = grp), linewidth = 0.9) +
  geom_point(data = fan_long, aes(week, value, colour = grp), size = 1.4) +
  scale_colour_manual(values = c('Biomarker 1 SD above the mean' = '#1a4f9c',
                                 'Biomarker 1 SD below the mean' = '#7fb2ea')) +
  scale_linetype_manual(values = c(TV = 'dotdash', PB = 'dashed'),
                        labels = c(TV = 'TV (grey, for reference)', PB = 'PB (grey)')) +
  facet_wrap(~path_f, nrow = 1) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20)) +
  labs(x = 'Week', y = 'Expected total improvement (points)',
       title = 'The interaction on the trajectory template: expected improvement by biomarker level',
       subtitle = paste('Published coupling rule, no carryover, c_bm = 0.25. The band is the',
                        'difference the biomarker makes: open on drug, closed off drug.\nThe',
                        'interaction is how much wider the band is on drug than off drug.')) +
  theme_fig()
save_fig(g4f, '40-fig4-interaction-fan.png', 11, 4.6)

## Figure 5: the estimand as the on-minus-off gap, by coupling rule (path 1)
f5 <- cond[path == 1 & week > 0]
f5 <- merge(f5, scen[, .(scenario = sid, lab)], by = 'scenario')
f5[, lab := factor(lab, levels = scen$lab)]
f5[, state := factor(ifelse(on, 'On drug', 'Off drug'), levels = c('On drug', 'Off drug'))]
seg <- f5[, .(m = mean(gap)), by = .(lab, state)]
e1 <- merge(est[path == 1], scen[, .(scenario = sid, lab)], by = 'scenario')
e1[, lab := factor(lab, levels = scen$lab)]
e1[, txt := sprintf('on - off = %.2f\nbeta_bm:D = %.3f', gap_on - gap_off,
                    ifelse(abs(beta) < 5e-4, 0, beta))]
## Zero-width visits carry no bar, so every visit also gets a marker.
g5f <- ggplot(f5, aes(week, gap, fill = state)) +
  geom_col(width = 0.7) +
  geom_point(aes(colour = state, shape = state), size = 2, show.legend = FALSE) +
  scale_shape_manual(values = c('On drug' = 16, 'Off drug' = 21)) +
  geom_hline(data = seg, aes(yintercept = m, colour = state, linetype = state),
             linewidth = 0.7) +
  geom_text(data = e1, aes(x = 0.3, y = 5.3, label = txt), hjust = 0, vjust = 1,
            size = 3.2, colour = ink, inherit.aes = FALSE) +
  facet_wrap(~lab, nrow = 1) +
  scale_fill_manual(values = c('On drug' = '#2a78d6', 'Off drug' = '#bdbcb7')) +
  scale_colour_manual(values = c('On drug' = '#1a4f9c', 'Off drug' = '#52514e')) +
  scale_linetype_manual(values = c('On drug' = 'solid', 'Off drug' = 'dashed')) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20)) +
  coord_cartesian(ylim = c(0, 5.5)) +
  labs(x = 'Week (path 1)', y = 'Width of the biomarker band (points)',
       title = 'The interaction estimand is a difference in band width, on drug minus off drug',
       subtitle = paste('Bars: improvement at +1 SD minus improvement at -1 SD of the biomarker,',
                        'at each visit; lines: their on-drug and off-drug means.\nThe E9',
                        'slope is -(on - off) / (2 sd_bm). Under the step rule, any carryover',
                        'makes the band equally wide off drug, so the estimand is zero.')) +
  theme_fig()
save_fig(g5f, '40-fig5-interaction-estimand.png', 12, 4.8)

## ---------------- Figure 4: what the time terms can absorb ----------------
## Population mean of the outcome, E(Y) = mu_BL - [TV + PB + BR], at every
## visit of every path, weighted by the published allocation. The
## nuisance part TV + PB is a function of the visit alone (all paths share
## the schedule and the expectancy pattern), so factor(t) absorbs it
## exactly. OLS on the population means gives the probability limit of
## the fixed effects under an independence working model; the random-
## intercept GLS of the mixed model differs, so these are illustrations.
Ns <- c(18, 18, 17, 17)
pm <- cm[t_half == 0]
pm[, `:=`(w = Ns[path], Db = as.numeric(on), De = e, Ey = -(TV + PB + BR))]
fits <- list(
  'Linear t' = Ey ~ Db + week,
  'Linear t + De' = Ey ~ Db + week + De,
  'Spline in t (3 df)' = Ey ~ Db + ns(week, df = 3),
  'Visit as factor' = Ey ~ Db + factor(week))
plim <- rbindlist(lapply(c(0, 1), function(th) {
  d <- cm[t_half == th]
  d[, `:=`(w = Ns[path], Db = as.numeric(on), De = e, Ey = -(TV + PB + BR))]
  rbindlist(lapply(names(fits), function(f) {
    m <- lm(fits[[f]], data = d, weights = w)
    data.table(t_half = th, model = f, Db = coef(m)[['Db']],
               resid_rms = sqrt(weighted.mean(resid(m)^2, d$w)))
  }))
}))
cat('\nPopulation-mean (OLS) Db coefficient by time adjustment\n')
print(plim, digits = 3)
fwrite(plim, file.path(out_dir, 'time-adjustment-plim.csv'))
## Fitted against true population means, path by path. The time terms
## absorb whatever the binary Db cannot represent at visits common to all
## paths, including the Gompertz rise of BR during open label, so the
## comparison is made on the whole mean, not on TV + PB alone.
path_pal <- c('Path 1' = '#2a78d6', 'Path 2' = '#1baf7a', 'Path 3' = '#eda100',
              'Path 4' = '#e87ba4')
mfit <- rbindlist(lapply(names(fits), function(f) {
  m <- lm(fits[[f]], data = pm, weights = w)
  data.table(model = f, path = pm$path, week = pm$week, truth = pm$Ey, fitted = fitted(m))
}))
rms <- plim[t_half == 0, .(model, lab = sprintf('%s (RMS misfit %.2f)', model, resid_rms))]
mfit <- merge(mfit, rms, by = 'model')
mfit[, lab := factor(lab, levels = rms$lab)]
mfit[, path_f := factor(paste('Path', path), levels = names(path_pal))]
g4 <- ggplot(mfit, aes(week, colour = path_f)) +
  geom_line(aes(y = fitted, linetype = path_f), linewidth = 0.7) +
  geom_point(aes(y = truth, shape = path_f), size = 1.9) +
  facet_wrap(~lab, ncol = 2) +
  scale_colour_manual(values = path_pal) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_linetype_manual(values = c('solid', 'dashed', 'dotted', 'dotdash')) +
  scale_x_continuous(breaks = c(0, 4, 8, 12, 16, 20)) +
  labs(x = 'Week', y = 'E(Y) relative to baseline (points)',
       title = 'How each time adjustment reproduces the mean of the published Hybrid data',
       subtitle = paste('Points: true population mean, -(TV + PB + BR), by path. Lines: the',
                        'allocation-weighted least-squares fit of each fixed-effects form with',
                        'binary Db,\nno carryover. Visit as a factor misfits only the',
                        'drug response that one Db coefficient cannot represent.')) +
  theme_fig()
save_fig(g4, '40-fig7-time-adjustment-fit.png', 10.5, 7.4)

## ---------------- Figure 5: the correlation matrix ----------------
pi1 <- dgp$path_info(hyb$t_wk, hyb$e, hyb$paths[[1]])
vlab <- c('B', 'BL', paste0('TV', 1:8), paste0('PB', 1:8), paste0('BR', 1:8))
cmat <- rbindlist(lapply(c(0, 0.1), function(th) {
  R <- dgp$build(pi1, 'A', 0.25, th)$R
  dimnames(R) <- list(vlab, vlab)
  d <- as.data.table(as.table(R))
  setnames(d, c('row', 'col', 'r'))
  d[, panel := sprintf('Path 1, c_bm = 0.25, half-life %s', ifelse(th == 0, '0', '0.1 week'))]
}))
cmat[, row := factor(row, levels = rev(vlab))]
cmat[, col := factor(col, levels = vlab)]
cmat[r == 0, r := NA]
g5 <- ggplot(cmat, aes(col, row, fill = r)) +
  geom_tile(colour = surface, linewidth = 0.4) +
  scale_fill_gradient(low = '#eaf1fb', high = '#1a4f9c', limits = c(0, 1),
                      na.value = '#f7f6f3', name = 'Correlation') +
  facet_wrap(~panel) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       title = 'The published correlation matrix (compound symmetry), one Hybrid path',
       subtitle = paste('Blocks: within-factor 0.8 at every lag; cross-factor 0.2 same visit,',
                        '0.1 otherwise; biomarker row: c_bm at BR visits\nwhose mean is',
                        'nonzero. Blank cells are 0. With any carryover every BR visit is',
                        'coupled (right).')) +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7), legend.position = 'right',
        panel.grid.major = element_blank())
save_fig(g5, '40-fig6-correlation-matrix.png', 12, 6.6)

## ---------------- collinearity of the analysis terms ----------------
## Per-observation design of the published long data (baseline row
## included, t = 0, Db = 0, De = 0), patients weighted by allocation.
X <- cm[t_half == 0, .(w = Ns[path], t = week, Db = as.numeric(on), De = e)]
wcor <- function(a, b, w) cov.wt(cbind(a, b), wt = w, cor = TRUE)$cor[1, 2]
cat(sprintf('\nHybrid design correlations (weighted, baseline included): Db-De %.2f, Db-t %.2f, De-t %.2f\n',
            wcor(X$Db, X$De, X$w), wcor(X$Db, X$t, X$w), wcor(X$De, X$t, X$w)))
r2 <- summary(lm(De ~ Db + t, data = X, weights = w))$r.squared
cat(sprintf('R^2 of De on (Db, t): %.2f; VIF %.1f\n', r2, 1 / (1 - r2)))

## ---------------- Figure 6: the time-adjustment Monte Carlo ----------------
f_alt <- file.path(ta_dir, 'summary-reps500.csv')
f_null <- file.path(ta_dir, 'summary-reps1000-null.csv')
ta <- fread(f_alt)
form_lab <- c(linear = 'Linear t (published)', spline = 'Spline in t (3 df)',
              visit = 'Visit as factor')
rule_lab <- c(step = 'Step coupling (published)', graded = 'Graded coupling',
              mean = 'Mean moderation')
ta[, form_f := factor(form_lab[form], levels = form_lab)]
ta[, rule_f := factor(rule_lab[rule], levels = rule_lab)]
ta[, design_f := factor(design, levels = c('Hybrid', 'CO'))]
form_pal <- c('#2a78d6', '#eda100', '#1baf7a'); names(form_pal) <- form_lab
form_shape <- c(16, 17, 15); names(form_shape) <- form_lab
g6 <- ggplot(ta, aes(factor(t_half), power_int, colour = form_f, shape = form_f)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.55)) +
  geom_linerange(aes(ymin = power_int - 2 * mcse_int, ymax = power_int + 2 * mcse_int),
                 position = position_dodge(width = 0.55), linewidth = 0.5) +
  facet_grid(design_f ~ rule_f) +
  scale_colour_manual(values = form_pal) + scale_shape_manual(values = form_shape) +
  scale_y_continuous(limits = c(0, 0.85)) +
  labs(x = 'Carryover half-life (weeks)', y = 'Power of the bm:Db test',
       title = 'Interaction power by the form of the time adjustment',
       subtitle = paste('Published construct (CS), N = 70, strength 0.25, 500 replicates per cell;',
                        'bars are 2 Monte Carlo SE. All three forms fitted to the same replicates.')) +
  theme_fig()
save_fig(g6, '40-fig8-time-adjustment-power.png', 10.5, 6.2)
if (file.exists(f_null)) {
  nl <- fread(f_null)
  cat('\nType I error by time adjustment\n')
  print(nl[, .(design, t_half, form, type1 = power_int,
               mcse = sqrt(power_int * (1 - power_int) / fits), fails)])
}
cat('\nInteraction power, alternative cells\n')
print(dcast(ta, design + rule + t_half ~ form, value.var = 'power_int'))
