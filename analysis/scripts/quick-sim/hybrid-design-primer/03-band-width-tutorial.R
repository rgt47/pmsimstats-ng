## Figures, tables and checks for docs/76 (the interaction estimand as a
## difference in band width: a tutorial). Expands Section 5 of docs/40,
## its Figure 5 and its Table 9 (per-path interaction slopes).
##
## Everything uses the published construct (compound symmetry,
## rho_cs = 0.8, 58b32a9 parameters) as reproduced in
## hendrickson-problems/04-power-simulation.R, with anchored carryover on
## the mean. Only the biomarker row of the correlation matrix changes
## between the correlation rules. Closed-form results follow docs/37,
## Section 4.4; a Monte Carlo run checks them.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/hybrid-design-primer/03-band-width-tutorial.R
suppressMessages({library(data.table); library(ggplot2); library(MASS)})
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
fig_dir <- 'docs/figures'
out_dir <- 'analysis/data/quick-sim/hybrid-design-primer/band-width'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd', 'gomp', 'designs',
          'path_info', 'build')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))

STR <- 0.25; SD_BR <- 8; SD_BM <- dgp$bm_sd; BM_M <- dgp$bm_m
NS <- c(18L, 18L, 17L, 17L); N <- sum(NS)
hyb <- dgp$designs$Hybrid
PI <- lapply(hyb$paths, function(od) dgp$path_info(hyb$t_wk, hyb$e, od))

ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'; grid_col <- '#e6e5e1'
on_col <- '#2a78d6'; on_dark <- '#1a4f9c'; off_col <- '#bdbcb7'
rule_pal <- c('Constant rule' = '#e0533d', 'Decaying rule' = '#2a78d6',
              'Mean moderation' = '#1baf7a')
path_pal <- c('Path 1' = '#2a78d6', 'Path 2' = '#1baf7a', 'Path 3' = '#eda100',
              'Path 4' = '#e87ba4', 'Pooled' = '#0b0b0b')
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
          axis.text = element_text(colour = ink2),
          axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
          plot.title = element_text(colour = ink, face = 'bold'),
          plot.subtitle = element_text(colour = ink2),
          plot.title.position = 'plot',
          legend.position = 'top', legend.justification = 'left',
          legend.title = element_blank(),
          legend.text = element_text(colour = ink))
}
save_fig <- function(g, name, w, h) {
  f <- file.path(fig_dir, name)
  ggsave(f, g, width = w, height = h, dpi = 200, bg = surface)
  cat('Wrote', f, '\n')
}

## ---------------- the correlation profile g_t ----------------
## constant: g = 1 wherever the drug-response mean is nonzero (the
##   published step rule); decaying: 1 on drug, 2^(-tsd / t_half) off
##   drug after exposure, 0 before; mean: the on-drug indicator, used
##   only for the band (the correlation row stays zero).
profile <- function(pi, rule, th) {
  br_mu <- dgp$build(pi, 'A', 0, th)$mu[(3 + 2 * length(pi$w)):(2 + 3 * length(pi$w))]
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, length(pi$w))
  phi[pi$on | pi$tsd <= 0] <- 0
  switch(rule,
         constant = as.numeric(br_mu != 0),
         decaying = ifelse(pi$on, 1, phi),
         mean = as.numeric(pi$on))
}

## Joint distribution of one path under a rule: mean vector and covariance
## of (B, BL, TV_1..8, PB_1..8, BR_1..8).
path_model <- function(k, rule, th, cbm = STR) {
  pi <- PI[[k]]
  b <- dgp$build(pi, 'A', 0, th)
  n <- b$n
  g <- profile(pi, rule, th)
  sds <- sqrt(diag(b$Sigma))
  R <- b$R
  idx <- (3 + 2 * n):(2 + 3 * n)
  if (rule != 'mean') R[1, idx] <- R[idx, 1] <- cbm * g
  S <- outer(sds, sds) * R
  list(pi = pi, mu = b$mu, S = S, n = n, g = g, idx = idx, rule = rule,
       cbm = cbm, min_eig = min(eigen(S, symmetric = TRUE, only.values = TRUE)$values))
}

## E9 for one path: Delta = mean of Y on drug minus mean of Y off drug over
## the eight post-baseline visits. Y = BL - (TV + PB + BR), so BL cancels.
path_e9 <- function(m) {
  pi <- m$pi
  a <- ifelse(pi$on, 1 / sum(pi$on), -1 / sum(!pi$on))
  L <- c(0, 0, -a, -a, -a)
  v <- drop(t(L) %*% m$S %*% L)
  cv <- sum(L * m$S[, 1])
  if (m$rule == 'mean') {
    v <- v + (m$cbm * SD_BR)^2
    cv <- -m$cbm * SD_BR * SD_BM
  }
  data.table(n_on = sum(pi$on), n_off = sum(!pi$on),
             g_on = mean(m$g[pi$on]), g_off = mean(m$g[!pi$on]),
             mu = sum(L * m$mu), beta = cv / SD_BM^2,
             s2 = v - cv^2 / SD_BM^2, min_eig = m$min_eig)
}

## Pooled E9 (common intercept) and its power: docs/37, Section 4.4.
pow_t <- function(beta, sres, n_tot = N, alpha = 0.05) {
  q <- qchisq((seq_len(400) - 0.5) / 400, n_tot - 1)
  ncp <- beta * sqrt(SD_BM^2 * q) / sres
  tc <- qt(1 - alpha / 2, n_tot - 2)
  mean(pt(tc, n_tot - 2, ncp, lower.tail = FALSE) + pt(-tc, n_tot - 2, ncp))
}
pooled <- function(pp, ns = NS) {
  w <- ns / sum(ns)
  bbar <- sum(w * pp$beta); mbar <- sum(w * pp$mu)
  sres2 <- sum(w * (pp$s2 + (pp$mu - mbar)^2 + SD_BM^2 * (pp$beta - bbar)^2))
  vhat <- sum(w * (pp$s2 + (pp$mu - mbar)^2)) / ((sum(ns) - 3) * SD_BM^2) +
    2 / sum(ns) * sum(w * (pp$beta - bbar)^2)
  se_norm <- sqrt(sres2 / (sum(ns) - 3)) / SD_BM
  data.table(G = sum(w * (pp$g_on - pp$g_off)), beta = bbar, sd = sqrt(vhat),
             sres = sqrt(sres2), power = pow_t(bbar, sqrt(sres2), sum(ns)),
             power_normal = pnorm(abs(bbar) / se_norm - qnorm(0.975)) +
               pnorm(-abs(bbar) / se_norm - qnorm(0.975)))
}
rule_cell <- function(rule, th) {
  pp <- rbindlist(lapply(1:4, function(k) cbind(path = k, path_e9(path_model(k, rule, th)))))
  list(paths = pp, pooled = pooled(pp))
}

## ---------------- the four scenarios of docs/40, Figure 5 ----------------
scen <- data.table(
  sid = c('none', 'const01', 'decay1', 'mean1'),
  rule = c('decaying', 'constant', 'decaying', 'mean'),
  th = c(0, 0.1, 1, 1),
  lab = c('A. No carryover', 'B. Constant rule, half-life 0.1 week',
          'C. Decaying rule, half-life 1 week', 'D. Mean moderation, half-life 1 week'))
cells <- lapply(seq_len(nrow(scen)), function(s) rule_cell(scen$rule[s], scen$th[s]))
tab9 <- rbindlist(lapply(seq_len(nrow(scen)), function(s)
  cbind(scenario = scen$sid[s], cells[[s]]$paths)))
tab9[, G := g_on - g_off]
pool4 <- rbindlist(lapply(seq_len(nrow(scen)), function(s)
  cbind(scenario = scen$sid[s], cells[[s]]$pooled)))
stopifnot(all(tab9$min_eig > 0))
## docs/40 Table 9 and docs/37 (closed-form E9 power and mean slope).
t9 <- dcast(tab9, scenario ~ path, value.var = 'beta')
stopifnot(abs(t9[scenario == 'none', `1`] + 0.130) < 5e-4,
          all(abs(unlist(t9[scenario == 'const01', -1])) < 1e-12),
          abs(t9[scenario == 'decay1', `1`] + 0.095) < 5e-4,
          abs(t9[scenario == 'decay1', `4`] + 0.101) < 5e-4,
          abs(pool4[scenario == 'none', power] - 0.721) < 5e-4,
          abs(pool4[scenario == 'decay1', power] - 0.473) < 5e-4,
          abs(pool4[scenario == 'decay1', beta] + 0.098) < 5e-4,
          abs(pool4[scenario == 'mean1', power] - 0.682) < 1e-3,
          abs(pool4[scenario == 'const01', power] - 0.05) < 1e-9)
cat('\nPer-path E9 quantities (docs/40 Table 9 expanded)\n')
print(tab9, digits = 4)
cat('\nPooled E9, closed form\n')
print(pool4, digits = 4)
fwrite(tab9, file.path(out_dir, 'per-path-e9.csv'))
fwrite(pool4, file.path(out_dir, 'pooled-e9.csv'))

## ---------------- the decaying rule across half-lives ----------------
th_tab <- c(0.1, 0.2, 0.5, 1, 2)
dec_tab <- rbindlist(lapply(th_tab, function(th) {
  cc <- rule_cell('decaying', th)
  rbind(cc$paths[, .(t_half = th, path = as.character(path), g_off, G = g_on - g_off, beta)],
        data.table(t_half = th, path = 'Pooled', g_off = NA_real_, G = cc$pooled$G,
                   beta = cc$pooled$beta))
}))
dec_tab[, n_factor := 1 / G^2]
cat('\nDecaying rule by half-life\n')
print(dcast(dec_tab, t_half ~ path, value.var = 'G'), digits = 4)
print(dcast(dec_tab, t_half ~ path, value.var = 'beta'), digits = 4)
fwrite(dec_tab, file.path(out_dir, 'decaying-by-half-life.csv'))

## Continuous curves for Figures 5 and 8.
th_grid <- sort(unique(c(10^seq(-2, 1, length.out = 61), th_tab)))
curve_G <- rbindlist(lapply(th_grid, function(th) {
  cc <- rule_cell('decaying', th)
  rbind(cc$paths[, .(t_half = th, path = paste('Path', path), G = g_on - g_off)],
        data.table(t_half = th, path = 'Pooled', G = cc$pooled$G))
}))
curve_pow <- rbindlist(lapply(th_grid, function(th) {
  rbindlist(lapply(c('constant', 'decaying', 'mean'), function(r) {
    cc <- rule_cell(r, th)
    data.table(t_half = th, rule = r, power = cc$pooled$power, beta = cc$pooled$beta)
  }))
}))
p0 <- pool4[scenario == 'none', power]
fwrite(curve_pow, file.path(out_dir, 'power-by-half-life.csv'))

## ---------------- Monte Carlo check ----------------
## 5,000 trials of N = 70 (18, 18, 17, 17) per scenario; E9 as specified
## (OLS of Delta on B with a common intercept, two-sided t-test).
sim_trials <- function(s, reps, seed) {
  set.seed(seed)
  ms <- lapply(1:4, function(k) path_model(k, scen$rule[s], scen$th[s]))
  ch <- lapply(ms, function(m) chol(m$S))
  out <- matrix(NA_real_, reps, 2)
  for (r in seq_len(reps)) {
    D <- B <- numeric(0)
    for (k in 1:4) {
      m <- ms[[k]]
      X <- matrix(rnorm(NS[k] * length(m$mu)), NS[k]) %*% ch[[k]]
      X <- sweep(X, 2, m$mu, '+')
      n <- m$n
      br <- X[, m$idx, drop = FALSE]
      if (m$rule == 'mean') {
        z <- (X[, 1] - BM_M) / SD_BM
        br[, m$pi$on] <- br[, m$pi$on] + STR * SD_BR * z
      }
      Y <- X[, 2] - (X[, 3:(2 + n)] + X[, (3 + n):(2 + 2 * n)] + br)
      a <- ifelse(m$pi$on, 1 / sum(m$pi$on), -1 / sum(!m$pi$on))
      D <- c(D, drop(Y %*% a)); B <- c(B, X[, 1])
    }
    sxx <- sum((B - mean(B))^2)
    bh <- sum((B - mean(B)) * (D - mean(D))) / sxx
    rss <- sum((D - mean(D) - bh * (B - mean(B)))^2)
    tval <- bh / sqrt(rss / (N - 2) / sxx)
    out[r, ] <- c(bh, 2 * pt(-abs(tval), N - 2))
  }
  data.table(scenario = scen$sid[s], slope = out[, 1], p = out[, 2])
}
mc_file <- file.path(out_dir, 'mc-e9-reps5000.csv')
if (file.exists(mc_file)) {
  mc <- fread(mc_file)
} else {
  mc <- rbindlist(lapply(seq_len(nrow(scen)), function(s) sim_trials(s, 5000, 7600 + s)))
  fwrite(mc, mc_file)
}
mc_sum <- mc[, .(mc_mean = mean(slope), mc_sd = sd(slope), mc_power = mean(p < 0.05),
                 reps = .N), by = scenario]
mc_sum <- merge(pool4[, .(scenario, beta, sd, power, power_normal)], mc_sum, by = 'scenario')
mc_sum[, `:=`(mcse_mean = mc_sd / sqrt(reps),
              mcse_power = sqrt(mc_power * (1 - mc_power) / reps))]
mc_sum[, `:=`(z_mean = (mc_mean - beta) / mcse_mean, z_power = (mc_power - power) / mcse_power)]
cat('\nClosed form against Monte Carlo\n')
print(mc_sum, digits = 4)
stopifnot(all(abs(mc_sum$z_mean) < 3), all(abs(mc_sum$z_power) < 3))
fwrite(mc_sum, file.path(out_dir, 'closed-form-vs-mc.csv'))

## ---------------- per-visit bands ----------------
band <- rbindlist(lapply(seq_len(nrow(scen)), function(s) {
  rbindlist(lapply(1:4, function(k) {
    pi <- PI[[k]]
    b <- dgp$build(pi, 'A', 0, scen$th[s])
    n <- b$n
    g <- profile(pi, scen$rule[s], scen$th[s])
    data.table(scenario = scen$sid[s], lab = scen$lab[s], path = k, week = pi$w,
               on = pi$on, tsd = pi$tsd, g = g,
               TV = b$mu[3:(2 + n)], PB = b$mu[(3 + n):(2 + 2 * n)],
               BR = b$mu[(3 + 2 * n):(2 + 3 * n)])
  }))
}))
band[, `:=`(mid = TV + PB + BR, half = STR * SD_BR * g)]
band[, `:=`(hi = mid + half, lo = mid - half, gap = 2 * half)]
band[, lab := factor(lab, levels = scen$lab)]
fwrite(band, file.path(out_dir, 'band-by-visit.csv'))
cat('\nPath 1, band width by visit\n')
print(dcast(band[path == 1], week + on + tsd ~ scenario, value.var = 'gap'), digits = 3)

## ---------------- Figure 1: one visit, a bivariate normal ----------------
set.seed(761)
m1 <- band[scenario == 'decay1' & path == 1]
pan <- data.table(week = c(8, 12),
                  lab = c('A. Week 8, on drug: Cor(B, BR) = 0.25',
                          'B. Week 12, two weeks after stopping: Cor(B, BR) = 0.0625'))
i1 <- match(pan$week, m1$week)
pan[, `:=`(mu = m1$BR[i1], r = STR * m1$g[i1])]
pts <- rbindlist(lapply(seq_len(nrow(pan)), function(i) {
  X <- mvrnorm(400, c(0, pan$mu[i]), matrix(c(1, pan$r[i] * SD_BR, pan$r[i] * SD_BR,
                                               SD_BR^2), 2))
  data.table(lab = pan$lab[i], z = X[, 1], br = X[, 2])
}))
lines1 <- pan[, .(lab, mu, slope = r * SD_BR)]
marks1 <- pan[, .(lab, z = c(-1, 1), y = mu + c(-1, 1) * r * SD_BR), by = seq_len(nrow(pan))]
txt1 <- pan[, .(lab, txt = sprintf('E(BR | b) = %.2f + %.2f b\nband width = 2 x %.2f = %.2f',
                                   mu, r * SD_BR, r * SD_BR, 2 * r * SD_BR))]
g1 <- ggplot(pts, aes(z, br)) +
  geom_point(colour = '#9ec5f0', size = 1, alpha = 0.7) +
  geom_abline(data = lines1, aes(intercept = mu, slope = slope), colour = on_dark,
              linewidth = 1) +
  geom_vline(xintercept = c(-1, 1), colour = ink2, linetype = 'dashed', linewidth = 0.4) +
  geom_point(data = marks1, aes(z, y), colour = ink, size = 2.6) +
  geom_text(data = txt1, aes(x = -3.2, y = 27, label = txt), hjust = 0, vjust = 1,
            size = 3.3, colour = ink) +
  facet_wrap(~lab) +
  coord_cartesian(xlim = c(-3.2, 3.2), ylim = c(-20, 28)) +
  labs(x = 'Standardized biomarker b', y = 'Drug response BR (points)',
       title = 'The interaction at one visit is a regression slope',
       subtitle = paste('400 simulated participants, path 1, decaying rule, half-life 1 week.',
                        'The line is the conditional mean E(BR | b);\nthe dots mark it at',
                        'b = -1 and b = +1. Their vertical distance is the band width at that visit.')) +
  theme_fig()
save_fig(g1, '76-fig1-one-visit-regression.png', 10.5, 4.8)

## ---------------- Figure 2: the band on the trajectory, path 1 ----------------
b2 <- band[path == 1]
b2 <- rbind(b2[, .(lab, week, on, mid, hi, lo)],
            unique(b2[, .(lab, week = 0, on = FALSE, mid = 0, hi = 0, lo = 0)]))
bars2 <- band[path == 1 & on == TRUE, .(lab, week)]
bars2[, start := c(0, 4, 8, 9, 10, 11, 12, 16)[match(week, c(4, 8, 9, 10, 11, 12, 16, 20))]]
lg2 <- melt(b2, id.vars = c('lab', 'week'), measure.vars = c('hi', 'mid', 'lo'),
            variable.name = 'grp')
lg2[, grp := factor(grp, levels = c('hi', 'mid', 'lo'),
                    labels = c('Biomarker 1 SD above the mean (b = +1)',
                               'Biomarker at the mean (b = 0)',
                               'Biomarker 1 SD below the mean (b = -1)'))]
g2 <- ggplot() +
  geom_rect(data = bars2, aes(xmin = start, xmax = week, ymin = -2.2, ymax = -1.2),
            fill = ink2) +
  geom_ribbon(data = b2, aes(week, ymin = lo, ymax = hi), fill = on_col, alpha = 0.18) +
  geom_line(data = lg2, aes(week, value, colour = grp, linetype = grp), linewidth = 0.8) +
  geom_point(data = lg2, aes(week, value, colour = grp), size = 1.3) +
  scale_colour_manual(values = c(on_dark, ink2, '#7fb2ea')) +
  scale_linetype_manual(values = c('solid', 'dotted', 'solid')) +
  facet_wrap(~lab, nrow = 2) +
  scale_x_continuous(breaks = c(0, 4, 8, 9, 10, 11, 12, 16, 20)) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  labs(x = 'Week (dark bars: on drug)', y = 'Expected improvement TV + PB + BR (points)',
       title = 'The biomarker band on the trajectory of path 1',
       subtitle = paste('c_bm = 0.25, published parameters. The shaded band is the difference',
                        'the biomarker makes at each visit.\nThe interaction asks how much',
                        'wider the band is at on-drug visits than at off-drug visits.')) +
  theme_fig() + theme(legend.text = element_text(size = 8.5))
save_fig(g2, '76-fig2-band-on-trajectory.png', 11, 7.4)

## ---------------- Figure 4: Figure 5 of docs/40 for every path ----------------
f3 <- copy(band)
f3[, state := factor(ifelse(on, 'On drug', 'Off drug'), levels = c('On drug', 'Off drug'))]
f3[, path_f := factor(paste('Path', path))]
seg3 <- f3[, .(m = mean(gap)), by = .(lab, path_f, state)]
ann3 <- merge(tab9, scen[, .(scenario = sid, lab)], by = 'scenario')
ann3[, `:=`(lab = factor(lab, levels = scen$lab), path_f = factor(paste('Path', path)),
            txt = sprintf('G = %.3f\nbeta = %.3f', G, ifelse(abs(beta) < 5e-5, 0, beta)))]
g3 <- ggplot(f3, aes(week, gap)) +
  geom_col(aes(fill = state), width = 0.75) +
  geom_point(aes(colour = state, shape = state), size = 1.6, show.legend = FALSE) +
  geom_hline(data = seg3, aes(yintercept = m, colour = state, linetype = state),
             linewidth = 0.6) +
  geom_text(data = ann3, aes(x = 0.2, y = 5.6, label = txt), hjust = 0, vjust = 1,
            size = 2.9, colour = ink, lineheight = 0.9) +
  facet_grid(lab ~ path_f, labeller = label_wrap_gen(width = 22)) +
  scale_fill_manual(values = c('On drug' = on_col, 'Off drug' = off_col)) +
  scale_colour_manual(values = c('On drug' = on_dark, 'Off drug' = ink2)) +
  scale_shape_manual(values = c('On drug' = 16, 'Off drug' = 21)) +
  scale_linetype_manual(values = c('On drug' = 'solid', 'Off drug' = 'dashed')) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20)) +
  coord_cartesian(ylim = c(0, 5.8)) +
  labs(x = 'Week', y = 'Band width at the visit (points)',
       title = 'Band width by visit, path and correlation rule',
       subtitle = paste('Bars: expected improvement at b = +1 minus at b = -1. Lines: on-drug',
                        '(solid) and off-drug (dashed) means.\nG = on mean minus off mean,',
                        'divided by the full width 2 c_bm sd_BR = 4; beta = -c_bm (sd_BR / sd_bm) G.')) +
  theme_fig() + theme(strip.text.y = element_text(angle = 0, size = 8.5))
save_fig(g3, '76-fig4-band-width-all-paths.png', 12, 9.6)

## ---------------- Figure 3: the decay of the correlation profile ----------------
tsd_grid <- seq(0, 8, by = 0.02)
dec4 <- rbindlist(lapply(c(0.1, 0.2, 0.5, 1, 2), function(th)
  data.table(t_half = sprintf('%s week', th), tsd = tsd_grid, g = 0.5^(tsd_grid / th))))
dec4[, t_half := factor(t_half, levels = sprintf('%s week', c(0.1, 0.2, 0.5, 1, 2)))]
vis4 <- unique(rbindlist(lapply(1:4, function(k)
  data.table(tsd = PI[[k]]$tsd[!PI[[k]]$on & PI[[k]]$tsd > 0]))))
ramp5 <- c('#c6dbf3', '#9ec5f0', '#5c9be0', '#2a78d6', '#1a4f9c')
g4 <- ggplot(dec4, aes(tsd, g, colour = t_half)) +
  geom_vline(data = vis4, aes(xintercept = tsd), colour = grid_col, linewidth = 2.2) +
  geom_line(linewidth = 0.9) +
  annotate('segment', x = 0.02, xend = 8, y = 1, yend = 1,
           colour = rule_pal[['Constant rule']], linetype = 'dashed', linewidth = 0.7) +
  annotate('text', x = 8, y = 1.05, hjust = 1, size = 3.2, colour = rule_pal[['Constant rule']],
           label = 'constant rule: g = 1 at every off-drug visit, for any half-life above 0') +
  scale_colour_manual(values = ramp5, name = 'Half-life') +
  scale_x_continuous(breaks = 0:8) +
  coord_cartesian(ylim = c(0, 1.08)) +
  labs(x = 'Weeks since the drug was stopped (t_sd)',
       y = 'Correlation profile g at an off-drug visit',
       title = 'How the decaying rule shrinks the correlation after the drug stops',
       subtitle = paste('g = 2^(-t_sd / t_half). Gray bands mark the values of t_sd at which',
                        'the Hybrid paths have off-drug visits (1, 2, 3, 4, 6 and 7 weeks).')) +
  theme_fig() + theme(legend.title = element_text(colour = ink2))
save_fig(g4, '76-fig3-decay-of-profile.png', 10, 5)

## ---------------- Figure 5: the attenuation factor G by half-life ----------------
pts5 <- dec_tab[, .(t_half, path = ifelse(path == 'Pooled', 'Pooled', paste('Path', path)), G)]
g5 <- ggplot(curve_G, aes(t_half, G, colour = path, linetype = path == 'Pooled')) +
  geom_hline(yintercept = 1, colour = rule_pal[['Mean moderation']], linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = rule_pal[['Constant rule']], linewidth = 0.6) +
  annotate('text', x = 0.01, y = 1.04, hjust = 0, size = 3.1,
           colour = rule_pal[['Mean moderation']], label = 'mean moderation: G = 1 at every half-life') +
  annotate('text', x = 0.01, y = 0.04, hjust = 0, size = 3.1,
           colour = rule_pal[['Constant rule']],
           label = 'constant rule: G = 0 at every half-life above 0 (G = 1 at exactly 0)') +
  geom_line(linewidth = 0.8) +
  geom_point(data = pts5, aes(t_half, G, colour = path), size = 1.8, inherit.aes = FALSE) +
  scale_colour_manual(values = path_pal) +
  scale_linetype_manual(values = c(`FALSE` = 'solid', `TRUE` = 'dashed'), guide = 'none') +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.2, 0.5, 1, 2, 5, 10),
                labels = c('0.01', '0.03', '0.1', '0.2', '0.5', '1', '2', '5', '10')) +
  coord_cartesian(ylim = c(0, 1.07)) +
  labs(x = 'Carryover half-life (weeks, log scale)',
       y = 'On/off correlation difference G',
       title = 'The decaying rule attenuates the interaction smoothly as the half-life grows',
       subtitle = paste('G = mean g on drug - mean g off drug, decaying rule, by path, and the',
                        'allocation-weighted pooled value (dashed).\nPoints: the tabulated',
                        'half-lives. The slope is beta = -0.130 x G,',
                        'and the sample size needed for a given power grows roughly as 1 / G^2.')) +
  theme_fig()
save_fig(g5, '76-fig5-attenuation-by-half-life.png', 10, 5.4)

## ---------------- Figure 6: what the band looks like in data ----------------
## 20,000 participants fix the lines close to their expectations; 1,000
## of them are drawn as points.
set.seed(766)
N6 <- 20000L
sim6 <- rbindlist(lapply(seq_len(nrow(scen)), function(s) {
  m <- path_model(1, scen$rule[s], scen$th[s])
  X <- mvrnorm(N6, m$mu, m$S)
  n <- m$n
  br <- X[, m$idx]
  if (m$rule == 'mean') br[, m$pi$on] <- br[, m$pi$on] + STR * SD_BR * (X[, 1] - BM_M) / SD_BM
  imp <- X[, 3:(2 + n)] + X[, (3 + n):(2 + 2 * n)] + br
  rbind(data.table(lab = scen$lab[s], B = X[, 1], state = 'On-drug mean',
                   y = rowMeans(imp[, m$pi$on])),
        data.table(lab = scen$lab[s], B = X[, 1], state = 'Off-drug mean',
                   y = rowMeans(imp[, !m$pi$on])))
}))
sim6[, `:=`(lab = factor(lab, levels = scen$lab),
            state = factor(state, levels = c('On-drug mean', 'Off-drug mean')))]
fit6 <- sim6[, as.list(coef(lm(y ~ B))), by = .(lab, state)]
setnames(fit6, c('lab', 'state', 'a', 'b'))
dif6 <- dcast(fit6, lab ~ state, value.var = 'b')
k6 <- STR * SD_BR / SD_BM
exp6 <- tab9[path == 1, .(lab = factor(scen$lab[match(scenario, scen$sid)], levels = scen$lab),
                          e_on = k6 * g_on, e_off = k6 * g_off)]
dif6 <- merge(dif6, exp6, by = 'lab')
dif6[, txt := sprintf(paste0('slope on = %.3f (expected %.3f)\nslope off = %.3f (expected %.3f)\n',
                             'difference = %.3f (expected %.3f)'),
                      `On-drug mean`, e_on, `Off-drug mean`, e_off,
                      `On-drug mean` - `Off-drug mean`, e_on - e_off)]
fwrite(dif6, file.path(out_dir, 'path1-on-off-slopes-vs-expected.csv'))
g6 <- ggplot(sim6[, .SD[seq_len(1000)], by = .(lab, state)], aes(B, y, colour = state)) +
  geom_point(size = 0.5, alpha = 0.25) +
  geom_abline(data = fit6, aes(intercept = a, slope = b, colour = state), linewidth = 1) +
  geom_text(data = dif6, aes(x = 75, y = 44, label = txt), hjust = 0, vjust = 1, size = 3.1,
            colour = ink, inherit.aes = FALSE) +
  facet_wrap(~lab, nrow = 2) +
  scale_colour_manual(values = c('On-drug mean' = on_dark, 'Off-drug mean' = '#8a8984')) +
  coord_cartesian(xlim = c(75, 175), ylim = c(-30, 45)) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = 'Biomarker B (standing systolic blood pressure, mmHg)',
       y = 'Mean improvement over the visits (points)',
       title = 'The band in data: two regression lines per participant group',
       subtitle = paste('20,000 simulated participants on path 1 (1,000 drawn). Each has an',
                        'on-drug mean and an off-drug mean improvement.\nThe interaction is the',
                        'difference between the two slopes (in improvement units; the E9 slope is its negative).')) +
  theme_fig()
save_fig(g6, '76-fig6-band-in-data.png', 11, 7.6)
fwrite(fit6, file.path(out_dir, 'path1-on-off-slopes.csv'))

## ---------------- Figure 7: the sampling distribution of the E9 slope ----------------
mc[, lab := factor(scen$lab[match(scenario, scen$sid)], levels = scen$lab)]
dens7 <- rbindlist(lapply(seq_len(nrow(scen)), function(s) {
  r <- pool4[scenario == scen$sid[s]]
  x <- seq(-0.35, 0.15, length.out = 300)
  data.table(lab = scen$lab[s], x = x, d = dnorm(x, r$beta, r$sd))
}))
dens7[, lab := factor(lab, levels = scen$lab)]
ann7 <- mc_sum[, .(lab = factor(scen$lab[match(scenario, scen$sid)], levels = scen$lab),
                   txt = sprintf('expected slope %.3f\nSD %.3f\npower: closed form %.3f, simulated %.3f',
                                 beta, sd, power, mc_power))]
g7 <- ggplot(mc, aes(slope)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.01, fill = '#9ec5f0',
                 colour = surface, boundary = 0) +
  geom_line(data = dens7, aes(x, d), colour = on_dark, linewidth = 0.8) +
  geom_vline(xintercept = 0, colour = ink2, linetype = 'dashed') +
  geom_text(data = ann7, aes(x = -0.35, y = 11.4, label = txt), hjust = 0, vjust = 1,
            size = 3, colour = ink) +
  facet_wrap(~lab, nrow = 2) +
  coord_cartesian(xlim = c(-0.35, 0.15), ylim = c(0, 11.5)) +
  labs(x = 'Estimated E9 slope (points per mmHg)', y = 'Density',
       title = 'From the expected slope to power: 5,000 simulated trials of N = 70',
       subtitle = paste('Histogram: simulated slopes; curve: the normal distribution with the',
                        'closed-form mean and SD. Power is the share of trials whose\nslope is',
                        'far enough from zero (two-sided t-test, alpha = 0.05).')) +
  theme_fig()
save_fig(g7, '76-fig7-sampling-distribution.png', 10.5, 7)

## ---------------- Figure 8: power by half-life ----------------
rl <- c(constant = 'Constant rule', decaying = 'Decaying rule', mean = 'Mean moderation')
curve_pow[, rule_f := factor(rl[rule], levels = rl)]
pts8 <- merge(mc_sum, scen[, .(scenario = sid, rule, th)], by = 'scenario')[th > 0]
pts8[, rule_f := factor(rl[rule], levels = rl)]
g8 <- ggplot(curve_pow, aes(t_half, power, colour = rule_f)) +
  geom_hline(yintercept = c(0.05, p0), colour = ink2, linetype = c('dotted', 'dashed'),
             linewidth = 0.4) +
  annotate('text', x = 10, y = p0 + 0.025, hjust = 1, size = 3.1, colour = ink2,
           label = sprintf('no carryover (covariance rules): %.3f', p0)) +
  annotate('text', x = 0.011, y = 0.085, hjust = 0, size = 3.1, colour = ink2,
           label = 'test size 0.05') +
  geom_line(linewidth = 0.9) +
  geom_pointrange(data = pts8, aes(th, mc_power, ymin = mc_power - 2 * mcse_power,
                                   ymax = mc_power + 2 * mcse_power), size = 0.35) +
  scale_colour_manual(values = setNames(rule_pal, rl)) +
  scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.2, 0.5, 1, 2, 5, 10),
                labels = c('0.01', '0.03', '0.1', '0.2', '0.5', '1', '2', '5', '10')) +
  coord_cartesian(ylim = c(0, 0.8)) +
  labs(x = 'Carryover half-life (weeks, log scale)', y = 'Power of the E9 test',
       title = 'Power follows the on/off correlation difference',
       subtitle = paste('Hybrid, N = 70, c_bm = 0.25, compound symmetry; closed form (lines)',
                        'and simulation (points, 2 Monte Carlo SE).')) +
  theme_fig()
save_fig(g8, '76-fig8-power-by-half-life.png', 10, 5)
cat('\nPower at the tabulated half-lives\n')
print(dcast(curve_pow[t_half %in% th_tab], t_half ~ rule, value.var = 'power'), digits = 3)
