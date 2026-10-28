## 17-rho-stage.R
##
## Closed-form inputs for paper 01's staged Results (revision of
## 2026-10-27): the serial-correlation stage, the planning section and
## the re-exposure trajectory figure. Everything is computed from the
## design definitions and published parameters of
## 04-power-simulation.R, as 01-e9-closed-form.R does, with rho varied.
##
##   1. Ceilings on c_bm (graded coupling, t1/2 = 0 and 0.2 week) for
##      rho in {0.5, 0.6, 0.7, 0.8} under compound symmetry, the
##      two-constant AR(1) form (c_1 at lag 0, c_x rho^|lag| otherwise)
##      and the separable form (c_1 rho^|lag|), by design.
##   2. Closed-form E9 and E9b power at c_bm = 0.25, N = 70, Hybrid, for
##      the same grid.
##   3. N for 80% power of E9 under the refined process (separable,
##      rho = 0.7) at t1/2 = 0.2 and 1 week, by design and c_bm.
##   4. Correlation matrices failing corpcor::is.positive.definite() on
##      the published and extended grids at rho = 0.7 (Table nonpd at
##      the recommended rho).
##   5. The non-separable form's validity condition on the paper's
##      schedules and on denser ones.
##   6. The drug-response mean by visit in the Hybrid design under
##      anchored and compounding carryover (the re-exposure step).
##   7. The stress test's option cells (measurement error, dropout,
##      Weibull decay, N, half-life) from the stored first-run summary.
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/17-rho-stage.R
suppressMessages({library(data.table); library(ggplot2); library(corpcor)})
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/rho-stage'
fig_dir <- 'docs/figures'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep <- c('RHO', 'C1', 'CX', 'N_TOTAL', 'bm_m', 'bm_sd', 'bl_m', 'bl_sd',
          'gomp', 'designs', 'path_info', 'build')
dgp <- new.env()
for (e in parse(here04)) {
  if (is.call(e) && identical(e[[1]], as.name('<-')) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% keep) eval(e, dgp)
}
stopifnot(all(keep %in% ls(dgp)))
dgp$designs$OL <- list(t_wk = rep(2.5, 8), e = rep(1, 8), paths = list(rep(1, 8)))

ALPHA <- 0.05
sd_B <- dgp$bm_sd
sd_BR <- 8
ARM <- c(CS = 'A', two = 'F', sep = 'C')
CON_LAB <- c(CS = 'Compound symmetry', two = 'AR(1), two-constant cross-component',
             sep = 'Separable AR(1)')
DESIGNS <- c('Hybrid', 'OL+BDC', 'CO')
RHOS <- c(0.5, 0.6, 0.7, 0.8)

design_paths <- function(name) {
  d <- dgp$designs[[name]]
  lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
}
alloc <- function(N, nP) N %/% nP + c(rep(1L, N %% nP), rep(0L, nP - N %% nP))

## Graded coupling: 1 on drug, 2^(-tsd / t1/2) off drug after exposure,
## 0 before exposure or without carryover.
graded <- function(pi, th) {
  n <- length(pi$w)
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  ifelse(pi$on, 1, phi)
}

## Response block and means at a given rho; c_bm enters the biomarker
## row only, so the block is built at c_bm = 0.
response <- function(pi, construct, rho, th) {
  dgp$RHO <- rho
  b <- dgp$build(pi, ARM[[construct]], 0, th)
  n <- length(pi$w)
  idx <- 3:(2 + 3 * n)
  list(M = b$R[idx, idx], Sigma = b$Sigma[idx, idx], mu = b$mu[idx], n = n)
}

ceiling_of <- function(design, construct, rho, th = 0) {
  min(vapply(design_paths(design), function(pi) {
    r <- response(pi, construct, rho, th)
    u <- c(rep(0, 2 * r$n), graded(pi, th))
    if (!any(u != 0)) return(Inf)
    1 / sqrt(drop(crossprod(u, solve(r$M, u))))
  }, numeric(1)))
}

path_moments <- function(pi, construct, rho, th, variant, cbm) {
  r <- response(pi, construct, rho, th)
  n <- r$n
  L <- cbind(diag(n), diag(n), diag(n))
  C_S <- L %*% r$Sigma %*% t(L)
  mu_S <- drop(L %*% r$mu)
  n_on <- sum(pi$on); n_off <- sum(!pi$on) + (variant == 'E9b')
  a <- ifelse(pi$on, 1 / n_on, -1 / n_off)
  g <- graded(pi, th)
  gamma <- -cbm * sd_BR / sd_B * sum(a * g)
  var_D <- drop(t(a) %*% C_S %*% a)
  data.table(gamma = gamma, tau2 = var_D - gamma^2 * sd_B^2, mu = -sum(a * mu_S))
}

power_of <- function(design, construct, rho, cbm, N = 70L, th = 0, variant = 'E9') {
  paths <- design_paths(design)
  nP <- length(paths)
  Ns <- alloc(N, nP)
  m <- rbindlist(lapply(paths, path_moments, construct, rho, th, variant, cbm))
  w <- Ns / N
  gbar <- sum(w * m$gamma); mubar <- sum(w * m$mu)
  resid <- sum(w * (m$tau2 + (m$mu - mubar)^2 + sd_B^2 * (m$gamma - gbar)^2))
  chi_q <- qchisq((seq_len(400) - 0.5) / 400, N - 1)
  tc <- qt(1 - ALPHA / 2, N - 2)
  ncp <- gbar * sd_B * sqrt(chi_q) / sqrt(resid)
  mean(pt(-tc, N - 2, ncp = ncp) + pt(tc, N - 2, ncp = ncp, lower.tail = FALSE))
}

## Checks against the values already in the paper.
chk <- c(cs08 = power_of('Hybrid', 'CS', 0.8, 0.25),
         cs08b = power_of('Hybrid', 'CS', 0.8, 0.25, variant = 'E9b'),
         sep08 = power_of('Hybrid', 'sep', 0.8, 0.25),
         cs08_t1 = power_of('Hybrid', 'CS', 0.8, 0.25, th = 1))
cat('Checks (paper: 0.721, 0.654, 0.335, 0.473):\n'); print(round(chk, 3))
stopifnot(abs(chk - c(0.721, 0.654, 0.335, 0.473)) < 0.002)
stopifnot(abs(ceiling_of('Hybrid', 'CS', 0.8) - 0.2555) < 1e-3,
          abs(ceiling_of('Hybrid', 'sep', 0.8) - 0.440) < 1e-3)

## ---------------- 1. ceilings by rho ----------------
ceil <- CJ(rho = RHOS, construct = names(ARM), design = DESIGNS, t_half = c(0, 0.2))
ceil[, ceiling := mapply(ceiling_of, design, construct, rho, t_half)]
ceil_min <- ceil[, .(min_ceiling = min(ceiling), binding = design[which.min(ceiling)]),
                 by = .(rho, construct, t_half)]
fwrite(ceil, file.path(out_dir, 'ceilings-by-rho.csv'))
fwrite(ceil_min, file.path(out_dir, 'ceilings-by-rho-min.csv'))
cat('\nCeiling on c_bm, graded coupling, t1/2 = 0 (min over designs; binding design)\n')
print(dcast(ceil_min[t_half == 0], rho ~ construct, value.var = 'min_ceiling'), digits = 3)
print(dcast(ceil_min[t_half == 0], rho ~ construct, value.var = 'binding'))
cat('\nCeiling by design, t1/2 = 0\n')
print(dcast(ceil[t_half == 0], rho + construct ~ design, value.var = 'ceiling'), digits = 3)

## ---------------- 2. power by rho ----------------
pw <- CJ(rho = RHOS, construct = names(ARM), t_half = c(0, 0.2), variant = c('E9', 'E9b'))
pw[, power := mapply(function(r, cn, th, v) power_of('Hybrid', cn, r, 0.25, th = th, variant = v),
                     rho, construct, t_half, variant)]
fwrite(pw, file.path(out_dir, 'power-by-rho.csv'))
cat('\nE9 power, Hybrid, N = 70, c_bm = 0.25, t1/2 = 0\n')
print(dcast(pw[variant == 'E9' & t_half == 0], rho ~ construct, value.var = 'power'), digits = 3)
cat('\nE9 power, t1/2 = 0.2\n')
print(dcast(pw[variant == 'E9' & t_half == 0.2], rho ~ construct, value.var = 'power'), digits = 3)

## ---------------- 3. N for 80% power, refined process ----------------
n80 <- function(design, cbm, th, construct = 'sep', rho = 0.7) {
  Ns <- seq(20L, 800L, by = 2L)
  for (N in Ns) if (power_of(design, construct, rho, cbm, N = N, th = th) >= 0.8) return(N)
  NA_integer_
}
plan <- CJ(design = DESIGNS, c_bm = c(0.25, 0.35, 0.45), t_half = c(0, 0.2, 1))
plan[, N80 := mapply(n80, design, c_bm, t_half)]
plan[, power_N70 := mapply(function(d, cb, th) power_of(d, 'sep', 0.7, cb, th = th),
                           design, c_bm, t_half)]
plan[, slope_per_sd := -c_bm * sd_BR]
fwrite(plan, file.path(out_dir, 'planning-n80.csv'))
cat('\nN for 80% power, E9, refined process (separable, rho = 0.7)\n')
print(dcast(plan, design + c_bm ~ t_half, value.var = 'N80'))
cat('\nPower at N = 70\n')
print(dcast(plan, design + c_bm ~ t_half, value.var = 'power_N70'), digits = 3)
cat('\nImplied slope per SD of biomarker (points):', paste(unique(plan$slope_per_sd)), '\n')
cat('Ceiling 0.474 implies', -0.474 * sd_BR, 'points per SD; biomarker SD', sd_B, 'mm Hg\n')

## ---------------- 4. non-PD counts at rho = 0.7 ----------------
count_arm <- function(construct, rho, dz, cbms, halves) {
  dgp$RHO <- rho
  rows <- list()
  for (dn in dz) for (cbm in cbms) for (th in halves) {
    d <- dgp$designs[[dn]]
    for (g in seq_along(d$paths)) {
      pi <- dgp$path_info(d$t_wk, d$e, d$paths[[g]])
      b <- dgp$build(pi, ARM[[construct]], cbm, th)
      rows[[length(rows) + 1]] <- data.table(construct = construct, rho = rho, design = dn,
                                             c.bm = cbm, t_half = th,
                                             repaired = !is.positive.definite(b$Sigma))
    }
  }
  rbindlist(rows)
}
pub_d <- c('OL', 'OL+BDC', 'CO', 'Hybrid')
grids <- rbindlist(lapply(c(0.8, 0.7), function(r) rbind(
  cbind(grid = 'published', count_arm('CS', r, pub_d, c(0, .3, .6), c(0, .1, .2))),
  cbind(grid = 'published', count_arm('sep', r, pub_d, c(0, .3, .6), c(0, .1, .2))),
  cbind(grid = 'extended', count_arm('CS', r, DESIGNS, c(0, .25, .3, .6), c(0, .1, .2, .5, 1))),
  cbind(grid = 'extended', count_arm('sep', r, DESIGNS, c(0, .25, .3, .6), c(0, .1, .2, .5, 1))))))
nonpd <- rbind(
  grids[grid == 'published', .(matrices = .N, repaired = sum(repaired)), by = .(grid, construct, rho)],
  grids[grid == 'extended' & c.bm > 0, .(matrices = .N, repaired = sum(repaired)),
        by = .(grid, construct, rho, c.bm)], fill = TRUE)
nonpd[, pct := round(100 * repaired / matrices)]
fwrite(nonpd, file.path(out_dir, 'nonpd-by-rho.csv'))
cat('\nMatrices failing is.positive.definite() (CS: step coupling; sep: graded coupling)\n')
print(nonpd)
## Check against the stored counts at rho = 0.8 (Table nonpd).
stopifnot(nonpd[grid == 'published' & rho == 0.8 & construct == 'CS', repaired] == 20,
          nonpd[grid == 'published' & rho == 0.8 & construct == 'sep', repaired] == 21)

## ---------------- 5. validity of the non-separable form ----------------
thr <- (dgp$C1 - dgp$CX) / (1 - dgp$CX)
eig_min <- function(w, rho) min(eigen(rho^abs(outer(w, w, '-')), TRUE, TRUE)$values)
sched <- list(Hybrid = cumsum(dgp$designs$Hybrid$t_wk), `OL+BDC` = cumsum(dgp$designs$`OL+BDC`$t_wk),
              CO = cumsum(dgp$designs$CO$t_wk), weekly20 = 1:20, daily56 = (1:56) / 7)
val <- CJ(schedule = names(sched), rho = c(0.7, 0.8))
val[, min_eig := mapply(function(s, r) eig_min(sched[[s]], r), schedule, rho)]
val[, valid := min_eig > thr]
fwrite(val, file.path(out_dir, 'nonseparable-validity.csv'))
cat(sprintf('\nNon-separable form: threshold (c_1 - c_x)/(1 - c_x) = %.3f\n', thr))
print(val, digits = 3)

## ---------------- 6. re-exposure trajectory ----------------
traj <- rbindlist(lapply(seq_along(dgp$designs$Hybrid$paths), function(p) {
  pi <- design_paths('Hybrid')[[p]]
  n <- length(pi$w)
  rbindlist(lapply(c(0.2, 1), function(th) rbindlist(lapply(c('anchored', 'recursive'), function(cr) {
    b <- dgp$build(pi, 'C', 0, th, carryover = cr)
    data.table(path = p, t_half = th, carryover = cr, week = pi$w, on = pi$on,
               br_mean = b$mu[(3 + 2 * n):(2 + 3 * n)])
  }))))
}))
fwrite(traj, file.path(out_dir, 'hybrid-br-mean-trajectory.csv'))
traj[, path_lab := sprintf('Path %d (%s)', path,
                            sapply(path, function(p) paste(ifelse(dgp$designs$Hybrid$paths[[p]] == 1, 'D', 'P'), collapse = '')))]
traj[, th_lab := factor(sprintf('t1/2 = %g week', t_half), levels = c('t1/2 = 0.2 week', 't1/2 = 1 week'))]
traj[, cr_lab := factor(ifelse(carryover == 'anchored', 'Anchored (refined)', 'Compounding (published)'),
                        levels = c('Anchored (refined)', 'Compounding (published)'))]
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
g <- ggplot(traj, aes(week, br_mean, colour = cr_lab, linetype = cr_lab)) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = on), size = 2.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = 'On drug', `FALSE` = 'Off drug'), name = NULL) +
  scale_colour_manual(values = c('#1baf7a', '#2a78d6'), name = NULL) +
  scale_linetype_manual(values = c('solid', 'dashed'), name = NULL) +
  scale_x_continuous(breaks = c(4, 8, 12, 16, 20)) +
  facet_grid(th_lab ~ path_lab) +
  labs(x = 'Week', y = 'Mean drug response (points)',
       title = 'Drug-response mean by visit in the Hybrid design',
       subtitle = 'Time on drug restarts at re-exposure (week 16 in paths 1 and 3); off drug the mean decays from its last on-drug value',
       caption = 'D: on drug; P: placebo. Visits at weeks 4, 8, 9, 10, 11, 12, 16 and 20.') +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = surface, colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
        strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
        plot.title = element_text(colour = ink, face = 'bold'),
        plot.subtitle = element_text(colour = ink2, size = 9),
        plot.caption = element_text(colour = ink2, hjust = 0),
        plot.title.position = 'plot', plot.caption.position = 'plot',
        legend.position = 'top', legend.justification = 'left')
out_fig <- file.path(fig_dir, '01-fig-hybrid-br-mean-reexposure.png')
ggsave(out_fig, g, width = 10, height = 5.4, dpi = 200, bg = surface)
cat('\nWrote', out_fig, '\n')
cat('Re-exposure step (anchored, path 1): week 12 and week 16 means\n')
print(traj[path == 1 & carryover == 'anchored' & week %in% c(12, 16), .(t_half, week, br_mean)])

## ---------------- 7. option cells of the stress test ----------------
ss <- fread('analysis/data/quick-sim/carryover-closed-form/strawman-stress/summary-full.csv')
opt_labels <- c('curve-C07', 'meas-error', 'dropout-MCAR', 'dropout-MAR', 'weibull0.5',
                'weibull2', 'thalf0', 'thalf1', 'N35', 'N150')
opt <- ss[label %in% opt_labels & cbm %in% c(0, 0.3) & analysis %in% c('R', 'S', 'F1')]
opt[, label := factor(label, levels = opt_labels)]
opt_w <- dcast(opt, label ~ cbm + analysis, value.var = 'reject')
fwrite(opt_w, file.path(out_dir, 'stress-option-cells.csv'))
cat('\nStress test, first run: null rejection (c_bm = 0) and power (c_bm = 0.30)\n')
print(opt_w, digits = 3)
