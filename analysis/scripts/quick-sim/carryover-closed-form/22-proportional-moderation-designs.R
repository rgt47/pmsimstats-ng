## Closed-form design comparison under proportional moderation
## (docs/38, Sections 4 and 6), against the decayed mean moderation of
## docs/39. No whitepaper yet; outputs in
## analysis/data/quick-sim/carryover-closed-form/proportional-moderation/.
##
## Rules for the moderation profile h_t (the fraction of the full
## biomarker moderation present at visit t):
##   current        1 on drug, 2^(-tsd / t1/2) off drug after exposure,
##                  0 before (docs/39, 'mean_decayed')
##   prop_zero      h_t = f_t, the normalized drug-response curve
##                  (modified Gompertz, d = 5, r = 0.42), restarting from
##                  zero at every re-exposure; off drug f decays from its
##                  last on-drug value with half-life t1/2
##   prop_resid     the same, but a re-exposure restarts the curve from
##                  the residual level at the restart (docs/38, Section 7)
## Patient i's drug response is (omega_bar zeta b_i + v_i) h_t + xi_it,
## with omega_bar zeta = c_bm sd_BR (matched to docs/39 at h = 1) and
## v_i ~ N(0, sigma_v^2) sensitivity variation that the biomarker does
## not explain, applied under every rule. xi is the BR block of the
## response covariance, as in docs/39.
##
## Analysis: E9 with sequence-specific intercepts (docs/39, stratified),
## N = 70, c_bm = 0.25. Response correlation: CS, separable AR(1) with
## 20% or 50% measurement error (docs/39).
##
## Run from anywhere:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/22-proportional-moderation-designs.R \
##     [--mc-reps R] [--cores K]
suppressMessages({library(data.table); library(MASS); library(parallel)})
root <- here::here()
setwd(root)
src <- readLines('analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R')
stop_at <- grep('^## ---------------- per-path table', src)[1] - 1
eval(parse(text = src[1:stop_at]))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, dflt) if (k %in% args) as.integer(args[which(args == k) + 1]) else dflt
mc_reps <- arg('--mc-reps', 5000L)
cores <- arg('--cores', max(1L, parallelly::availableCores() - 1L))
out_dir <- 'analysis/data/quick-sim/carryover-closed-form/proportional-moderation'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

STRUCTS <- c('CS', 'nug0.2', 'nug0.5')
HALVES <- c(0, 0.2, 0.5, 1, 2)
RULES_P <- c('current', 'prop_zero', 'prop_resid')
SIGMA_V <- c(0, 2)
MOD_FULL <- CBM * sd_BR
K3 <- (1 - dgp$C1) * diag(3) + dgp$C1

## Normalized drug-response curve and its inverse.
f_curve <- function(t) dgp$gomp(t, 1, 5, 0.42)
# Residual levels below 1e-10 are treated as zero: f(0) carries
# rounding error of that order, which would leave no sign change.
f_inv <- function(y) {
  if (y <= 1e-10) return(0)
  uniroot(function(t) f_curve(t) - y, c(0, 200), tol = 1e-10)$root
}
stopifnot(abs(f_curve(0)) < 1e-12, abs(f_curve(100) - 1) < 1e-6,
          abs(f_curve(f_inv(0.39)) - 0.39) < 1e-8)

## Moderation profile of one path. path_info() places a start or stop
## at the visit before the first visit of the new state, so the stop
## time is the week of the last on-drug visit and the restart time is
## w_t - tod_t.
profile <- function(pi, rule, th) {
  if (rule == 'current') return(decay_profile(pi, th))
  h <- numeric(n)
  f_stop <- NA_real_; w_stop <- NA_real_; t0_run <- 0
  for (t in seq_len(n)) {
    if (pi$on[t]) {
      # The residual is fixed once per on-drug run, at its first visit;
      # later visits of the run continue the same curve.
      if (t == 1 || !pi$on[t - 1]) {
        w_start <- pi$w[t] - pi$tod[t]
        r0 <- if (rule == 'prop_resid' && !is.na(f_stop) && th > 0)
          f_stop * 0.5^((w_start - w_stop) / th) else 0
        stopifnot(r0 >= 0, r0 < 1)
        t0_run <- f_inv(r0)
      }
      h[t] <- f_curve(t0_run + pi$tod[t])
      f_stop <- h[t]; w_stop <- pi$w[t]
    } else if (!is.na(f_stop) && th > 0) {
      h[t] <- f_stop * 0.5^((pi$w[t] - w_stop) / th)
    }
  }
  h
}

resp_cor <- function(pi, structure) {
  if (structure == 'CS') return(dgp$build(pi, 'A', 0, 0)$R[-(1:2), -(1:2)])
  nu <- as.numeric(sub('nug', '', structure))
  kronecker(K3, (1 - nu) * dgp$RHO^abs(outer(pi$w, pi$w, '-')) + nu * diag(n))
}
resp_cov <- function(pi, structure) {
  sds <- c(rep(10, n), 10 * pi$e, rep(sd_BR, n))
  L <- cbind(diag(n), diag(n), diag(n))
  L %*% (outer(sds, sds) * resp_cor(pi, structure)) %*% t(L)
}

seq_moments <- function(pi, rule, th, structure, sigma_v) {
  a <- weights(pi, 'E9')
  h <- profile(pi, rule, th)
  Fx <- sum(a * h)
  aCa <- drop(t(a) %*% resp_cov(pi, structure) %*% a)
  tau2 <- aCa + sigma_v^2 * Fx^2
  data.table(gamma = -MOD_FULL / sd_B * Fx, var_D = tau2 + (MOD_FULL * Fx)^2,
             tau2 = tau2, F = Fx, aCa = aCa)
}

strat_power <- function(m, Ns) {
  P <- length(Ns); Nt <- sum(Ns)
  w <- (Ns - 1) / sum(Ns - 1)
  gbar <- sum(w * m$gamma)
  resid <- sum(w * (m$tau2 + sd_B^2 * (m$gamma - gbar)^2))
  df <- Nt - P - 1
  q <- qchisq((seq_len(400) - 0.5) / 400, Nt - P)
  ncp <- gbar * sd_B * sqrt(q) / sqrt(resid)
  tc <- qt(1 - ALPHA / 2, df)
  mean(pt(-tc, df, ncp = ncp) + pt(tc, df, ncp = ncp, lower.tail = FALSE))
}

## Designs: the three published designs, and blinded counterbalanced
## ABAB and alternating designs at 2.5-week spacing (20 weeks) and the
## alternating design at 4-week spacing (32 weeks; the robust design of
## docs/39, Section 6).
blinded <- function(pattern, d) {
  list(t_wk = rep(d, n), e = rep(0.5, n), paths = list(pattern, 1 - pattern))
}
DES <- list(
  CO = dgp$designs$CO,
  Hybrid = dgp$designs$Hybrid,
  `OL+BDC` = dgp$designs$`OL+BDC`,
  ABAB = blinded(c(1, 1, 0, 0, 1, 1, 0, 0), 2.5),
  Alternating = blinded(c(1, 0, 1, 0, 1, 0, 1, 0), 2.5),
  `Alternating 4 wk` = blinded(c(1, 0, 1, 0, 1, 0, 1, 0), 4))
alloc <- function(P) N %/% P + c(rep(1L, N %% P), rep(0L, P - N %% P))
paths_of <- function(des) lapply(des$paths, function(od) dgp$path_info(des$t_wk, des$e, od))

eval_cell <- function(dname, rule, th, structure, sigma_v) {
  ps <- paths_of(DES[[dname]])
  m <- rbindlist(lapply(ps, seq_moments, rule = rule, th = th,
                        structure = structure, sigma_v = sigma_v))
  Ns <- alloc(length(ps))
  data.table(design = dname, rule = rule, t_half = th, structure = structure,
             sigma_v = sigma_v, F_mean = sum(Ns * m$F) / N,
             aCa_mean = sum(Ns * m$aCa) / N,
             power = strat_power(m, Ns))
}

grid <- CJ(design = names(DES), rule = RULES_P, t_half = HALVES,
           structure = STRUCTS, sigma_v = SIGMA_V, sorted = FALSE)
## prop_resid differs from prop_zero only when a re-exposure follows a
## decay; at t1/2 = 0 the two are identical by construction.
res <- rbindlist(mclapply(seq_len(nrow(grid)), function(k)
  eval_cell(grid$design[k], grid$rule[k], grid$t_half[k], grid$structure[k],
            grid$sigma_v[k]), mc.cores = cores))
stopifnot(nrow(res) == nrow(grid), all(res$power > 0 & res$power < 1))
fwrite(res, file.path(out_dir, 'closed-form-power.csv'))

## Validation 1: the current rule reproduces docs/39 (CS, t1/2 = 0,
## Hybrid 0.663, CO 0.765, both at the published schedules).
chk <- res[rule == 'current' & t_half == 0 & structure == 'CS' & sigma_v == 0]
cat('\nCurrent rule, CS, no carryover (docs/39: Hybrid 0.663, CO 0.765)\n')
print(chk[, .(design, power = round(power, 3))])

## Profiles for the record.
prof_rows <- rbindlist(lapply(names(DES), function(dn) {
  ps <- paths_of(DES[[dn]])
  rbindlist(lapply(seq_along(ps), function(p) rbindlist(lapply(RULES_P, function(r)
    rbindlist(lapply(c(0.2, 1), function(th)
      data.table(design = dn, path = p, rule = r, t_half = th,
                 week = ps[[p]]$w, on = as.integer(ps[[p]]$on),
                 h = round(profile(ps[[p]], r, th), 4))))))))
}))
fwrite(prof_rows, file.path(out_dir, 'moderation-profiles.csv'))

## Validation 2: visit-level Monte Carlo of the stratified E9 for a
## subset of cells. Each replicate draws the response factors from their
## joint covariance, adds the moderated drug response at every visit,
## forms each patient's on-minus-off contrast from the visit data and
## fits the stratified regression; the closed form is not used.
mc_cell <- function(dname, rule, th, structure, sigma_v, reps, seed) {
  set.seed(seed)
  ps <- paths_of(DES[[dname]]); Ns <- alloc(length(ps))
  pre <- lapply(ps, function(pi) list(
    a = weights(pi, 'E9'), h = profile(pi, rule, th),
    R = chol(resp_cov(pi, structure) + diag(1e-10, n))))
  rej <- 0L
  for (r in seq_len(reps)) {
    D <- B <- S <- NULL
    for (p in seq_along(ps)) {
      np <- Ns[p]; pp <- pre[[p]]
      b <- rnorm(np)
      eps <- matrix(rnorm(np * n), np) %*% pp$R
      br <- outer(MOD_FULL * b + rnorm(np, 0, sigma_v), pp$h)
      Y <- -(eps + br)
      D <- c(D, drop(Y %*% pp$a)); B <- c(B, b * sd_B)
      S <- c(S, rep(p, np))
    }
    fit <- lm(D ~ factor(S) + B)
    rej <- rej + (summary(fit)$coefficients['B', 'Pr(>|t|)'] < ALPHA)
  }
  rej / reps
}
mc_grid <- CJ(design = c('CO', 'Hybrid', 'Alternating'),
              rule = c('current', 'prop_zero', 'prop_resid'),
              t_half = c(0.2, 1), structure = c('CS', 'nug0.2'),
              sigma_v = c(0, 2), sorted = FALSE)
mc_grid[, seed := 20261030L + .I]
mc <- rbindlist(mclapply(seq_len(nrow(mc_grid)), function(k) {
  g <- mc_grid[k]
  data.table(g, mc_power = mc_cell(g$design, g$rule, g$t_half, g$structure,
                                   g$sigma_v, mc_reps, g$seed))
}, mc.cores = cores))
mc <- merge(mc, res, by = c('design', 'rule', 't_half', 'structure', 'sigma_v'))
mc[, mcse := sqrt(mc_power * (1 - mc_power) / mc_reps)]
mc[, z := (mc_power - power) / mcse]
fwrite(mc, file.path(out_dir, 'monte-carlo-check.csv'))
cat(sprintf('\nMonte Carlo check: %d cells, %d reps each; |z| max %.2f, mean z %+.2f, %d of %d beyond 2\n',
            nrow(mc), mc_reps, max(abs(mc$z)), mean(mc$z), sum(abs(mc$z) > 2), nrow(mc)))

## Summaries.
show <- function(st, sv) {
  x <- res[structure == st & sigma_v == sv]
  dcast(x[, .(design, rule, t_half, p = round(power, 3))],
        design + rule ~ t_half, value.var = 'p')
}
for (st in STRUCTS) for (sv in SIGMA_V) {
  cat(sprintf('\nPower, structure %s, sigma_v %g (columns: t1/2 in weeks)\n', st, sv))
  print(show(st, sv))
}
rank_tab <- res[, .(best = design[which.max(power)],
                    order = paste(design[order(-power)], collapse = ' > ')),
                by = .(rule, structure, sigma_v, t_half)]
fwrite(rank_tab, file.path(out_dir, 'design-ranking.csv'))
cat('\nBest design by condition\n')
print(dcast(rank_tab, structure + sigma_v + t_half ~ rule, value.var = 'best'))
