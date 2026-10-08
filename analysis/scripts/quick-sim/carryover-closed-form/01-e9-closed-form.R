## Closed-form moments and power of the paired-difference interaction
## statistic (E9, compendium-mathematics Section 'Analysis specifications
## for carryover') for the Hybrid design, N = 70, as a function of
## carryover half-life. Whitepaper: docs/37-carryover-power-closed-form.md.
##
## E9: Delta_i = (on-drug mean) - (off-drug mean) of patient i's outcome;
## OLS of Delta_i on bm_i over all patients; test the slope.
## E9b: the same with the baseline visit counted as an off-drug
## observation, which is how the published lme_analysis() codes Db.
##
## For each path the script computes, exactly,
##   gamma_p = Cov(Delta, B) / Var(B)     (regression slope)
##   tau2_p  = Var(Delta) - gamma_p^2 Var(B)   (residual variance)
##   mu_p    = E(Delta)
## and from these the mean and variance of the pooled slope and its power.
##
## Response structures (from 04-power-simulation.R, which reproduces the
## 58b32a9 matrices): CS (the published construct, configuration A),
## covar (AR(1), decaying cross-factor, configuration F) and sep
## (separable AR(1), configuration C).
## Rules for the biomarker-treatment interaction:
##   step          coupling 1 wherever the mean drug response is nonzero
##                 (published)
##   graded        coupling 1 on drug, 2^(-tsd / t1/2) off drug
##   mean          on-drug BR shifted by beta_bm * sd_BR * b, b the
##                 standardized biomarker (configuration D); beta_bm is
##                 matched to c_bm (both CBM) as in NOTATION.md rule 2
##   mean_decayed  the same shift, scaled by 2^(-tsd / t1/2) off drug
##
## Validation: (1) a direct Monte Carlo of E9 and E9b, cached per cell;
## (2) the lmer results of 04 (configurations A and D, published
## analysis).
##
## Run from the repository root:
##   Rscript analysis/scripts/quick-sim/carryover-closed-form/01-e9-closed-form.R \
##     [--mc-reps R] [--cores K]
suppressMessages({library(data.table); library(MASS); library(ggplot2); library(parallel)})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(k, d) if (k %in% args) as.integer(args[which(args == k) + 1]) else d
mc_reps <- arg('--mc-reps', 5000L)
cores <- arg('--cores', 6L)
here04 <- 'analysis/scripts/quick-sim/hendrickson-problems/04-power-simulation.R'
out_dir <- 'analysis/data/quick-sim/carryover-closed-form'
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

CBM <- 0.25
ALPHA <- 0.05
N <- 70L
DESIGNS <- c('Hybrid', 'OL+BDC', 'CO')
## The design in use: its paths and the allocation of the N patients to
## them (as in the published code). Every design has 8 visits.
use_design <- function(name) {
  d <- dgp$designs[[name]]
  nP <<- length(d$paths)
  Ns <<- N %/% nP + c(rep(1L, N %% nP), rep(0L, nP - N %% nP))
  paths <<- lapply(d$paths, function(od) dgp$path_info(d$t_wk, d$e, od))
  invisible(name)
}
use_design('Hybrid')
n <- length(dgp$designs$Hybrid$t_wk)
stopifnot(all(vapply(DESIGNS, function(x) length(dgp$designs[[x]]$t_wk), 1L) == n))
i_tv <- 3:(2 + n); i_pb <- (3 + n):(2 + 2 * n); i_br <- (3 + 2 * n):(2 + 3 * n)
sd_B <- dgp$bm_sd
sd_BR <- 8
ARM <- c(CS = 'A', covar = 'F', sep = 'C')
RULES <- c('step', 'graded', 'mean', 'mean_decayed')
is_cov <- function(rule) rule %in% c('step', 'graded')

## Anchored decay profile: 1 on drug, 2^(-tsd / t1/2) off drug after an
## exposure, 0 before any exposure or without carryover.
decay_profile <- function(pi, th) {
  phi <- if (th > 0) 0.5^(pi$tsd / th) else rep(0, n)
  phi[pi$on | pi$tsd <= 0] <- 0
  ifelse(pi$on, 1, phi)
}

## Coupling pattern g_t (correlation of B with BR_t is c_bm * g_t).
## The step rule reads the mean drug response, which is the same in
## every response structure.
coupling <- function(pi, rule, th) {
  switch(rule,
    step = as.numeric(dgp$build(pi, 'A', 0, th)$mu[i_br] != 0),
    graded = decay_profile(pi, th),
    rep(0, n))
}

## Moderation profile h_t of the mean rules (shift beta_bm * sd_BR * b * h_t).
mod_profile <- function(pi, rule, th) {
  switch(rule,
    mean = as.numeric(pi$on),
    mean_decayed = decay_profile(pi, th),
    rep(0, n))
}

## Contrast weights on the 8 post-baseline visits. E9b counts the
## baseline (outcome = BL, no response factors) as one more off visit;
## BL cancels in both because the weights over all rows sum to zero.
weights <- function(pi, variant) {
  n_on <- sum(pi$on); n_off <- sum(!pi$on) + (variant == 'E9b')
  ifelse(pi$on, 1 / n_on, -1 / n_off)
}

## Exact per-path moments. Delta = -a' S, S_t = TV_t + PB_t + BR_t.
## g, if given, overrides the rule's coupling pattern (used for the
## power-versus-gap curve).
path_moments <- function(pi, rule, th, variant, cbm = CBM, g = NULL, construct = 'CS') {
  b <- dgp$build(pi, ARM[[construct]], 0, th)
  L <- cbind(diag(n), diag(n), diag(n))
  idx <- c(i_tv, i_pb, i_br)
  C_S <- L %*% b$Sigma[idx, idx] %*% t(L)
  mu_S <- drop(L %*% b$mu[idx])
  a <- weights(pi, variant)
  aCa <- drop(t(a) %*% C_S %*% a)
  if (is_cov(rule) || !is.null(g)) {
    if (is.null(g)) g <- coupling(pi, rule, th)
    gamma <- -cbm * sd_BR / sd_B * sum(a * g)
    var_D <- aCa
    prof <- g
  } else {
    ## The shift -cbm*sd_BR*z*sum(a*h) adds variance rather than
    ## explaining it.
    h <- mod_profile(pi, rule, th)
    gamma <- -cbm * sd_BR / sd_B * sum(a * h)
    var_D <- aCa + (cbm * sd_BR * sum(a * h))^2
    prof <- h
  }
  data.table(gamma = gamma, var_D = var_D, tau2 = var_D - gamma^2 * sd_B^2,
             mu = -sum(a * mu_S), gbar_off = mean(prof[!pi$on]),
             wsum = sum(a * prof))
}

## Pooled E9 slope (common intercept, as specified). Path mean
## differences (mu) and slope differences (gamma) act as extra residual
## variance. Conditional on S = sum (B - Bbar)^2 the t-statistic is
## noncentral t with N - 2 df and noncentrality gbar sqrt(S) / sd_resid;
## S ~ sd_B^2 chi^2_{N-1}, so power is that tail probability averaged
## over the chi-square law (quadrature on its quantiles), and
## Var(slope) uses E(1/S) = 1 / (sd_B^2 (N - 3)).
chi_q <- qchisq((seq_len(400) - 0.5) / 400, N - 1)
pooled <- function(m) {
  w <- Ns / N
  gbar <- sum(w * m$gamma)
  mubar <- sum(w * m$mu)
  resid <- sum(w * (m$tau2 + (m$mu - mubar)^2 + sd_B^2 * (m$gamma - gbar)^2))
  V <- sum(w * (m$tau2 + (m$mu - mubar)^2)) / ((N - 3) * sd_B^2) +
    2 * sum(w * (m$gamma - gbar)^2) / N
  tc <- qt(1 - ALPHA / 2, N - 2)
  ncp_S <- gbar * sd_B * sqrt(chi_q) / sqrt(resid)
  power <- mean(pt(-tc, N - 2, ncp = ncp_S) + pt(tc, N - 2, ncp = ncp_S,
                                                  lower.tail = FALSE))
  data.table(mean = gbar, sd = sqrt(V), resid_var = resid,
             ncp = gbar * sd_B * sqrt(N - 1) / sqrt(resid), power = power)
}

closed <- function(rule, th, variant, cbm = CBM, construct = 'CS', design = 'Hybrid') {
  use_design(design)
  m <- rbindlist(lapply(paths, path_moments, rule = rule, th = th, variant = variant,
                        cbm = cbm, construct = construct))
  cbind(data.table(design = design, construct = construct, rule = rule, t_half = th,
                   variant = variant, c.bm = cbm), pooled(m))
}

## c.bm ceiling for a coupling rule in a response structure: the largest
## c.bm with every path's correlation matrix positive definite,
## min over paths (u' M^-1 u)^(-1/2), u the coupling pattern. The mean
## rules put nothing in the matrix and have no ceiling.
ceiling_of <- function(rule, th, construct = 'CS', design = 'Hybrid') {
  if (!is_cov(rule)) return(Inf)
  use_design(design)
  min(vapply(paths, function(pi) {
    R <- dgp$build(pi, ARM[[construct]], 0, th)$R
    M <- R[-(1:2), -(1:2)]
    u <- c(rep(0, 2 * n), coupling(pi, rule, th))
    if (!any(u != 0)) return(Inf)
    1 / sqrt(drop(crossprod(u, solve(M, u))))
  }, numeric(1)))
}

grid_closed <- function(gr) {
  if (!'design' %in% names(gr)) gr[, design := 'Hybrid']
  out <- rbindlist(lapply(seq_len(nrow(gr)), function(k)
    closed(gr$rule[k], gr$t_half[k], gr$variant[k], gr$c.bm[k], gr$construct[k],
           gr$design[k])))
  use_design('Hybrid')
  out
}

## ---------------- per-path table (the structure of the result) --------
pp <- rbindlist(lapply(DESIGNS, function(dn) {
  use_design(dn)
  pp_grid <- CJ(construct = names(ARM), rule = RULES, t_half = c(0, 0.1, 0.5, 1),
                variant = c('E9', 'E9b'), path = seq_len(nP), sorted = FALSE)
  rbindlist(lapply(seq_len(nrow(pp_grid)), function(k) {
    g <- pp_grid[k]
    cbind(design = dn, g, n_on = sum(paths[[g$path]]$on), n_off = sum(!paths[[g$path]]$on),
          path_moments(paths[[g$path]], g$rule, g$t_half, g$variant,
                       construct = g$construct))
  }))
}))
fwrite(pp, file.path(out_dir, 'e9-per-path.csv'))
for (dn in DESIGNS) {
  use_design(dn)
  cat(sprintf('%s paths (n per path %s): on-drug pattern and time since discontinuation\n',
              dn, paste(Ns, collapse = '/')))
  for (p in seq_len(nP))
    cat(sprintf('  path %d: on %s; tsd %s\n', p, paste(as.integer(paths[[p]]$on), collapse = ''),
                paste(paths[[p]]$tsd, collapse = ' ')))
}
cat('\nPer-path slope gamma by design, CS, E9 (no-carryover value',
    round(-CBM * sd_BR / sd_B, 4), ')\n')
print(dcast(pp[construct == 'CS' & variant == 'E9', .(design, rule, path, t_half,
                                                     g = round(gamma, 4))],
            design + rule + path ~ t_half, value.var = 'g'))
pp <- pp[design == 'Hybrid' & construct == 'CS']
use_design('Hybrid')
cat('\nVar(Delta) by path and t1/2 (covariance rules; must not depend on t1/2)\n')
print(dcast(pp[rule == 'step', .(variant, path, t_half, v = round(var_D, 4))],
            variant + path ~ t_half, value.var = 'v'))
cat('\nPer-path slope gamma (expected at no carryover: -c sd_BR/sd_B =',
    round(-CBM * sd_BR / sd_B, 4), ')\n')
print(dcast(pp[, .(rule, variant, path, t_half, g = round(gamma, 4))],
            rule + variant + path ~ t_half, value.var = 'g'))

## ---------------- continuous curves ----------------
## At very short half-lives the recursive off-drug mean underflows to
## exactly 0 at the latest off-drug visits (Hybrid: 2^-1300 at t1/2 =
## 0.01, below about 0.012 weeks; CO, whose off-drug visits are 2.5-10
## weeks after stopping: 2^(-25 / t1/2), below about 0.023 weeks), so the
## step rule, in this code as in 58b32a9, drops their coupling; the
## curves start at 0.03 to keep that floating-point effect out.
th_grid <- c(0, exp(seq(log(0.03), log(4), length.out = 120)))
curves_all <- grid_closed(CJ(design = DESIGNS, construct = names(ARM), rule = RULES,
                             variant = c('E9', 'E9b'), t_half = th_grid, c.bm = CBM,
                             sorted = FALSE))
fwrite(curves_all, file.path(out_dir, 'e9-closed-form-curves-designs.csv'))
curves <- curves_all[design == 'Hybrid']
fwrite(curves, file.path(out_dir, 'e9-closed-form-curves.csv'))
mono <- curves_all[rule %in% c('graded', 'mean_decayed') & t_half > 0][
  order(design, construct, rule, variant, t_half)][
  , .(power_nonincreasing = all(diff(power) <= 1e-12),
      slope_shrinking = all(diff(abs(mean)) <= 1e-12)),
  by = .(design, construct, rule, variant)]
cat('\nDecaying rules, monotonicity over the t1/2 grid\n'); print(mono)

## ---------------- power against the coupling gap ----------------
## Coupling 1 on drug and 1 - gap at every post-baseline off-drug visit,
## in every path, under the CS construct with no carryover in the means
## (t1/2 = 0; the means enter only the negligible path-mean term). The
## gap is defined on the post-baseline visits; under E9b the baseline
## keeps zero coupling, so its effective gap is larger.
closed_gap <- function(gap, variant) {
  m <- rbindlist(lapply(paths, function(pi) {
    g <- ifelse(pi$on, 1, 1 - gap)
    path_moments(pi, 'step', 0, variant, g = g)
  }))
  cbind(data.table(gap = gap, variant = variant), pooled(m))
}
gap_curve <- rbindlist(lapply(c('E9', 'E9b'), function(v)
  rbindlist(lapply(seq(0, 1, by = 0.01), closed_gap, variant = v))))
## Where the coupling rules sit: the patient-weighted mean of
## 1 - gbar_off over paths (post-baseline visits).
rule_pts <- rbindlist(lapply(c(0, 0.2, 0.5, 1, 2), function(th) {
  gp <- vapply(paths, function(pi) mean(coupling(pi, 'graded', th)[!pi$on]), numeric(1))
  data.table(rule = 'graded', t_half = th, gap = 1 - sum(Ns / N * gp))
}))
rule_pts <- rbind(rule_pts, data.table(rule = 'step', t_half = 0.1, gap = 0))
fwrite(gap_curve, file.path(out_dir, 'e9-power-vs-gap.csv'))
cat('\nPower against the coupling gap (closed form)\n')
print(dcast(gap_curve[gap %in% c(0, 0.25, 0.5, 0.75, 1)],
            gap ~ variant, value.var = 'power')[, lapply(.SD, round, 3)])
cat('\nWhere the rules sit (gap)\n'); print(rule_pts[, gap := round(gap, 3)])

## ---------------- heatmap grid: c.bm x t1/2 per rule (CS) ----------------
## Cells above the CS ceiling have no joint distribution; they are
## reported as infeasible, not evaluated.
hm_th <- c(0, 0.1, 0.2, 0.5, 1, 2)
hm_grid <- CJ(rule = RULES, variant = c('E9b', 'E9'), c.bm = c(0, 0.1, 0.2, 0.3),
              t_half = hm_th, sorted = FALSE)
hm <- rbindlist(lapply(seq_len(nrow(hm_grid)), function(k) {
  g <- hm_grid[k]
  ceil <- ceiling_of(g$rule, g$t_half)
  if (g$c.bm >= ceil) return(cbind(g, ceiling = ceil, feasible = FALSE,
                                    power = NA_real_, mean = NA_real_))
  r <- closed(g$rule, g$t_half, g$variant, cbm = g$c.bm)
  cbind(g, ceiling = ceil, feasible = TRUE, power = r$power, mean = r$mean)
}))
fwrite(hm, file.path(out_dir, 'e9-power-heatmap.csv'))
cat('\nCS ceiling by rule and t1/2\n')
print(dcast(unique(hm[, .(rule, t_half, ceiling = round(ceiling, 3))]),
            rule ~ t_half, value.var = 'ceiling'))
cat('\nPower by c.bm and t1/2 (NA = above the CS ceiling)\n')
print(dcast(hm[, .(rule, variant, c.bm, t_half, p = round(power, 3))],
            rule + variant + c.bm ~ t_half, value.var = 'p'))
ceil_all <- CJ(design = DESIGNS, construct = names(ARM), rule = c('step', 'graded'),
                t_half = hm_th)[
  , .(ceiling = ceiling_of(rule, t_half, construct, design)),
  by = .(design, construct, rule, t_half)]
use_design('Hybrid')
fwrite(ceil_all, file.path(out_dir, 'e9-ceilings-designs.csv'))
cat('\nCeilings by design, structure and covariance rule\n')
print(dcast(ceil_all[, .(design, construct, rule, t_half, c = round(ceiling, 3))],
            design + construct + rule ~ t_half, value.var = 'c'))

## ---------------- factor x time factorization ----------------
## With s_t = (sd_TV, sd_PB,t, sd_BR) the factor SDs at visit t and the
## response correlation written as a sum of Kronecker terms (docs/36,
## Appendix A.7),
##   sep    K (x) A                      K = (1-c1) I + c1 J
##   covar  K_a (x) A + K_b (x) I        K_a = (1-cx) I + cx J,
##                                       K_b = (c1-cx)(J - I)
##   CS     K_o (x) I + K_p (x) J        K_o = (1-rho) I + (c1-cx)(J - I),
##                                       K_p = rho I + cx (J - I)
## the contrast variance is Var(Delta) = sum_{t,s} a_t a_s T_ts (s_t' K s_s)
## summed over the terms, T the time kernel (A, I or J). Checked here
## against the full-matrix a' C_S a, and the separable ceiling
## ([K^-1]_33 u' A^-1 u)^(-1/2) against the direct one.
J3 <- matrix(1, 3, 3); I3 <- diag(3)
K_terms <- function(construct) {
  rho <- dgp$RHO; c1 <- dgp$C1; cx <- dgp$CX
  switch(construct,
    sep = list(list(K = (1 - c1) * I3 + c1 * J3, T = 'A')),
    covar = list(list(K = (1 - cx) * I3 + cx * J3, T = 'A'),
                 list(K = (c1 - cx) * (J3 - I3), T = 'I')),
    CS = list(list(K = (1 - rho) * I3 + (c1 - cx) * (J3 - I3), T = 'I'),
              list(K = rho * I3 + cx * (J3 - I3), T = 'J')))
}
time_kernel <- function(pi, type) {
  switch(type, A = dgp$RHO^abs(outer(pi$w, pi$w, '-')), I = diag(n), J = matrix(1, n, n))
}
factored_var <- function(pi, variant, construct) {
  a <- weights(pi, variant)
  s <- cbind(10, 10 * pi$e, sd_BR)
  sum(vapply(K_terms(construct), function(tm) {
    sum(outer(a, a) * time_kernel(pi, tm$T) * (s %*% tm$K %*% t(s)))
  }, numeric(1)))
}
fz <- rbindlist(lapply(DESIGNS, function(dn) {
  use_design(dn)
  rbindlist(lapply(names(ARM), function(cn) rbindlist(lapply(c('E9', 'E9b'), function(v)
    rbindlist(lapply(seq_len(nP), function(p) {
      pi <- paths[[p]]
      data.table(design = dn, construct = cn, variant = v, path = p,
                 full = path_moments(pi, 'step', 0, v, construct = cn)$var_D,
                 factored = factored_var(pi, v, cn),
                 time_part = drop(t(weights(pi, v)) %*% time_kernel(pi, 'A') %*%
                                  weights(pi, v)))
    }))))))
}))
use_design('Hybrid')
fwrite(fz, file.path(out_dir, 'e9-factorization.csv'))
cat(sprintf('\nFactorization: Var(Delta) factored against full matrix, %d path cases: max |diff| %.1e\n',
            nrow(fz), max(abs(fz$full - fz$factored))))
cat('Var(Delta), E9, mean over paths, by design and structure\n')
print(dcast(fz[variant == 'E9', .(v = round(mean(full), 2)), by = .(design, construct)],
            design ~ construct, value.var = 'v'))
Kinv33 <- solve((1 - dgp$C1) * I3 + dgp$C1 * J3)[3, 3]
sep_ceil <- rbindlist(lapply(DESIGNS, function(dn) {
  use_design(dn)
  rbindlist(lapply(c('step', 'graded'), function(rule)
    rbindlist(lapply(hm_th, function(th) {
      q <- vapply(paths, function(pi) {
        u <- coupling(pi, rule, th)
        if (!any(u != 0)) return(0)
        drop(t(u) %*% solve(time_kernel(pi, 'A'), u))
      }, numeric(1))
      data.table(design = dn, rule = rule, t_half = th,
                 factored = 1 / sqrt(Kinv33 * max(q)),
                 direct = ceiling_of(rule, th, 'sep', dn))
    }))))
}))
use_design('Hybrid')
cat(sprintf('Separable ceiling, factored against direct, %d cases: max |diff| %.1e\n',
            nrow(sep_ceil), max(abs(sep_ceil$factored - sep_ceil$direct))))

## ---------------- Monte Carlo validation ----------------
## One cell: a design, a construct, a rule, a c.bm and a t1/2; both
## variants from the same simulated trials. Draws come from the full joint
## distribution.
mc_cell <- function(design, construct, rule, cbm, th, reps, seed) {
  use_design(design)
  set.seed(seed)
  Bs <- Ds9 <- Ds9b <- NULL
  for (p in seq_len(nP)) {
    pi <- paths[[p]]
    b <- dgp$build(pi, ARM[[construct]], 0, th)
    R <- b$R
    if (is_cov(rule)) R[1, i_br] <- R[i_br, 1] <- cbm * coupling(pi, rule, th)
    if (min(eigen(R, TRUE, TRUE)$values) <= 0) stop('matrix not PD')
    sds <- sqrt(diag(b$Sigma))
    X <- mvrnorm(Ns[p] * reps, b$mu, outer(sds, sds) * R)
    if (!is_cov(rule)) {
      h <- mod_profile(pi, rule, th)
      z <- (X[, 1] - b$mu[1]) / sds[1]
      X[, i_br] <- X[, i_br] + cbm * sd_BR * outer(z, h)
    }
    St <- X[, i_tv] + X[, i_pb] + X[, i_br]
    m <- function(v) matrix(v, nrow = reps, byrow = TRUE)
    Bs <- cbind(Bs, m(X[, 1]))
    Ds9 <- cbind(Ds9, m(-drop(St %*% weights(pi, 'E9'))))
    Ds9b <- cbind(Ds9b, m(-drop(St %*% weights(pi, 'E9b'))))
  }
  ols <- function(Bm, Dm, v) {
    bc <- Bm - rowMeans(Bm); dc <- Dm - rowMeans(Dm)
    sxx <- rowSums(bc^2); slope <- rowSums(bc * dc) / sxx
    rv <- (rowSums(dc^2) - slope^2 * sxx) / (N - 2)
    pval <- 2 * pt(-abs(slope / sqrt(rv / sxx)), N - 2)
    data.table(design = design, construct = construct, rule = rule, c.bm = cbm,
               t_half = th, variant = v, mc_mean = mean(slope), mc_sd = sd(slope),
               mc_power = mean(pval < ALPHA))
  }
  rbind(ols(Bs, Ds9, 'E9'), ols(Bs, Ds9b, 'E9b'))
}

## Cells: Hybrid (a) CS at c_bm = 0.25, all rules; (b) the CS heatmap
## rows at c_bm = 0.1, 0.2, 0.3, feasible cells only; (c) the two AR(1)
## structures at c_bm = 0.25; and OL+BDC and CO, all three structures at
## c_bm = 0.25. Each cell is cached; the first run's Hybrid CS cells for
## the three original rules are taken from the earlier cache. Hybrid keys
## carry no design prefix, so the cells from before the design extension
## remain valid.
mc_th <- c(0, 0.1, 0.2, 0.5, 1, 2)
jobs <- unique(rbind(
  CJ(design = 'Hybrid', construct = 'CS', rule = RULES, c.bm = c(0.1, 0.2, 0.25, 0.3),
     t_half = mc_th),
  CJ(design = 'Hybrid', construct = c('covar', 'sep'), rule = RULES, c.bm = CBM,
     t_half = mc_th),
  CJ(design = c('OL+BDC', 'CO'), construct = names(ARM), rule = RULES, c.bm = CBM,
     t_half = mc_th)))
jobs[, ceiling := ceiling_of(rule, t_half, construct, design),
     by = .(design, construct, rule, t_half)]
use_design('Hybrid')
jobs <- jobs[c.bm < ceiling]
jobs[, key := paste0(ifelse(design == 'Hybrid', '', paste0(gsub('[^A-Za-z]', '', design), '_')),
                     sprintf('%s_%s_c%.2f_t%.2f', construct, rule, c.bm, t_half))]
cell_dir <- file.path(out_dir, sprintf('mc-cells-reps%d', mc_reps))
dir.create(cell_dir, showWarnings = FALSE)
legacy <- file.path(out_dir, sprintf('e9-mc-reps%d.csv', mc_reps))
if (file.exists(legacy)) {
  lg <- fread(legacy)
  for (k in which(jobs$design == 'Hybrid' & jobs$construct == 'CS' & jobs$c.bm == CBM &
                  jobs$rule %in% c('step', 'graded', 'mean'))) {
    f <- file.path(cell_dir, paste0(jobs$key[k], '.csv'))
    if (!file.exists(f))
      fwrite(cbind(data.table(construct = 'CS', c.bm = CBM),
                   lg[rule == jobs$rule[k] & t_half == jobs$t_half[k]]), f)
  }
}
todo <- jobs[!file.exists(file.path(cell_dir, paste0(key, '.csv')))]
cat(sprintf('\nMonte Carlo: %d cells, %d cached, %d to run (%d reps, %d cores)\n',
            nrow(jobs), nrow(jobs) - nrow(todo), nrow(todo), mc_reps, cores))
if (nrow(todo)) {
  done <- mclapply(seq_len(nrow(todo)), function(k) {
    j <- todo[k]
    seed <- sum(utf8ToInt(j$key) * seq_along(utf8ToInt(j$key))) %% 1e6 + 7
    r <- tryCatch(mc_cell(j$design, j$construct, j$rule, j$c.bm, j$t_half, mc_reps, seed),
                  error = function(e) conditionMessage(e))
    if (is.character(r)) return(paste(j$key, r))
    fwrite(r, file.path(cell_dir, paste0(j$key, '.csv')))
    NA_character_
  }, mc.cores = cores, mc.preschedule = FALSE)
  errs <- Filter(function(x) !is.na(x), unlist(done))
  if (length(errs)) stop('Monte Carlo cells failed: ', paste(errs, collapse = '; '))
}
## Cells cached before the design extension have no design column.
mcs <- rbindlist(lapply(file.path(cell_dir, paste0(jobs$key, '.csv')), fread),
                 use.names = TRUE, fill = TRUE)
mcs[is.na(design), design := 'Hybrid']
cf <- grid_closed(unique(jobs[, .(design, construct, rule, c.bm, t_half)])[
  , .(variant = c('E9', 'E9b')), by = .(design, construct, rule, c.bm, t_half)])
val <- merge(cf, mcs, by = c('design', 'construct', 'rule', 'c.bm', 't_half', 'variant'))
stopifnot(nrow(val) == 2 * nrow(jobs))
val[, power_z := (mc_power - power) / sqrt(power * (1 - power) / mc_reps)]
val[, block := fifelse(design != 'Hybrid', paste0(design, ', all structures, c_bm 0.25'),
               fifelse(construct != 'CS', 'Hybrid, AR(1) structures, c_bm 0.25',
               fifelse(c.bm == CBM, 'Hybrid, CS, c_bm 0.25', 'Hybrid, CS heatmap rows')))]
fwrite(val, file.path(out_dir, 'e9-closed-form-vs-mc.csv'))
cat(sprintf('\nClosed form against Monte Carlo (%d reps per cell)\n', mc_reps))
print(val[, .(cells = .N, max_abs_z = round(max(abs(power_z)), 2),
              mean_z = round(mean(power_z), 2), beyond_2 = sum(abs(power_z) > 2),
              max_mean_diff = signif(max(abs(mean - mc_mean)), 2),
              max_sd_diff = signif(max(abs(sd - mc_sd)), 2)), by = block])
cat('\nE9 power by design, structure and rule (closed form; MC in parentheses)\n')
print(dcast(val[variant == 'E9' & c.bm == CBM,
                .(design, construct, rule, t_half,
                  p = sprintf('%.3f (%.3f)', power, mc_power))],
            design + construct + rule ~ t_half, value.var = 'p'))

## ---------------- against the lmer simulation (04, published analysis) ----
## Configurations A (step) and D (mean moderation) only: the published
## analysis is not calibrated on the AR(1) configurations (docs/36,
## Section 6.8).
cs25 <- val[construct == 'CS' & c.bm == CBM]
lm <- fread('analysis/data/quick-sim/hendrickson-problems/power-sim/power-summary-reps500.csv')
lm <- lm[design %in% DESIGNS & c.bm == CBM & arm %in% c('A', 'D'),
         .(design, rule = ifelse(arm == 'A', 'step', 'mean'), t_half,
           lmer_power = power, lmer_mean = mean_beta, lmer_sd = emp_sd)]
cmp <- merge(lm, cs25[variant == 'E9b', .(design, rule, t_half, e9b_mean = mean,
                                           e9b_sd = sd, e9b_power = power)],
             by = c('design', 'rule', 't_half'))
cmp <- merge(cmp, cs25[variant == 'E9', .(design, rule, t_half, e9_power = power)],
             by = c('design', 'rule', 't_half'))
fwrite(cmp, file.path(out_dir, 'e9-vs-lmer.csv'))
cat('\nE9b closed form against lmer (04, 500 reps, published analysis)\n')
print(cmp[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

## ---------------- figures ----------------
ink <- '#0b0b0b'; ink2 <- '#52514e'; surface <- '#fcfcfb'
rule_lab <- c(step = 'Step coupling (published)', graded = 'Graded coupling',
              mean = 'Mean moderation', mean_decayed = 'Decayed mean moderation')
## Validated categorical slots (CVD and normal-vision checks pass; three
## are below 3:1 contrast, so each rule also has its own line type and the
## whitepaper carries every value in tables).
pal <- setNames(c('#2a78d6', '#1baf7a', '#eda100', '#e87ba4'), rule_lab)
lty <- setNames(c('solid', 'solid', 'solid', 'dashed'), rule_lab)
var_lev <- c(E9b = 'E9b: baseline counted off drug (as in the published lmer)',
             E9 = 'E9: post-baseline visits only')
con_lab <- c(CS = 'Compound symmetry (published)', covar = 'AR(1), covar cross-factor',
             sep = 'AR(1), separable')
theme_fig <- function() {
  theme_minimal(base_size = 11) +
    theme(plot.background = element_rect(fill = surface, colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = '#e6e5e1', linewidth = 0.3),
          axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
          strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
          plot.title = element_text(colour = ink, face = 'bold'),
          plot.subtitle = element_text(colour = ink2),
          plot.caption = element_text(colour = ink2, hjust = 0),
          plot.title.position = 'plot', plot.caption.position = 'plot',
          legend.position = 'top', legend.justification = 'left',
          legend.title = element_blank(), legend.text = element_text(colour = ink))
}
## Half-life 0 plotted at the left edge of the log axis.
x0 <- 0.018
prep <- function(dt) {
  dt <- copy(dt)
  dt[, rule := factor(rule_lab[rule], levels = rule_lab)]
  dt[, variant := factor(var_lev[variant], levels = var_lev)]
  dt[, construct := factor(con_lab[construct], levels = con_lab)]
  dt[, x := ifelse(t_half == 0, x0, t_half)]
  dt
}
fd <- prep(curves)
pts <- prep(val[c.bm == CBM & design == 'Hybrid'])
power_plot <- function(fdd, ptd, facet, title, subtitle) {
  ggplot(fdd[t_half > 0], aes(x, power, colour = rule, linetype = rule)) +
    geom_hline(yintercept = ALPHA, colour = ink2, linewidth = 0.4, linetype = 'dashed') +
    geom_line(linewidth = 0.9) +
    geom_point(data = fdd[t_half == 0], size = 2.6, show.legend = FALSE) +
    geom_point(data = ptd, aes(y = mc_power), shape = 21, fill = 'white', size = 2.2,
               stroke = 0.9, show.legend = FALSE) +
    scale_colour_manual(values = pal) +
    scale_linetype_manual(values = lty) +
    scale_x_log10(breaks = c(x0, 0.03, 0.1, 1, 4), labels = c('0', '0.03', '0.1', '1', '4')) +
    scale_y_continuous(limits = c(0, 1)) +
    facet +
    labs(x = 'Carryover half-life (weeks, log scale; 0 at the left edge)',
         y = 'Power (alpha = 0.05)', title = title, subtitle = subtitle,
         caption = paste('Lines: closed form. Filled points: closed form at t1/2 = 0.',
                         sprintf('Open points: Monte Carlo, %d replicates. Dashed line: alpha.',
                                 mc_reps), sep = '\n')) +
    theme_fig() +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
}
g <- power_plot(fd[construct == con_lab[['CS']]], pts[construct == con_lab[['CS']]],
                facet_wrap(~ variant, nrow = 1),
                'Power of the paired-difference statistic against carryover half-life',
                sprintf('Hybrid design, N = 70, c_bm = beta_bm = %.2f, compound-symmetry construct',
                        CBM))
ggsave(file.path(fig_dir, '37-fig1-power-closed-form.png'), g, width = 10, height = 5.1,
       dpi = 200, bg = surface)
g4 <- power_plot(fd, pts, facet_grid(construct ~ variant),
                 'Power against carryover half-life, by response correlation structure',
                 sprintf('Hybrid design, N = 70, c_bm = beta_bm = %.2f', CBM))
ggsave(file.path(fig_dir, '37-fig4-power-by-structure.png'), g4, width = 10, height = 9,
       dpi = 200, bg = surface)
## Figure 5: E9 power by design (columns) and response structure (rows).
fd5 <- prep(curves_all[variant == 'E9'])
fd5[, design := factor(design, levels = DESIGNS)]
pts5 <- prep(val[c.bm == CBM & variant == 'E9'])
pts5[, design := factor(design, levels = DESIGNS)]
g5 <- power_plot(fd5, pts5, facet_grid(construct ~ design),
                 'Power against carryover half-life, by design and response structure',
                 sprintf('E9, N = 70, c_bm = beta_bm = %.2f', CBM))
ggsave(file.path(fig_dir, '37-fig5-power-by-design.png'), g5, width = 11, height = 9,
       dpi = 200, bg = surface)
cat('\nWrote Figures 1, 4 and 5\n')

## Figure 2: power against the coupling gap, with the rules marked.
var_lab <- c(E9 = 'E9: post-baseline visits only',
             E9b = 'E9b: baseline counted off drug (as in the published lmer)')
var_pal <- setNames(c('#2a78d6', '#eb6834'), var_lab)
gd <- copy(gap_curve)[, variant := factor(var_lab[variant], levels = var_lab)]
## t1/2 = 0.2 sits at gap 0.991, on top of the no-carryover point.
mk <- rule_pts[!(rule == 'graded' & t_half == 0.2)][, .(rule, t_half, gap,
                   power = approx(gap_curve[variant == 'E9']$gap,
                                  gap_curve[variant == 'E9']$power, xout = gap)$y)]
mk[, label := ifelse(rule == 'step', 'step rule, any t1/2 > 0',
                     ifelse(t_half == 0, 'no carryover',
                            sprintf('graded, t1/2 = %s', format(t_half))))]
g2 <- ggplot(gd, aes(gap, power, colour = variant)) +
  geom_hline(yintercept = ALPHA, colour = ink2, linewidth = 0.4, linetype = 'dashed') +
  geom_line(linewidth = 0.9) +
  geom_point(data = mk, aes(gap, power), inherit.aes = FALSE, colour = ink, size = 2.2) +
  ## Labels below and to the right where the E9b line runs just above
  ## the point.
  geom_text(data = mk, aes(gap, power, label = label,
                           hjust = ifelse(rule == 'step' | t_half == 2, -0.08, 1.08),
                           vjust = ifelse(rule == 'step' | t_half == 2, 1.6, -0.6)),
            inherit.aes = FALSE, colour = ink2, size = 3) +
  scale_colour_manual(values = var_pal) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Coupling gap: mean on-drug coupling minus mean off-drug coupling',
       y = 'Power (alpha = 0.05)',
       title = 'Power depends on carryover only through the coupling gap',
       subtitle = sprintf('Hybrid design, N = 70, c_bm = %.2f, compound-symmetry construct', CBM),
       caption = paste('Lines: closed form. Points: where the coupling rules sit (E9).',
                       'Gap on the post-baseline visits; under E9b the baseline keeps zero coupling.',
                       sep = '\n')) +
  theme_fig() +
  theme(strip.text = element_blank())
ggsave(file.path(fig_dir, '37-fig2-power-vs-gap.png'), g2, width = 8, height = 5.6,
       dpi = 200, bg = surface)
cat('Wrote', file.path(fig_dir, '37-fig2-power-vs-gap.png'), '\n')

## Figure 3: power heatmap per rule, c.bm x t1/2 (CS). Power is a
## magnitude, so one hue light to dark; infeasible cells are neutral
## gray and labeled, so they never read as low power.
hd <- copy(hm)
hd[, rule := factor(rule_lab[rule], levels = rule_lab)]
hd[, variant := factor(variant, levels = c('E9b', 'E9'),
                       labels = c('E9b (as published)', 'E9'))]
hd[, t_lab := factor(format(t_half, drop0trailing = TRUE),
                     levels = c('0', '0.1', '0.2', '0.5', '1', '2'))]
hd[, c_lab := factor(format(c.bm, nsmall = 1), levels = c('0.0', '0.1', '0.2', '0.3'))]
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
       y = expression(paste('Strength: ', c[bm], ' (covariance rules) or ', beta[bm],
                            ' (mean rules)')),
       title = 'Power by interaction strength and carryover half-life, per rule',
       subtitle = 'Hybrid design, N = 70, compound-symmetry construct, closed form',
       caption = paste('Gray cells: c_bm at or above the CS ceiling, so the correlation matrix',
                       'is not positive definite and no distribution exists.',
                       'Zero strength gives the test size.')) +
  theme_minimal(base_size = 11) +
  theme(plot.background = element_rect(fill = surface, colour = NA),
        panel.grid = element_blank(),
        axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
        strip.text = element_text(colour = ink, face = 'bold', hjust = 0),
        plot.title = element_text(colour = ink, face = 'bold'),
        plot.subtitle = element_text(colour = ink2),
        plot.caption = element_text(colour = ink2, hjust = 0),
        plot.title.position = 'plot', plot.caption.position = 'plot',
        legend.title = element_text(colour = ink2), legend.text = element_text(colour = ink2))
ggsave(file.path(fig_dir, '37-fig3-power-heatmap.png'), g3, width = 13, height = 5.6,
       dpi = 200, bg = surface)
cat('Wrote', file.path(fig_dir, '37-fig3-power-heatmap.png'), '\n')
