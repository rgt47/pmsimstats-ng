## analysis/scripts/carryover-sensitivity/diagnose-kappa-analytic.R
##
## Population (pseudo-true) calibration ratio kappa for paper 10
## (analysis/report/10-interaction-test-calibration/), computed
## exactly from the data-generating covariance and mean of Sx on each
## Hybrid path. No Monte Carlo error. Adopted from the 2026-10-28
## referee's p10-kappa-lib.R / p10-kappa-analytic.R, with a mean-model
## factor added.
##
## Method. At c_bm = 0 the biomarker b is N(0, 1) and independent of
## Sx, so the fixed-effect information is block diagonal: the block
## for the biomarker terms (bm, bm:D, and bm:De when present) is
## E(b^2) Z' V0^{-1} Z = Z' V0^{-1} Z with Z the drug-term columns.
## The variance parameters of the working covariance V0 are set to
## their pseudo-true values, the minimizers of the Kullback-Leibler
## criterion log|V0| + tr(V0^{-1} S), with S = V + r r' the second
## moment of Sx about the pseudo-true working mean and r the
## path-specific mean misfit. The model-based variance of the
## interaction is v0 = [A^{-1}]_{22} with A = sum_p w_p Z_p' V0^{-1}
## Z_p, and its true variance is the sandwich v = [A^{-1} B A^{-1}]_22
## with B = sum_p w_p Z_p' V0^{-1} S_p V0^{-1} Z_p. kappa = sqrt(v0/v).
##
## Factors.
##   DGP     stored  : c.cfct = 0.1 and the original (compounding)
##                     carryover rule, as the 2026-05-20 runs used
##           refined : c.cfct = 0.2 = c.cf1t and carryover that decays
##                     from the last on-drug mean (current package)
##   mean    fitted  : the E-specification mean, beta_0 + beta_T t + D
##                     (+ L for E3); the misfit r is carried into S
##           expect  : fitted mean plus the design expectancy De and
##                     bm:De (the F1 mean structure at the E drug term)
##           visit   : a free mean for each of the eight visits
##                     (common to all paths) plus the drug term, the
##                     usual MMRM mean
##           correct : r = 0 (a saturated path-by-visit mean, or the
##                     true mean; identical in the population)
##   cov     M0 to M3 (independence, random intercept, CAR(1),
##           random intercept plus CAR(1))
##   spec    E1, E2, E3
##
## Equal path weights (1/4) are used; the simulation allocates
## 18, 18, 17, 17 at N = 70.
##
## Usage (from the repository root):
##   Rscript analysis/scripts/carryover-sensitivity/diagnose-kappa-analytic.R
##
## Writes output/paper10-refined/diag-kappa-analytic.rds

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

repo_root <- here::here()
source(file.path(repo_root, 'implementations/tidyverse/R/functions.R'))
source(file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/simulation-core.R'))

rp <- default_resp_param()
bp <- default_baseline_param()
ds <- build_design_set('Hybrid')

mk_mp <- function(cfct, cbm = 0) list(
  N = 70, c.bm = cbm, carryover_t1half = 1,
  carryover_form = 'exponential', weibull_shape = 1,
  dgp_architecture = 'mvn',
  c.tv = 0.7, c.pb = 0.7, c.br = 0.7, c.cf1t = 0.2, c.cfct = cfct)

## The original carryover rule recurses on the already-carried mean;
## the current package decays from the last on-drug mean.
recursive_br <- function(td, halflife = 1) {
  trd <- prepare_trial_data(td)
  base <- mod_gompertz(trd$tod, rp$max[3], rp$disp[3], rp$rate[3])
  mu <- base
  for (k in seq_along(mu)) {
    if (!trd$on_drug[k] && trd$tsd[k] > 0 && k > 1)
      mu[k] <- base[k] + mu[k - 1] * carryover_decay(trd$tsd[k], halflife)
  }
  mu
}

## Mean and covariance of Sx = BL - TV - PB - BR on one path.
path_moments <- function(td, mp, carry = c('anchored', 'recursive')) {
  carry <- match.arg(carry)
  s <- build_sigma_matrix(mp, rp, bp, td, dgp_architecture = 'mvn',
                          stop_on_nonpd = TRUE)
  lab <- s$labels
  nt <- nrow(td)
  tpn <- td$timepoint_name
  A <- matrix(0, nt, length(lab))
  for (k in seq_len(nt)) {
    A[k, lab == 'BL'] <- 1
    for (f in c('tv', 'pb', 'br'))
      A[k, lab == paste(tpn[k], f, sep = '.')] <- -1
  }
  mu_all <- s$means
  if (carry == 'recursive') {
    idx <- match(paste(tpn, 'br', sep = '.'), lab)
    mu_all[idx] <- recursive_br(td)
  }
  list(V = A %*% s$sigma %*% t(A), mu = drop(A %*% mu_all))
}

## Long-form design columns per path, in visit order.
long_p <- map(seq_along(ds), function(i) {
  f <- tibble(ptID = 1L, path = 1L, bm = 0)
  for (nm in ds[[i]]$timepoint_name) f[[nm]] <- 0
  dl <- prepare_long_data(f, ds[i], carryover_t1half = 1)
  stopifnot(identical(dl$timepoint_name, ds[[i]]$timepoint_name))
  dl
})

V0_fun <- function(par, tt, model) {
  n <- length(tt)
  R <- outer(tt, tt, function(a, b) par$rho_a^abs(a - b))
  switch(model,
    M0 = par$s2 * diag(n),
    M1 = par$su2 * matrix(1, n, n) + par$s2 * diag(n),
    M2 = par$s2 * R,
    M3 = par$su2 * matrix(1, n, n) + par$s2 * R)
}

unpack <- function(th, model) {
  out <- list(s2 = exp(th[1]), su2 = 0, rho_a = 0)
  if (model %in% c('M1', 'M3')) out$su2 <- exp(th[2])
  if (model %in% c('M2', 'M3')) out$rho_a <- plogis(th[length(th)])
  out
}

kappa_pop <- function(mom, spec, model, mean_model) {
  w <- rep(1 / length(ds), length(ds))
  D <- function(dl) if (spec == 'E2') dl$Dbc else dl$Db
  use_r <- mean_model != 'correct'
  Mlist <- map(long_p, function(dl) {
    cols <- if (mean_model == 'visit')
      cbind(diag(nrow(dl)), D(dl)) else cbind(1, dl$t, D(dl))
    if (spec == 'E3') cols <- cbind(cols, dl$L)
    if (mean_model == 'expect') cols <- cbind(cols, dl$De)
    cols
  })
  Zlist <- map(long_p, function(dl) {
    z <- cbind(1, D(dl))
    if (mean_model == 'expect') z <- cbind(z, dl$De)
    z
  })
  ttl <- map(long_p, ~ .x$t)
  pt_mean <- function(Vi) {
    G <- Reduce(`+`, map(seq_along(ds), ~ w[.x] * t(Mlist[[.x]]) %*%
      Vi[[.x]] %*% Mlist[[.x]]))
    h <- Reduce(`+`, map(seq_along(ds), ~ w[.x] * t(Mlist[[.x]]) %*%
      Vi[[.x]] %*% mom[[.x]]$mu))
    solve(G, h)
  }
  resid <- function(p, b) {
    if (use_r) mom[[p]]$mu - Mlist[[p]] %*% b
    else matrix(0, length(mom[[p]]$mu), 1)
  }
  obj <- function(th) {
    par <- unpack(th, model)
    V0s <- map(ttl, ~ V0_fun(par, .x, model))
    Vi <- map(V0s, solve)
    b <- pt_mean(Vi)
    val <- 0
    for (p in seq_along(ds)) {
      r <- resid(p, b)
      S <- mom[[p]]$V + r %*% t(r)
      val <- val + w[p] * (determinant(V0s[[p]])$modulus +
        sum(diag(Vi[[p]] %*% S)))
    }
    as.numeric(val)
  }
  th0 <- switch(model, M0 = log(200), M1 = c(log(150), log(100)),
                M2 = c(log(250), qlogis(0.8)),
                M3 = c(log(150), log(100), qlogis(0.7)))
  opt <- if (length(th0) == 1)
    optimize(function(x) obj(x), c(0, 12)) else optim(th0, obj,
      method = 'BFGS', control = list(maxit = 500, reltol = 1e-12))
  th <- if (length(th0) == 1) opt$minimum else opt$par
  par <- unpack(th, model)
  V0s <- map(ttl, ~ V0_fun(par, .x, model))
  Vi <- map(V0s, solve)
  b <- pt_mean(Vi)
  A <- Reduce(`+`, map(seq_along(ds), ~ w[.x] * t(Zlist[[.x]]) %*%
    Vi[[.x]] %*% Zlist[[.x]]))
  B <- Reduce(`+`, map(seq_along(ds), function(p) {
    r <- resid(p, b)
    w[p] * t(Zlist[[p]]) %*% Vi[[p]] %*% (mom[[p]]$V + r %*% t(r)) %*%
      Vi[[p]] %*% Zlist[[p]]
  }))
  Ai <- solve(A)
  v0 <- Ai[2, 2]
  v <- (Ai %*% B %*% Ai)[2, 2]
  rms_r <- if (use_r) sqrt(mean(unlist(map(seq_along(ds),
    ~ resid(.x, b)^2)))) else 0
  tibble(spec = spec, model = model, mean_model = mean_model,
         kappa = sqrt(v0 / v), se_true_N70 = sqrt(v / 70),
         se_model_N70 = sqrt(v0 / 70),
         sigma_u = sqrt(par$su2), sigma_e = sqrt(par$s2),
         rho_a = par$rho_a, rms_misfit = rms_r,
         size_pred = 2 * pnorm(-sqrt(v0 / v) * qnorm(0.975)))
}

run_scenario <- function(cfct, carry, label) {
  mom <- map(ds, ~ path_moments(.x, mk_mp(cfct), carry))
  grid <- expand_grid(spec = c('E1', 'E2', 'E3'),
                      model = c('M0', 'M1', 'M2', 'M3'),
                      mean_model = c('fitted', 'expect', 'visit',
                                     'correct'))
  out <- pmap_dfr(grid, function(spec, model, mean_model)
    kappa_pop(mom, spec, model, mean_model))
  out$dgp <- label
  out
}

## Positive definiteness of every Hybrid path under the current
## package (stop_on_nonpd = TRUE refuses a non-PD matrix).
pd_check <- expand_grid(c_bm = c(0, 0.45), c_cfct = c(0.1, 0.2)) |>
  mutate(min_eig = map2_dbl(c_bm, c_cfct, function(cb, cf) {
    min(map_dbl(ds, function(td) {
      s <- build_sigma_matrix(mk_mp(cf, cb), rp, bp, td,
                              dgp_architecture = 'mvn',
                              stop_on_nonpd = TRUE)$sigma
      pos <- diag(s) > 0
      min(eigen(s[pos, pos], only.values = TRUE)$values)
    }))
  }))

## Correctly specified check: homoscedastic placebo component, so the
## covariance of Sx is exactly random intercept plus a separable CAR(1)
## sum; with the correct mean M3 should return kappa = 1.
homo_check <- local({
  pm_h <- function(td, mp) {
    td2 <- td
    td2$e <- 1
    s <- build_sigma_matrix(mp, rp, bp, td2, dgp_architecture = 'mvn',
                            stop_on_nonpd = TRUE)
    lab <- s$labels
    nt <- nrow(td)
    tpn <- td$timepoint_name
    A <- matrix(0, nt, length(lab))
    for (k in seq_len(nt)) {
      A[k, lab == 'BL'] <- 1
      for (f in c('tv', 'pb', 'br'))
        A[k, lab == paste(tpn[k], f, sep = '.')] <- -1
    }
    list(V = A %*% s$sigma %*% t(A), mu = rep(0, nt))
  }
  momh <- map(ds, ~ pm_h(.x, mk_mp(0.2)))
  map_dfr(c('E1', 'E2'), ~ kappa_pop(momh, .x, 'M3', 'correct'))
})

t0 <- Sys.time()
res <- bind_rows(
  run_scenario(0.1, 'recursive', 'stored'),
  run_scenario(0.2, 'anchored', 'refined'))
elapsed <- as.numeric(Sys.time() - t0, units = 'secs')

misfit <- map_dfr(c('stored', 'refined'), function(lab) {
  cf <- if (lab == 'stored') 0.1 else 0.2
  ca <- if (lab == 'stored') 'recursive' else 'anchored'
  mom <- map(ds, ~ path_moments(.x, mk_mp(cf), ca))
  tibble(dgp = lab, path = seq_along(ds),
         sd_of_path_mean = map_dbl(mom, ~ sqrt(mean((.x$mu -
           mean(.x$mu))^2))))
})

print(as.data.frame(res |> filter(model %in% c('M1', 'M3')) |>
  mutate(across(where(is.numeric), ~ round(.x, 4)))))
cat('\nCorrectly specified check (expect kappa = 1, sigma_u = 10,',
    'rho_a = 0.7):\n')
print(as.data.frame(homo_check), digits = 4)
cat('\nPD check:\n')
print(as.data.frame(pd_check), digits = 4)

out_dir <- file.path(repo_root,
  'analysis/scripts/carryover-sensitivity/output/paper10-refined')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(list(kappa = res, homo_check = homo_check, pd_check = pd_check,
             misfit = misfit, elapsed_secs = elapsed,
             date = Sys.Date()),
        file.path(out_dir, 'diag-kappa-analytic.rds'))
cat(sprintf('Wrote %s (%.0f s)\n',
            file.path(out_dir, 'diag-kappa-analytic.rds'), elapsed))
