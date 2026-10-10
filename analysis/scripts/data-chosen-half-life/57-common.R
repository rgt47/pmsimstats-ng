## Shared code for docs/57-data-chosen-carryover-half-life.md.
##
## Loads the data-generating process, the trial builder and the fitting
## helpers of the stress-test driver
## analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R
## by evaluating its top-level definitions (everything before the
## pre-run checks, skipping the lines that create output directories),
## so that the trials generated here are those of the driver, seed for
## seed. Then defines the half-life procedures compared in docs/57.
## Source from the repository root.
suppressMessages({
  library(data.table); library(MASS); library(nlme); library(parallel)
  library(clubSandwich)
})

DRIVER <- 'analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R'

load_driver <- function(carryover, allocation) {
  ex <- parse(DRIVER)
  heads <- vapply(ex, function(e) paste(deparse(e), collapse = ' '), character(1))
  stop_at <- which(grepl('^setups <- ', heads))[1]
  for (i in seq_len(stop_at - 1)) {
    if (grepl('dir.create|^out_dir <-|^run_tag <-|^cell_dir <-', heads[i])) next
    eval(ex[[i]], globalenv())
  }
  assign('CARRYOVER', carryover, globalenv())
  assign('ALLOCATION', allocation, globalenv())
  invisible(TRUE)
}

G5 <- c(0, 0.25, 0.5, 1, 2)
FM1 <- Sx ~ bmc + t + X + bmc:X + De + bmc:De
FM0 <- Sx ~ bmc + t + X + De + bmc:De
TARGET <- 'bmc:X'
na5 <- c(est = NA_real_, se = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_)

## CR2 standard error with Satterthwaite degrees of freedom, as in the
## driver's cr2(), returning the t statistic and df as well.
cr2t <- function(fit) {
  if (is.null(fit)) return(na5)
  ct <- tryCatch(clubSandwich::coef_test(fit, vcov = 'CR2', test = 'Satterthwaite'),
                 error = function(e) NULL)
  if (is.null(ct)) return(na5)
  k <- which(ct$Coef == TARGET)
  c(est = ct$beta[k], se = ct$SE[k], t = ct$beta[k] / ct$SE[k],
    df = ct$df_Satt[k], p = ct$p_Satt[k])
}
ll_or_na <- function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f))

## Per analysis half-life h: ML log-likelihoods with (ll1) and without
## (ll0) the interaction, and the REML + CR2 test of the interaction.
## All candidates have the same parameter count, so AIC differences are
## -2 times log-likelihood differences and the AIC rule is argmax ll.
grid_fits <- function(d, grid = G5, test = TRUE) {
  rbindlist(lapply(grid, function(h) {
    dd <- with_x(d, h)
    out <- list(h = h, ll1 = ll_or_na(fit_lme(FM1, dd, 'ML')),
                ll0 = ll_or_na(fit_lme(FM0, dd, 'ML')))
    if (test) out <- c(out, as.list(cr2t(fit_lme(FM1, dd))))
    out
  }))
}

pick <- function(g, k, name) {
  if (length(k) == 0 || is.na(k)) return(data.table(method = name, h = NA_real_,
    est = NA_real_, se = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_))
  data.table(method = name, h = g$h[k], est = g$est[k], se = g$se[k],
             t = g$t[k], df = g$df[k], p = g$p[k])
}
wmax <- function(x) { if (all(is.na(x))) NA_integer_ else which.max(x) }

## Methods that use only the full-sample grid fits.
grid_methods <- function(g, h_fixed = 0.5) {
  out <- rbind(pick(g, which(g$h == h_fixed), 'F1'),
               pick(g, wmax(g$ll1), 'F1A'),
               pick(g, wmax(g$ll0), 'F1A0'))
  ## Smoothed-AIC model averaging (Buckland et al. 1997), flat weights
  ## over the grid: a frequentist stand-in for marginalizing over t1/2.
  ok <- is.finite(g$ll1) & is.finite(g$est)
  if (any(ok)) {
    w <- exp(g$ll1[ok] - max(g$ll1[ok])); w <- w / sum(w)
    eb <- sum(w * g$est[ok])
    sb <- sum(w * sqrt(g$se[ok]^2 + (g$est[ok] - eb)^2))
    dfb <- sum(w * g$df[ok])
    out <- rbind(out, data.table(method = 'MA', h = sum(w * g$h[ok]), est = eb, se = sb,
                                 t = eb / sb, df = dfb, p = 2 * pt(-abs(eb / sb), dfb)))
  }
  ## Bonferroni over the grid, and the raw max-|t| for bootstrap calibration.
  km <- if (all(is.na(g$t))) NA_integer_ else which.max(abs(g$t))
  bon <- pick(g, km, 'MaxT')
  bon$p <- min(1, length(g$h) * min(g$p, na.rm = TRUE))
  bon$t <- abs(bon$t)
  out <- rbind(out, bon)
  out
}

## Path of each participant from the driver's ids, (k - 1) * 1000 + j.
path_of <- function(d) floor((as.integer(as.character(d$ptID)) - 1) / 1000) + 1

## Sample splitting and cross-fitting. Participants are split into two
## halves at random within path; the half-life is chosen on one half by
## the F1A0 criterion (ML, no interaction term) and the interaction is
## tested on the other half at that value, then the roles are swapped.
split_methods <- function(d, grid = G5) {
  ids <- unique(data.table(ptID = d$ptID, path = path_of(d)))
  ids[, half := sample(rep(1:2, length.out = .N)), by = path]
  dA <- d[ptID %in% ids[half == 1, ptID]]
  dB <- d[ptID %in% ids[half == 2, ptID]]
  choose <- function(dd) {
    ll <- vapply(grid, function(h) ll_or_na(fit_lme(FM0, with_x(dd, h), 'ML')), numeric(1))
    if (all(is.na(ll))) NA_real_ else grid[which.max(ll)]
  }
  hA <- choose(dA); hB <- choose(dB)
  tB <- if (is.na(hA)) na5 else cr2t(fit_lme(FM1, with_x(dB, hA)))
  tA <- if (is.na(hB)) na5 else cr2t(fit_lme(FM1, with_x(dA, hB)))
  single <- data.table(method = 'Split', h = hA, est = tB['est'], se = tB['se'],
                       t = tB['t'], df = tB['df'], p = tB['p'])
  est <- (tA['est'] + tB['est']) / 2
  se <- sqrt(tA['se']^2 + tB['se']^2) / 2
  df <- (tA['se']^2 + tB['se']^2)^2 / (tA['se']^4 / tA['df'] + tB['se']^4 / tB['df'])
  avg <- data.table(method = 'CrossAvg', h = (hA + hB) / 2, est = est, se = se,
                    t = est / se, df = df, p = 2 * pt(-abs(est / se), df))
  cct <- 0.5 * tan((0.5 - tA['p']) * pi) + 0.5 * tan((0.5 - tB['p']) * pi)
  cau <- data.table(method = 'CrossCauchy', h = (hA + hB) / 2, est = est, se = se,
                    t = NA_real_, df = NA_real_, p = 0.5 - atan(cct) / pi)
  rbind(single, avg, cau)
}

## Parametric bootstrap from the null model (no interaction) refitted by
## REML at the half-life chosen by its own ML criterion (F1A0's choice):
## Y* = fixed part + random intercept + CAR(1) errors, with the
## estimated variance parameters, on the same participants, visits,
## biomarker values and paths.
null_fit <- function(d, h0) fit_lme(FM0, with_x(d, h0))
sim_null <- function(fit, d, h0) {
  dd <- with_x(d, h0)
  mu <- drop(model.matrix(FM0[-2], dd) %*% nlme::fixef(fit))
  tau <- as.numeric(nlme::VarCorr(fit)['(Intercept)', 'StdDev'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  sig <- fit$sigma
  y <- mu
  for (idx in split(seq_len(nrow(dd)), dd$ptID)) {
    tt <- dd$t[idx]
    S <- tau^2 + sig^2 * phi^abs(outer(tt, tt, '-'))
    y[idx] <- mu[idx] + drop(crossprod(chol(S), rnorm(length(idx))))
  }
  ds <- copy(d); ds[, Sx := y]; ds
}
## Bootstrap statistics for one resample: |t| for F1, F1A, F1A0 and the
## max over the grid, each rerunning its own selection on Y*.
boot_stats <- function(ds, grid = G5) {
  g <- grid_fits(ds, grid)
  m <- grid_methods(g)
  setNames(abs(m[match(c('F1', 'F1A', 'F1A0', 'MaxT'), method), t]),
           c('F1', 'F1A', 'F1A0', 'MaxT'))
}

## One trial, all methods. aux_seed seeds the split and the bootstrap
## draw in a separate stream, restoring the trial stream afterwards so
## that trial r here is trial r of the driver.
analyze_trial <- function(d, aux_seed, n_boot = 1L) {
  old <- .Random.seed
  on.exit(assign('.Random.seed', old, globalenv()))
  set.seed(aux_seed)
  g <- grid_fits(d)
  m <- rbind(grid_methods(g), split_methods(d))
  h0 <- m[method == 'F1A0', h]
  bs <- NULL
  if (n_boot > 0 && is.finite(h0)) {
    nf <- null_fit(d, h0)
    if (!is.null(nf)) bs <- rbindlist(lapply(seq_len(n_boot), function(b)
      as.list(boot_stats(sim_null(nf, d, h0)))))
  }
  list(grid = g, methods = m, boot = bs)
}
