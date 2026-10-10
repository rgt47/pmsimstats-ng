## Shared code for docs/58-main-effect-testing-and-carryover-half-life.md.
##
## Loads the data-generating functions of the stress-test driver
## analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R
## by evaluating its top-level definitions (everything before the
## pre-run checks, skipping the lines that create output directories),
## as analysis/scripts/data-chosen-half-life/57-common.R does. The
## driver's build() gives the refined process (separable AR(1),
## rho = 0.7, c1 = 0.2, anchored carryover by default); the drug-response
## ceiling is hard-coded there, so the drug effect is scaled here by the
## multiplier f_m after build(), which is exact because anchored and
## recursive carryover are both linear in the drug-response mean.
##
## Two groups of helpers follow:
##   closed form: path moments of the outcome, analysis design matrices,
##     the GLS probability limit and sandwich variance of the drug
##     coefficient under a working covariance, estimand values;
##   simulation: one trial with a drug effect and optional
##     participant-by-treatment heterogeneity, main-effect fits with CR2.
## Source from the repository root.
suppressMessages({
  library(data.table); library(MASS); library(nlme); library(parallel)
  library(clubSandwich); library(mvtnorm)
})

DRIVER <- 'analysis/scripts/quick-sim/carryover-closed-form/11-strawman-stress-test.R'
OUT <- 'analysis/data/quick-sim/main-effect-half-life'
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
M_BR <- 10.98604
## The analysis grid: docs/57's grid with 0.25 replaced by the
## protocol's default analysis half-life, 0.2 week.
G5 <- c(0, 0.2, 0.5, 1, 2)

load_driver <- function(carryover = 'anchored', allocation = 'random') {
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

## A cell of the refined process with no biomarker coupling.
mk_cl <- function(t_half, rho = 0.7, decay_k = 1)
  list(construct = 'C07', rho = rho, t_half = t_half, decay_k = decay_k, cbm = 0,
       c_tv = 0, b_tv = 0, c_pb = 0, c_bl = 0, meas = FALSE, dropout = 'none')

design_pis <- function(design) {
  D <- dgp$designs[[design]]
  lapply(D$paths, function(od) dgp$path_info(D$t_wk, D$e, od))
}

## ---------------- closed form ----------------
## Outcome moments on one path: rows are baseline (t = 0) and the visits.
## q is the unit drug-response profile (BR mean / M_BR, anchored or
## recursive per `carry`); the drug effect is f_m * M_BR * q; tau adds
## tau^2 q q' (participant-specific ceiling deviation v_i).
path_moments <- function(pi, cl, f_m = 1, tau = 0, carry = 'anchored') {
  old <- CARRYOVER; on.exit(assign('CARRYOVER', old, globalenv()))
  assign('CARRYOVER', carry, globalenv())
  b <- build(pi, cl); n <- b$n
  i_tv <- 3:(2 + n); i_pb <- (3 + n):(2 + 2 * n); i_br <- (3 + 2 * n):(2 + 3 * n)
  mu <- b$mu
  q <- c(0, mu[i_br] / M_BR)
  mu[i_br] <- f_m * mu[i_br]
  L <- matrix(0, n + 1, length(mu)); L[, 2] <- 1
  for (t in seq_len(n)) L[t + 1, c(i_tv[t], i_pb[t], i_br[t])] <- -1
  list(mu = drop(L %*% mu), S = L %*% b$Sigma %*% t(L) + tau^2 * tcrossprod(q),
       q = q, w = c(0, pi$w), e = c(0, pi$e), on = c(FALSE, pi$on),
       tsd = c(0, pi$tsd), n = n + 1)
}
design_moments <- function(design, t_half, f_m = 1, tau = 0, rho = 0.7,
                           carry = 'anchored', decay_k = 1) {
  cl <- mk_cl(t_half, rho, decay_k)
  lapply(design_pis(design), path_moments, cl = cl, f_m = f_m, tau = tau, carry = carry)
}

## Exposure-decayed indicator with the pre-exposure guard (docs/56, 3.2);
## h = Inf gives 1 at every visit after first exposure.
dbc_h <- function(on, tsd, h) {
  if (h <= 0) return(as.numeric(on))
  if (!is.finite(h)) return(ifelse(on, 1, ifelse(tsd > 0, 1, 0)))
  ifelse(on, 1, ifelse(tsd > 0, 0.5^(tsd / h), 0))
}

## Analysis design matrix on one path; the drug column is named X and
## is last. F1 = t + X + eta; S = t + X; VF = visit factor + X.
design_X <- function(pm, h, type = 'F1') {
  X <- dbc_h(pm$on, pm$tsd, h)
  M <- switch(type,
    F1 = cbind(int = 1, t = pm$w, eta = pm$e, X = X),
    S  = cbind(int = 1, t = pm$w, X = X),
    VF = cbind(diag(pm$n), X = X))
  M
}

## Working covariance: random intercept plus CAR(1) on weeks.
V_rc <- function(par, w) {
  exp(par[1]) + exp(par[2]) * plogis(par[3])^abs(outer(w, w, '-'))
}

## GLS probability limit and sandwich variance of every coefficient, at
## expected path shares wts (1/P) and total N. The meat includes the
## mean misfit r_p r_p', as it does for participant contributions under
## randomized allocation.
gls_plim <- function(Xs, V, pms, wts = rep(1 / length(Xs), length(Xs)), N = 70) {
  Vi <- solve(V)
  B <- Reduce('+', Map(function(X, w) w * crossprod(X, Vi %*% X), Xs, wts))
  Bi <- solve(B)
  beta <- drop(Bi %*% Reduce('+', Map(function(X, pm, w) w * crossprod(X, Vi %*% pm$mu),
                                      Xs, pms, wts)))
  res <- Map(function(X, pm) pm$mu - drop(X %*% beta), Xs, pms)
  meat <- Reduce('+', Map(function(X, pm, r, w)
    w * crossprod(X, Vi %*% (pm$S + tcrossprod(r)) %*% Vi %*% X), Xs, pms, res, wts))
  kD <- which(colnames(Xs[[1]]) == 'X')
  om <- lapply(Xs, function(X) drop((Bi %*% t(X) %*% Vi)[kD, ]))
  list(beta = beta, b = beta[kD], var = (Bi %*% meat %*% Bi)[kD, kD] / N,
       omega = om, res = res, V = V)
}

## ML probability limit of the random-intercept + CAR(1) working
## covariance: minimizes the expected negative log-likelihood per
## participant, with the mean at its GLS limit for each candidate V.
rc_plim <- function(Xs, pms, wts = rep(1 / length(Xs), length(Xs)),
                    start = c(log(300), log(300), qlogis(0.6))) {
  w <- pms[[1]]$w
  obj <- function(par) {
    V <- V_rc(par, w)
    g <- tryCatch(gls_plim(Xs, V, pms, wts), error = function(e) NULL)
    if (is.null(g)) return(1e10)
    ld <- as.numeric(determinant(V)$modulus); Vi <- solve(V)
    sum(mapply(function(pm, r, wt) wt * (ld + sum(diag(Vi %*% (pm$S + tcrossprod(r))))),
               pms, g$res, wts))
  }
  o <- optim(start, obj, method = 'BFGS', control = list(maxit = 500))
  list(par = o$par, V = V_rc(o$par, w), su2 = exp(o$par[1]), se2 = exp(o$par[2]),
       phi = plogis(o$par[3]), conv = o$convergence)
}

## Probability limit and SE of the drug coefficient for one analysis at
## half-life h under working covariance 'indep', 'true' or 'rc'.
analysis_plim <- function(pms, h, type = 'F1', cov = 'rc', N = 70, rc = NULL) {
  Xs <- lapply(pms, design_X, h = h, type = type)
  V <- switch(cov,
    indep = diag(pms[[1]]$n),
    true = Reduce('+', lapply(pms, `[[`, 'S')) / length(pms),
    rc = { if (is.null(rc)) rc <- rc_plim(Xs, pms); rc$V })
  g <- gls_plim(Xs, V, pms, N = N)
  g$Xs <- Xs; g
}

## Two-sided rejection probability of a z test with mean b / se.
rej <- function(b, se, crit = qnorm(0.975)) pnorm(b / se - crit) + pnorm(-b / se - crit)

## Estimands (Section 4.1 of the protocol), equal path weights. 'exp':
## on-minus-off mean drug response over the post-baseline visits;
## 'bl': the same over blinded visits; 'full': -f_m * M_BR. Returned on
## the beta_D scale (negative when the drug helps).
estimand <- function(pms, f_m, kind = 'exp') {
  if (kind == 'full') return(-f_m * M_BR)
  v <- sapply(pms, function(pm) {
    sel <- seq_len(pm$n) > 1
    if (kind == 'bl') sel <- sel & pm$e == 0.5
    br <- f_m * M_BR * pm$q
    on <- sel & pm$on; off <- sel & !pm$on
    if (!any(on) || !any(off)) return(NA_real_)
    mean(br[on]) - mean(br[off])
  })
  -mean(v, na.rm = TRUE)
}
## On-minus-off mean of D_bc(h) over post-baseline visits, equal path
## weights: the factor that converts the coefficient of D_bc(h) to the
## model-implied exposure contrast (the beta_D^exp scale).
dbar <- function(pms, h) mean(sapply(pms, function(pm) {
  X <- dbc_h(pm$on, pm$tsd, h); sel <- seq_len(pm$n) > 1
  mean(X[sel & pm$on]) - mean(X[sel & !pm$on])
}))

## ---------------- simulation ----------------
## One Hybrid trial of the refined process (PIS from the driver) with
## drug-effect multiplier f_m and heterogeneity tau (v_i ~ N(0, tau^2)
## scaling the unit profile q). Returns long data with a path column.
make_trial_me <- function(cl, f_m, tau, Ns) {
  long <- rbindlist(lapply(seq_along(PIS), function(k) {
    pi <- PIS[[k]]
    b <- build(pi, cl); n <- b$n
    i_br <- (3 + 2 * n):(2 + 3 * n)
    q <- b$mu[i_br] / M_BR
    b$mu[i_br] <- f_m * b$mu[i_br]
    dr <- draw(b, pi, cl, Ns[k])
    v <- rnorm(Ns[k], 0, tau)
    Y <- dr$Y - outer(v, q)
    id <- (k - 1) * 1000 + seq_len(Ns[k])
    post <- data.table(ptID = rep(id, times = n), path = k, Sx = as.vector(Y),
                       vis = rep(seq_len(n), each = Ns[k]), t = rep(pi$w, each = Ns[k]),
                       Db = rep(pi$on, each = Ns[k]), tsd = rep(pi$tsd, each = Ns[k]),
                       De = rep(pi$e, each = Ns[k]))
    rbind(data.table(ptID = id, path = k, Sx = dr$BL, vis = 0L, t = 0, Db = FALSE,
                     tsd = 0, De = 0), post)
  }))
  setorder(long, ptID, vis)
  long[, `:=`(ptID = factor(ptID),
              visit = factor(as.character(t), levels = as.character(c(0, PIS[[1]]$w))))]
  long
}

## Relabel participants to new paths (vector over participants in the
## order of levels(d$ptID)); outcome vectors are kept, design rows change.
relabel <- function(d, newpath) {
  d <- copy(d)
  on_m <- sapply(PIS, function(pi) c(FALSE, pi$on))
  tsd_m <- sapply(PIS, function(pi) c(0, pi$tsd))
  pp <- newpath[as.integer(d$ptID)]
  d[, `:=`(path = pp, Db = on_m[cbind(vis + 1L, pp)], tsd = tsd_m[cbind(vis + 1L, pp)])]
  d
}

FMS <- list(F1 = Sx ~ t + X + De, S = Sx ~ t + X, VF = Sx ~ visit + X)
F1_0 <- Sx ~ t + De
with_xh <- function(d, h) { d <- copy(d); d[, X := dbc_h(Db, tsd, h)]; d }
na5 <- c(est = NA_real_, se = NA_real_, t = NA_real_, df = NA_real_, p = NA_real_)
cr2x <- function(fit, target = 'X') {
  if (is.null(fit)) return(na5)
  ct <- tryCatch(clubSandwich::coef_test(fit, vcov = 'CR2', test = 'Satterthwaite'),
                 error = function(e) NULL)
  if (is.null(ct)) return(na5)
  k <- which(ct$Coef == target)
  c(est = ct$beta[k], se = ct$SE[k], t = ct$beta[k] / ct$SE[k], df = ct$df_Satt[k],
    p = ct$p_Satt[k])
}
ll_or_na <- function(f) if (is.null(f)) NA_real_ else as.numeric(logLik(f))

## Participant-level influence of the drug coefficient for an lme fit
## with random intercept + CAR(1), used for the correlation of the T(h).
infl_X <- function(fit, d) {
  su2 <- as.numeric(VarCorr(fit)['(Intercept)', 'Variance'])
  phi <- as.numeric(coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  s2 <- fit$sigma^2
  Xm <- model.matrix(formula(fit)[-2], d)
  r <- d$Sx - drop(Xm %*% fixef(fit))
  idx <- split(seq_len(nrow(d)), d$ptID)
  parts <- lapply(idx, function(ii) {
    tt <- d$t[ii]; V <- su2 + s2 * phi^abs(outer(tt, tt, '-')); Vi <- solve(V)
    list(B = crossprod(Xm[ii, , drop = FALSE], Vi %*% Xm[ii, , drop = FALSE]),
         u = drop(crossprod(Xm[ii, , drop = FALSE], Vi %*% r[ii])))
  })
  Bi <- solve(Reduce('+', lapply(parts, `[[`, 'B')))
  kD <- which(colnames(Xm) == 'X')
  sapply(parts, function(p) drop(Bi %*% p$u)[kD])
}

## The grid of F1 fits on one trial: REML + CR2 test of X, ML
## log-likelihood with X, the model-implied exposure contrast factor,
## and the influence vectors for the correlation of the T(h).
grid_me <- function(d, grid = G5, form = FMS$F1, infl = TRUE) {
  pp <- unique(d[, .(ptID, path)])$path
  res <- lapply(grid, function(h) {
    dd <- with_xh(d, h)
    fr <- fit_lme(form, dd)
    ct <- cr2x(fr)
    ## observed on-minus-off mean of D_bc(h), post-baseline, averaged
    ## over participants
    db <- dd[vis > 0, .(m = mean(X[Db]) - mean(X[!Db])), by = ptID][, mean(m)]
    list(row = c(h = h, ct, ll1 = ll_or_na(fit_lme(form, dd, 'ML')), dbar = db),
         u = if (infl && !is.null(fr)) infl_X(fr, dd) else NULL)
  })
  g <- as.data.table(do.call(rbind, lapply(res, `[[`, 'row')))
  U <- do.call(cbind, lapply(res, `[[`, 'u'))
  list(g = g, U = U)
}

## Procedures computed from one grid: fixed h, AIC with the tested
## term (F1A, equivalent to argmax LR(h)), AIC-weight averaging on the
## beta_D^exp scale, and max-|t| with Bonferroni and with the
## multivariate-normal calibration from the influence correlation.
procedures <- function(G, h_fixed = 0.2) {
  g <- G$g
  out <- list()
  k0 <- which(g$h == h_fixed)
  out$fixed <- c(h = h_fixed, est_exp = g$est[k0] * g$dbar[k0], t = g$t[k0], p = g$p[k0])
  ka <- which.max(g$ll1)
  out$aic <- c(h = g$h[ka], est_exp = g$est[ka] * g$dbar[ka], t = g$t[ka], p = g$p[ka])
  w <- exp(g$ll1 - max(g$ll1)); w <- w / sum(w)
  ee <- g$est * g$dbar; se <- g$se * g$dbar
  em <- sum(w * ee); sm <- sum(w * sqrt(se^2 + (ee - em)^2)); dfm <- sum(w * g$df)
  out$avg <- c(h = sum(w * g$h), est_exp = em, t = em / sm, p = 2 * pt(-abs(em / sm), dfm))
  km <- which.max(abs(g$t)); tmax <- abs(g$t[km])
  out$bonf <- c(h = g$h[km], est_exp = g$est[km] * g$dbar[km], t = tmax,
                p = min(1, length(g$h) * min(g$p)))
  pm <- NA_real_
  if (!is.null(G$U) && ncol(G$U) == nrow(g)) {
    R <- cov2cor(crossprod(G$U))
    pm <- 1 - as.numeric(pmvnorm(lower = rep(-tmax, nrow(g)), upper = rep(tmax, nrow(g)),
                                 corr = R))
  }
  out$maxt_mvn <- c(h = g$h[km], est_exp = g$est[km] * g$dbar[km], t = tmax, p = pm)
  rbindlist(lapply(names(out), function(k) as.data.table(as.list(out[[k]]))[, method := k]))
}

## ---------------- randomization tests with label-invariant weights ----------------
## The working covariance is estimated once, from the null model
## Sx ~ t + De, whose design is the same on every path; its estimate is
## therefore invariant to relabeling, and a GLS statistic computed with
## it is a fixed function of (Y, labels). Relabeling the paths and
## recomputing the statistic gives an exact test under the sharp null.
## Each relabeling costs only linear algebra: B depends on the path
## counts, and the score on path sums of Y. Complete data (9 rows per
## participant) are assumed, as in the trials simulated here.
null_V <- function(d) {
  f0 <- fit_lme(F1_0, d)
  su2 <- as.numeric(VarCorr(f0)['(Intercept)', 'Variance'])
  phi <- as.numeric(coef(f0$modelStruct$corStruct, unconstrained = FALSE))
  w <- c(0, PIS[[1]]$w)
  list(V = su2 + f0$sigma^2 * phi^abs(outer(w, w, '-')), su2 = su2, s2 = f0$sigma^2, phi = phi)
}
perm_setup <- function(d, V, hs = G5, types = c('F1', 'VF')) {
  Y <- matrix(d$Sx, ncol = 9, byrow = TRUE)
  path <- unique(d[, .(ptID, path)])$path
  Vi <- solve(V)
  pps <- lapply(PIS, function(pi) list(w = c(0, pi$w), e = c(0, pi$e), on = c(FALSE, pi$on),
                                        tsd = c(0, pi$tsd), n = 9))
  specs <- CJ(type = types, h = hs, sorted = FALSE)
  parts <- lapply(seq_len(nrow(specs)), function(j) {
    Xs <- lapply(pps, design_X, h = specs$h[j], type = specs$type[j])
    list(Xs = Xs, As = lapply(Xs, function(X) crossprod(X, Vi)),
         AX = lapply(Xs, function(X) crossprod(X, Vi %*% X)), kD = ncol(Xs[[1]]))
  })
  ## per-path on-minus-off mean of D_bc(h) over post-baseline visits
  for (j in seq_along(parts)) parts[[j]]$dbar_p <- sapply(seq_along(pps), function(k) {
    x <- parts[[j]]$Xs[[k]][-1, 'X']; on <- pps[[k]]$on[-1]; mean(x[on]) - mean(x[!on]) })
  list(Y = Y, path = path, specs = specs, parts = parts, N = nrow(Y))
}
## Statistics for one labeling: per (type, h) the GLS estimate, its
## CR0 sandwich t, its model-based t (whose square is the fixed-V
## likelihood-ratio statistic, so its argmax over h is the AIC choice),
## and the estimate on the beta_D^exp scale.
perm_stats <- function(S, path) {
  P <- length(PIS)
  idx <- lapply(seq_len(P), function(k) which(path == k))
  Ysum <- lapply(idx, function(ii) colSums(S$Y[ii, , drop = FALSE]))
  Np <- sapply(idx, length)
  out <- vapply(S$parts, function(pt) {
    B <- Reduce('+', Map(function(ax, n) n * ax, pt$AX, Np))
    Bi <- solve(B)
    beta <- drop(Bi %*% Reduce('+', Map(function(A, ys) drop(A %*% ys), pt$As, Ysum)))
    meat <- Reduce('+', lapply(seq_len(P), function(k) {
      ii <- idx[[k]]; if (!length(ii)) return(0)
      E <- sweep(S$Y[ii, , drop = FALSE], 2, drop(pt$Xs[[k]] %*% beta))
      crossprod(E %*% t(pt$As[[k]]))
    }))
    kD <- pt$kD
    vc <- (Bi %*% meat %*% Bi)[kD, kD] * S$N / (S$N - ncol(B))
    db <- sum(Np * pt$dbar_p) / sum(Np)
    b <- unname(beta[kD])
    c(est = b, t = b / sqrt(vc), tmb = b / sqrt(Bi[kD, kD]), est_exp = b * db)
  }, numeric(4))
  out
}
## Summary statistics from perm_stats output: F1 and VF at h = 0.2
## (|t| and |estimate|), max |t| over the grid, and |t| at the AIC
## choice (argmax |tmb|).
perm_summary <- function(S, st, h_fixed = 0.2) {
  sp <- S$specs
  res <- c()
  for (ty in unique(sp$type)) {
    j <- which(sp$type == ty)
    k0 <- j[sp$h[j] == h_fixed]
    if (!length(k0)) k0 <- j[1]
    ka <- j[which.max(abs(st['tmb', j]))]
    res[paste0(ty, '_fix_t')] <- abs(st['t', k0])
    res[paste0(ty, '_fix_est')] <- abs(st['est', k0])
    res[paste0(ty, '_max_t')] <- max(abs(st['t', j]))
    res[paste0(ty, '_aic_t')] <- abs(st['t', ka])
    res[paste0(ty, '_aic_h')] <- sp$h[ka]
  }
  res
}
## Exact randomization p-values (fixed path counts: permutations of the
## labels; or re-randomization by the trial's allocation rule).
rand_test <- function(S, M = 199, scheme = 'permute', min_n = 2L) {
  obs <- perm_summary(S, perm_stats(S, S$path))
  sims <- vapply(seq_len(M), function(m) {
    if (scheme == 'permute') pp <- sample(S$path) else repeat {
      pp <- sample.int(length(PIS), S$N, replace = TRUE)
      if (all(tabulate(pp, length(PIS)) >= min_n)) break
    }
    perm_summary(S, perm_stats(S, pp))
  }, numeric(length(obs)))
  keys <- names(obs)[!grepl('_h$', names(obs))]
  p <- sapply(keys, function(k) (1 + sum(sims[k, ] >= obs[k] - 1e-12)) / (M + 1))
  list(obs = obs, p = p, sims = sims)
}
