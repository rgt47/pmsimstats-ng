# pilot-core.R
#
# Shared core for the paper 04 linear-DGP pilot: data-generating
# mechanisms, analysis rules, and the closed-form (GLS) performance
# calculations. Sourced by production-run.R, check-pilot.R and
# report.Rmd. No random-number seeding happens in this file.
#
# Units: weeks. One period is one week. The carryover half-life
# t_half is in weeks (NOTATION.md, Units).

PILOT <- list(
  n_patients = 35L,
  n_periods = 8L,
  mu0 = 8.0,
  sigma_u = 1.5,
  sigma_w = 2.3,
  rho = 0.3,
  alpha = 0.05,
  t_half_a = 0.2,
  r_bl = 0.8,
  oracle_cut = 0.05
)

#' Mean carryover fraction in an off-drug period
#'
#' The outcome of a period is the weekly mean, so the carryover is the
#' average of 2^(-s / t_half) over the period. `k` is the number of
#' off-drug periods since the last on-drug period (0 for the first
#' off-drug period). A half-life of 0 means no carryover.
#'
#' @param k Integer vector, completed off-drug periods since last drug.
#' @param t_half Half-life in weeks.
#' @param period Period length in weeks.
carry_fraction <- function(k, t_half, period = 1) {
  if (t_half <= 0) {
    return(rep(0, length(k)))
  }
  t_half / (period * log(2)) *
    (2^(-k * period / t_half) - 2^(-(k + 1) * period / t_half))
}

#' Exposure structure of one drug sequence
#'
#' @param d Integer vector of drug indicators, one per period.
#' @return Data frame with k (periods since last drug, NA before first
#'   exposure and on drug), first_off (first off-drug period after a
#'   discontinuation) and exposed (any earlier drug period).
exposure_structure <- function(d) {
  n <- length(d)
  k <- rep(NA_integer_, n)
  exposed <- rep(FALSE, n)
  run <- NA_integer_
  seen <- FALSE
  for (j in seq_len(n)) {
    exposed[j] <- seen
    if (d[j] == 1L) {
      run <- NA_integer_
      seen <- TRUE
    } else if (seen) {
      run <- if (is.na(run)) 0L else run + 1L
      k[j] <- run
    }
  }
  data.frame(k = k, exposed = exposed, first_off = !is.na(k) & k == 0L)
}

#' Regressors and retained rows for each N-of-1 analysis rule
#'
#' ALL: binary D, every period (no carryover handling).
#' XS: binary D, first off-drug period after each discontinuation
#'   dropped by schedule (no information about the truth used).
#' DBC: Dbc built at the analysis half-life t_half_a, every period.
#' DBCL: DBC plus an explicit lag regressor L1 (first off-drug period).
#' XO: ORACLE. Binary D; drops off-drug periods whose true carryover
#'   fraction exceeds `oracle_cut`. Uses the true half-life and is not
#'   available to a trial.
#'
#' @param d Drug indicators for one participant.
#' @param t_half True half-life (weeks), used only by XO.
#' @param t_half_a Analysis half-life (weeks).
analysis_design <- function(d, t_half, t_half_a = PILOT$t_half_a,
                            cut = PILOT$oracle_cut) {
  es <- exposure_structure(d)
  off_exposed <- d == 0L & es$exposed
  f_a <- ifelse(off_exposed, carry_fraction(es$k, t_half_a), 0)
  f_true <- ifelse(off_exposed, carry_fraction(es$k, t_half), 0)
  l1 <- as.numeric(es$first_off)
  list(
    ALL = list(x = d, l1 = NULL, keep = rep(TRUE, length(d))),
    XS = list(x = d, l1 = NULL, keep = !es$first_off),
    DBC = list(x = d + f_a, l1 = NULL, keep = rep(TRUE, length(d))),
    DBCL = list(x = d + f_a, l1 = l1, keep = rep(TRUE, length(d))),
    XO = list(x = d, l1 = NULL, keep = !(off_exposed & f_true > cut))
  )
}

ANALYSES <- c('ALL', 'XS', 'DBC', 'DBCL', 'XO')

#' Long-format regressors for every rule over a participant panel
#'
#' @param dat Data frame from `generate_nof1()` over all participants.
#' @param t_half True half-life, used only by the oracle rule XO.
rules_long <- function(dat, t_half) {
  parts <- split(dat, dat$ptID)
  per <- lapply(parts, function(p) analysis_design(p$D, t_half))
  out <- lapply(ANALYSES, function(a) {
    list(
      x_long = unlist(lapply(per, function(z) z[[a]]$x), use.names = FALSE),
      l1_long = if (is.null(per[[1L]][[a]]$l1)) NULL else
        unlist(lapply(per, function(z) z[[a]]$l1), use.names = FALSE),
      keep_long = unlist(lapply(per, function(z) z[[a]]$keep),
                         use.names = FALSE)
    )
  })
  names(out) <- ANALYSES
  out
}

#' Within-participant covariance of the DGP
#'
#' Random intercept plus continuous-time AR(1) residual with lag-one
#' correlation rho per week.
v_matrix <- function(t_wk, sigma_u = PILOT$sigma_u,
                     sigma_w = PILOT$sigma_w, rho = PILOT$rho) {
  sigma_u^2 + sigma_w^2 * outer(t_wk, t_wk, function(a, b) rho^abs(a - b))
}

#' Generate one N-of-1 participant
#'
#' @param pid Participant id.
#' @param u Participant intercept deviation.
#' @param beta_d True on-drug effect (nightmares per week).
#' @param t_half True carryover half-life (weeks).
#' @param chol_w Upper Cholesky factor of the residual covariance.
generate_nof1 <- function(pid, u, beta_d, t_half, chol_w,
                          n_periods = PILOT$n_periods,
                          mu0 = PILOT$mu0) {
  d <- sample(rep(c(1L, 0L), each = n_periods %/% 2L))
  es <- exposure_structure(d)
  carry <- ifelse(d == 0L & es$exposed,
                  beta_d * carry_fraction(es$k, t_half), 0)
  eps <- as.numeric(crossprod(chol_w, stats::rnorm(n_periods)))
  data.frame(
    ptID = pid, t_wk = seq_len(n_periods), D = d,
    Sx = mu0 + u + beta_d * d + carry + eps
  )
}

#' Fit one N-of-1 analysis rule with two inference methods
#'
#' Returns the coefficient of the drug regressor with (a) the `nlme`
#' containment (between-within) degrees of freedom and (b) CR2
#' cluster-robust standard errors with Satterthwaite degrees of
#' freedom from clubSandwich.
fit_nof1 <- function(dat, rule, true_effect) {
  na_out <- list(
    beta = NA_real_, se_c = NA_real_, df_c = NA_real_, p_c = NA_real_,
    lo_c = NA_real_, hi_c = NA_real_, se_r = NA_real_, df_r = NA_real_,
    p_r = NA_real_, lo_r = NA_real_, hi_r = NA_real_, n_obs = NA_integer_,
    ok = FALSE
  )
  rhs <- if (is.null(rule$l1)) 'Dx + t_wk' else 'Dx + L1 + t_wk'
  d <- dat
  d$Dx <- rule$x_long
  if (!is.null(rule$l1)) d$L1 <- rule$l1_long
  d <- d[rule$keep_long, ]
  fit <- tryCatch(
    nlme::lme(
      stats::as.formula(paste('Sx ~', rhs)),
      random = ~ 1 | ptID,
      correlation = nlme::corCAR1(form = ~ t_wk | ptID),
      data = d, method = 'REML',
      control = nlme::lmeControl(opt = 'optim', maxIter = 200L,
                                 msMaxIter = 200L)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(na_out)
  st <- summary(fit)$tTable
  beta <- st['Dx', 'Value']
  se_c <- st['Dx', 'Std.Error']
  df_c <- st['Dx', 'DF']
  hw <- stats::qt(0.975, df_c) * se_c
  cr <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', cluster = d$ptID,
                            test = 'Satterthwaite'),
    error = function(e) NULL
  )
  se_r <- df_r <- p_r <- NA_real_
  if (!is.null(cr)) {
    row <- cr[rownames(cr) == 'Dx' | cr$Coef == 'Dx', , drop = FALSE]
    if (nrow(row) == 1L) {
      se_r <- row$SE
      df_r <- row$df_Satt
      p_r <- row$p_Satt
    }
  }
  hw_r <- stats::qt(0.975, df_r) * se_r
  list(
    beta = beta, se_c = se_c, df_c = df_c,
    p_c = 2 * stats::pt(abs(beta / se_c), df_c, lower.tail = FALSE),
    lo_c = beta - hw, hi_c = beta + hw,
    se_r = se_r, df_r = df_r, p_r = p_r,
    lo_r = beta - hw_r, hi_r = beta + hw_r,
    n_obs = nrow(d), ok = TRUE
  )
}

#' Simulated covariate with reliability r
#'
#' @param u True participant intercept deviations.
#' @param r Reliability, Var(u) / Var(observed covariate).
noisy_baseline <- function(u, r, sigma_u = PILOT$sigma_u) {
  if (r >= 1) return(u)
  u + stats::rnorm(length(u), 0, sigma_u * sqrt((1 - r) / r))
}

#' Parallel-group RCT arms on the shared participant pool
#'
#' RCT1: one post-baseline measurement, ANCOVA on the baseline
#'   covariate (reliability `r_bl`; `r_bl = 1` is the error-free
#'   covariate of the 2026-05 run).
#' RCT8: eight post-baseline weekly visits (matched observation
#'   burden), ANCOVA on the participant mean of the eight visits with
#'   the same baseline covariate. For balanced complete data this
#'   equals the GLS estimator of the mixed model.
generate_rct <- function(u, beta_d, chol_w, r_bl, n_visits) {
  n <- length(u)
  trt <- sample(c(rep(1L, n %/% 2L), rep(0L, n - n %/% 2L)))
  x <- noisy_baseline(u, r_bl)
  cw <- chol_w[seq_len(n_visits), seq_len(n_visits), drop = FALSE]
  eps <- t(crossprod(cw, matrix(stats::rnorm(n_visits * n),
                                    n_visits, n)))
  y <- PILOT$mu0 + u + beta_d * trt + rowMeans(eps)
  data.frame(y = y, trt = trt, x = x)
}

fit_rct <- function(dat) {
  fit <- stats::lm(y ~ trt + x, data = dat)
  cf <- summary(fit)$coefficients
  ci <- stats::confint(fit, 'trt')
  list(beta = cf['trt', 'Estimate'], se = cf['trt', 'Std. Error'],
       p = cf['trt', 'Pr(>|t|)'], lo = ci[1L], hi = ci[2L])
}

# ---- Closed form -------------------------------------------------

#' All balanced drug sequences (70 for 4 of 8)
all_sequences <- function(n_periods = PILOT$n_periods) {
  idx <- utils::combn(n_periods, n_periods %/% 2L)
  lapply(seq_len(ncol(idx)), function(j) {
    d <- integer(n_periods)
    d[idx[, j]] <- 1L
    d
  })
}

#' Closed-form GLS performance of one N-of-1 analysis rule
#'
#' Averages the GLS information and the bias numerator over the
#' equally likely balanced sequences. The covariance model is the
#' DGP's (random intercept plus AR(1)), so the model-based standard
#' error is the correct one up to the sampling error of the estimated
#' variance components.
#'
#' @param rule One of ANALYSES.
#' @param beta_d True effect.
#' @param t_half True carryover half-life (weeks).
#' @return Named list: se, bias, mean_obs, power, coverage, df.
cf_nof1 <- function(rule, beta_d, t_half, n = PILOT$n_patients,
                    t_half_a = PILOT$t_half_a, alpha = PILOT$alpha) {
  seqs <- all_sequences()
  tw <- seq_len(PILOT$n_periods)
  info <- 0
  num <- 0
  nobs <- 0
  w <- 1 / length(seqs)
  for (d in seqs) {
    ad <- analysis_design(d, t_half, t_half_a)[[rule]]
    keep <- ad$keep
    es <- exposure_structure(d)
    cols <- list(rep(1, length(d)), ad$x)
    if (!is.null(ad$l1)) cols <- c(cols, list(ad$l1))
    cols <- c(cols, list(tw))
    x <- do.call(cbind, cols)[keep, , drop = FALSE]
    vi <- solve(v_matrix(tw[keep]))
    off_exposed <- d == 0L & es$exposed
    mu_carry <- ifelse(off_exposed,
                       beta_d * carry_fraction(es$k, t_half), 0)
    mu <- (beta_d * d + mu_carry)[keep]
    info <- info + w * crossprod(x, vi %*% x)
    num <- num + w * as.numeric(crossprod(x, vi %*% mu))
    nobs <- nobs + w * sum(keep)
  }
  ev <- solve(n * info)
  se <- sqrt(ev[2L, 2L])
  coef_mean <- as.numeric(ev %*% (n * num))
  bias <- coef_mean[2L] - beta_d
  df <- n * nobs - n - ncol(info) + 1
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- coef_mean[2L] / se
  power <- stats::pt(crit, df, ncp, lower.tail = FALSE) +
    stats::pt(-crit, df, ncp)
  ncp_cov <- bias / se
  cover <- stats::pt(crit, df, ncp_cov) - stats::pt(-crit, df, ncp_cov)
  list(se = se, bias = bias, mean_obs = nobs, df = df,
       power = power, coverage = cover)
}

#' Closed-form performance of the RCT comparators
#'
#' @param r_bl Baseline reliability (1 = error-free covariate).
#' @param n_visits Post-baseline weekly visits (1 or 8).
cf_rct <- function(beta_d, r_bl, n_visits, n = PILOT$n_patients,
                   alpha = PILOT$alpha) {
  tw <- seq_len(n_visits)
  r_mat <- outer(tw, tw, function(a, b) PILOT$rho^abs(a - b))
  s <- sum(r_mat) / n_visits^2
  resid_var <- PILOT$sigma_w^2 * s + PILOT$sigma_u^2 * (1 - r_bl)
  n1 <- n %/% 2L
  n0 <- n - n1
  se <- sqrt(resid_var * (1 / n1 + 1 / n0))
  df <- n - 3
  crit <- stats::qt(1 - alpha / 2, df)
  ncp <- beta_d / se
  power <- stats::pt(crit, df, ncp, lower.tail = FALSE) +
    stats::pt(-crit, df, ncp)
  list(se = se, bias = 0, df = df, power = power,
       coverage = 1 - alpha)
}

#' Closed-form table over a grid of effects and half-lives
cf_table <- function(effects, t_halves) {
  grid <- expand.grid(beta_d = effects, t_half = t_halves)
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    b <- grid$beta_d[i]
    h <- grid$t_half[i]
    n1 <- do.call(rbind, lapply(ANALYSES, function(a) {
      r <- cf_nof1(a, b, h)
      data.frame(beta_d = b, t_half = h, arm = a, se = r$se,
                 bias = r$bias, mean_obs = r$mean_obs,
                 power = r$power, coverage = r$coverage)
    }))
    rc <- rbind(
      cbind(arm = 'RCT1_exact', as.data.frame(cf_rct(b, 1, 1))),
      cbind(arm = 'RCT1_r08', as.data.frame(cf_rct(b, PILOT$r_bl, 1))),
      cbind(arm = 'RCT8_r08', as.data.frame(cf_rct(b, PILOT$r_bl, 8)))
    )
    rc$beta_d <- b
    rc$t_half <- h
    rc$mean_obs <- c(1, 1, 8)
    rbind(n1, rc[, names(n1)])
  }))
  rownames(out) <- NULL
  out
}

# ---- Simulation: one replicate and its summaries ------------------

ARMS <- c(ANALYSES, 'RCT1_exact', 'RCT1_r08', 'RCT8_r08')

#' One replicate: shared participant pool, all arms
#'
#' @param stream A `.Random.seed` vector for this replicate (L'Ecuyer).
run_replicate <- function(rep_id, stream, true_effect, t_half,
                          n = PILOT$n_patients) {
  assign('.Random.seed', stream, envir = globalenv())
  u <- stats::rnorm(n, 0, PILOT$sigma_u)
  tw <- seq_len(PILOT$n_periods)
  chol_w <- chol(PILOT$sigma_w^2 *
    outer(tw, tw, function(a, b) PILOT$rho^abs(a - b)))
  dat <- do.call(rbind, lapply(seq_len(n), function(i) {
    generate_nof1(i, u[i], true_effect, t_half, chol_w)
  }))
  rl <- rules_long(dat, t_half)
  rows <- lapply(ANALYSES, function(a) {
    f <- fit_nof1(dat, rl[[a]], true_effect)
    data.frame(rep = rep_id, arm = a, f, stringsAsFactors = FALSE)
  })
  rct_arm <- function(arm, r_bl, v) {
    f <- fit_rct(generate_rct(u, true_effect, chol_w, r_bl, v))
    data.frame(
      rep = rep_id, arm = arm, beta = f$beta, se_c = f$se,
      df_c = n - 3, p_c = f$p, lo_c = f$lo, hi_c = f$hi,
      se_r = NA_real_, df_r = NA_real_, p_r = NA_real_,
      lo_r = NA_real_, hi_r = NA_real_, n_obs = v * n, ok = TRUE
    )
  }
  rows <- c(rows, list(
    rct_arm('RCT1_exact', 1, 1),
    rct_arm('RCT1_r08', PILOT$r_bl, 1),
    rct_arm('RCT8_r08', PILOT$r_bl, 8)
  ))
  do.call(rbind, rows)
}

#' Independent L'Ecuyer streams, one per replicate, for all cells
#'
#' `parallel::mclapply(mc.set.seed = TRUE)` does not advance the
#' master stream between calls, so two calls repeat the same draws.
#' This function advances one chain with `nextRNGStream()` so that
#' every (cell, replicate) pair has its own stream, independent of
#' the number of workers.
make_streams <- function(n_total, seed) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  s <- .Random.seed
  out <- vector('list', n_total)
  for (i in seq_len(n_total)) {
    out[[i]] <- s
    s <- parallel::nextRNGStream(s)
  }
  out
}

#' Cell summary from replicate-level rows
#'
#' @param reps Data frame of replicate rows for one cell.
#' @param true_effect True effect in the cell.
summarize_cell <- function(reps, true_effect, alpha = PILOT$alpha) {
  do.call(rbind, lapply(split(reps, reps$arm), function(d) {
    d <- d[d$ok, ]
    n <- nrow(d)
    pw <- function(p) mean(p < alpha, na.rm = TRUE)
    mc <- function(p) sqrt(p * (1 - p) / n)
    cov_c <- mean(d$lo_c <= true_effect & d$hi_c >= true_effect)
    cov_r <- mean(d$lo_r <= true_effect & d$hi_r >= true_effect,
                  na.rm = TRUE)
    data.frame(
      arm = d$arm[1L], n_ok = n,
      power_c = pw(d$p_c), power_c_mcse = mc(pw(d$p_c)),
      power_r = pw(d$p_r), power_r_mcse = mc(pw(d$p_r)),
      bias = mean(d$beta) - true_effect,
      bias_mcse = stats::sd(d$beta) / sqrt(n),
      emp_se = stats::sd(d$beta), mean_se_c = mean(d$se_c),
      mean_se_r = mean(d$se_r, na.rm = TRUE),
      cov_c = cov_c, cov_c_mcse = mc(cov_c),
      cov_r = cov_r, cov_r_mcse = mc(cov_r),
      mean_obs = mean(d$n_obs) / PILOT$n_patients
    )
  }))
}

#' Paired power difference with MCSE from replicate-level indicators
paired_diff <- function(reps, arm_a, arm_b, col = 'p_c',
                        alpha = PILOT$alpha) {
  a <- reps[reps$arm == arm_a, c('rep', col)]
  b <- reps[reps$arm == arm_b, c('rep', col)]
  m <- merge(a, b, by = 'rep', suffixes = c('_a', '_b'))
  ia <- m[[paste0(col, '_a')]] < alpha
  ib <- m[[paste0(col, '_b')]] < alpha
  dd <- as.numeric(ia) - as.numeric(ib)
  c(diff = mean(dd), mcse = stats::sd(dd) / sqrt(length(dd)),
    mcse_independent = sqrt(mean(ia) * (1 - mean(ia)) / length(ia) +
      mean(ib) * (1 - mean(ib)) / length(ib)))
}
