## 20-paired-dropout.R
##
## Paired (common-random-numbers) rerun for paper 09,
## analysis/report/09-informative-dropout-by-design/report.Rmd.
## Written 2026-10-29 in response to the 2026-10-28 referee report
## (M1, M2, M4 to M7, M12, M13).
##
## What it does. For each design (OL+BDC, CO, Hybrid) and replicate,
## one complete dataset is generated at a total N = 70 under the
## mean-shift encoding (beta_bm = 0.45, carryover half-life 0.5 week).
## Every dropout condition is then applied to that same dataset with
## the same uniform draws, so the conditions differ only in the
## hazard. Losses and shifts are paired differences.
##
## The OL design is not run. With every visit on drug it has no
## off-drug contrast, lme_analysis() therefore tests bm:t instead of
## bm:Dbc, and the OL row of the earlier run did not estimate the
## interaction (referee M4).
##
## Dropout conditions (per-visit hazard, monotone; alpha solved on
## each realized path dataset so that the expected fraction of
## participants who ever drop equals the target; see report Sec. 2):
##   none            no dropout
##   MCAR_25         constant hazard
##   magnitude_25/40 beta0 + beta1 * z^2, z the change from baseline
##                   divided by its cross-sectional SD at that visit
##                   (not centered): the prototype's 'biased' hazard,
##                   relabeled because it rises with the size of the
##                   change in either direction
##   worsening_25    exp(gamma * zc), zc the centered and scaled
##                   change from baseline (positive = worse)
##   improvement_25  exp(-gamma * zc)
##   return_25       exp(gamma * zr), zr the centered and scaled
##                   change since the previous visit (symptom return)
## gamma = log(3): the hazard triples per SD of the driver.
##
## Analyses (each fitted to every post-dropout dataset):
##   ta0.5  historical model, Dbc at t_a = 0.5 week (the DGP half-life
##          of the mean drug response; reproduces lme_analysis())
##   ta0    the same model with t_a = 0, i.e. binary drug state, which
##          matches the generating interaction: the mean shift is
##          applied at on-drug visits only and has no carryover
## Each reports the model-based test (lme t, nlme df) and the CR2
## cluster-robust test with Satterthwaite df (clubSandwich), as in
## analysis/scripts/carryover-sensitivity/simulation-core.R.
## De enters for OL+BDC and Hybrid, as in lme_analysis() (CO has a
## constant post-baseline expectancy).
##
## Regimes: 'alt' (beta_bm = 0.45) and, when P09_NULL_REPS > 0,
## 'null' (beta_bm = 0; same seeds, so identical draws before the
## shift).
##
## Seeds: set.seed(P09_SEED + 100000 * design_index + rep) per
## (design, replicate); the null regime adds 50000.
##
## Environment variables:
##   P09_REPS       replicates per design, alt regime (default 1000)
##   P09_NULL_REPS  replicates per design, null regime (default 1000)
##   P09_CORES      worker count (default 8)
##   P09_OUT        output directory
##                  (default analysis/data/informative-dropout-by-design)
##   P09_TAG        file tag (default 'paired')
##
## Run from the repository root through the wrapper
## analysis/scripts/informative-dropout-by-design/20-run-paired.sh,
## which pins BLAS threads to 1.
##
## Outputs:
##   <out>/<tag>-replicates.rds   one row per design x regime x rep x
##                                condition x analysis
##   <out>/<tag>-dropout.rds      dropout diagnostics per design x
##                                regime x rep x condition
##   <out>/<tag>-meta.rds         settings, timing, session info

suppressPackageStartupMessages({
  pkgload::load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
  loadNamespace('clubSandwich')
})

P09_SEED <- 20261029L
n_reps <- as.integer(Sys.getenv('P09_REPS', '1000'))
n_null <- as.integer(Sys.getenv('P09_NULL_REPS', '1000'))
n_cores <- as.integer(Sys.getenv('P09_CORES', '8'))
out_dir <- Sys.getenv('P09_OUT', 'analysis/data/informative-dropout-by-design')
tag <- Sys.getenv('P09_TAG', 'paired')

N_TOTAL <- 70L
BETA_BM <- 0.45
T_HALF <- 0.5
SIGMA_BR <- 5
TRUE_BETA <- -BETA_BM * SIGMA_BR / 1
GAMMA <- log(3)
ANALYSIS_HL <- c(ta0.5 = 0.5, ta0 = 0)

td_OLBDC <- buildtrialdesign(
  name_longform = 'open label + blinded discontinuation',
  name_shortform = 'OLBDC',
  timepoints = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames = c('OL1', 'OL2', 'OL3', 'OL4',
                  'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

td_CO <- buildtrialdesign(
  name_longform = 'traditional crossover',
  name_shortform = 'CO',
  timepoints = cumulative(rep(2.5, 8)),
  timeptnames = c(paste0('COa', 1:4), paste0('COb', 1:4)),
  expectancies = rep(0.5, 8),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 0, 0),
    pathB = c(0, 0, 0, 0, 1, 1, 1, 1)
  )
)

td_Hybrid <- buildtrialdesign(
  name_longform = 'hybrid N-of-1',
  name_shortform = 'Hybrid',
  timepoints = c(4, 8, 9, 10, 11, 12, 16, 20),
  timeptnames = c('OL1', 'OL2', 'BD1', 'BD2',
                  'BD3', 'BD4', 'COd', 'COp'),
  expectancies = c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 0, 0, 1, 0),
    pathB = c(1, 1, 1, 1, 0, 0, 0, 1),
    pathC = c(1, 1, 1, 0, 0, 0, 1, 0),
    pathD = c(1, 1, 1, 0, 0, 0, 0, 1)
  )
)

designs <- list(OLBDC = td_OLBDC, CO = td_CO, Hybrid = td_Hybrid)

bp <- data.table(cat = c('bm', 'BL'), m = c(0, 70), sd = c(1, 10))
rp <- data.table(
  cat = c('tv', 'pb', 'br'),
  max = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd = c(5, 2, SIGMA_BR)
)

make_mp <- function(beta_bm) {
  data.table(N = NA_integer_, c.bm = beta_bm,
             carryover_t1half = T_HALF,
             c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
             c.cf1t = 0.1, c.cfct = 0.05)
}

allocate_across_paths <- function(n_total, n_paths) {
  base <- n_total %/% n_paths
  rep(base, n_paths) +
    c(rep(1L, n_total %% n_paths),
      rep(0L, n_paths - n_total %% n_paths))
}

## The mean-shift encoding leaves Sigma independent of beta_bm, so one
## factorization per path serves both regimes.
sigma_cache <- lapply(designs, function(td) {
  lapply(td$trialpaths, function(path) {
    buildSigma(make_mp(BETA_BM), rp, bp, path,
               makePositiveDefinite = FALSE,
               dgp_architecture = 'mean_moderation',
               stopOnNonPD = TRUE)
  })
})

CONDITIONS <- list(
  none           = list(mech = 'none',        rate = 0),
  MCAR_25        = list(mech = 'MCAR',        rate = 0.25),
  magnitude_25   = list(mech = 'magnitude',   rate = 0.25),
  magnitude_40   = list(mech = 'magnitude',   rate = 0.40),
  worsening_25   = list(mech = 'worsening',   rate = 0.25),
  improvement_25 = list(mech = 'improvement', rate = 0.25),
  return_25      = list(mech = 'return',      rate = 0.25)
)

scale_cols <- function(m, center) {
  apply(m, 2, function(x) {
    s <- sd(x)
    if (!is.finite(s) || s < 1e-6) s <- 1
    (x - if (center) mean(x) else 0) / s
  })
}

hazard_raw <- function(y, bl, mech, beta0 = 0.05, beta1 = 0.5) {
  n <- nrow(y)
  v <- ncol(y)
  if (mech == 'MCAR') return(matrix(1, n, v))
  delta <- y - bl
  h <- switch(
    mech,
    magnitude = beta0 + beta1 * scale_cols(delta, FALSE)^2,
    worsening = exp(GAMMA * scale_cols(delta, TRUE)),
    improvement = exp(-GAMMA * scale_cols(delta, TRUE)),
    return = exp(GAMMA * scale_cols(y - cbind(bl, y[, -v, drop = FALSE]),
                                    TRUE)),
    stop('unknown mechanism ', mech)
  )
  h[!is.finite(h)] <- 1
  if (mech == 'magnitude') h <- pmin(pmax(h, 0), 1)
  h
}

## alpha solves mean(1 - prod_t (1 - alpha * h_it)) = target on the
## realized path dataset, as in the prototype driver.
calibrate_alpha <- function(h, target) {
  obj <- function(a) {
    mean(1 - exp(rowSums(log(pmax(1 - pmin(a * h, 1), 1e-9))))) -
      target
  }
  if (obj(1e-6) > 0) return(1e-6)
  if (obj(1) < 0) return(1)
  uniroot(obj, c(1e-6, 1), tol = 1e-6)$root
}

## Returns the first dropout visit per participant (NA if none).
first_dropout <- function(y, bl, cond, u) {
  if (cond$mech == 'none') return(rep(NA_integer_, nrow(y)))
  h <- hazard_raw(y, bl, cond$mech)
  a <- calibrate_alpha(h, cond$rate)
  ev <- u < pmin(a * h, 1)
  fd <- max.col(ev, ties.method = 'first')
  fd[rowSums(ev) == 0] <- NA_integer_
  fd
}

path_long_template <- function(path) {
  td <- as.data.table(path)[, .(timeptnames, t_wk, e, tod, tsd)]
  td <- rbind(data.table(timeptnames = 'BL', t_wk = 0, e = 0,
                         tod = 0, tsd = 0), td)
  td[, t := cumsum(t_wk)]
  td
}

build_long <- function(y_list, bl_list, bm_list, fd_list, templates) {
  out <- vector('list', length(y_list))
  off <- 0L
  for (g in seq_along(y_list)) {
    tp <- templates[[g]]
    y <- cbind(bl_list[[g]], y_list[[g]])
    fd <- fd_list[[g]]
    keep <- matrix(TRUE, nrow(y), ncol(y))
    for (i in which(!is.na(fd))) keep[i, (fd[i] + 1L):ncol(y)] <- FALSE
    n <- nrow(y)
    ids <- off + seq_len(n)
    off <- off + n
    l <- data.table(
      ptID = rep(ids, times = ncol(y)),
      bm = rep(bm_list[[g]], times = ncol(y)),
      Sx = as.vector(y),
      t = rep(tp$t, each = n),
      De = rep(tp$e, each = n),
      Db = rep(tp$tod > 0, each = n),
      tsd = rep(tp$tsd, each = n)
    )
    out[[g]] <- l[as.vector(keep)]
  }
  rbindlist(out)
}

fit_analysis <- function(long, t_a, use_de) {
  long <- copy(long)
  long[, Dbc := fifelse(Db, 1,
                        fifelse(tsd > 0 & t_a > 0, 0.5^(tsd / t_a), 0))]
  form <- if (use_de) {
    Sx ~ bm + t + Dbc + bm:Dbc + De
  } else {
    Sx ~ bm + t + Dbc + bm:Dbc
  }
  car1 <- TRUE
  fit <- tryCatch(
    nlme::lme(form, random = ~1 | ptID,
              correlation = nlme::corCAR1(form = ~t | ptID),
              data = long,
              control = nlme::lmeControl(opt = 'optim', maxIter = 200,
                                         msMaxIter = 200)),
    error = function(e) NULL)
  if (is.null(fit)) {
    car1 <- FALSE
    fit <- tryCatch(
      nlme::lme(form, random = ~1 | ptID, data = long,
                control = nlme::lmeControl(opt = 'optim')),
      error = function(e) NULL)
  }
  na_row <- list(est = NA_real_, se_mb = NA_real_, df_mb = NA_real_,
                 p_mb = NA_real_, se_cr2 = NA_real_, df_cr2 = NA_real_,
                 p_cr2 = NA_real_, car1 = car1, fit_ok = FALSE)
  if (is.null(fit)) return(na_row)
  tt <- summary(fit)$tTable
  tgt <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(tt))[1]
  if (is.na(tgt)) return(na_row)
  ct <- tryCatch(
    clubSandwich::coef_test(fit, vcov = 'CR2', test = 'Satterthwaite'),
    error = function(e) NULL)
  cr <- if (is.null(ct)) NULL else ct[ct$Coef == tgt, , drop = FALSE]
  list(est = tt[tgt, 'Value'], se_mb = tt[tgt, 'Std.Error'],
       df_mb = tt[tgt, 'DF'], p_mb = tt[tgt, 'p-value'],
       se_cr2 = if (is.null(cr) || !nrow(cr)) NA_real_ else cr$SE,
       df_cr2 = if (is.null(cr) || !nrow(cr)) NA_real_ else cr$df_Satt,
       p_cr2 = if (is.null(cr) || !nrow(cr)) NA_real_ else cr$p_Satt,
       car1 = car1, fit_ok = TRUE)
}

run_unit <- function(design, regime, rep) {
  d_idx <- match(design, names(designs))
  seed <- P09_SEED + 100000L * d_idx + rep +
    if (regime == 'null') 50000L else 0L
  set.seed(seed)
  td <- designs[[design]]
  paths <- td$trialpaths
  n_pp <- allocate_across_paths(N_TOTAL, length(paths))
  mp <- make_mp(if (regime == 'alt') BETA_BM else 0)
  templates <- lapply(paths, path_long_template)
  use_de <- length(unique(paths[[length(paths)]]$e)) > 1

  y_list <- bl_list <- bm_list <- u_list <- ond_list <- impr_list <-
    vector('list', length(paths))
  for (g in seq_along(paths)) {
    mp$N <- n_pp[g]
    di <- generateData(mp, rp, bp, paths[[g]], empirical = FALSE,
                       makePositiveDefinite = FALSE,
                       cached_sigma = sigma_cache[[design]][[g]],
                       dgp_architecture = 'mean_moderation',
                       stopOnNonPD = TRUE)
    tn <- paths[[g]]$timeptnames
    y_list[[g]] <- as.matrix(di[, ..tn])
    bl_list[[g]] <- di$BL
    bm_list[[g]] <- di$bm
    u_list[[g]] <- matrix(runif(length(y_list[[g]])),
                          nrow = n_pp[g])
    ond_list[[g]] <- paths[[g]]$tod > 0
    impr_list[[g]] <- rowMeans(di$BL - y_list[[g]])
  }

  res <- vector('list', length(CONDITIONS))
  diag <- vector('list', length(CONDITIONS))
  for (k in seq_along(CONDITIONS)) {
    cond <- CONDITIONS[[k]]
    fd_list <- lapply(seq_along(paths), function(g) {
      first_dropout(y_list[[g]], bl_list[[g]], cond, u_list[[g]])
    })
    long <- build_long(y_list, bl_list, bm_list, fd_list, templates)

    fd_all <- unlist(fd_list)
    dropped <- !is.na(fd_all)
    off_at_drop <- unlist(lapply(seq_along(paths), function(g) {
      f <- fd_list[[g]]
      f <- f[!is.na(f)]
      !ond_list[[g]][f]
    }))
    n_post <- sum(vapply(y_list, length, integer(1)))
    n_lost <- sum(vapply(seq_along(paths), function(g) {
      f <- fd_list[[g]]
      v <- ncol(y_list[[g]])
      sum(ifelse(is.na(f), 0L, v - f + 1L))
    }, numeric(1)))
    impr <- unlist(impr_list)
    path_frac <- vapply(fd_list, function(f) mean(!is.na(f)), numeric(1))
    diag[[k]] <- data.table(
      design = design, regime = regime, rep = rep,
      condition = names(CONDITIONS)[k],
      ever_drop = mean(dropped),
      obs_missing = n_lost / n_post,
      frac_drop_offdrug = if (any(dropped)) mean(off_at_drop) else
        NA_real_,
      cor_drop_improve = if (any(dropped) && !all(dropped))
        cor(as.numeric(dropped), impr) else NA_real_,
      cor_drop_bm = if (any(dropped) && !all(dropped))
        cor(as.numeric(dropped), unlist(bm_list)) else NA_real_,
      drop_p1 = path_frac[1], drop_p2 = path_frac[2],
      drop_p3 = if (length(path_frac) >= 3) path_frac[3] else NA_real_,
      drop_p4 = if (length(path_frac) >= 4) path_frac[4] else NA_real_
    )

    res[[k]] <- rbindlist(lapply(names(ANALYSIS_HL), function(an) {
      r <- fit_analysis(long, ANALYSIS_HL[[an]], use_de)
      data.table(design = design, regime = regime, rep = rep,
                 condition = names(CONDITIONS)[k], analysis = an,
                 as.data.table(r))
    }))
  }
  list(res = rbindlist(res), diag = rbindlist(diag))
}

units <- rbind(
  CJ(design = names(designs), regime = 'alt', rep = seq_len(n_reps)),
  if (n_null > 0) {
    CJ(design = names(designs), regime = 'null', rep = seq_len(n_null))
  }
)

cat(sprintf('Units: %d (alt %d, null %d per design); cores %d\n',
            nrow(units), n_reps, n_null, n_cores))
t0 <- Sys.time()
out <- mclapply(seq_len(nrow(units)), function(j) {
  tryCatch(run_unit(units$design[j], units$regime[j], units$rep[j]),
           error = function(e) {
             list(res = data.table(design = units$design[j],
                                   regime = units$regime[j],
                                   rep = units$rep[j],
                                   error = conditionMessage(e)),
                  diag = NULL)
           })
}, mc.cores = n_cores, mc.preschedule = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = 'secs'))

reps_dt <- rbindlist(lapply(out, `[[`, 'res'), fill = TRUE)
diag_dt <- rbindlist(lapply(out, `[[`, 'diag'), fill = TRUE)
reps_dt[, true_beta := fifelse(regime == 'alt', TRUE_BETA, 0)]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(reps_dt, file.path(out_dir, paste0(tag, '-replicates.rds')))
saveRDS(diag_dt, file.path(out_dir, paste0(tag, '-dropout.rds')))
saveRDS(list(seed = P09_SEED, n_reps = n_reps, n_null = n_null,
             n_cores = n_cores, N_total = N_TOTAL, beta_bm = BETA_BM,
             t_half = T_HALF, true_beta = TRUE_BETA, gamma = GAMMA,
             analysis_hl = ANALYSIS_HL, conditions = CONDITIONS,
             elapsed_secs = elapsed, finished = Sys.time(),
             session = sessionInfo()),
        file.path(out_dir, paste0(tag, '-meta.rds')))

cat(sprintf('Done: %d rows, %.1f s\n', nrow(reps_dt), elapsed))
if ('error' %in% names(reps_dt)) {
  print(reps_dt[!is.na(error), .N, by = .(design, regime, error)])
}
