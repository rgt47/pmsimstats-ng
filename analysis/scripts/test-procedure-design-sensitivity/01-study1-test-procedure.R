## 01-study1-test-procedure.R
##
## Study 1 (revised) for paper 08-test-procedure-design-sensitivity.
##
## Supersedes analysis/scripts/quick-sim/08-test-procedure-prototype.R.
## Changes addressing the 2026-06-17 referee report:
##   M2  add a continuous-biomarker within-subject contrast arm so the
##       dichotomisation cost is separated from the test-procedure cost.
##   M4  restore the sample-size grid to N in {25, 35, 70, 100}.
##   M6  compute and report bias, empirical SE, mean model SE, SE ratio,
##       and 95% CI coverage for every coefficient-bearing arm.
##   M7  run under parallel::mclapply with L'Ecuyer-CMRG streams so the
##       reproducibility statement (parallel, per-stream seeding) is true.
##   Suggestion 4: benchmark the inline Mancl-DeRouen SE against
##       geesmv::GEE.var.md on a probe cell.
##
## Five procedures per replicate:
##   1. rmanova       strict split-plot F, median-dichotomised bm,
##                    within-phase pre-averaged (dichotomised estimand).
##   2. contrast_cont within-subject (on - off) contrast regressed on the
##                    CONTINUOUS bm; isolates the test from dichotomisation.
##   3. lme           nlme::lme with corCAR1 (lme_analysis()).
##   4. gee           geepack AR(1), naive sandwich.
##   5. gee_mancl     GEE AR(1) with inline Mancl-DeRouen (2001) sandwich.
##
## Arms 2-5 share the interaction estimand on the Sx scale
## (true bm:Dbc ~= c.bm * sd_br); arm 1 targets a dichotomised estimand
## and is excluded from bias/coverage.
##
## Outputs (overwrite):
##   analysis/data/quick-sim/08-test-replicates.rds
##   analysis/data/quick-sim/08-test-summary.txt
##   analysis/data/quick-sim/08-test-betatrue.txt
##   analysis/data/quick-sim/08-test-mancl-benchmark.txt
##   analysis/figures/quick-sim/08-test-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(nlme)
  library(parallel)
})

have_geepack <- requireNamespace('geepack', quietly = TRUE)
if (have_geepack) suppressPackageStartupMessages(library(geepack))
have_geesmv <- requireNamespace('geesmv', quietly = TRUE)
cat('geepack available:', have_geepack, '\n')
cat('geesmv available:', have_geesmv, '\n')

## --------------------------------------------------------------
## Run configuration (env-overridable for probe vs production).
## --------------------------------------------------------------

reps_null  <- as.integer(Sys.getenv('REPS_NULL',  '2000'))
reps_power <- as.integer(Sys.getenv('REPS_POWER', '1000'))
n_cores    <- as.integer(Sys.getenv('N_CORES',
                  as.character(max(1L, parallel::detectCores() - 1L))))
N_levels   <- as.integer(strsplit(Sys.getenv('N_LEVELS', '25,35,70,100'),
                                   ',')[[1]])
c_bm_levels <- as.numeric(strsplit(Sys.getenv('CBM_LEVELS', '0,0.45'),
                                    ',')[[1]])
master_seed <- 20260507L

RNGkind("L'Ecuyer-CMRG")
set.seed(master_seed)

t_start <- Sys.time()

## --------------------------------------------------------------
## Trial design: OL+BDC (open-label run-in then blinded
## discontinuation). N is split across two paths.
## --------------------------------------------------------------

td <- buildtrialdesign(
  name_longform = 'OL+BDC',
  name_shortform = 'OLBDC',
  timepoints   = c(4, 8, 12, 16, 17, 18, 19, 20),
  timeptnames  = c('OL1', 'OL2', 'OL3', 'OL4',
                   'BD1', 'BD2', 'BD3', 'BD4'),
  expectancies = c(1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5),
  ondrug = list(
    pathA = c(1, 1, 1, 1, 1, 1, 0, 0),
    pathB = c(1, 1, 1, 1, 1, 0, 0, 0)
  )
)

resp_param <- data.table(
  cat  = c('tv', 'pb', 'br'),
  max  = c(10.98604, 6.50647, 10.98604),
  disp = c(5, 5, 5),
  rate = c(0.42, 0.35, 0.42),
  sd   = c(5, 2, 5)
)

baseline_param <- data.table(
  cat = c('bm', 'BL'),
  m   = c(0, 70),
  sd  = c(1, 10)
)

sd_br <- resp_param[cat == 'br', sd]

split_N <- function(N_total) c(ceiling(N_total / 2), floor(N_total / 2))

simulate_one_dataset <- function(c_bm, N_total, t1half = 0) {
  npp <- split_N(N_total)
  mp1 <- data.table(N = npp[1], c.bm = c_bm,
                    carryover_t1half = t1half,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  mp2 <- data.table(N = npp[2], c.bm = c_bm,
                    carryover_t1half = t1half,
                    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
                    c.cf1t = 0.1, c.cfct = 0.05)
  d1 <- generateData(mp1, resp_param, baseline_param,
                     td$trialpaths[[1]], empirical = FALSE,
                     makePositiveDefinite = TRUE,
                     dgp_architecture = 'mean_moderation')
  d1[, path := 1]
  d2 <- generateData(mp2, resp_param, baseline_param,
                     td$trialpaths[[2]], empirical = FALSE,
                     makePositiveDefinite = TRUE,
                     dgp_architecture = 'mean_moderation')
  d2[, path := 2]
  d2[, ptID := ptID + max(d1$ptID)]
  rbind(d1, d2)
}

build_long_data <- function(dat) {
  groups <- sort(unique(dat$path))
  long_list <- vector('list', length(groups))
  for (g in seq_along(groups)) {
    tdg <- copy(td$trialpaths[[g]])
    tdg <- rbind(list(timeptnames = 'BL', t_wk = 0, e = 0,
                      tod = 0, tsd = 0, tpb = 0), tdg)
    tdg[, t := cumulative(t_wk)]
    timeptnames <- union(tdg$timeptnames, 'BL')
    datg <- dat[path == groups[g]]
    eval_str <- paste0('long <- datg[,.(ptID,bm,',
                       paste(timeptnames, collapse = ','), ')]')
    eval(parse(text = eval_str))
    long_m <- melt(long, id.vars = c('ptID', 'bm'),
                   measure.vars = tdg$timeptnames,
                   variable.name = 'timeptnames',
                   value.name = 'Sx', na.rm = FALSE)
    long_m <- merge(long_m,
                    tdg[, .(timeptnames, t,
                            Db = (tod > 0), tsd)],
                    by = 'timeptnames', all = TRUE)
    long_m[, Dbc := as.numeric(Db)]
    long_m[, path := groups[g]]
    long_list[[g]] <- long_m
  }
  out <- rbindlist(long_list, fill = TRUE)
  out[!is.na(Sx) & !is.na(bm)]
}

na_arm <- function() list(p = NA_real_, beta = NA_real_, se = NA_real_,
                          ci_lo = NA_real_, ci_hi = NA_real_,
                          converged = FALSE)

## --------------------------------------------------------------
## Procedure 1: strict RM-ANOVA (dichotomised, pre-averaged).
## --------------------------------------------------------------

fit_rmanova <- function(long) {
  med_bm <- median(unique(long[, .(ptID, bm)])$bm)
  long2 <- copy(long)
  long2[, bm_high := factor(bm > med_bm, levels = c(FALSE, TRUE))]
  long2[, Db_f := factor(Db, levels = c(FALSE, TRUE))]
  agg <- long2[, .(Sx = mean(Sx, na.rm = TRUE)),
               by = .(ptID, bm_high, Db_f)]
  ok_pts <- agg[, .N, by = ptID][N == 2L, ptID]
  agg <- agg[ptID %in% ok_pts]
  if (nrow(agg) < 8L || length(unique(agg$bm_high)) < 2L) return(na_arm())
  agg[, ptID := factor(ptID)]
  fit <- tryCatch(aov(Sx ~ bm_high * Db_f + Error(ptID / Db_f), data = agg),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(fit)) return(na_arm())
  s <- summary(fit)
  p_int <- NA_real_
  for (stratum in s) {
    tab <- stratum[[1]]
    rn <- trimws(rownames(tab))
    hit <- which(rn %in% c('bm_high:Db_f', 'Db_f:bm_high'))
    if (length(hit) > 0) { p_int <- tab[hit[1], 'Pr(>F)']; break }
  }
  list(p = as.numeric(p_int), beta = NA_real_, se = NA_real_,
       ci_lo = NA_real_, ci_hi = NA_real_, converged = !is.na(p_int))
}

## --------------------------------------------------------------
## Procedure 2: continuous-biomarker within-subject contrast.
## Delta_l = mean(Sx | on-drug) - mean(Sx | off-drug); regress on bm.
## --------------------------------------------------------------

fit_contrast <- function(long) {
  sub <- long[, .(on  = mean(Sx[Db],  na.rm = TRUE),
                  off = mean(Sx[!Db], na.rm = TRUE),
                  bm  = bm[1]), by = ptID]
  sub <- sub[is.finite(on) & is.finite(off)]
  sub[, Delta := on - off]
  if (nrow(sub) < 5L || length(unique(sub$bm)) < 3L) return(na_arm())
  fit <- tryCatch(lm(Delta ~ bm, data = sub),
                  error = function(e) NULL)
  if (is.null(fit)) return(na_arm())
  sm <- summary(fit)$coefficients
  if (!('bm' %in% rownames(sm))) return(na_arm())
  est <- sm['bm', 'Estimate']; se <- sm['bm', 'Std. Error']
  p   <- sm['bm', 'Pr(>|t|)']; df <- fit$df.residual
  tc  <- qt(0.975, df)
  list(p = as.numeric(p), beta = as.numeric(est), se = as.numeric(se),
       ci_lo = est - tc * se, ci_hi = est + tc * se, converged = TRUE)
}

## --------------------------------------------------------------
## Procedure 3: linear mixed-effects with corCAR1.
## --------------------------------------------------------------

fit_lme <- function(dat) {
  op <- list(useDE = TRUE, t_random_slope = FALSE,
             full_model_out = FALSE, carryover_t1half = 0,
             simplecarryover = FALSE, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res) || is.na(res$p) || is.na(res$beta)) return(na_arm())
  se <- as.numeric(res$betaSE)
  beta <- as.numeric(res$beta)
  ci_lo <- if (is.finite(se)) beta - 1.96 * se else NA_real_
  ci_hi <- if (is.finite(se)) beta + 1.96 * se else NA_real_
  list(p = as.numeric(res$p), beta = beta, se = se,
       ci_lo = ci_lo, ci_hi = ci_hi, converged = TRUE)
}

## --------------------------------------------------------------
## Procedure 4: GEE with naive sandwich variance, AR(1).
## --------------------------------------------------------------

fit_gee <- function(long) {
  if (!have_geepack) return(c(na_arm(), list(fit = NULL, long2 = NULL)))
  long2 <- copy(long); setorder(long2, ptID, t)
  long2[, ptID_f := as.integer(factor(ptID))]
  fit <- tryCatch(geepack::geeglm(Sx ~ bm + t + Dbc + bm:Dbc,
                    id = ptID_f, data = long2,
                    family = gaussian(), corstr = 'ar1'),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(fit)) return(c(na_arm(), list(fit = NULL, long2 = NULL)))
  cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(cf)) return(c(na_arm(), list(fit = NULL, long2 = NULL)))
  rn <- rownames(cf); hit <- which(rn %in% c('bm:Dbc', 'Dbc:bm'))
  if (length(hit) == 0L) return(c(na_arm(), list(fit = NULL, long2 = NULL)))
  p_col <- intersect(c('Pr(>|W|)', 'Pr(>W)'), colnames(cf))[1]
  est <- as.numeric(cf[hit[1], 'Estimate'])
  se  <- as.numeric(cf[hit[1], 'Std.err'])
  list(p = as.numeric(cf[hit[1], p_col]), beta = est, se = se,
       ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       converged = TRUE, fit = fit, long2 = long2)
}

## --------------------------------------------------------------
## Procedure 5: Mancl-DeRouen (2001) small-sample sandwich.
## --------------------------------------------------------------

mancl_derouen_se <- function(gee_fit, dat_long) {
  X <- model.matrix(gee_fit$formula, data = dat_long)
  if (!setequal(rownames(X), as.character(seq_len(nrow(dat_long))))) {
    keep <- complete.cases(dat_long[, all.vars(gee_fit$formula),
                                    with = FALSE])
    dat_long <- dat_long[keep]
    X <- model.matrix(gee_fit$formula, data = dat_long)
  }
  y <- model.response(model.frame(gee_fit$formula, data = dat_long))
  beta <- coef(gee_fit); p <- length(beta)
  res <- as.numeric(y - X %*% beta)
  ids <- dat_long$ptID_f; uniq <- unique(ids); K <- length(uniq)
  alpha <- tryCatch({
    a <- gee_fit$geese$alpha
    if (is.null(a) || length(a) == 0L) a <- summary(gee_fit)$corr[1, 'Estimate']
    as.numeric(a)
  }, error = function(e) 0)
  if (is.na(alpha) || abs(alpha) >= 0.999) alpha <- sign(alpha) * 0.95
  phi <- as.numeric(summary(gee_fit)$dispersion[1, 'Estimate'])
  if (is.na(phi) || phi <= 0) phi <- var(res, na.rm = TRUE)
  bread <- matrix(0, p, p)
  for (k in seq_len(K)) {
    rows <- which(ids == uniq[k]); n_i <- length(rows)
    Xi <- X[rows, , drop = FALSE]
    R_i <- alpha ^ abs(outer(seq_len(n_i), seq_len(n_i), '-'))
    V_inv <- tryCatch(solve(phi * R_i), error = function(e) NULL)
    if (is.null(V_inv)) next
    bread <- bread + crossprod(Xi, V_inv) %*% Xi
  }
  bread_inv <- tryCatch(solve(bread), error = function(e) NULL)
  if (is.null(bread_inv)) return(list(se = NA_real_, p = NA_real_,
                                      beta = NA_real_))
  meat <- matrix(0, p, p)
  for (k in seq_len(K)) {
    rows <- which(ids == uniq[k]); n_i <- length(rows)
    Xi <- X[rows, , drop = FALSE]; ri <- res[rows]
    R_i <- alpha ^ abs(outer(seq_len(n_i), seq_len(n_i), '-'))
    V_inv <- tryCatch(solve(phi * R_i), error = function(e) NULL)
    if (is.null(V_inv)) next
    H_ii <- Xi %*% bread_inv %*% crossprod(Xi, V_inv)
    IminusH_inv <- tryCatch(solve(diag(n_i) - H_ii),
                            error = function(e) NULL)
    if (is.null(IminusH_inv)) next
    u_i <- crossprod(Xi, V_inv) %*% (IminusH_inv %*% ri)
    meat <- meat + tcrossprod(u_i)
  }
  vcov_md <- bread_inv %*% meat %*% bread_inv
  se <- sqrt(diag(vcov_md)); names(se) <- names(beta)
  hit_name <- intersect(c('bm:Dbc', 'Dbc:bm'), names(beta))
  if (length(hit_name) == 0L) return(list(se = NA_real_, p = NA_real_,
                                          beta = NA_real_))
  z <- as.numeric(beta[hit_name[1]] / se[hit_name[1]])
  list(se = as.numeric(se[hit_name[1]]), p = 2 * pnorm(-abs(z)),
       beta = as.numeric(beta[hit_name[1]]))
}

fit_gee_mancl <- function(gee_res) {
  if (is.null(gee_res$fit) || !gee_res$converged) return(na_arm())
  out <- tryCatch(mancl_derouen_se(gee_res$fit, gee_res$long2),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out) || is.na(out$p)) return(na_arm())
  est <- as.numeric(out$beta); se <- as.numeric(out$se)
  list(p = as.numeric(out$p), beta = est, se = se,
       ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       converged = TRUE)
}

## --------------------------------------------------------------
## One replicate across all procedures.
## --------------------------------------------------------------

proc_levels <- c('rmanova', 'contrast_cont', 'lme', 'gee', 'gee_mancl')
if (!have_geepack) proc_levels <- c('rmanova', 'contrast_cont', 'lme')

run_one_rep <- function(c_bm, N_total) {
  dat <- tryCatch(simulate_one_dataset(c_bm, N_total),
                  error = function(e) NULL)
  if (is.null(dat)) {
    z <- na_arm()
    return(setNames(replicate(length(proc_levels), z, simplify = FALSE),
                    proc_levels))
  }
  long <- build_long_data(dat)
  out <- list(rmanova = fit_rmanova(long),
              contrast_cont = fit_contrast(long),
              lme = fit_lme(dat))
  if (have_geepack) {
    gee_res <- fit_gee(long)
    out$gee <- list(p = gee_res$p, beta = gee_res$beta, se = gee_res$se,
                    ci_lo = gee_res$ci_lo, ci_hi = gee_res$ci_hi,
                    converged = gee_res$converged)
    out$gee_mancl <- fit_gee_mancl(gee_res)
  }
  out
}

## --------------------------------------------------------------
## Reference run to pin the interaction estimand (beta_true).
## Large-N lme fit at c.bm = 0.45; averaged over a few draws.
## Analytic expectation: c.bm * sd_br = 0.45 * 5 = 2.25.
## --------------------------------------------------------------

cat('\nEstimating beta_true (reference, N=2000 x 6 draws)...\n')
bt_draws <- unlist(lapply(seq_len(6L), function(i) {
  d <- tryCatch(simulate_one_dataset(0.45, 2000L), error = function(e) NULL)
  if (is.null(d)) return(NA_real_)
  fit_lme(d)$beta
}))
beta_true <- mean(bt_draws, na.rm = TRUE)
cat(sprintf('beta_true = %.4f (analytic c.bm*sd_br = %.2f); draws: %s\n',
            beta_true, 0.45 * sd_br,
            paste(sprintf('%.3f', bt_draws), collapse = ', ')))
writeLines(c(sprintf('beta_true_mc\t%.6f', beta_true),
             sprintf('beta_true_analytic\t%.6f', 0.45 * sd_br),
             sprintf('n_draws\t%d', sum(is.finite(bt_draws)))),
           'analysis/data/quick-sim/08-test-betatrue.txt')

## --------------------------------------------------------------
## Build the job list (cell x rep) and run in parallel.
## --------------------------------------------------------------

cells <- CJ(c.bm = c_bm_levels, N = N_levels, sorted = FALSE)
setorder(cells, c.bm, N)
cells[, nrep := ifelse(c.bm == 0, reps_null, reps_power)]

jobs <- rbindlist(lapply(seq_len(nrow(cells)), function(i)
  data.table(cell_i = i, c.bm = cells$c.bm[i], N = cells$N[i],
             rep_idx = seq_len(cells$nrep[i]))))

cat(sprintf('\nJobs: %d (cells: %d; null reps %d, power reps %d).\n',
            nrow(jobs), nrow(cells), reps_null, reps_power))
cat(sprintf('Cores: %d. Procedures: %s.\n', n_cores,
            paste(proc_levels, collapse = ', ')))

run_job <- function(j) {
  pr <- run_one_rep(jobs$c.bm[j], jobs$N[j])
  rbindlist(lapply(proc_levels, function(proc) {
    a <- pr[[proc]]
    data.table(procedure = proc, c.bm = jobs$c.bm[j], N = jobs$N[j],
               rep_idx = jobs$rep_idx[j], beta = a$beta, se = a$se,
               ci_lo = a$ci_lo, ci_hi = a$ci_hi,
               p_interaction = a$p, converged = a$converged)
  }))
}

res_list <- mclapply(seq_len(nrow(jobs)), run_job,
                     mc.cores = n_cores, mc.set.seed = TRUE,
                     mc.preschedule = TRUE)

ok_jobs <- vapply(res_list, is.data.table, logical(1))
if (any(!ok_jobs)) cat(sprintf('WARNING: %d job(s) failed and were dropped.\n',
                               sum(!ok_jobs)))
replicates <- rbindlist(res_list[ok_jobs])

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('\nDone in %.1f s. Total fits: %d.\n',
            elapsed_total, nrow(replicates)))

dir.create('analysis/data/quick-sim', showWarnings = FALSE, recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE, recursive = TRUE)
saveRDS(replicates, 'analysis/data/quick-sim/08-test-replicates.rds')

## --------------------------------------------------------------
## Extended ADEMP summary: rejection rate + bias/SE/coverage.
## --------------------------------------------------------------

summary_dt <- replicates[, {
  n_total <- .N
  ok <- converged & !is.na(p_interaction)
  n_conv <- sum(ok)
  pwr <- if (n_conv > 0L) mean(p_interaction[ok] < 0.05) else NA_real_
  mcse_pwr <- if (n_conv > 0L) sqrt(pwr * (1 - pwr) / n_conv) else NA_real_
  ## Per-cell estimand: 0 under the null, beta_true under c.bm = 0.45.
  bt <- if (c.bm[1] == 0) 0 else beta_true
  has_beta <- any(is.finite(beta[ok]))
  if (has_beta) {
    b   <- beta[ok]; s <- se[ok]; lo <- ci_lo[ok]; hi <- ci_hi[ok]
    bias    <- mean(b, na.rm = TRUE) - bt
    emp_se  <- sd(b, na.rm = TRUE)
    mod_se  <- mean(s, na.rm = TRUE)
    cover   <- mean(lo <= bt & bt <= hi, na.rm = TRUE)
  } else {
    bias <- NA_real_; emp_se <- NA_real_; mod_se <- NA_real_; cover <- NA_real_
  }
  list(n_reps = n_total, rejection_rate = pwr, mcse_rate = mcse_pwr,
       conv_rate = n_conv / n_total, bias = bias, emp_se = emp_se,
       mean_model_se = mod_se,
       se_ratio = if (is.finite(emp_se) && emp_se > 0) mod_se / emp_se else NA_real_,
       coverage = cover)
}, by = .(procedure, c.bm, N)]

setorder(summary_dt, procedure, c.bm, N)
write.table(summary_dt, file = 'analysis/data/quick-sim/08-test-summary.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)
cat('\nExtended Morris summary:\n'); print(summary_dt)

## --------------------------------------------------------------
## geesmv benchmark for the inline Mancl-DeRouen SE (probe cell).
## --------------------------------------------------------------

if (have_geepack && have_geesmv) {
  cat('\nBenchmarking inline Mancl-DeRouen vs geesmv::GEE.var.md ...\n')
  bench <- vector('list', 0); set.seed(master_seed + 99L)
  for (b in seq_len(25L)) {
    d <- tryCatch(simulate_one_dataset(0.45, 35L), error = function(e) NULL)
    if (is.null(d)) next
    long <- build_long_data(d)
    gr <- fit_gee(long)
    if (is.null(gr$fit) || !gr$converged) next
    inline <- tryCatch(mancl_derouen_se(gr$fit, gr$long2),
                       error = function(e) NULL)
    ref <- tryCatch({
      l2 <- as.data.frame(gr$long2)
      gv <- suppressMessages(utils::capture.output(
        gvr <- geesmv::GEE.var.md(Sx ~ bm + t + Dbc + bm:Dbc, id = 'ptID_f',
                                  family = gaussian, data = l2,
                                  corstr = 'AR-M')))
      nm <- names(coef(gr$fit))
      idx <- which(nm %in% c('bm:Dbc', 'Dbc:bm'))[1]
      sqrt(gvr$cov.beta[idx])
    }, error = function(e) NULL)
    if (!is.null(inline) && !is.null(ref) && is.finite(inline$se) &&
        is.finite(ref)) {
      bench[[length(bench) + 1L]] <- data.table(
        rep = b, inline_se = inline$se, geesmv_se = ref,
        ratio = inline$se / ref)
    }
  }
  if (length(bench) > 0L) {
    bench_dt <- rbindlist(bench)
    cat(sprintf('  n=%d; mean inline/geesmv SE ratio = %.4f (sd %.4f)\n',
                nrow(bench_dt), mean(bench_dt$ratio), sd(bench_dt$ratio)))
    write.table(bench_dt,
                'analysis/data/quick-sim/08-test-mancl-benchmark.txt',
                sep = '\t', row.names = FALSE, quote = FALSE)
  } else {
    cat('  benchmark produced no usable pairs.\n')
  }
}

## --------------------------------------------------------------
## Power/Type-I figure.
## --------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, role := ifelse(c.bm == 0, 'Type-I (c.bm = 0)', 'Power (c.bm = 0.45)')]
plot_dt[, lo := pmax(0, rejection_rate - 1.96 * mcse_rate)]
plot_dt[, hi := pmin(1, rejection_rate + 1.96 * mcse_rate)]
proc_factor_levels <- c('rmanova', 'contrast_cont', 'lme', 'gee', 'gee_mancl')
proc_factor_labels <- c('RM-ANOVA (dich.)', 'Contrast (cont. bm)',
                        'LME (corCAR1)', 'GEE (AR1)', 'GEE Mancl-DeRouen')
plot_dt[, proc_lbl := factor(procedure, levels = proc_factor_levels,
                             labels = proc_factor_labels)]
plot_dt[, N_lbl := factor(sprintf('N = %d', N),
                          levels = sprintf('N = %d', sort(N_levels)))]

p_fig <- ggplot(plot_dt, aes(x = proc_lbl, y = rejection_rate,
                             ymin = lo, ymax = hi)) +
  geom_col(fill = 'grey75', colour = 'grey30') +
  geom_errorbar(width = 0.2) +
  geom_text(aes(label = sprintf('%.2f', rejection_rate)), vjust = -0.6,
            size = 2.5) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
  facet_grid(N_lbl ~ role) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = 'Test procedure', y = 'Rejection rate (p < 0.05)',
       title = 'Test-procedure comparison for the biomarker x treatment interaction',
       subtitle = sprintf('OL+BDC, Architecture A; null reps=%d, power reps=%d',
                          reps_null, reps_power),
       caption = 'Error bars: 1.96 x MCSE (binomial); dashed line at 0.05.') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave('analysis/figures/quick-sim/08-test-power.pdf', p_fig,
       width = 11.0, height = 6.5)

cat('\nWrote replicates, summary, beta_true, benchmark, figure.\n')
invisible(NULL)
