## 08-test-procedure-prototype.R
##
## 30-min prototype simulation for paper
## 08-test-procedure-design-sensitivity.
##
## Compares four test procedures for the biomarker x treatment
## interaction, across two sample sizes:
##   1. Strict repeated-measures ANOVA (split-plot, between-subject
##      median-dichotomised biomarker, two within-Db levels).
##   2. Linear mixed-effects with corCAR1 residual correlation
##      (existing lme_analysis()).
##   3. GEE with sandwich variance (geepack::geeglm), AR(1) working
##      correlation, model-based Wald test from summary().
##   4. GEE with the Mancl-DeRouen (2001) small-sample bias-corrected
##      sandwich variance (hand-rolled; same point estimate as #3 but
##      with a small-sample-corrected SE and Wald-z p-value).
##
## Grid: procedure {rmanova, lme, gee, gee_mancl_smallN} x
##   c.bm {0, 0.45} x N {25, 35} = 16 cells, target 800 reps.
##
## Outputs (overwrite):
##   analysis/data/quick-sim/08-test-replicates.rds
##   analysis/data/quick-sim/08-test-summary.txt
##   analysis/figures/quick-sim/08-test-power.pdf

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(ggplot2)
  library(nlme)
})

have_geepack <- requireNamespace('geepack', quietly = TRUE)
if (have_geepack) {
  suppressPackageStartupMessages(library(geepack))
}
cat('geepack available:', have_geepack, '\n')

## Mancl-DeRouen is implemented inline below using only geepack
## internals; flag separately for clarity.
have_mancl <- have_geepack
cat('gee_mancl_smallN implemented inline:', have_mancl, '\n')

t_start <- Sys.time()
## Internal budget: target 1700 s overall. The harness runs each
## Bash call up to 10 min; the orchestrator can chain multiple
## passes by appending to the replicates RDS. For a single
## chained pass, the effective budget is set to BUDGET_SECS env
## var if present, defaulting to 1700.
budget_secs <- as.numeric(Sys.getenv('BUDGET_SECS', '1700'))
n_reps_target <- 800L

set.seed(20260507)

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

split_N <- function(N_total) c(ceiling(N_total / 2), floor(N_total / 2))

## --------------------------------------------------------------
## Simulate one trial dataset (wide form).
## --------------------------------------------------------------

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
                     td$trialpaths[[1]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE,
                     dgp_architecture = 'mean_moderation')
  d1[, path := 1]
  d2 <- generateData(mp2, resp_param, baseline_param,
                     td$trialpaths[[2]],
                     empirical = FALSE,
                     makePositiveDefinite = TRUE,
                     dgp_architecture = 'mean_moderation')
  d2[, path := 2]
  d2[, ptID := ptID + max(d1$ptID)]
  rbind(d1, d2)
}

## --------------------------------------------------------------
## Build long-form data with bm, t, Dbc columns.
## --------------------------------------------------------------

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

## --------------------------------------------------------------
## Procedure 1: strict RM-ANOVA on the bm_high x Db interaction.
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
  if (nrow(agg) < 8L || length(unique(agg$bm_high)) < 2L) {
    return(list(p = NA_real_, beta = NA_real_, converged = FALSE))
  }
  agg[, ptID := factor(ptID)]
  fit <- tryCatch(
    aov(Sx ~ bm_high * Db_f + Error(ptID / Db_f), data = agg),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) {
    return(list(p = NA_real_, beta = NA_real_, converged = FALSE))
  }
  s <- summary(fit)
  p_int <- NA_real_
  for (stratum in s) {
    tab <- stratum[[1]]
    rn <- trimws(rownames(tab))
    hit <- which(rn %in% c('bm_high:Db_f', 'Db_f:bm_high'))
    if (length(hit) > 0) {
      p_int <- tab[hit[1], 'Pr(>F)']
      break
    }
  }
  list(p = as.numeric(p_int), beta = NA_real_,
       converged = !is.na(p_int))
}

## --------------------------------------------------------------
## Procedure 2: linear mixed-effects with corCAR1.
## --------------------------------------------------------------

fit_lme <- function(dat) {
  op <- list(useDE = TRUE, t_random_slope = FALSE,
             full_model_out = FALSE, carryover_t1half = 0,
             simplecarryover = FALSE, carryover_scalefactor = 1)
  res <- tryCatch(
    lme_analysis(td$trialpaths, dat, op),
    error = function(e) NULL
  )
  if (is.null(res) || is.na(res$p)) {
    return(list(p = NA_real_, beta = NA_real_, converged = FALSE))
  }
  list(p = as.numeric(res$p), beta = as.numeric(res$beta),
       converged = TRUE)
}

## --------------------------------------------------------------
## Procedure 3: GEE with model-based sandwich variance, AR(1).
## --------------------------------------------------------------

fit_gee <- function(long) {
  if (!have_geepack) return(list(p = NA_real_, beta = NA_real_,
                                 converged = FALSE,
                                 fit = NULL))
  long2 <- copy(long)
  setorder(long2, ptID, t)
  long2[, ptID_f := as.integer(factor(ptID))]
  fit <- tryCatch(
    geepack::geeglm(Sx ~ bm + t + Dbc + bm:Dbc,
                    id = ptID_f, data = long2,
                    family = gaussian(), corstr = 'ar1'),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(fit)) return(list(p = NA_real_, beta = NA_real_,
                                 converged = FALSE, fit = NULL))
  cf <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(cf)) return(list(p = NA_real_, beta = NA_real_,
                               converged = FALSE, fit = NULL))
  rn <- rownames(cf)
  hit <- which(rn %in% c('bm:Dbc', 'Dbc:bm'))
  if (length(hit) == 0L) return(list(p = NA_real_, beta = NA_real_,
                                     converged = FALSE,
                                     fit = NULL))
  p_col <- intersect(c('Pr(>|W|)', 'Pr(>W)'), colnames(cf))[1]
  list(p = as.numeric(cf[hit[1], p_col]),
       beta = as.numeric(cf[hit[1], 'Estimate']),
       converged = TRUE,
       fit = fit, long2 = long2)
}

## --------------------------------------------------------------
## Procedure 4: GEE with Mancl-DeRouen (2001) small-sample
## bias-corrected sandwich SE.
##
## Mancl & DeRouen (2001) Biometrics 57:126-134.
## For cluster i, replace the residual r_i in the sandwich meat
## with a leverage-corrected (I - H_ii)^{-1} r_i, where
##   H_ii = X_i (X' V^{-1} X)^{-1} X_i' V_i^{-1}
## with V_i the working covariance (working AR(1) correlation
## inflated by the GEE-estimated dispersion).
##
## This is implemented from the geeglm fit returned by fit_gee().
## We reuse the same beta estimate; only the SE, Wald z, and the
## resulting p-value are replaced.
## --------------------------------------------------------------

mancl_derouen_se <- function(gee_fit, dat_long) {
  X <- model.matrix(gee_fit$formula, data = dat_long)
  if (!setequal(rownames(X), as.character(seq_len(nrow(dat_long))))) {
    ## Re-align rows used in the fit (NA-handling safe).
    keep <- complete.cases(dat_long[, all.vars(gee_fit$formula),
                                    with = FALSE])
    dat_long <- dat_long[keep]
    X <- model.matrix(gee_fit$formula, data = dat_long)
  }
  y <- model.response(model.frame(gee_fit$formula, data = dat_long))
  beta <- coef(gee_fit)
  p <- length(beta)
  res <- as.numeric(y - X %*% beta)

  ids <- dat_long$ptID_f
  uniq <- unique(ids)
  K <- length(uniq)

  ## Working AR(1) correlation; geepack stores alpha in
  ## summary(fit)$corr$Estimate (or fit$geese$alpha for some
  ## versions). Use a robust extraction.
  alpha <- tryCatch(
    {
      a <- gee_fit$geese$alpha
      if (is.null(a) || length(a) == 0L) {
        a <- summary(gee_fit)$corr[1, 'Estimate']
      }
      as.numeric(a)
    },
    error = function(e) 0
  )
  if (is.na(alpha) || abs(alpha) >= 0.999) alpha <- sign(alpha) * 0.95
  phi <- as.numeric(summary(gee_fit)$dispersion[1, 'Estimate'])
  if (is.na(phi) || phi <= 0) phi <- var(res, na.rm = TRUE)

  ## First pass: bread = (X' V^{-1} X)^{-1}
  bread <- matrix(0, p, p)
  for (k in seq_len(K)) {
    rows <- which(ids == uniq[k])
    n_i <- length(rows)
    Xi <- X[rows, , drop = FALSE]
    R_i <- alpha ^ abs(outer(seq_len(n_i), seq_len(n_i), '-'))
    V_i <- phi * R_i
    V_inv <- tryCatch(solve(V_i), error = function(e) NULL)
    if (is.null(V_inv)) next
    bread <- bread + crossprod(Xi, V_inv) %*% Xi
  }
  bread_inv <- tryCatch(solve(bread), error = function(e) NULL)
  if (is.null(bread_inv)) return(list(se = NA_real_,
                                       p = NA_real_))

  ## Second pass: meat with Mancl-DeRouen (I - H_ii)^{-1} adjustment.
  meat <- matrix(0, p, p)
  for (k in seq_len(K)) {
    rows <- which(ids == uniq[k])
    n_i <- length(rows)
    Xi <- X[rows, , drop = FALSE]
    ri <- res[rows]
    R_i <- alpha ^ abs(outer(seq_len(n_i), seq_len(n_i), '-'))
    V_i <- phi * R_i
    V_inv <- tryCatch(solve(V_i), error = function(e) NULL)
    if (is.null(V_inv)) next
    H_ii <- Xi %*% bread_inv %*% crossprod(Xi, V_inv)
    IminusH <- diag(n_i) - H_ii
    IminusH_inv <- tryCatch(solve(IminusH),
                            error = function(e) NULL)
    if (is.null(IminusH_inv)) next
    r_corr <- IminusH_inv %*% ri
    u_i <- crossprod(Xi, V_inv) %*% r_corr
    meat <- meat + tcrossprod(u_i)
  }

  vcov_md <- bread_inv %*% meat %*% bread_inv
  se <- sqrt(diag(vcov_md))
  names(se) <- names(beta)

  hit_name <- intersect(c('bm:Dbc', 'Dbc:bm'), names(beta))
  if (length(hit_name) == 0L) return(list(se = NA_real_,
                                           p = NA_real_))
  z <- as.numeric(beta[hit_name[1]] / se[hit_name[1]])
  pval <- 2 * pnorm(-abs(z))
  list(se = as.numeric(se[hit_name[1]]),
       p = pval,
       beta = as.numeric(beta[hit_name[1]]))
}

fit_gee_mancl <- function(gee_res) {
  if (is.null(gee_res$fit) || !gee_res$converged) {
    return(list(p = NA_real_, beta = NA_real_, converged = FALSE))
  }
  out <- tryCatch(
    mancl_derouen_se(gee_res$fit, gee_res$long2),
    error = function(e) NULL,
    warning = function(w) NULL
  )
  if (is.null(out) || is.na(out$p)) {
    return(list(p = NA_real_, beta = gee_res$beta,
                converged = FALSE))
  }
  list(p = as.numeric(out$p),
       beta = as.numeric(out$beta),
       converged = TRUE)
}

## --------------------------------------------------------------
## One full rep across all procedures for a (c.bm, N) cell.
## --------------------------------------------------------------

run_one_rep <- function(c_bm, N_total) {
  dat <- tryCatch(simulate_one_dataset(c_bm, N_total),
                  error = function(e) NULL)
  if (is.null(dat)) {
    na_row <- list(p = NA_real_, beta = NA_real_, converged = FALSE)
    return(list(rmanova = na_row, lme = na_row,
                gee = na_row, gee_mancl_smallN = na_row))
  }
  long <- build_long_data(dat)
  rm_res <- fit_rmanova(long)
  lme_res <- fit_lme(dat)
  gee_res <- fit_gee(long)
  gee_md  <- fit_gee_mancl(gee_res)
  gee_pub <- list(p = gee_res$p, beta = gee_res$beta,
                  converged = gee_res$converged)
  list(rmanova = rm_res, lme = lme_res,
       gee = gee_pub, gee_mancl_smallN = gee_md)
}

## --------------------------------------------------------------
## Main loop over the 16-cell grid.
## --------------------------------------------------------------

c_bm_levels <- c(0, 0.45)
N_levels <- c(25L, 35L)
proc_levels <- c('rmanova', 'lme', 'gee', 'gee_mancl_smallN')
if (!have_mancl) proc_levels <- c('rmanova', 'lme', 'gee')

cat(sprintf(
  '\nMain loop: %d c.bm x %d N x %d procedures, target %d reps/cell.\n',
  length(c_bm_levels), length(N_levels), length(proc_levels),
  n_reps_target))
cat(sprintf('Budget: %d s.\n', budget_secs))

results <- vector('list', 0)
row_i <- 0L
fits_done <- 0L
aborted <- FALSE
reps_done <- list()

cells <- CJ(c.bm = c_bm_levels, N = N_levels, sorted = FALSE)
setorder(cells, c.bm, N)

for (cell_i in seq_len(nrow(cells))) {
  cb <- cells$c.bm[cell_i]
  Nv <- cells$N[cell_i]
  cell_key <- sprintf('c.bm=%.2f_N=%d', cb, Nv)
  reps_done[[cell_key]] <- 0L
  for (k in seq_len(n_reps_target)) {
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
    if (elapsed > budget_secs) {
      cat(sprintf(
        'Budget hit at rep %d cell %s (elapsed %.0fs); stopping.\n',
        k, cell_key, elapsed))
      aborted <- TRUE
      break
    }
    pr <- run_one_rep(cb, Nv)
    for (proc in proc_levels) {
      row_i <- row_i + 1L
      results[[row_i]] <- data.table(
        procedure = proc, c.bm = cb, N = Nv, rep_idx = k,
        beta = pr[[proc]]$beta,
        p_interaction = pr[[proc]]$p,
        converged = pr[[proc]]$converged
      )
      fits_done <- fits_done + 1L
      if (fits_done %% 200L == 0L) {
        cat(sprintf(
          '  fits=%d cell=%s rep=%d elapsed=%.0fs\n',
          fits_done, cell_key, k, elapsed))
      }
    }
    reps_done[[cell_key]] <- k
  }
  if (aborted) break
}

replicates <- rbindlist(results)

elapsed_total <- as.numeric(difftime(Sys.time(), t_start, units = 'secs'))
cat(sprintf('\nDone in %.1f s. Aborted: %s. Total fits: %d.\n',
            elapsed_total, aborted, fits_done))
cat('Reps per cell:\n')
print(reps_done)

## --------------------------------------------------------------
## Persist replicate-level data
## --------------------------------------------------------------

dir.create('analysis/data/quick-sim', showWarnings = FALSE,
           recursive = TRUE)
dir.create('analysis/figures/quick-sim', showWarnings = FALSE,
           recursive = TRUE)

saveRDS(replicates,
        'analysis/data/quick-sim/08-test-replicates.rds')

## --------------------------------------------------------------
## Morris ADEMP summary (rejection rate, MCSE, convergence rate)
## --------------------------------------------------------------

summary_dt <- replicates[, {
  n_total <- .N
  ok <- converged & !is.na(p_interaction)
  n_conv <- sum(ok)
  pwr <- if (n_conv > 0L) mean(p_interaction[ok] < 0.05) else NA_real_
  mcse_pwr <- if (n_conv > 0L) sqrt(pwr * (1 - pwr) / n_conv) else NA_real_
  list(n_reps = n_total,
       rejection_rate = pwr,
       mcse_rate = mcse_pwr,
       conv_rate = n_conv / n_total)
}, by = .(procedure, c.bm, N)]

setorder(summary_dt, procedure, c.bm, N)

write.table(summary_dt,
            file = 'analysis/data/quick-sim/08-test-summary.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

cat('\nMorris summary table:\n')
print(summary_dt)

## --------------------------------------------------------------
## Power figure: facet by N x c.bm role
## --------------------------------------------------------------

plot_dt <- copy(summary_dt)
plot_dt[, role := ifelse(c.bm == 0,
                          'Type-I (c.bm = 0)',
                          'Power (c.bm = 0.45)')]
plot_dt[, lo := pmax(0, rejection_rate - 1.96 * mcse_rate)]
plot_dt[, hi := pmin(1, rejection_rate + 1.96 * mcse_rate)]
proc_factor_levels <- c('rmanova', 'lme', 'gee', 'gee_mancl_smallN')
proc_factor_labels <- c('RM-ANOVA', 'LME (corCAR1)',
                        'GEE (AR1)', 'GEE Mancl-DeRouen')
plot_dt[, proc_lbl := factor(procedure,
                             levels = proc_factor_levels,
                             labels = proc_factor_labels)]
plot_dt[, N_lbl := factor(sprintf('N = %d', N),
                          levels = sprintf('N = %d', sort(N_levels)))]

p_fig <- ggplot(plot_dt,
                aes(x = proc_lbl, y = rejection_rate,
                    ymin = lo, ymax = hi)) +
  geom_col(fill = 'grey75', colour = 'grey30') +
  geom_errorbar(width = 0.2) +
  geom_text(aes(label = sprintf('%.2f', rejection_rate)),
            vjust = -0.6, size = 2.9) +
  geom_hline(yintercept = 0.05, linetype = 'dashed', alpha = 0.4) +
  facet_grid(N_lbl ~ role) +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  labs(x = 'Test procedure',
       y = 'Rejection rate (p < 0.05)',
       title = paste('Test-procedure comparison for the',
                     'biomarker x treatment interaction'),
       subtitle = sprintf(
         paste0('OL+BDC, Architecture A (mean moderation), ',
                'reps target = %d/cell, total fits = %d'),
         n_reps_target, fits_done),
       caption = 'Error bars: 1.96 x MCSE (binomial); dashed line at 0.05.') +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave('analysis/figures/quick-sim/08-test-power.pdf',
       p_fig, width = 9.0, height = 6.0)

cat('\nWrote:\n')
cat('  analysis/data/quick-sim/08-test-replicates.rds\n')
cat('  analysis/data/quick-sim/08-test-summary.txt\n')
cat('  analysis/figures/quick-sim/08-test-power.pdf\n')
if (aborted) cat('NOTE: budget aborted - some cells short of target.\n')

invisible(NULL)
