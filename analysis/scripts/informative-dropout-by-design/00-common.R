#' 00-common.R
#'
#' Shared infrastructure for paper 09 (informative-dropout-by-design)
#' production drivers. Sourced by per-study scripts (01, 02, 03, 04)
#' and by run-all.R.
#'
#' Defines:
#'   - Four trial designs (OL, OL+BDC, CO, hybrid) at 20 weeks each.
#'   - Five dropout-pattern parameter sets (β0, β1) per pre-reg.
#'   - apply_hazard_dropout(): quadratic-hazard dropout module
#'     lifted from analysis/scripts/quick-sim/09-dropout-prototype.R.
#'   - Path-tracking helper for Study 3.
#'   - Replicate executor and chunked cell driver.
#'   - Morris ADEMP cell-summary helper with MCSE on every measure.
#'
#' Outputs land under
#' analysis/data/derived_data/informative-dropout-by-design/.
#'
#' Reconciliation note. The pmsimstats package contains R/censordata.R
#' which exposes a (β0, β1, eb1) censoring function with
#' multiplicative shape eb1; the shape exponent and intercept-vs-
#' slope semantics differ from this paper's pre-registered hazard
#' Pr(drop at t) = β0 + β1 · Δ_Sx(t)^2. We use the prototype's
#' apply_hazard_dropout() directly so the simulation matches the
#' published mechanism exactly.

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(parallel)
})

#=====================================================================
# Paths
#=====================================================================

P9_DATA_DIR <- file.path('analysis', 'data', 'derived_data',
                         'informative-dropout-by-design')
P9_FIG_DIR  <- file.path('analysis', 'figures',
                         'informative-dropout-by-design')
if (!dir.exists(P9_DATA_DIR)) dir.create(P9_DATA_DIR, recursive = TRUE)
if (!dir.exists(P9_FIG_DIR))  dir.create(P9_FIG_DIR,  recursive = TRUE)

#=====================================================================
# Trial designs (all 20 weeks total, weekly observations)
#=====================================================================
#
# Per pre-reg: OL (20w on); OL+BDC (16w OL + 4w BDC); CO (10w P1 +
# 10w P2 with within-participant randomisation); hybrid (8w OL + 4w
# BDC + 4w XO1 + 4w XO2 = Hendrickson 2020 primary).

design_OL <- function() {
  buildtrialdesign(
    name_longform  = 'open label 20w',
    name_shortform = 'OL',
    timepoints     = seq(1, 20, by = 1),
    timeptnames    = paste0('OL', 1:20),
    expectancies   = rep(1, 20),
    ondrug = list(pathA = rep(1, 20))
  )
}

design_OL_BDC <- function() {
  buildtrialdesign(
    name_longform  = 'open label + 4w blinded discontinuation',
    name_shortform = 'OLBDC',
    timepoints     = seq(1, 20, by = 1),
    timeptnames    = c(paste0('OL',  1:16),
                       paste0('BDC', 1:4)),
    expectancies   = c(rep(1, 16), rep(0.5, 4)),
    ondrug = list(pathA = c(rep(1, 16), rep(0, 4)))
  )
}

design_CO <- function() {
  ## Two-period crossover, 10w + 10w, two paths (AB and BA).
  buildtrialdesign(
    name_longform  = 'two-period crossover',
    name_shortform = 'CO',
    timepoints     = seq(1, 20, by = 1),
    timeptnames    = c(paste0('P1_', 1:10), paste0('P2_', 1:10)),
    expectancies   = rep(1, 20),
    ondrug = list(
      AB = c(rep(1, 10), rep(0, 10)),
      BA = c(rep(0, 10), rep(1, 10))
    )
  )
}

design_hybrid <- function() {
  ## 8w OL + 4w BDC + 4w XO1 + 4w XO2; two crossover paths.
  buildtrialdesign(
    name_longform  = 'hybrid OL+BDC+crossover (Hendrickson 2020)',
    name_shortform = 'Hybrid',
    timepoints     = seq(1, 20, by = 1),
    timeptnames    = c(paste0('OL',  1:8),
                       paste0('BDC', 1:4),
                       paste0('XO1_', 1:4),
                       paste0('XO2_', 1:4)),
    expectancies   = c(rep(1, 8), rep(0.5, 12)),
    ondrug = list(
      AB = c(rep(1, 8), rep(0, 4), rep(1, 4), rep(0, 4)),
      BA = c(rep(1, 8), rep(0, 4), rep(0, 4), rep(1, 4))
    )
  )
}

ALL_DESIGNS <- function() {
  list(OL    = design_OL(),
       OLBDC = design_OL_BDC(),
       CO    = design_CO(),
       Hybrid = design_hybrid())
}

#=====================================================================
# Dropout patterns (per pre-registration)
#=====================================================================

DROPOUT_PATTERNS <- list(
  none           = list(beta0 = 0,    beta1 = 0  ),
  balanced       = list(beta0 = 0.05, beta1 = 0.5),
  more_of_flat   = list(beta0 = 0.05, beta1 = 0.2),
  more_of_biased = list(beta0 = 0.05, beta1 = 0.8),
  high_dropout   = list(beta0 = 0.15, beta1 = 0.5)
)

#=====================================================================
# Hazard-based dropout (lifted from prototype)
#=====================================================================
#
# Per-visit hazard h(t) = beta0 + beta1 · (Δ_Sx_z(t))^2 where
# Δ_Sx(t) = Sx(t) − BL is standardised across participants by the
# per-visit SD. Once a participant drops, all subsequent
# measurements become NA. No alpha-rescaling: the patterns'
# (β0, β1) values determine the natural dropout rate.

apply_hazard_dropout <- function(dat, trialdesign,
                                  beta0 = 0.05, beta1 = 0.5) {
  if (beta0 == 0 && beta1 == 0) return(dat)
  d <- data.table(trialdesign)
  tnames <- d$timeptnames
  n_pt <- nrow(dat)

  bl <- dat$BL
  sx_mat <- as.matrix(dat[, mget(tnames)])
  delta_mat <- sx_mat - bl
  v_sd <- apply(delta_mat, 2, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-6) 1 else s
  })
  delta_z <- sweep(delta_mat, 2, v_sd, '/')

  haz <- beta0 + beta1 * (delta_z^2)
  haz[haz < 0] <- 0
  haz[haz > 1] <- 1

  draws <- matrix(runif(n_pt * length(tnames)),
                  nrow = n_pt, ncol = length(tnames))
  drop_event <- draws < haz
  first_drop <- apply(drop_event, 1, function(x) {
    w <- which(x)
    if (length(w) == 0) NA_integer_ else w[1]
  })

  for (i in seq_len(n_pt)) {
    fi <- first_drop[i]
    if (!is.na(fi)) {
      for (j in fi:length(tnames)) {
        set(dat, i = i, j = tnames[j], value = NA_real_)
      }
    }
  }

  attr(dat, 'dropout_fraction') <- mean(!is.na(first_drop))
  attr(dat, 'first_drop')       <- first_drop
  dat
}

#=====================================================================
# Parameter constructors (prazosin-PTSD baseline)
#=====================================================================

make_resp_params_09 <- function() {
  data.table(
    cat  = c('tv', 'pb', 'br'),
    max  = c(10.98604, 6.50647, 10.98604),
    disp = c(5, 5, 5),
    rate = c(0.42, 0.35, 0.42),
    sd   = c(5, 2, 5)
  )
}

make_bl_params_09 <- function() {
  data.table(cat = c('bm', 'BL'),
             m   = c(0, 70),
             sd  = c(1, 10))
}

#' @param t1half Carryover half-life in weeks. Pre-registration
#'   baseline: 3 days = 3/7 weeks.
make_model_params_09 <- function(N = 35, c_bm = 0.45,
                                  t1half = 3 / 7) {
  data.table(
    N = N, c.bm = c_bm,
    carryover_t1half = t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
}

#=====================================================================
# Analysis wrapper
#=====================================================================

lme_fit_09 <- function(td, dat, t1half = 3 / 7) {
  op <- list(useDE = TRUE, t_random_slope = FALSE,
             full_model_out = FALSE,
             carryover_t1half = t1half,
             simplecarryover = FALSE,
             carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) {
                    data.table(beta = NA_real_, betaSE = NA_real_,
                               p = NA_real_, issingular = NA,
                               warning = conditionMessage(e))
                  })
  if (!('beta' %in% names(res))) {
    res <- data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, issingular = NA, warning = NA)
  }
  res
}

#=====================================================================
# True-effect helper
#=====================================================================

TRUE_BETA_09 <- function(c_bm, sigma_br = 5, sigma_bm = 1) {
  if (c_bm == 0) 0 else -c_bm * sigma_br / sigma_bm
}

#=====================================================================
# One replicate
#=====================================================================

#' @param cell list with: design_name, dropout_name, N, c_bm,
#'   optionally analysis_name (for Study 4).
one_rep_09 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- ALL_DESIGNS()[[cell$design_name]]
  paths <- td$trialpaths
  mp <- make_model_params_09(N = cell$N, c_bm = cell$c_bm)
  rp <- make_resp_params_09()
  bp <- make_bl_params_09()
  drop <- DROPOUT_PATTERNS[[cell$dropout_name]]

  dat_list <- vector('list', length(paths))
  drop_fracs <- numeric(length(paths))
  for (g in seq_along(paths)) {
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]],
                   empirical = FALSE, makePositiveDefinite = TRUE,
                   dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    di_dropped <- apply_hazard_dropout(di, paths[[g]],
                                        beta0 = drop$beta0,
                                        beta1 = drop$beta1)
    drop_fracs[g] <- attr(di_dropped, 'dropout_fraction') %||% 0
    dat_list[[g]] <- di_dropped
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, dropout_fraction = NA_real_,
                      issingular = NA, warning = 'gen_failed',
                      rep_idx = rep_idx))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)

  res <- lme_fit_09(td, dat)
  res[, dropout_fraction := mean(drop_fracs[good])]
  res[, rep_idx := rep_idx]
  res
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

#=====================================================================
# Cell-level execution with chunked checkpoint loop
#=====================================================================

run_cell_09 <- function(cell, n_reps, study_seed, cell_id,
                        n_cores = max(1L,
                                      parallel::detectCores() - 1L),
                        chunk_size = 50L,
                        progress_every = 100L,
                        time_budget_sec = Inf) {
  t0 <- Sys.time()
  elapsed <- function()
    as.numeric(difftime(Sys.time(), t0, units = 'secs'))

  cell_str <- sprintf(
    'cell %d (design=%s drop=%s N=%d c_bm=%g)', cell_id,
    cell$design_name, cell$dropout_name, cell$N, cell$c_bm)
  cat(sprintf('[%6.1fs] Starting %s; n_reps=%d\n',
              elapsed(), cell_str, n_reps))

  results <- vector('list', 0L)
  done <- 0L
  while (done < n_reps) {
    if (elapsed() > time_budget_sec) {
      cat(sprintf('[%6.1fs] Time budget hit; aborting cell.\n',
                  elapsed()))
      break
    }
    take <- min(chunk_size, n_reps - done)
    chunk <- parallel::mclapply(
      seq.int(done + 1L, done + take),
      function(r) one_rep_09(r, cell, study_seed, cell_id),
      mc.cores = n_cores)
    good <- !vapply(chunk, function(x)
      inherits(x, 'try-error') || is.null(x), logical(1))
    results[[length(results) + 1L]] <- rbindlist(chunk[good],
                                                 fill = TRUE)
    done <- done + take
    if (done %% progress_every < chunk_size) {
      cat(sprintf('[%6.1fs] %s: %d/%d reps\n',
                  elapsed(), cell_str, done, n_reps))
    }
  }
  out <- rbindlist(results, fill = TRUE)
  out[, `:=`(design_name = cell$design_name,
             dropout_name = cell$dropout_name,
             N = cell$N, c_bm = cell$c_bm,
             cell_id = cell_id)]
  setcolorder(out, c('cell_id', 'design_name', 'dropout_name',
                     'N', 'c_bm', 'rep_idx',
                     'beta', 'betaSE', 'p',
                     'dropout_fraction'))
  cat(sprintf('[%6.1fs] %s: complete (%d reps).\n',
              elapsed(), cell_str, done))
  out
}

#=====================================================================
# Morris ADEMP cell summary
#=====================================================================

summarise_cell_09 <- function(d, true_beta) {
  conv <- d[!is.na(beta) & !is.na(p)]
  n <- nrow(conv)
  n_total <- nrow(d)
  if (n == 0) {
    return(data.table(
      n_reps = n_total, n_converged = 0L, conv_rate = 0,
      power = NA_real_, mcse_power = NA_real_,
      mean_beta = NA_real_, sd_beta = NA_real_,
      mcse_mean_beta = NA_real_,
      bias = NA_real_, mcse_bias = NA_real_,
      coverage = NA_real_, mcse_coverage = NA_real_,
      mean_dropout = NA_real_))
  }
  power <- mean(conv$p < 0.05, na.rm = TRUE)
  mcse_power <- sqrt(power * (1 - power) / n)
  mean_beta <- mean(conv$beta, na.rm = TRUE)
  sd_beta   <- sd(conv$beta,   na.rm = TRUE)
  mcse_mean_beta <- sd_beta / sqrt(n)
  bias <- mean_beta - true_beta
  mcse_bias <- mcse_mean_beta
  lo <- conv$beta - 1.96 * conv$betaSE
  hi <- conv$beta + 1.96 * conv$betaSE
  cov_ <- mean(true_beta >= lo & true_beta <= hi, na.rm = TRUE)
  mcse_cov <- sqrt(cov_ * (1 - cov_) / n)
  data.table(
    n_reps = n_total, n_converged = n,
    conv_rate = n / n_total,
    power = power, mcse_power = mcse_power,
    mean_beta = mean_beta, sd_beta = sd_beta,
    mcse_mean_beta = mcse_mean_beta,
    bias = bias, mcse_bias = mcse_bias,
    coverage = cov_, mcse_coverage = mcse_cov,
    mean_dropout = mean(d$dropout_fraction, na.rm = TRUE))
}

#=====================================================================
# Save helpers
#=====================================================================

save_cell_09 <- function(cell_dt, study_id, cell_id) {
  fn <- file.path(P9_DATA_DIR,
                  sprintf('study-%s-cell-%03d-replicates.rds',
                          study_id, cell_id))
  saveRDS(cell_dt, fn)
  invisible(fn)
}

save_summary_09 <- function(summary_dt, study_id) {
  fn_rds <- file.path(P9_DATA_DIR,
                      sprintf('study-%s-summary.rds', study_id))
  fn_txt <- file.path(P9_DATA_DIR,
                      sprintf('study-%s-summary.txt', study_id))
  saveRDS(summary_dt, fn_rds)
  write.table(summary_dt, file = fn_txt,
              sep = '\t', quote = FALSE, row.names = FALSE)
  cat(sprintf('Wrote %s and %s.\n', fn_rds, fn_txt))
  invisible(c(fn_rds, fn_txt))
}
