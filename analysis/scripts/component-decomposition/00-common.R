#' 00-common.R
#'
#' Shared infrastructure for paper 06 (component-decomposition)
#' production drivers. Sourced by the per-study scripts (02, 03, 04)
#' and by run-all.R. Lifts the analysis wrappers, parameter
#' constructors, and chunked-checkpoint run loop from the
#' 01-prototype-quick-sim.R pilot, generalises parameters over the
#' Study A grid (m_PB x m_TV x N), and adds a stub for the full
#' nonlinear three-component fit.
#'
#' Outputs land under analysis/data/derived_data/component-decomposition/
#' (one rds per cell plus a study-level Morris ADEMP summary txt).
#'
#' Naming: pmsimstats short-form variables (cat, max, disp, rate,
#' sd; bm, BL); analysis short-form options (useDE, t_random_slope,
#' simplecarryover, carryover_t1half, carryover_scalefactor).

suppressPackageStartupMessages({
  library(devtools)
  load_all('.', quiet = TRUE)
  library(data.table)
  library(nlme)
  library(parallel)
})

#=====================================================================
# Paths
#=====================================================================

P6_DATA_DIR <- file.path('analysis', 'data', 'derived_data',
                         'component-decomposition')
P6_FIG_DIR  <- file.path('analysis', 'figures',
                         'component-decomposition')
if (!dir.exists(P6_DATA_DIR)) dir.create(P6_DATA_DIR, recursive = TRUE)
if (!dir.exists(P6_FIG_DIR))  dir.create(P6_FIG_DIR,  recursive = TRUE)

#=====================================================================
# Trial design: Hendrickson hybrid OL+BDC
#=====================================================================
#
# Pre-registration calls for an 8-week OL + 4-week BDC + 4-week
# crossover (16 weeks). The 01-prototype uses a simplified 8-visit
# OL+BDC at weekly density (no crossover phase). For Study A's bias
# estimand the simplified design is sufficient; Study B's
# identifiability question and Study C's subadditivity question need
# the full 16-week design with the crossover phase to be defensible.
#
# We expose both via design_hybrid_simple() and design_hybrid_full().
# Study A defaults to the simple design (matching the prototype);
# Studies B and C use the full design.

design_hybrid_simple <- function() {
  buildtrialdesign(
    name_longform  = 'OL + BDC (simple)',
    name_shortform = 'OLBDC',
    timepoints   = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20),
    timeptnames  = c(paste0('OL', 1:5), paste0('BDC', 1:3)),
    expectancies = c(rep(1, 5), rep(0.5, 3)),
    ondrug = list(pathA = c(rep(1, 5), 1, 0, 0))
  )
}

design_hybrid_full <- function() {
  ## 8 OL weekly + 4 BDC weekly + 4 crossover-1 weekly + 4 crossover-2
  ## weekly = 16 visits. Two paths: AB (on first then off second) and
  ## BA (off first then on second).
  buildtrialdesign(
    name_longform  = 'OL + BDC + crossover (full)',
    name_shortform = 'OLBDCxo',
    timepoints   = seq(1, 20, by = 1.25),
    timeptnames  = c(paste0('OL', 1:8), paste0('BDC', 1:4),
                     paste0('XO1_', 1:2), paste0('XO2_', 1:2)),
    expectancies = c(rep(1, 8), rep(0.5, 4), rep(0.5, 4)),
    ondrug = list(
      AB = c(rep(1, 8), rep(0, 4), rep(1, 2), rep(0, 2)),
      BA = c(rep(1, 8), rep(0, 4), rep(0, 2), rep(1, 2))
    )
  )
}

#=====================================================================
# Parameter constructors
#=====================================================================

#' Response-parameter table parameterised over (m_PB, m_TV).
#'
#' @param m_PB Population mean of the placebo-belief modGompertz
#'   amplitude; pre-registration grid {0, 1, 3, 6, 10}. Reproduced
#'   here on the same scale as the prototype's pb_max axis (the
#'   prototype used 0, 6.5, 13).
#' @param m_TV Population mean of the natural-history modGompertz
#'   amplitude; pre-registration grid {-1, 0, 1, 2}. Negative values
#'   represent placebo-arm worsening of natural history.
make_resp_params <- function(m_PB = 6.5, m_TV = 0,
                              m_BR = 10.98604) {
  data.table(
    cat  = c('tv', 'pb', 'br'),
    max  = c(m_TV, m_PB, m_BR),
    disp = c(5, 5, 5),
    rate = c(0.42, 0.35, 0.42),
    sd   = c(5, 2, 5)
  )
}

make_bl_params <- function() {
  data.table(
    cat = c('bm', 'BL'),
    m   = c(0, 70),
    sd  = c(1, 10)
  )
}

make_model_params <- function(N = 70, c_bm = 0.45, t1half = 1) {
  data.table(
    N = N, c.bm = c_bm,
    carryover_t1half = t1half,
    c.tv = 0.7, c.pb = 0.7, c.br = 0.7,
    c.cf1t = 0.1, c.cfct = 0.05
  )
}

#=====================================================================
# True-effect helper (Architecture B, MVN)
#=====================================================================
# bm:Dbc on Sx = BL - (BR + PB + TV) is -c.bm * sigma_BR / sigma_bm.
# At c.bm = 0.45, sigma_BR = 5, sigma_bm = 1: TRUE_BETA = -2.25.
TRUE_BETA <- function(c_bm = 0.45, sigma_br = 5, sigma_bm = 1) {
  -c_bm * sigma_br / sigma_bm
}

#=====================================================================
# Analysis wrappers (ported from 01-prototype-quick-sim.R)
#=====================================================================

one_component_fit <- function(td, dat) {
  op <- list(useDE = FALSE, t_random_slope = FALSE,
             full_model_out = FALSE, simplecarryover = FALSE,
             carryover_t1half = 1, carryover_scalefactor = 1)
  res <- tryCatch(lme_analysis(td$trialpaths, dat, op),
                  error = function(e) NULL)
  if (is.null(res)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'fit_failed'))
  }
  data.table(beta            = res$beta,
             betaSE          = res$betaSE,
             p               = res$p,
             converged       = !is.na(res$beta),
             formula_dropped = 'none')
}

phase_augmented_fit <- function(td, dat, N) {
  op_prep <- list(useDE = FALSE, t_random_slope = FALSE,
                  full_model_out = TRUE, simplecarryover = FALSE,
                  carryover_t1half = 1, carryover_scalefactor = 1)
  prep <- tryCatch(lme_analysis(td$trialpaths, dat, op_prep),
                   error = function(e) NULL)
  if (is.null(prep) || is.null(prep$datamerged)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'prep_failed'))
  }
  d <- copy(prep$datamerged)
  d <- d[!is.na(Sx) & !is.na(bm) & !is.na(t) & !is.na(Dbc) & !is.na(De)]
  d <- d[t > 0]
  d[, phase := factor(ifelse(De >= 0.99, 'OL', 'BDC'),
                      levels = c('OL', 'BDC'))]
  if (length(unique(d$phase)) < 2) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'single_phase'))
  }

  ctl <- nlme::lmeControl(opt = 'optim', maxIter = 200, msMaxIter = 200)
  formula_dropped <- 'none'
  fit <- NULL

  if (N >= 70) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + Dbc:phase
                     + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL)
    if (is.null(fit)) formula_dropped <- 'Dbc:phase'
  } else {
    formula_dropped <- 'Dbc:phase'
  }

  if (is.null(fit)) {
    fit <- tryCatch(
      nlme::lme(Sx ~ bm + t + Dbc + phase + bm:Dbc + bm:Dbc:phase,
                random = ~1 | ptID,
                correlation = nlme::corCAR1(form = ~t | ptID),
                data = d, control = ctl),
      error = function(e) NULL)
  }
  if (is.null(fit)) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+all')))
  }
  ct <- summary(fit)$tTable
  target <- intersect(c('bm:Dbc', 'Dbc:bm'), rownames(ct))
  if (length(target) == 0) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = paste0(formula_dropped, '+bmDbc')))
  }
  data.table(beta            = unname(ct[target[1], 'Value']),
             betaSE          = unname(ct[target[1], 'Std.Error']),
             p               = unname(ct[target[1], 'p-value']),
             converged       = TRUE,
             formula_dropped = formula_dropped)
}

#' Full nonlinear three-component fit (STUB).
#'
#' Per the pre-registration this fits an nlme::nlme model of
#' BL - (modGompertz_BR + modGompertz_PB + modGompertz_TV) with
#' participant-specific random effects on each modGompertz amplitude.
#' Convergence is fragile at small N; the pre-registration restricts
#' this analysis to N >= 100. The closed-form gradient and starting-
#' value strategy are non-trivial and not yet fully implemented. This
#' stub returns NA values and marks formula_dropped = 'nlme_stub' so
#' downstream summaries flag it.
#'
#' TODO: implement the nlme::nlme call with self-starting modGompertz
#' (pmsimstats::utilities exposes modGompertz; gradient via deriv()),
#' Pinheiro-Bates pdSymm/pdDiag random-effect specification, and
#' starting values from the participant-specific OL-phase end-of-run
#' Sx changes. Track convergence-rate target >= 80% per pre-reg.
full_nonlinear_fit <- function(td, dat, N) {
  if (N < 100) {
    return(data.table(beta = NA_real_, betaSE = NA_real_,
                      p = NA_real_, converged = FALSE,
                      formula_dropped = 'N_below_100'))
  }
  data.table(beta = NA_real_, betaSE = NA_real_,
             p = NA_real_, converged = FALSE,
             formula_dropped = 'nlme_stub')
}

#=====================================================================
# One replicate
#=====================================================================

#' @param rep_idx integer replicate index.
#' @param cell named list with at least m_PB, m_TV, N, c_bm, analyses
#'   (character vector subset of c('one_component', 'phase_augmented',
#'   'full_nonlinear')); optional design (default 'simple').
#' @param study_seed master seed for the study.
#' @param cell_id integer cell identifier (used in seed derivation).
one_rep_06 <- function(rep_idx, cell, study_seed, cell_id) {
  set.seed(study_seed + 1000L * cell_id + rep_idx)
  td <- if (!is.null(cell$design) && cell$design == 'full')
          design_hybrid_full() else design_hybrid_simple()
  mp <- make_model_params(N = cell$N, c_bm = cell$c_bm,
                          t1half = 1)
  rp <- make_resp_params(m_PB = cell$m_PB, m_TV = cell$m_TV)
  bp <- make_bl_params()

  ## Generate data: union of paths if multiple. Use first path
  ## (or AB) for simplicity; pmsimstats handles multi-path data.
  paths <- td$trialpaths
  dat_list <- vector('list', length(paths))
  for (g in seq_along(paths)) {
    di <- tryCatch(
      generateData(mp, rp, bp, paths[[g]],
                   empirical = FALSE, makePositiveDefinite = TRUE,
                   dgp_architecture = 'mvn'),
      error = function(e) NULL)
    if (is.null(di)) next
    di[, path := g]
    dat_list[[g]] <- di
  }
  good <- !vapply(dat_list, is.null, logical(1))
  if (!any(good)) {
    na_row <- data.table(beta = NA_real_, betaSE = NA_real_,
                         p = NA_real_, converged = FALSE,
                         formula_dropped = 'gen_failed')
    return(rbindlist(lapply(cell$analyses, function(an)
      data.table(analysis = an, rep_idx = rep_idx, na_row))))
  }
  dat <- rbindlist(dat_list[good], fill = TRUE)

  out <- vector('list', 0L)
  for (an in cell$analyses) {
    fit <- switch(an,
      one_component   = one_component_fit(td, dat),
      phase_augmented = phase_augmented_fit(td, dat, cell$N),
      full_nonlinear  = full_nonlinear_fit(td, dat, cell$N),
      data.table(beta = NA_real_, betaSE = NA_real_,
                 p = NA_real_, converged = FALSE,
                 formula_dropped = 'unknown_analysis'))
    out[[length(out) + 1L]] <- cbind(analysis = an,
                                     rep_idx = rep_idx, fit)
  }
  rbindlist(out, fill = TRUE)
}

#=====================================================================
# Cell-level execution with chunked checkpoint loop
#=====================================================================

#' Run all replicates for one cell with parallel mclapply and a
#' chunked progress/checkpoint pattern. Returns the per-replicate
#' data.table for the cell.
run_cell_06 <- function(cell, n_reps, study_seed, cell_id,
                         n_cores = max(1L,
                                       parallel::detectCores() - 1L),
                         chunk_size = 50L,
                         progress_every = 100L,
                         time_budget_sec = Inf) {
  t0 <- Sys.time()
  elapsed <- function()
    as.numeric(difftime(Sys.time(), t0, units = 'secs'))

  cell_str <- paste0('cell ', cell_id, ' (m_PB=', cell$m_PB,
                     ', m_TV=', cell$m_TV, ', N=', cell$N,
                     ', c_bm=', cell$c_bm, ')')
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
      function(r) one_rep_06(r, cell, study_seed, cell_id),
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
  out[, `:=`(m_PB = cell$m_PB, m_TV = cell$m_TV,
             N = cell$N, c_bm = cell$c_bm,
             cell_id = cell_id)]
  setcolorder(out, c('cell_id', 'm_PB', 'm_TV', 'N', 'c_bm',
                     'analysis', 'rep_idx',
                     'beta', 'betaSE', 'p',
                     'converged', 'formula_dropped'))
  cat(sprintf('[%6.1fs] %s: complete (%d reps).\n',
              elapsed(), cell_str, done))
  out
}

#=====================================================================
# Morris ADEMP cell summary
#=====================================================================

summarise_cell_06 <- function(d, true_beta) {
  conv <- d[converged == TRUE]
  n <- nrow(conv)
  if (n == 0) {
    return(data.table(
      n_reps = nrow(d), n_converged = 0L, conv_rate = 0,
      power = NA_real_, mcse_power = NA_real_,
      mean_beta = NA_real_, sd_beta = NA_real_,
      mcse_mean_beta = NA_real_,
      bias = NA_real_, mcse_bias = NA_real_,
      coverage = NA_real_, mcse_coverage = NA_real_))
  }
  power <- mean(conv$p < 0.05, na.rm = TRUE)
  mcse_power <- sqrt(power * (1 - power) / n)
  mean_beta <- mean(conv$beta, na.rm = TRUE)
  sd_beta   <- sd(conv$beta,   na.rm = TRUE)
  mcse_mean_beta <- sd_beta / sqrt(n)
  bias <- mean_beta - true_beta
  mcse_bias <- mcse_mean_beta
  ## Wald 95% CI coverage
  lo <- conv$beta - 1.96 * conv$betaSE
  hi <- conv$beta + 1.96 * conv$betaSE
  cov <- mean(true_beta >= lo & true_beta <= hi, na.rm = TRUE)
  mcse_cov <- sqrt(cov * (1 - cov) / n)
  data.table(n_reps = nrow(d), n_converged = n,
             conv_rate = mean(d$converged, na.rm = TRUE),
             power = power, mcse_power = mcse_power,
             mean_beta = mean_beta, sd_beta = sd_beta,
             mcse_mean_beta = mcse_mean_beta,
             bias = bias, mcse_bias = mcse_bias,
             coverage = cov, mcse_coverage = mcse_cov)
}

#=====================================================================
# Save helpers
#=====================================================================

save_cell_06 <- function(cell_dt, study_id, cell_id) {
  fn <- file.path(P6_DATA_DIR,
                  sprintf('study-%s-cell-%03d-replicates.rds',
                          study_id, cell_id))
  saveRDS(cell_dt, fn)
  invisible(fn)
}

save_summary_06 <- function(summary_dt, study_id) {
  fn_rds <- file.path(P6_DATA_DIR,
                      sprintf('study-%s-summary.rds', study_id))
  fn_txt <- file.path(P6_DATA_DIR,
                      sprintf('study-%s-summary.txt', study_id))
  saveRDS(summary_dt, fn_rds)
  write.table(summary_dt, file = fn_txt,
              sep = '\t', quote = FALSE, row.names = FALSE)
  cat(sprintf('Wrote %s and %s.\n', fn_rds, fn_txt))
  invisible(c(fn_rds, fn_txt))
}
